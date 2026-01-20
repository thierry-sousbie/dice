#ifndef __DICE_FFTW_INTERFACE_HXX__
#define __DICE_FFTW_INTERFACE_HXX__

#include <stdio.h>

#include <complex>
#include <complex.h>
#include <fftw3.h>

#include "./fftw3-mpi_dummy.h"

#include "../dice_globals.hxx"

#include "FFTWTools.hxx"
#include "FFTWKernelFunctorInterface.hxx"

#include "./internal/regularGridSlicer_FFTW.hxx"

/**
 * @file
 * @brief  An interface to FFTW for MPI shared regular grids that can be used to
 *  apply convolution kernels in fourier space
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup FFT
 *   \{
 */

/**
 * \class FFTWConvolverT
 * \brief An interface to FFTW for MPI shared regular grids that takes care of local grids
 * allocation, FFT plans computation and applying convolution kernels in fourier space
 * \tparam G An MPI shared regular grid, typically RegularGridT ?
 */
template <class G>
class FFTWConvolverT
{
public:
  typedef FFTWConvolverT<G> MyType;
  typedef G Grid;
  typedef typename G::Data Data;
  typedef typename G::Params GridParams;

  typedef std::complex<Data> Complex;

  static const long BOUNDARY_TYPE = G::BOUNDARY_TYPE;
  static const int NDIM = G::NDIM;

  /** regular grid slicer for FFTW */
  friend class internal::RegularGridSlicerFFTWT<G>;
  /** A regular grid slicer that load balances as expected by FFTW*/
  typedef internal::RegularGridSlicerFFTWT<G> Slicer;
  typedef FFTWToolsT<Data, NDIM> FFTWTools;
  typedef FFTWKernelFunctorInterfaceT<MyType> KernelFunctorInterface;

  //! constructor
  FFTWConvolverT() : grid(NULL),
                     gridFacade(NULL),
                     kernelFunctorInterface(NULL),
                     nForwardPlans(0),
                     nBackwardPlans(0),
                     // kernel(NULL),
                     temp(NULL),
                     ownTemp(false)
  {
    STATIC_ASSERT_ERROR_TEMPLATE_GRID_FIELD_LAYOUT_MUST_BE_CONSECUTIVE(typename hlp::IsTrueT<G::IS_INTERLEAVED>::Result());
  }

  ~FFTWConvolverT()
  {
    freePlans();
    freeKernels();
    if (gridFacade != grid)
      delete grid;
  }

  /** \brief returns the FFTW flag corresponding to the wisdom level passed as argument.
   *  values in [0..3], in growing order of optimization, correspond to a level
   *  of FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT and FFTW_EXHAUSTIVE respectively.
   */
  static int getFFTW_Wisdom(int level)
  {
    if (level == 0)
      return FFTW_ESTIMATE;
    else if (level == 1)
      return FFTW_MEASURE;
    else if (level == 2)
      return FFTW_PATIENT;
    else if (level == 3)
      return FFTW_EXHAUSTIVE;
    else
      return FFTW_MEASURE;
  }

  /** \brief initialize FFTW routines. This needs to be called once before any FFTW
   *  function is called. Calling more than once has no effect, and initializeGridAndPlans
   *  does take care of initializing FFTW for you. However, you are responsible for making
   *  sure it was called at least once if you plan on using FFTWConvolver::Slicer directly.
   */
  static void initializeFFTW()
  {
#ifdef HAVE_FFTW3_THREADS

    fftw_init_threads();
#ifdef HAVE_FFTW3_MPI
    fftw_mpi_init();
#endif // HAVE_FFTW3_MPI

#else // HAVE_FFTW3_THREADS

#ifdef HAVE_FFTW3_MPI
    fftw_mpi_init();
#else  // HAVE_FFTW3_MPI
    fftw_init();
#endif // HAVE_FFTW3_MPI

#endif // HAVE_FFTW3_THREADS
  }

  /** \brief Export acquired FFTW wisdom from file 'filename'
   */
  int exportWisdom(const char *filename, bool quiet = false)
  {
    if (!quiet)
      glb::console->print<LOG_STD>("Exporting FFTW wisdom to file '%s'.\n", filename);
    return fftw_export_wisdom_to_filename(filename);
  }

  /** \brief Import FFTW wisdom from file 'filename'. File 'filename' may not exist,
   *   in which case a simple warning is issued and 0 is returned.
   *   \note calling this function will initialize FFTW
   */
  int importWisdom(const char *filename, bool quiet = false)
  {
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
      if (!quiet)
      {
        // glb::console->print<LOG_WARNING>
        //("Unable to open file '%s' for reading.",filename);
        glb::console->print<LOG_STD>("Cannot import FFTW wisdom: unable to open file '%s' for reading.\n", filename);
      }
      return 0;
    }
    else
    {
      fclose(f);
      if (!quiet)
        glb::console->print<LOG_STD>("Importing FFTW wisdom from file '%s'.\n", filename);
    }

    do_initializeFFTW();
    return fftw_import_wisdom_from_filename(filename);
  }

  /** \brief Initialize the convolver and the grid: create the FFT plans, load-balance
   *  and allocate the grid and build the kernel(s). Note that \a g should not have been
   *  initialized before (or after) calling this function (this fucntion takes care of it)
   * \param g The regular shared grid to intialize. All convolutions will be
   *  performed on that grid for every call to execute(). It should not be intialized
   *  before (or after) calling this function.
   * \param gp The grid parameters of type G::Params of the global grid.
   * \param kernelFunctor a functor that defines the kernel. Its class must implement the
   * interface defined by KernelFunctorInterface (see FFTWKernelFunctorInterfaceT)
   * \param wisdomLevel specifies how thorough the optimization of FFTW plan creation
   *  should be. Values belong to [0..3], in growing order of optimization, corresponding
   *  to a level of FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT and FFTW_EXHAUSTIVE
   *  respectively.
   * \param storeKernels if storeKernels is true, the fourier kernel is always stored even
   * when this is not necessary (i.e. when the functor is purely defined in fourier space).
   * Storing the kernel results in a larger memory footprint (one full grid must be
   * stored for each kernel output field) but as a result convolution is faster ...
   * \param com the MPI communicator to use for load balancing
   * \param nThreads The number of openMP threads to use for FFTs
   * \param quiet Do not write anything to console if true
   * \tparam K the class of the kernel functor, which must publicly inherit from
   * KernelFunctorInterface (see FFTWKernelFunctorInterfaceT)
   * \warning This function takes care of initializing the grid \a g, so the grid itself
   * does not need to be initialized beforehand (i.e. no need to call Grid::initialize)
   */
  template <class K>
  void initializeGridAndPlans(Grid *g,
                              const GridParams &gp,
                              const K &kernelFunctor,
                              int wisdomLevel = 0,
                              bool storeKernels = true,
                              MpiCommunication *com = glb::mpiComWorld,
                              int nThreads = glb::num_omp_threads,
                              bool quiet = false)
  {
    int wisdom = getFFTW_Wisdom(wisdomLevel);

    if (!quiet)
    {
      if (glb::console->willPrint<LOG_INFO>())
      {
        char wisdomStr[64];
        if (wisdom == FFTW_ESTIMATE)
          sprintf(wisdomStr, "FFTW_ESTIMATE");
        else if (wisdom == FFTW_MEASURE)
          sprintf(wisdomStr, "FFTW_MEASURE");
        else if (wisdom == FFTW_PATIENT)
          sprintf(wisdomStr, "FFTW_PATIENT");
        else if (wisdom == FFTW_EXHAUSTIVE)
          sprintf(wisdomStr, "FFTW_EXHAUSTIVE");
        glb::console->printFlush<LOG_INFO>("Initializing FFTW convolver (%s):\n", wisdomStr);
        glb::console->indent();
      }
      else
        glb::console->printFlush<LOG_STD>("Initializing FFTW convolver ... ");
    }

    if (!quiet)
      glb::console->printFlush<LOG_INFO>(" * Initializing FFTW ... ");

    if ((grid != NULL) && (grid != gridFacade))
      delete grid;
    gridFacade = g;
    grid = g;

    mpiCom = com;
    numThreads = nThreads;
    wisdomFlags = wisdom;
    do_initializeFFTW();

#ifdef HAVE_FFTW3_THREADS
    fftw_plan_with_nthreads(numThreads);
#else
    if (numThreads > 1)
    {
      PRINT_SRC_INFO(LOG_WARNING);
      glb::console->print<LOG_WARNING>("%d threads are required, but FFTW3-threads is not available: only 1 thread will be used.\n", numThreads);
      glb::console->print<LOG_WARNING>("Recompile with FFTW3-threads enabled to fix this !\n");
    }
#endif

    if (!quiet)
      glb::console->printFlush<LOG_INFO>("done.\n");
    if (!quiet)
      glb::console->printFlush<LOG_INFO>(" * Initializing the grid ... ");

    Slicer slicer(gp, mpiCom, numThreads, periodic,
                  kernelFunctor.getOutputFieldsCount());

    if (slicer.needFacade())
      grid = new Grid();

    grid->initialize(slicer);

    if (slicer.toFacade(gp))
      gridFacade->initialize(slicer, grid->getLocalGrid()->getDataPtr());
    /*
    glb::console->print<LOG_INFO_ALL>("Grid :(periodic=%d)\n",periodic);
    grid->getGlobalGridParams().template print<LOG_INFO>();
    grid->getLocalGridParams().template print<LOG_INFO_ALL>();

    glb::console->print<LOG_INFO_ALL>("gridFacade :\n");
    //if (mpiCom->rank()==1) gridFacade->getGlobalGridParams().template print<LOG_INFO>();
    gridFacade->getLocalGridParams().template print<LOG_INFO_ALL>();
    */
    if (!quiet)
      glb::console->printFlush<LOG_INFO>("done.\n");

    createPlansAndKernels(kernelFunctor, storeKernels, quiet);

    if (!quiet)
    {
      if (glb::console->willPrint<LOG_INFO>())
      {
        glb::console->unIndent();
        glb::console->print<LOG_INFO>("All done.\n");
      }
      else
        glb::console->print<LOG_STD>("done.\n");
    }
  }

  /** \brief apply the kernel(s) in fourier space and transform back.
   *  \note Depending on the number of output fields in the convolution kernel, the number
   *  of fields in the grid may be changed after a call to execute. The number of fields
   *  in the grid can then be reset to its initial state by calling
   *  grid->reinterpretNFields().
   */
  void execute(bool quiet = false)
  {

    if (!quiet)
    {
      if (glb::console->willPrint<LOG_INFO>())
      {
        glb::console->print<LOG_INFO>("\n");
        glb::console->indent();
      }
      // else
      // glb::console->printFlush<LOG_STD>("Applying %s operator ... ");
    }

    if (grid != gridFacade)
    {
      if (!quiet)
        glb::console->printFlush<LOG_INFO>("Rearranging arrays ... ");
      Slicer::rearrangeFromFacadeToGrid(grid, gridFacade, mpiCom);
      if (!quiet)
        glb::console->print<LOG_INFO>("done.\n");
    }

    if (!quiet)
    {
      if (nForwardPlans > 1)
        glb::console->printFlush<LOG_INFO>("* Forward FFT (%d plans) ... ", nForwardPlans);
      else
        glb::console->printFlush<LOG_INFO>("* Forward FFT ... ");
    }

    // forward fourier transform first
    for (int i = 0; i < nForwardPlans; ++i)
    {
      // add padding for r2c transform
      FFTWTools::rearrangeToFFTW(grid->getDataPtr(), grid->getDataPtr(),
                                 info.localDims, forwardPlanIsPadded[i], periodic);
      fftw_execute(forwardPlan[i]);
    }

    if (!quiet)
      glb::console->printFlush<LOG_INFO>("done.\n");

    // Then apply operator and backward transform
    for (int i = 0; i < nBackwardPlans; ++i)
    {
      if (!quiet)
        glb::console->printFlush<LOG_INFO>("* Applying kernel + backward FFT (%d/%d) ... ", i + 1, nBackwardPlans);

      Data *destPtr = grid->getDataPtr(0, i); // + info.nLocalElements*i;
      fftw_complex *srcPtr;

      if (i < nBackwardPlans - 1)
      {
        // The NDIM-1 first backward transforms are done inplace
        memcpy(destPtr, temp, sizeof(fftw_complex) * info.fftAlloc); // nFourierElements);
        srcPtr = (fftw_complex *)destPtr;
      }
      else
        srcPtr = temp;

      applyKernel(srcPtr, i);

      fftw_execute(backwardPlan[i]);

      FFTWTools::rearrangeFromFFTW(destPtr, destPtr, info.localDims, backwardPlanIsPadded[i], periodic);

      if (!quiet)
        glb::console->printFlush<LOG_INFO>("done.\n");
    }

    if (grid != gridFacade)
    {
      if (!quiet)
        glb::console->printFlush<LOG_INFO>("Rearranging arrays ... ");
      Slicer::rearrangeFromGridToFacade(grid, gridFacade, mpiCom);
      if (!quiet)
        glb::console->print<LOG_INFO>("done.\n");
    }

    grid->reinterpretNFields(outputFieldsCount);

    if (!quiet)
    {
      if (glb::console->willPrint<LOG_INFO>())
      {
        glb::console->unIndent();
        // glb::console->print<LOG_INFO>("Done.\n");
      }
      // else
      // glb::console->print<LOG_STD>("done.\n");
    }
  }

private:
  static const int periodic = (BOUNDARY_TYPE == BoundaryType::PERIODIC);

  int initialized(int setVal = -1)
  {
    static int FFTW_initialized = 0;
    if (setVal >= 0)
      FFTW_initialized = setVal;

    return FFTW_initialized;
  }

  void do_initializeFFTW()
  {
    if (initialized())
      return;

    initialized(true);
    initializeFFTW();
  }

  void freePlans()
  {
    if ((temp != NULL) && (ownTemp))
    {
      fftw_free(temp);
      temp = NULL;
    }

    for (int i = 0; i < nForwardPlans; ++i)
      fftw_destroy_plan(forwardPlan[i]);
    for (int i = 0; i < nBackwardPlans; ++i)
      fftw_destroy_plan(backwardPlan[i]);

    nForwardPlans = 0;
    nBackwardPlans = 0;
  }

  void freeKernels()
  {
    for (long i = 0; i < kernel.size(); ++i)
      fftw_free(kernel[i]);
    kernel.clear();
    if (kernelFunctorInterface != NULL)
      delete kernelFunctorInterface;
    kernelFunctorInterface = NULL;
  }

  template <class K>
  void createPlansAndKernels(const K &kernelFunctor, bool storeKernels, bool quiet)
  {
    info.set(grid, mpiCom);
    outputFieldsCount = kernelFunctor.getOutputFieldsCount();
    createPlans(kernelFunctor.getOutputFieldsCount(), quiet);
    createKernels(kernelFunctor, storeKernels, quiet);
  }

  void createPlans(int count, bool quiet)
  {
    Data *rPtr;
    fftw_complex *cPtr;

    freePlans();

    if (!quiet)
      glb::console->printFlush<LOG_INFO>(" * Setting up FFT plans ... (F)");
    nForwardPlans = 1;

    // We need a temporary array if transforming to more than 1 field
    if (count > 1)
    {
      temp = (fftw_complex *)fftw_alloc_complex(info.fftAlloc);
      ownTemp = true;
    }
    else
    {
      temp = (fftw_complex *)grid->getDataPtr();
      ownTemp = false;
    }

    // Create the forward transform, from the grid to temp
    rPtr = grid->getDataPtr();
    cPtr = temp;
    forwardPlan[0] = fftw_mpi_plan_dft_r2c(NDIM, info.dimsR, rPtr, cPtr, mpiCom->getCom(),
                                           wisdomFlags | FFTW_DESTROY_INPUT |
                                               FFTW_MPI_TRANSPOSED_OUT);

    forwardPlanIsPadded[0] = fftwNeedsPadding(rPtr, cPtr);

    nBackwardPlans = count;
    // The derivative along the first dimensions is computed inplace, in the grid
    for (int i = 0; i < count - 1; ++i)
    {
      if (!quiet)
        glb::console->printFlush<LOG_INFO>("(B%d)", i + 1);
      rPtr = grid->getDataPtr(0, i); //+(info.nLocalElements*i);
      cPtr = (fftw_complex *)rPtr;
      backwardPlan[i] = fftw_mpi_plan_dft_c2r(NDIM, info.dimsR, cPtr, rPtr, mpiCom->getCom(),
                                              wisdomFlags | FFTW_DESTROY_INPUT |
                                                  FFTW_MPI_TRANSPOSED_IN);
      backwardPlanIsPadded[i] = fftwNeedsPadding(rPtr, cPtr);
    }

    if (!quiet)
      glb::console->printFlush<LOG_INFO>("(B%d) ", count);
    // The last derivative is computed out of place from temp to the grid
    rPtr = grid->getDataPtr(0, count - 1); //+(info.nLocalElements*(count-1));
    cPtr = temp;
    backwardPlan[count - 1] = fftw_mpi_plan_dft_c2r(NDIM, info.dimsR, cPtr, rPtr, mpiCom->getCom(),
                                                    wisdomFlags | FFTW_DESTROY_INPUT |
                                                        FFTW_MPI_TRANSPOSED_IN);
    backwardPlanIsPadded[count - 1] = fftwNeedsPadding(rPtr, cPtr);

    if (!quiet)
      glb::console->printFlush<LOG_INFO>("done.\n");
  }

  template <class K>
  void createKernels(const K &kernelFunctor, bool storeKernels, bool quiet = false)
  {
    static const double twopi = 6.283185307179586232L;
    // static const double pi =    3.141592653589793238L;
    int nKernels = kernelFunctor.getOutputFieldsCount();
    freeKernels();

    // We need to precompute kernels only if they have a real space part, in which
    // case we have to store their fourier transforms. If this is not the case, we
    // simply store the functor and regenerate the kernel for every convolution.
    if ((!storeKernels) && (!(K::SPACE & K::Interface::REAL_SPACE)))
    {
      if (!quiet)
        glb::console->print<LOG_INFO>(" * Creating %d kernels (%s) ... skipped (pure fourier kernel).\n",
                                      nKernels,
                                      kernelFunctor.getKernelName().c_str());
      K *tmp = new K;
      (*tmp) = kernelFunctor;
      kernelFunctorInterface = static_cast<KernelFunctorInterface *>(tmp);
      kernel.assign(outputFieldsCount, NULL);
      return;
    }

    long p0[NDIM] = {0};
    long k0[NDIM] = {0};
    long p[NDIM] = {0};
    long q[NDIM] = {0};
    fftw_plan plan;

    if (!quiet)
      glb::console->printFlush<LOG_INFO>(" * Creating %d kernels (%s) ... ",
                                         nKernels,
                                         kernelFunctor.getKernelName().c_str());

    for (int i = 0; i < NDIM; ++i)
      p0[i] = info.localPosition[i];
    for (int i = 0; i < NDIM; ++i)
      k0[i] = info.localPositionFT[i];

    for (int n = 0; n < nKernels; ++n)
    {
      if (!quiet)
        glb::console->printFlush<LOG_INFO>("(%d)", n);
      kernel.push_back((fftw_complex *)fftw_alloc_complex(info.fftAlloc));
      Complex *kernelF = reinterpret_cast<Complex *>(kernel.back());
      Data *kernelR = reinterpret_cast<Data *>(kernel.back());

      if (K::SPACE & K::Interface::REAL_SPACE)
      {
        double norm = pow(info.nElements, -2.0);
        double xNorm[NDIM];

        for (int i = 0; i < NDIM; ++i)
        {
          p[i] = 0;
          xNorm[i] = info.size[i] / info.dims[i];
        }

        plan = fftw_mpi_plan_dft_r2c(NDIM, info.dimsR, kernelR, kernel.back(),
                                     mpiCom->getCom(), FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);

        bool needPadding = fftwNeedsPadding(kernelR, kernel.back());

        double x[NDIM];
        for (long i = 0; i < info.localStride[NDIM]; ++i)
        {
          for (int j = 0; j < NDIM; ++j)
          {
            q[j] = p0[j] + p[j];
            x[j] = xNorm[j] * FFTWTools::indGen(q[j], info.dims[j]);
          }

          kernelR[i] = kernelFunctor.getReal(n, x, xNorm) * norm;
          hlp::getNext<NDIM>(p, info.localDims);
        }

        if (!quiet)
          glb::console->printFlush<LOG_INFO>("(fft)");
        FFTWTools::rearrangeToFFTW(kernelR, kernelR, info.localDims, needPadding, periodic);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        if (K::SPACE & K::Interface::FOURIER_SPACE)
        {
          double k[NDIM];
          double k1[NDIM];
          double kNorm[NDIM];
          double k1Norm[NDIM];

          for (int i = 0; i < NDIM; ++i)
          {
            p[i] = 0;
            kNorm[i] = twopi / (info.sizeT[i] * info.dimsT[i]);
            k1Norm[i] = 1.0L / (info.dimsT[i]);
          }

          for (long i = 0; i < info.localStrideFT[NDIM]; ++i)
          {
            for (int j = 0; j < NDIM; ++j)
              q[j] = k0[j] + p[j];

            fftwTransposedSwap(q[NDIM - 2], q[NDIM - 1]);

            for (int j = 0; j < NDIM; ++j)
            {
              k1[j] = FFTWTools::indGen(q[j], info.dimsT[j]);
              k[j] = k1[j] * kNorm[j];
              k1[j] = k1[j] * k1Norm[j];
            }

            kernelF[i] *= kernelFunctor.getFourier(n, k, k1);
            hlp::getNext<NDIM>(p, info.localDimsFT);
          }
        }
      }
      else if (K::SPACE & K::Interface::FOURIER_SPACE)
      {
        double norm = 1.0 / info.nElements;
        double kNorm[NDIM];
        double k1Norm[NDIM];
        double k[NDIM];
        double k1[NDIM];

        for (int i = 0; i < NDIM; ++i)
        {
          p[i] = 0;
          kNorm[i] = twopi / (info.sizeT[i]);
          k1Norm[i] = 1.0L / (info.dimsT[i]);
        }

        for (long i = 0; i < info.localStrideFT[NDIM]; ++i)
        {
          for (int j = 0; j < NDIM; ++j)
            q[j] = k0[j] + p[j];

          fftwTransposedSwap(q[NDIM - 2], q[NDIM - 1]);

          for (int j = 0; j < NDIM; ++j)
          {
            k1[j] = FFTWTools::indGen(q[j], info.dimsT[j]);
            k[j] = k1[j] * kNorm[j];
            k1[j] = k1[j] * k1Norm[j];
          }

          kernelF[i] = kernelFunctor.getFourier(n, k, k1) * Complex(norm);

          hlp::getNext<NDIM>(p, info.localDimsFT);
        }
      }
    }

    if (!quiet)
      glb::console->printFlush<LOG_INFO>(" done.\n");
  }

  void applyKernel(fftw_complex *data, int which)
  {
    Complex *d = reinterpret_cast<Complex *>(data);
    Complex *k = reinterpret_cast<Complex *>(kernel[which]);

    if (k != NULL)
    {
      for (long i = 0; i < info.localStrideFT[NDIM]; ++i)
        d[i] *= k[i];
    }
    else
      applyFunctorKernel(data, which);
  }

  void applyFunctorKernel(fftw_complex *data, int which)
  {
    static const double pi = acos(-1.0);
    static const double twopi = pi * 2; // 6.283185307179586232L;

    Complex *d = reinterpret_cast<Complex *>(data);
    double norm = 1.0L / info.nElements;
    long k0[NDIM] = {0};
    long p[NDIM] = {0};
    long q[NDIM] = {0};
    double kNorm[NDIM];
    double k1Norm[NDIM];
    double k[NDIM];
    double k1[NDIM];

    for (int i = 0; i < NDIM; ++i)
    {
      k0[i] = info.localPositionFT[i];
      p[i] = 0;
      kNorm[i] = twopi / (info.sizeT[i]);
      k1Norm[i] = 1.0L / (info.dimsT[i]);
    }

    for (long i = 0; i < info.localStrideFT[NDIM]; ++i)
    {
      for (int j = 0; j < NDIM; ++j)
        q[j] = k0[j] + p[j];

      fftwTransposedSwap(q[NDIM - 2], q[NDIM - 1]);

      for (int j = 0; j < NDIM; ++j)
      {
        k1[j] = FFTWTools::indGen(q[j], info.dimsT[j]);
        k[j] = k1[j] * kNorm[j];
        k1[j] = k1[j] * k1Norm[j];
      }

      d[i] *= kernelFunctorInterface->getFourier(which, k, k1) * norm;
      hlp::getNext<NDIM>(p, info.localDimsFT);
    }
  }

  struct GridInfo
  {
    ptrdiff_t fftAlloc;
    ptrdiff_t nElements;
    ptrdiff_t nLocalElements;
    ptrdiff_t dims[NDIM];
    ptrdiff_t dimsR[NDIM];
    ptrdiff_t dimsF[NDIM];
    ptrdiff_t dimsFR[NDIM];
    ptrdiff_t localDims[NDIM];
    ptrdiff_t localDimsF[NDIM];
    double size[NDIM];
    double sizeT[NDIM];
    double localSize[NDIM];
    ptrdiff_t localPosition[NDIM];
    ptrdiff_t stride[NDIM + 1];
    ptrdiff_t strideF[NDIM + 1];
    ptrdiff_t localStride[NDIM + 1];
    ptrdiff_t localStrideF[NDIM + 1];

    ptrdiff_t dimsT[NDIM];
    ptrdiff_t localDimsT[NDIM];
    ptrdiff_t localPositionT[NDIM];

    ptrdiff_t dimsFT[NDIM];
    ptrdiff_t localDimsFT[NDIM];
    ptrdiff_t localPositionFT[NDIM];
    ptrdiff_t localStrideFT[NDIM + 1];

    void set(Grid *g, MpiCommunication *com)
    {
      nElements = 1;
      nLocalElements = 1;
      stride[0] = 1;
      localStride[0] = 1;
      strideF[0] = 1;
      localStrideF[0] = 1;

      for (int i = 0; i < NDIM; ++i)
      {
        dims[i] = g->getResolution(i);
        dimsR[NDIM - i - 1] = dims[i];
        localDims[i] = g->getLocalResolution(i);
        nElements *= dims[i];
        nLocalElements *= localDims[i];
        size[i] = g->getSize(i);
        localSize[i] = g->getLocalSize(i);
        dimsF[i] = (i == 0) ? (dims[i] / 2 + 1) : dims[i];
        dimsFR[NDIM - i - 1] = dimsF[i];
        localDimsF[i] = (i == 0) ? (localDims[i] / 2 + 1) : localDims[i];

        stride[i + 1] = stride[i] * dims[i];
        localStride[i + 1] = localStride[i] * localDims[i];
        strideF[i + 1] = strideF[i] * dimsF[i];
        localStrideF[i + 1] = localStrideF[i] * localDimsF[i];
      }

      ptrdiff_t local_n0, local_0_start;
      ptrdiff_t local_n1, local_1_start;

      fftAlloc = fftw_mpi_local_size(NDIM, dimsFR,
                                     com->getCom(),
                                     &local_n0, &local_0_start);

      fftw_mpi_local_size_transposed(NDIM, dimsFR,
                                     com->getCom(),
                                     &local_n0, &local_0_start,
                                     &local_n1, &local_1_start);

      for (int i = 0; i < NDIM; ++i)
      {
        dimsT[i] = dims[i];
        localDimsT[i] = (i == NDIM - 1) ? (local_n0) : dims[i];
        localPositionT[i] = (i == NDIM - 1) ? (local_0_start) : 0; // g->getLocalPosition(i);

        dimsFT[i] = dimsF[i];
        localDimsFT[i] = (i == NDIM - 2) ? (local_n1) : dimsF[i];
        localPositionFT[i] = (i == NDIM - 2) ? (local_1_start) : 0;

        sizeT[i] = size[i];
      }

      fftwTransposedSwap(dimsT[NDIM - 1], dimsT[NDIM - 2]);
      fftwTransposedSwap(localDimsT[NDIM - 1], localDimsT[NDIM - 2]);
      fftwTransposedSwap(localPositionT[NDIM - 1], localPositionT[NDIM - 2]);

      fftwTransposedSwap(dimsFT[NDIM - 1], dimsFT[NDIM - 2]);
      fftwTransposedSwap(localDimsFT[NDIM - 1], localDimsFT[NDIM - 2]);
      fftwTransposedSwap(localPositionFT[NDIM - 1], localPositionFT[NDIM - 2]);

      fftwTransposedSwap(sizeT[NDIM - 1], sizeT[NDIM - 2]);

      localStrideFT[0] = 1;
      for (int i = 0; i < NDIM; ++i)
        localStrideFT[i + 1] = localStrideFT[i] * localDimsFT[i];
    }
  };

  // gridFacade is what the user will see. It may differ from grid for non-periodic
  // boundaries as we are padding ...
  Grid *grid;
  Grid *gridFacade;

  MpiCommunication *mpiCom;
  int numThreads;
  int wisdomFlags;
  int outputFieldsCount;

  KernelFunctorInterface *kernelFunctorInterface;

  fftw_plan forwardPlan[NDIM];
  int forwardPlanIsPadded[NDIM];
  int nForwardPlans;
  fftw_plan backwardPlan[NDIM];
  int backwardPlanIsPadded[NDIM];
  int nBackwardPlans;

  GridInfo info;

  std::vector<fftw_complex *> kernel;
  fftw_complex *temp;
  bool ownTemp;

private:
  // STATIC ASSERTS
  void STATIC_ASSERT_ERROR_TEMPLATE_GRID_FIELD_LAYOUT_MUST_BE_CONSECUTIVE(typename hlp::IsFalse isInterleaved) {}
};

/** \}*/
#include "../internal/namespace.footer"
#endif

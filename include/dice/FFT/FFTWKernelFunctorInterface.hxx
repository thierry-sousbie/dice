#ifndef __FFTW_KERNEL_FUNCTOR_INTERFACE_HXX__
#define __FFTW_KERNEL_FUNCTOR_INTERFACE_HXX__

/**
 * @file
 * @brief  defines a virtual interface for FFTW kernels (see FFTWConvolverT)
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup FFT
 *   \{
 */

/**
 * \class FFTWKernelFunctorInterfaceT
 * This class defines a virtual interface for FFTW kernels (see FFTWConvolverT). Any kernel
 * to be used with a convolver must implement it by publicly inheriting from it and also
 * defining a static const int variable named SPACE equal to a combination of SpaceType
 * enum values that defines in which space the kernel is defined (the appropriate functions
 * getFourier and getReal will be called depending on SPACE. For instance, if
 * SPACE=REAL_SPACE|FOURIER_SPACE, then the real part of the kernel will be created
 * in real space, fourier transformed, and then multiplied by its fourier part.
 */

template <class FFTC>
class FFTWKernelFunctorInterfaceT
{
public:
  typedef FFTC FFTWConvolver;
  typedef typename FFTWConvolver::FFTWTools FFTWTools;
  typedef typename FFTWConvolver::Complex Complex;
  static const int NDIM = FFTWConvolver::NDIM;

  enum SpaceType
  {
    REAL_SPACE = (1 << 0),
    FOURIER_SPACE = (1 << 1)
  };

  virtual ~FFTWKernelFunctorInterfaceT() {};

  /** \brief A function that returns the kernel value in fourier space at wave vectors \a k
   *  \param[in] index the index of the kernel to compute
   *  \param[in] k the wave vectors coordinates (correctly normalized to the physical box size)
   *  \param[in] k1 the wave vectors normalized to 1 along each direction (i.e. -0.5<=k1[i]<0.5)
   */
  virtual Complex getFourier(int index, const double k[NDIM], const double k1[NDIM]) const = 0;

  /** \brief A function that returns the value of the function to convolve in real space
   *  at position \a p
   *  \param[in] index the index of the kernel to compute
   *  \param[in] p the periodized coordinates within the grid
   *  \param[in] spacing the grid spacing
   */
  virtual double getReal(int index, const double p[NDIM], const double spacing[NDIM]) const = 0;
  /** \brief returns the number of fields this kernel outputs (e.g. this should be =NDIM
   * for a kernel that computes derivatives in fourier space)
   */
  virtual int getOutputFieldsCount() const = 0;
  /** \brief returns the name of the kernel as a string
   */
  virtual const std::string getKernelName() const = 0;
};

/** \}*/
#include "../internal/namespace.footer"
#endif

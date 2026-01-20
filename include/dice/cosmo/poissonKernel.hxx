#ifndef __POISSON_KERNEL_HXX__
#define __POISSON_KERNEL_HXX__

#include "../internal/namespace.header"

namespace cosmo
{

  enum PoissonKernelMode
  {
    pkPotential = 0,
    pkDisplacement = 1
  };

  template <class C, PoissonKernelMode Mode, bool Periodic = true> //, int FD_ORDER=-1>
  class PoissonKernelT : public C::KernelFunctorInterface
  {
  public:
    typedef typename C::KernelFunctorInterface Interface;

    // enum SpaceType {REAL_SPACE = (1<<0), FOURIER_SPACE=(1<<1)};

    typedef PoissonKernelT<C, Mode> MyType;
    typedef typename Interface::FFTWTools FFTWTools;
    typedef typename Interface::Complex Complex;

    static const int NDIM = Interface::NDIM;
    // static const int SPACE = Interface::REAL_SPACE;
    static const int SPACE = Interface::FOURIER_SPACE;
    // static const int SPACE = (Periodic)?(Interface::FOURIER_SPACE):(Interface::REAL_SPACE);
    static const PoissonKernelMode PK_MODE = Mode;

    PoissonKernelT(double factor_ = 1.0, double G_ = 1.0) : G(G_),
                                                            factor(factor_)
    {
      // const double pi=4.0*atan(1.0);
      // factor = -4.0*pi*G;
    }

    ~PoissonKernelT()
    {
    }

    // template <class T, class D>
    // T getReal(int index, D x[NDIM]) const
    double getReal(int index, const double x[NDIM], const double spacing[NDIM]) const
    {
      return getReal(index, x, spacing, hlp::ConstantValue<Mode>());
    }

    // template <class T, class D>
    // T getFourier(int index, D k[NDIM], D k1[NDIM]) const
    Complex getFourier(int index, const double k[NDIM], const double k1[NDIM]) const
    {
      return getFourier(index, k, k1, hlp::ConstantValue<Mode>());
    }

    int getOutputFieldsCount() const { return (Mode == pkDisplacement) ? NDIM : 1; }
    const std::string getKernelName() const { return std::string("Poisson kernel"); }

  private:
    double G;
    double factor;

    double getReal(int index, const double x[NDIM], const double spacing[NDIM],
                   hlp::ConstantValue<pkDisplacement>) const
    {
      printf("Error in poissonKernel: Cannot compute displacement for non periodic boundary conditions ! (NOT IMPLEMENTED)");
      exit(-1);
      /*
      double result;
      double d=0;
      for (int i=0;i<NDIM;++i)
  d+=x[i]*x[i];

      if (d==0)
  result = 1./sqrt(2.0/256.0);
      else
  result = -factor/sqrt(d);

      return result;
      */
    }

    Complex getFourier(int index, const double k[NDIM], const double k1[NDIM],
                       hlp::ConstantValue<pkDisplacement>) const
    {
      double result = 0;
      double laplacian = 0;
      double k1Norm2 = 0;

      for (int i = 0; i < NDIM; ++i)
      {
        laplacian += k[i] * k[i];
        k1Norm2 += k1[i] * k1[i];
      }

      if (laplacian != 0)
        result = -(k[index] * factor * FFTWTools::hanning(sqrt(k1Norm2))) / laplacian;

      // FiniteDifferenceFilter<FD_ORDER>(sqrt(k1Norm2)))/laplacian;
      // FFTWTools::hanning(sqrt(k1Norm2)))/laplacian;

      // result *= FFTWTools::hanning(sqrt(k1Norm2));

      // double *tmpd = reinterpret_cast<double*>(&tmp);
      // printf("result = (%lg,%lg) == (%lg,%lg) \n",0.0,result*k[index],tmpd[0],tmpd[1]);

      return Complex(0, result);
    }

    double getReal(int index, const double x[NDIM], const double spacing[NDIM],
                   hlp::ConstantValue<pkPotential>) const
    {
      double result;
      double d = 0;
      for (int i = 0; i < NDIM; ++i)
        d += x[i] * x[i];

      if (d == 0)
      {
        double h = spacing[0];
        for (int i = 1; i < NDIM; ++i)
          if (spacing[i] < h)
            h = spacing[i];
        result = -factor / h;
      }
      else
        result = -factor / sqrt(d);

      return result;
    }

    Complex getFourier(int index, const double k[NDIM], const double k1[NDIM],
                       hlp::ConstantValue<pkPotential>) const
    {
      double result = 0;
      double laplacian = 0;
      // double k1Norm2=0;
      for (int i = 0; i < NDIM; ++i)
      {
        laplacian += k[i] * k[i];
        // k1Norm2+=k1[i]*k1[i];
      }

      if (laplacian != 0)
        result = -factor / laplacian;

      return Complex(result, 0);
    }

    /*
    template <int ORDER>
    static double FiniteDifferenceFilter(double k1Norm)
    {
      return FiniteDifferenceFilter(k1Norm, hlp::ConstantValue<ORDER>());
    }

    // (Hanning)
    static double FiniteDifferenceFilter(double k1Norm, hlp::ConstantValue<-1>)
    {
      static const double twopi=atan(1.0L)*8.0L;
      return (k1Norm>0.5)?0:(0.5*(1.0+cos(twopi*k1Norm)));
    }

    // NO filter
    static double FiniteDifferenceFilter(double k1Norm, hlp::ConstantValue<0>)
    {
      return 1;
    }

    // FIXME : COMPUTE THE OPERATORS (Hockney & Eastwood, p284)
    // (2^d+1) points FD (Hanning)
    static double FiniteDifferenceFilter(double k1Norm, hlp::ConstantValue<1>)
    {

    }

    // (4^d+1) points FD
    static double FiniteDifferenceFilter(double k1Norm, hlp::ConstantValue<2>)
    {

    }
    */
  };
} // namespace cosmo

#include "../internal/namespace.footer"
#endif

#ifndef __SPECIALIZED_GAUSS_QUADRATURE_CUBOID1D_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_CUBOID1D_HXX__

#include "../../internal/namespace.header"

namespace internal
{
  namespace gaussQuadrature
  {

    template <>
    struct SpecializedGQT<1, 0, fETypeE::cuboid>
    {
      static const int NP = 1;

      template <class F>
      static double get(F &functor)
      {
        return functor(0);
      }
    };

    template <>
    struct SpecializedGQT<1, 1, fETypeE::cuboid>
    {
      static const int NP = 1;

      template <class F>
      static double get(F &functor)
      {
        return functor(0.5);
      }
    };

    template <>
    struct SpecializedGQT<1, 3, fETypeE::cuboid>
    {
      static const int NP = 2;

      template <class F>
      static double get(F &functor)
      {
        static const double A = 1.0 / sqrt(3.0);
        static const double W = 0.5;
        return W * (functor(0.5 * (1.0 - A)) + functor(0.5 * (1.0 + A)));
      }
    };
    template <>
    struct SpecializedGQT<1, 2, fETypeE::cuboid> : public SpecializedGQT<1, 3, fETypeE::cuboid>
    {
    };

    template <>
    struct SpecializedGQT<1, 5, fETypeE::cuboid>
    {
      static const int NP = 3;

      template <class F>
      static double get(F &functor)
      {
        static const double A = sqrt(3.0 / 5.0);
        static const double W1 = 5.0 / 18.0;
        static const double W2 = 8.0 / 18.0;
        return W1 * (functor(0.5 * (1.0 - A)) + functor(0.5 * (1.0 + A))) +
               W2 * functor(0.5);
      }
    };
    template <>
    struct SpecializedGQT<1, 4, fETypeE::cuboid> : public SpecializedGQT<1, 5, fETypeE::cuboid>
    {
    };

    template <>
    struct SpecializedGQT<1, 7, fETypeE::cuboid>
    {
      static const int NP = 4;

      template <class F>
      static double get(F &functor)
      {
        static const double A = sqrt(3.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
        static const double B = sqrt(3.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
        static const double W1 = 0.25 - 1.0 / 12.0 * sqrt(5.0 / 6.0);
        static const double W2 = 0.25 + 1.0 / 12.0 * sqrt(5.0 / 6.0);

        return W1 * (functor(0.5 * (1.0 - A)) + functor(0.5 * (1.0 + A))) +
               W2 * (functor(0.5 * (1.0 - B)) + functor(0.5 * (1.0 + B)));
      }
    };
    template <>
    struct SpecializedGQT<1, 6, fETypeE::cuboid> : public SpecializedGQT<1, 7, fETypeE::cuboid>
    {
    };

    template <>
    struct SpecializedGQT<1, 9, fETypeE::cuboid>
    {
      static const int NP = 5;

      template <class F>
      static double get(F &functor)
      {
        static const double A = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
        static const double B = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
        static const double C = 0;
        static const double W1 = (322.0 - 13.0 * sqrt(70.0)) / 1800.0;
        static const double W2 = (322.0 + 13.0 * sqrt(70.0)) / 1800.0;
        static const double W3 = 512.0 / 1800.0;

        return W1 * (functor(0.5 * (1.0 - A)) + functor(0.5 * (1.0 + A))) +
               W2 * (functor(0.5 * (1.0 - B)) + functor(0.5 * (1.0 + B))) +
               W3 * functor(0.5 * (1.0 + C));
      }
    };
    template <>
    struct SpecializedGQT<1, 8, fETypeE::cuboid> : public SpecializedGQT<1, 9, fETypeE::cuboid>
    {
    };

  }
}

#include "../../internal/namespace.footer"
#endif

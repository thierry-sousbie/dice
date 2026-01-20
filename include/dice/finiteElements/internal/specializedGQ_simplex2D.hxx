#ifndef __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX2D_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX2D_HXX__

#include "../../internal/namespace.header"

// From Williams, Shunn & Jameson, 2014, JCAM

namespace internal
{
  namespace gaussQuadrature
  {

    template <>
    struct SpecializedGQT<2, 0, fETypeE::simplex>
    {
      static const int NP = 1;

      template <class F>
      static double get(F &functor)
      {
        return functor(0, 0, 0);
      }
    };

    template <>
    struct SpecializedGQT<2, 1, fETypeE::simplex>
    {
      static const int NP = 1;

      template <class F>
      static double get(F &functor)
      {
        static const double A = 1.0 / 3.0;
        return functor(A, A, A);
      }
    };

    template <>
    struct SpecializedGQT<2, 2, fETypeE::simplex>
    {
      static const int NP = 3;

      template <class F>
      static double get(F &functor)
      {
        static const double A = 1.0 / 6.0;
        static const double B = 2.0 / 3.0;
        static const double W = 1.0 / 3.0;

        return Permutation<3>::template get(functor, W, A, B);
        // return W*(functor(A,A,B)+functor(A,B,A)+functor(B,A,A));
      }
    };

    // No rule for 3, use 4 ...
    // template <>
    // struct SpecializedGQT<2,3,fETypeE::simplex> :
    //   public SpecializedGQT<2,4,fETypeE::simplex>
    // {
    // };

    template <>
    struct SpecializedGQT<2, 4, fETypeE::simplex>
    {
      static const int NP = 6;

      template <class F>
      static double get(F &functor)
      {
        static const double AA = 0.091576213509780;
        static const double AB = 0.816847572980440;

        static const double BA = 0.445948490915964;
        static const double BB = 0.108103018168071;

        static const double WA = 0.109951743655333;
        static const double WB = 0.223381589678000;

        return Permutation<3>::template get(functor, WA, AA, AB) +
               Permutation<3>::template get(functor, WB, BA, BB);
        /*
        return
          WA*(functor(AA,AA,AB)+functor(AA,AB,AA)+functor(AB,AA,AA))+
          WB*(functor(BA,BA,BB)+functor(BA,BB,BA)+functor(BB,BA,BA));
        */
      }
    };

    template <>
    struct SpecializedGQT<2, 3, fETypeE::simplex> : public SpecializedGQT<2, 4, fETypeE::simplex>
    {
    };

    template <>
    struct SpecializedGQT<2, 5, fETypeE::simplex>
    {
      static const int NP = 10;

      template <class F>
      static double get(F &functor)
      {
        static const double AA = 0.055564052669793;
        static const double AB = 0.888871894660413;

        static const double BA = 0.295533711735893;
        static const double BB = 0.634210747745723;
        static const double BC = 0.070255540518384;

        static const double CA = 0.333333333333333;

        static const double WA = 0.041955512996649;
        static const double WB = 0.112098412070887;
        static const double WC = 0.201542988584730;

        return Permutation<3>::template get(functor, WA, AA, AB) +
               Permutation<3>::template get(functor, WB, BA, BB, BC) +
               Permutation<3>::template get(functor, WC, CA);
        /*
        return
          WA*(functor(AA,AA,AB)+functor(AA,AB,AA)+functor(AB,AA,AA))+
          WB*(functor(BA,BB,BC)+functor(BC,BA,BB)+functor(BB,BC,BA)+
              functor(BB,BA,BC)+functor(BC,BB,BA)+functor(BA,BC,BB))+
          WC*funcotr(CA,CA,CA);
        */
      }
    };

    // No rule for 6, use 7 ...
    // template <>
    // struct SpecializedGQT<2,6,fETypeE::simplex> :
    //   public SpecializedGQT<2,7,fETypeE::simplex>
    // {
    // };

    template <>
    struct SpecializedGQT<2, 7, fETypeE::simplex>
    {
      static const int NP = 15;

      template <class F>
      static double get(F &functor)
      {
        static const double AA = 0.035870877695734;
        static const double AB = 0.928258244608533;

        static const double BA = 0.241729395767967;
        static const double BB = 0.516541208464066;

        static const double CA = 0.474308787777079;
        static const double CB = 0.051382424445843;

        static const double DA = 0.201503881881800;
        static const double DB = 0.751183631106484;
        static const double DC = 0.047312487011716;

        static const double WA = 0.017915455012303;
        static const double WB = 0.127712195881265;
        static const double WC = 0.076206062385535;
        static const double WD = 0.055749810027115;

        return Permutation<3>::template get(functor, WA, AA, AB) +
               Permutation<3>::template get(functor, WB, BA, BB) +
               Permutation<3>::template get(functor, WC, CA, CB) +
               Permutation<3>::template get(functor, WD, DA, DB, DC);
      }
    };

    template <>
    struct SpecializedGQT<2, 6, fETypeE::simplex> : public SpecializedGQT<2, 7, fETypeE::simplex>
    {
    };

    template <>
    struct SpecializedGQT<2, 8, fETypeE::simplex>
    {
      static const int NP = 21;

      template <class F>
      static double get(F &functor)
      {
        static const double AA = 0.028112952182664;
        static const double AB = 0.943774095634672;

        static const double BA = 0.177139098469317;
        static const double BB = 0.645721803061365;

        static const double CA = 0.405508595867433;
        static const double CB = 0.188982808265134;

        static const double DA = 0.148565812270887;
        static const double DB = 0.817900980028499;
        static const double DC = 0.033533207700614;

        static const double EA = 0.357196298615681;
        static const double EB = 0.604978911775132;
        static const double EC = 0.037824789609186;

        static const double WA = 0.010359374696538;
        static const double WB = 0.075394884326738;
        static const double WC = 0.097547802373242;
        static const double WD = 0.028969269372473;
        static const double WE = 0.046046366595935;

        return Permutation<3>::template get(functor, WA, AA, AB) +
               Permutation<3>::template get(functor, WB, BA, BB) +
               Permutation<3>::template get(functor, WC, CA, CB) +
               Permutation<3>::template get(functor, WD, DA, DB, DC) +
               Permutation<3>::template get(functor, WE, EA, EB, EC);
      }
    };

  }
}

#include "../../internal/namespace.footer"
#endif

#ifndef __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX3D_HXX__
#define __SPECIALIZED_GAUSS_QUADRATURE_SIMPLEX3D_HXX__

#include "../../internal/namespace.header"

namespace internal {
  namespace gaussQuadrature {

    template <>
    struct SpecializedGQT<3,0,fETypeE::simplex>
    {
      static const int NP  = 1;
      
      template <class F>
      static double get(F& functor)
      {
	return functor(0,0,0,0);
      }      
    };

    template <>
    struct SpecializedGQT<3,1,fETypeE::simplex>
    {
      static const int NP  = 1;
      
      template <class F>
      static double get(F& functor)
      {
	static const double A= 1.0/4.0;
	return functor(A,A,A,A);
      }      
    };

    template <>
    struct SpecializedGQT<3,2,fETypeE::simplex>
    {
      static const int NP = 4;      
   
      template <class F>
      static double get(F &functor)
      {
	static const double A=(5.0-sqrt(5.0))/20.0;
	static const double B=(5.0+3.0*sqrt(5.0))/20.0;
	static const double W=1.0/4.0;

	return Permutation<4>::template get31(functor,W,A,B);
      }      
    };   

    template <>
    struct SpecializedGQT<3,3,fETypeE::simplex>
    {
      static const int NP = 8;      
   
      template <class F>
      static double get(F &functor)
      {
	static const double AA=(55.0-3.0*sqrt(17.0)+sqrt(1022.0-134.0*sqrt(17.0)))/196.0;
	static const double AB=1.0-3.0*AA;

	static const double BA=(55.0-3.0*sqrt(17.0)-sqrt(1022.0-134.0*sqrt(17.0)))/196.0;
	static const double BB=1.0-3.0*BA;

	static const double WA=1.0/8.0+
	  sqrt((1715161837.0-406006699.0*sqrt(17.0))/23101.0)/3120.0;
	static const double WB=1.0/8.0-
	  sqrt((1715161837.0-406006699.0*sqrt(17.0))/23101.0)/3120.0;

	return 
	  Permutation<4>::template get31(functor,WA,AA,AB)+
	  Permutation<4>::template get31(functor,WB,BA,BB);
      }      
    };   

    template <>
    struct SpecializedGQT<3,4,fETypeE::simplex>
    {
      static const int NP = 14;      
   
      template <class F>
      static double get(F &functor)
      {
	static const double AA=0.09273525031089122640232391373703060;
	static const double AB=1.0-3.0*AA;

	static const double BA=0.31088591926330060979734573376345783;
	static const double BB=1.0-3.0*BA;

	static const double CA=0.45449629587435035050811947372066056;
	static const double CB=0.5-CA;
	
	static const double WA=7.34930431163613451e-02;
	static const double WB=1.12687925718015294e-01;
	static const double WC=4.25460207770822357e-02;

	return 
	  Permutation<4>::template get31(functor,WA,AA,AB)+
	  Permutation<4>::template get31(functor,WB,BA,BB)+
	  Permutation<4>::template get22(functor,WC,CA,CB);
      }      
    };  

    template <>
    struct SpecializedGQT<3,5,fETypeE::simplex>
    {
      static const int NP = 15;      
   
      template <class F>
      static double get(F &functor)
      {
	static const double AA=(7.0-sqrt(15.0))/34.0;
	static const double AB=1.0-3.0*AA;

	static const double BA=7.0/17.0-AA;
	static const double BB=1.0-3.0*BA;

	static const double CA=(10.0-2.0*sqrt(15.0))/40.0;
	static const double CB=0.5-CA;

	static const double DA=1.0/4.0;
	
	static const double WA=(2665.0+14.0*sqrt(15.0))/37800.0;
	static const double WB=(2665.0-14.0*sqrt(15.0))/37800.0;
	static const double WC=10.0/189.0;
	static const double WD=16.0/135.0;

	return 
	  Permutation<4>::template get31(functor,WA,AA,AB)+
	  Permutation<4>::template get31(functor,WB,BA,BB)+
	  Permutation<4>::template get22(functor,WC,CA,CB)+
	  Permutation<4>::template get22(functor,WD,DA);
      }      
    };  

    template <>
    struct SpecializedGQT<3,6,fETypeE::simplex>
    {
      static const int NP = 24;      
   
      template <class F>
      static double get(F &functor)
      {
	static const double AA=0.214602871259152029288839219386284991;
	static const double AB=1.0-3.0*AA;

	static const double BA=0.040673958534611353115579448956410059;
	static const double BB=1.0-3.0*BA;

	static const double CA=0.322337890142275510343994470762492125;
	static const double CB=1.0-3.0*CA;

	static const double DA=(3.0-sqrt(5.0))/12.0;
	static const double DB=(5.0+sqrt(5.0))/12.0;
	static const double DC=(1.0+sqrt(5.0))/12.0;
	
	static const double WA=3.99227502581672861e-02;
	static const double WB=1.00772110553204019e-02;
	static const double WC=5.53571815436545017e-02;
	static const double WD=27.0/560.0;

	return 
	  Permutation<4>::template get31(functor,WA,AA,AB)+
	  Permutation<4>::template get31(functor,WB,BA,BB)+
	  Permutation<4>::template get31(functor,WC,CA,CB)+
	  Permutation<4>::template get211(functor,WD,DA,DB,DC);
      }      
    };  

    template <>
    struct SpecializedGQT<3,9,fETypeE::simplex>
    {
      static const int NP = 84;      
   
      template <class F>
      static double get(F &functor)
      {
	// From Williams, Shunn & Jameson, 2014, JCAM
	
	static const double AA=0.026878474414817;
	static const double AB=0.919364576755549;
	
	static const double BA=0.187140675803470;
	static const double BB=0.438577972589591;
	
	static const double CA=0.473575835127937;
	static const double CB=0.026424164872063;
	
	static const double DA=0.352045262027356;
	static const double DB=0.147954737972644;

	static const double EA=0.020953442220056;
	static const double EB=0.732309909692947;
	static const double EC=0.225783205866940;
	  
	static const double FA=0.096989733123466;
	static const double FB=0.647557594086975;
	static const double FC=0.158462939666092;

	static const double GA=0.322111431830857;
	static const double GB=0.033665704507429;

	static const double HA=0.097608162890442;
	static const double HB=0.011844417749498;
	static const double HC=0.792939256469618;

	static const double IA=0.133558160703568;
	static const double IB=0.296501020543124;
	static const double IC=0.028756405953071;
	static const double ID=0.541184412800237;

	static const double WA=0.002144935144316;
	static const double WB=0.020826641690769;
	static const double WC=0.007210136064455;
	static const double WD=0.030798919159712;
	static const double WE=0.004357844813864;
	static const double WF=0.008593530677833;
	static const double WG=0.023000681669286;
	static const double WH=0.004863063904912;
	static const double WI=0.015595140078259;

	return 
	  Permutation<4>::template get31(functor,WA,AA,AB)+
	  Permutation<4>::template get31(functor,WB,BA,BB)+
	  Permutation<4>::template get22(functor,WC,CA,CB)+
	  Permutation<4>::template get22(functor,WD,DA,DB)+
	  Permutation<4>::template get211(functor,WE,EA,EB,EC)+
	  Permutation<4>::template get211(functor,WF,FA,FB,FC)+
	  Permutation<4>::template get31(functor,WG,GA,GB)+
	  Permutation<4>::template get211(functor,WH,HA,HB,HC)+
	  Permutation<4>::template get(functor,WI,IA,IB,IC,ID);	  
      }
    };

    // No rule for 7 and 8 => use 9, 84 points
    template <>
    struct SpecializedGQT<3,7,fETypeE::simplex>
      : public SpecializedGQT<3,9,fETypeE::simplex>
    {};

    template <>
    struct SpecializedGQT<3,8,fETypeE::simplex>
      : public SpecializedGQT<3,9,fETypeE::simplex>
    {};

  }
}

#include "../../internal/namespace.footer"
#endif

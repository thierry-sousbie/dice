#ifndef __COSMO_PARAMETERS_HXX__
#define __COSMO_PARAMETERS_HXX__

#include <algorithm>

#include "../internal/namespace.header" 


namespace cosmo {

  struct CosmoParameters
  {
    typedef CosmoParameters MyType;

    static std::string parserCategory() {return "cosmology";}
    static std::string classHeader() {return "cosmo_parameters";}
    static float classVersion() {return 0.10;}
    static float compatibleSinceClassVersion() {return 0.10;}

    double oM;	//!< baryon+dark matter density    
    double oL;	//!< dark energy density
    double oB;	//!< baryon matter density
    double oK;  //!< curvature density
    double h0;	//!< Hubble constant
    double w;   //!< dark energy

    CosmoParameters(double oM_p=0.3L, double oL_p=0.7L, double oB_p=0.0L, double h0_p=0.72L, 
		    double w_p=-1.0L, double oK_p=0)
    {
      set(oM_p,oL_p,oB_p,h0_p,w_p,oK_p);
    }

    void set(double oM_p, double oL_p, double oB_p, double h0_p, double w_p=-1.0L, double oK_p=0)
    {
      oM=oM_p;
      oL=oL_p;
      oB=oB_p;
      h0=h0_p;
      w=std::min(w_p,double(0));
      //if (oK_p<=-1.0) oK= 1.0L-oM-oL;
      //else oK=oK_p;
      oK=oK_p;
    }

    template <class R, class PM>
    void parse(R *reader, PM &paramsManager)
    {
      oM=paramsManager.get("oM",parserCategory(),oM,reader,PM::FILE_FIRST,
			   "Matter fraction");
      oL=paramsManager.get("oL",parserCategory(),oL,reader,PM::FILE_FIRST,
			   "Cosmological constant fraction");  
      oB=paramsManager.get("oB",parserCategory(),oB,reader,PM::FILE_FIRST,
			   "Baryonic matter fraction (already included in oM)");
      h0=paramsManager.get("h0",parserCategory(),h0,reader,PM::FILE_FIRST,
			   "Hubble parameter (h0 = H0 / 100 kms-1 Mpc-1)");   
      w =paramsManager.get("w",parserCategory(),w,reader,PM::FILE_FIRST,
			   "Equation of state for dark energy");   
      oK=1.0L-oM-oL;
      // oK=paramsManager.get("oK",parserCategory(),oK,reader,PM::FILE_FIRST,
      // 			   "Fraction of curvature");  
      
    }

    template <class PP>
    void parse(PP &paramsParser)
    {
      oM=paramsParser.get("oM",parserCategory(),oM,
			  "Matter fraction");			
      oL=paramsParser.get("oL",parserCategory(),oL,
			  "Cosmological constant fraction");   
      oB=paramsParser.get("oB",parserCategory(),oB,
			  "Baryonic matter fraction (already included in oM)");
      h0=paramsParser.get("h0",parserCategory(),h0,
			  "Hubble parameter (h0 = H0 / 100 kms-1 Mpc-1)");
      w=paramsParser.get("w",parserCategory(),w,
			 "Equation of state for dark energy");
      oK = 1.0L-oM-oL;
      // oK=paramsParser.get("oK",parserCategory(),oK,
      // 			  "Fraction of curvature");
    }    
  };

} // namespace cosmo

#include "../internal/namespace.footer"
#endif

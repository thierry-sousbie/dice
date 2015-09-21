#ifndef __SIMPLICIAL_GRID_PARAMS_HXX__
#define __SIMPLICIAL_GRID_PARAMS_HXX__

/**
 * @file 
 * @brief  Parameters for SimplicialGridT
 * @author Thierry Sousbie
 */

#include "./tesselationType.hxx"

#include "../../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

template <int ND, int NDW ,class CT>
struct SimplicialGridParamsT {
  typedef CT Cell;  
  static const int NDIM=ND;
  static const int NDIM_W=NDW;

  double x0[NDIM_W];
  double delta[NDIM_W];
  int resolution[NDIM];
  TesselationType t;
    
  SimplicialGridParamsT()
  {
    std::fill_n(x0,NDIM,-1);
    std::fill_n(x0+NDIM,NDIM_W-NDIM,0);
    std::fill_n(delta,NDIM,2);
    std::fill_n(delta+NDIM,NDIM_W-NDIM,0);
    std::fill_n(resolution,NDIM,64);
    t=TesselationTypeV::ANY;
  }

  template <class BR, class PM>
  void parse(BR *reader, PM &manager, 
	     const std::string parserCategory="simplicalGrid",
	     int which=0)
  { 
    for (int i=0;i<NDIM_W;i++)
      {	
	x0[i]=manager.
	  get("x0",parserCategory,x0[i],i+which*NDIM_W,
	      reader,PM::FILE_FIRST,
	      "Initial coordinates of the mesh lower left corner");

	delta[i]=manager.
	  get("delta",parserCategory,delta[i],i+which*NDIM_W,
	      reader,PM::FILE_FIRST,
	      "Size of the initial mesh");			 

	if (i<NDIM) 
	  {
	    if (i>0) resolution[i]=resolution[i-1];
	    resolution[i]=manager.
	      get("resolution",parserCategory,resolution[i],i+which*NDIM,
		  reader,PM::FILE_FIRST,
		  "Resolution of the initial mesh (in pixels before tesselation)");
	  }
      }  

    std::string initTesselationTypeStr = manager.template 
      get<std::string>("initTesselationType",parserCategory,
		       TesselationTypeSelect().getString(t,true),which,
		       reader,PM::FILE_FIRST,
		       TesselationTypeSelect().getAllString
		       ("Topology of the initial grid simplicial tesselation (%s)"));
    t = TesselationTypeSelect().getVal(initTesselationTypeStr,true);
  }

  template <class PP>
  void parse(PP &parser, 
	     const std::string parserCategory="simplicalGrid",
	     int which=0)
  {
    for (int i=0;i<NDIM_W;i++)
      {
	x0[i]=parser.
	  get("x0",parserCategory,x0[i],i+which*NDIM_W,
	      "Initial coordinates of the mesh lower left corner");
	
	delta[i]=parser.
	  get("delta",parserCategory,delta[i],i+which*NDIM_W,
	      "Size of the initial mesh");			 

	if (i<NDIM) 
	  resolution[i]=parser.
	    get("resolution",parserCategory,resolution[i],i+which*NDIM,
		"Resolution of the initial mesh (in pixels before tesselation)");
      }  

    std::string initTesselationTypeStr = parser.template 
      get<std::string>("initTesselationType",parserCategory,
		       TesselationTypeSelect().getString(t,true),which,
		       TesselationTypeSelect().getAllString
		       ("Topology of the initial grid simplicial tesselation (%s)"));
    t = TesselationTypeSelect().getVal(initTesselationTypeStr,true);
  }
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

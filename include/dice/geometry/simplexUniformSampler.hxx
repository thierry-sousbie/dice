#ifndef __SIMPLEX_UNIFORMSAMPLER_HXX__
#define __SIMPLEX_UNIFORMSAMPLER_HXX__

#include <stdlib.h>
#include "../tools/wrappers/standardDrand48ReentrantWrapper.hxx"
#include "./internal/simplexVolume_implementation.hxx"
#include <cmath>

/**
 * @file 
 * @brief  A class to generate a uniform distribution of points within a simplex
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Geometry
 *   \{
 */

/**
 * \class SimplexUniformSamplerT
 * \brief A class to generate a uniform distribution of points within a simplex
 *
 * \tparam ND  Number of the dimensions of the simplex
 * \tparam NDW Number of dimensions of the embedding space
 * \tparam T   Data type
 */
template <int ND, int NDW> class SimplexUniformSamplerT;

template <int NDW>
class SimplexUniformSamplerT<2,NDW>
{    
public:
  typedef SimplexUniformSamplerT<2,NDW> MyType;
  static const int NDIM = 2;
  static const int NDIM_W = NDW;

  template <class T1, class IT>
  static long generateRandom(const T1 p[NDIM_W],
			     const T1 base[NDIM][NDIM_W], 
			     int N, IT coordsOut, 
			     typename DRand48_rWapper::RandData *randSeed)
  {
    for (int i=0;i<N;++i)
      {
	double x,y;
	DRand48_rWapper::drand48_r(randSeed,&x);
	DRand48_rWapper::drand48_r(randSeed,&y);
	generate(x,y,p,base,coordsOut);
      }
    return N;
  }

  template <class T1, class T2, class IT1, class IT2>
  static long generateRandom(const T1 p[NDIM_W],
				     const T1 base[NDIM][NDIM_W], 
				     const T2 value[NDIM+1],
				     int N, 
				     IT1 coordsOut,
				     IT2 valueOut, 
				     typename DRand48_rWapper::RandData *randSeed)
  {
    for (int i=0;i<N;++i)
      {
	double x,y;
	DRand48_rWapper::drand48_r(randSeed,&x);
	DRand48_rWapper::drand48_r(randSeed,&y);
	generate(x,y,p,base,value,coordsOut,valueOut);
      }
    return N;
  }

  
  template <class T1, class IT>
  static long generateRegular(const T1 p[NDIM_W],
			      const T1 base[NDIM][NDIM_W], 
			      int NPerDim, 
			      IT coordsOut)
  {
    const double delta = 1.0/(NPerDim+1);
    const double max   = 1.0-delta/4;
    long count=0;

    for (double x=delta/2; x<max; x+=delta)
      for (double y=delta/2; y<max; y+=delta)
	{ 
	  generate(x,y,p,base,coordsOut);
	  ++count;	  	  
	}
    return count;

  }

  template <class T1, class T2, class IT1, class IT2>
  static long generateRegular(const T1 p[NDIM_W],
			      const T1 base[NDIM][NDIM_W], 
			      const T2 value[NDIM+1],
			      int NPerDim, 
			      IT1 coordsOut,
			      IT2 valueOut)
  {
    const double delta = 1.0/(NPerDim+1);    
    const double max   = 1.0-delta/4;   
    long count=0;

    for (double x=delta/2; x<max; x+=delta)
      for (double y=delta/2; y<max; y+=delta)
	{
	  generate(x,y,p,base,value,coordsOut,valueOut);
	  ++count;	  
	}
    return count;
  }

private:
  template <class T1, class T2, class IT1, class IT2>
  static void generate(double x, double y,
		       const T1 p[NDIM_W],
		       const T1 base[NDIM][NDIM_W], 
		       const T2 value[NDIM+1],
		       IT1 &coordsOut,
		       IT2 &valueOut)
  {
    double r,s;
    fold(x,y,r,s);

    for (int d=0;d<NDIM_W;++d)
      {
	(*coordsOut)=p[d]+r*base[0][d]+s*base[1][d];
	++coordsOut;
      }
    (*valueOut)=(1.0-r-s)*value[0]+r*value[1]+s*value[2];
    ++valueOut;   
  }
  
  template <class T1, class IT1>
  static void generate(double x, double y,
		       const T1 p[NDIM_W],
		       const T1 base[NDIM][NDIM_W], 		
		       IT1 &coordsOut)
  {
    double r,s;
    fold(x,y,r,s);

    for (int d=0;d<NDIM_W;++d)
      {
	(*coordsOut)=p[d]+r*base[0][d]+s*base[1][d];
	++coordsOut;
      }   
  }

  static void fold(double x, double y, double &r, double &s)
  {
    if (x+y<=1.0)
      {
	r=x;
	s=y;
      }
    else
      {
	r=1.0-x;
	s=1.0-y;
      }
  }  
};




template <int NDW>
class SimplexUniformSamplerT<3,NDW>
{    
public:
  typedef SimplexUniformSamplerT<3,NDW> MyType;
  static const int NDIM = 3;
  static const int NDIM_W = NDW;

  template <class T1, class IT>
  static long generateRandom(const T1 p[NDIM_W],
			     const T1 base[NDIM][NDIM_W], 
			     int N, IT coordsOut, 
			     typename DRand48_rWapper::RandData *randSeed)
  {
    for (int i=0;i<N;++i)
      {
	double x,y,z;
	DRand48_rWapper::drand48_r(randSeed,&x);
	DRand48_rWapper::drand48_r(randSeed,&y);
	DRand48_rWapper::drand48_r(randSeed,&z);
	generate(x,y,z,p,base,coordsOut);
      }
    return N;
  }

  template <class T1, class T2, class IT1, class IT2>
  static long generateRandom(const T1 p[NDIM_W],
				     const T1 base[NDIM][NDIM_W], 
				     const T2 value[NDIM+1],
				     int N, 
				     IT1 coordsOut,
				     IT2 valueOut, 
				     typename DRand48_rWapper::RandData *randSeed)
  {
    for (int i=0;i<N;++i)
      {
	double x,y,z;
	DRand48_rWapper::drand48_r(randSeed,&x);
	DRand48_rWapper::drand48_r(randSeed,&y);
	DRand48_rWapper::drand48_r(randSeed,&z);
	generate(x,y,z,p,base,value,coordsOut,valueOut);
      }
    return N;
  }

  
  template <class T1, class IT>
  static long generateRegular(const T1 p[NDIM_W],
			      const T1 base[NDIM][NDIM_W], 
			      int NPerDim, 
			      IT coordsOut)
  {
    const double delta = 1.0/(NPerDim+1);
    const double max   = 1.0-delta/4;
    long count=0;
    // This is not very fast ...
    for (double x=delta/4; x<max; x+=delta)
      for (double y=delta/4; y<max; y+=delta)
	for (double z=delta/4; z<max; z+=delta)
	  {
	    generate(x,y,z,p,base,coordsOut);
	    ++count;	  	  
	  }
    return count;
  }

  template <class T1, class T2, class IT1, class IT2>
  static long generateRegular(const T1 p[NDIM_W],
			      const T1 base[NDIM][NDIM_W], 
			      const T2 value[NDIM+1],
			      int NPerDim, 
			      IT1 coordsOut,
			      IT2 valueOut)
  {
    const double delta = 1.0/(NPerDim+1);
    const double max   = 1.0-delta/4;
    long count=0;
    // This is not very fast ...
    for (double x=delta/4; x<max; x+=delta)
      for (double y=delta/4; y<max; y+=delta)
	for (double z=delta/4; z<max; z+=delta)
	  {
	    generate(x,y,p,base,value,coordsOut,valueOut);
	    ++count;	  	  
	  }
    return count;
  }

private:
  template <class T1, class T2, class IT1, class IT2>
  static void generate(double x, double y, double z,
		       const T1 p[NDIM_W],
		       const T1 base[NDIM][NDIM_W], 
		       const T2 value[NDIM+1],
		       IT1 &coordsOut,
		       IT2 &valueOut)
  {
    double r,s,t;
    fold(x,y,z,r,s,t);

    for (int d=0;d<NDIM_W;++d)
      {
	(*coordsOut)=p[d]+r*base[0][d]+s*base[1][d]+t*base[2][d];
	++coordsOut;
      }
    (*valueOut)=(1.0-r-s-t)*value[0]+r*value[1]+s*value[2]+t*value[3];
    ++valueOut;
  }
  
  template <class T1, class IT1>
  static void generate(double x, double y, double z,
		       const T1 p[NDIM_W],
		       const T1 base[NDIM][NDIM_W], 		
		       IT1 &coordsOut)
  {
    double r,s,t;
    fold(x,y,z,r,s,t);

    for (int d=0;d<NDIM_W;++d)
      {
	(*coordsOut)=p[d]+r*base[0][d]+s*base[1][d]+t*base[2][d];
	++coordsOut;
      }   
  }

  static void fold(double x, double y, double z, double &r, double &s, double &t)
  {
    r=x;s=y;t=z;
	    
    if (r+s>1.0)
      {
	r = 1.0 - r;
	s = 1.0 - s;
      }

    if(s+t>1.0) 
      {
	t = 1.0 - r - s;
	s = 1.0 - z;
      } 
    else if (r+s+t>1.0)
      {		
	t = r + s + t - 1.0;
	r = 1.0 - s - z;
      }	  
  }  
};

  /*
template <int NDW>
struct SimplexUniformSamplerT<3,NDW>
{    
  typedef SimplexUniformSamplerT<3,NDW> MyType;
  static const int NDIM = 3;
  static const int NDIM_W = NDW;

  template <class T, class IT>
  static long generate(const T p[NDIM_W],
		       const T base[NDIM][NDIM_W], 
		       int N, 
		       IT coordsOut)
  {
    int NN[NDIM]={N,N,N};
    return generate(p,base,NN,coordsOut);
  }

  template <class T1, class T2, class IT1, class IT2>
  static int generate(const T1 p[NDIM_W],
		      const T1 base[NDIM][NDIM_W], 
		      const T2 value[NDIM+1],
		      int N, 
		      IT1 coordsOut,
		      IT2 valueOut)
  {
    int NN[NDIM]={N,N};
    return generate(p,base,value,NN,coordsOut,valueOut);
  }

  template <class T1, class T2, class IT>
  static long generate(const T1 p[NDIM_W],
		       const T2 base[NDIM][NDIM_W], 
		       int N[NDIM], 
		       IT coordsOut)
  {
    const double delta[NDIM]= {1.0/(N[0]+1),1.0/(N[1]+1),1.0/(N[2]+1)};
    const double max[NDIM]= {1.0-delta[0]/4,1.0-delta[1]/4,1.0-delta[2]/4};
    
    // This is not very fast ...
    for (double x=delta[0]/2; x<max[0]; x+=delta[0])
      for (double y=delta[1]/2; y<max[1]; y+=delta[1])
	for (double z=delta[2]/2; z<max[2]; z+=delta[2])
	  {
	    double s=x;
	    double t=y;
	    double u=z;
	    
	    if (s+t>1.0)
	      {
		s = 1.0 - s;
		t = 1.0 - t;
	      }

	    if(t+u>1.0) 
	      {
		u = 1.0 - s - t;
		t = 1.0 - z;
	      } 
	    else if (s+t+u>1.0)
	      {		
		u = s + t + u - 1.0;
		s = 1.0 - t - z;
	      }
	    
	    for (int d=0;d<NDIM_W;++d)
	      {
		(*coordsOut)=p[d]+s*base[0][d]+t*base[1][d]+u*base[2][d];
		++coordsOut;
	      }
	  } 

    return N[0]*N[1]*N[2];
  }

  template <class T1, class T2, class T3, class IT1, class IT2>
  static int generate(const T1 p[NDIM_W],
		      const T1 base[NDIM][NDIM_W], 
		      const T2 value[NDIM+1],
		      T3 N[NDIM], 
		      IT1 coordsOut,
		      IT2 valueOut)
  {
    const double delta[NDIM]= {1.0/(N[0]+1),1.0/(N[1]+1),1.0/(N[2]+1)};
    const double max[NDIM]= {1.0-delta[0]/4,1.0-delta[1]/4,1.0-delta[2]/4};
    
    // This is not very fast ...
    for (double x=delta[0]/2; x<max[0]; x+=delta[0])
      for (double y=delta[1]/2; y<max[1]; y+=delta[1])
	for (double z=delta[2]/2; z<max[2]; z+=delta[2])
	  {
	    double s=x;
	    double t=y;
	    double u=z;
	    
	    if (s+t>1.0)
	      {
		s = 1.0 - s;
		t = 1.0 - t;
	      }

	    if(t+u>1.0) 
	      {
		u = 1.0 - s - t;
		t = 1.0 - z;
	      } 
	    else if (s+t+u>1.0)
	      {		
		u = s + t + u - 1.0;
		s = 1.0 - t - z;
	      }
	    
	    for (int d=0;d<NDIM_W;++d)
	      {
		(*coordsOut)=p[d]+s*base[0][d]+t*base[1][d]+u*base[2][d];
		++coordsOut;
	      }
	    (*valueOut)=(1.0-s-t-u)*value[0]+s*value[1]+t*value[2]+u*value[3];
	  } 

    return N[0]*N[1]*N[2];
  }
};

  */
/** \}*/
#include "../internal/namespace.footer"
#endif

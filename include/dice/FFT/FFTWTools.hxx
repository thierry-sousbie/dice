#ifndef __FFTW_TOOLS_HXX__
#define __FFTW_TOOLS_HXX__

/**
 * @file 
 * @brief defines a few helper functions to facilitate FFTW usage
 * @author Thierry Sousbie
 */

#include <fftw3.h>
#include <string.h>

#include "../tools/helpers/helpers.hxx"

#include "../internal/namespace.header"
/** \addtogroup FFT
 *   \{
 */

/**
 * \class FFTWToolsT
 * A class that contains static member functions that can be used as helpers for FFTW
 */
template <typename DT, int ND>
struct FFTWToolsT
{
public:
  typedef DT Data;
  static const int NDIM = ND;
 
  template <typename T1,typename T2>
  static long actualSize(const T1 dim[NDIM], T2 dimP[NDIM], 
			 T2 deltaP[NDIM+1], T2 strideP[NDIM+1], 
			 int FFTWPadding, int periodic)
  {
    T2 cDim[NDIM];
    for (int i=0;i<NDIM;++i) cDim[i]=dim[i];
    /*
    if (!periodic)
      for (int i=0;i<NDIM;++i) cDim[i]=dim[i]*2;
    else
      for (int i=0;i<NDIM;++i) cDim[i]=dim[i];
    */

    deltaP[0]=0;    
    strideP[0]=1;   

    if (FFTWPadding)
      dimP[0]=(cDim[0]/2+1)*2;
    else
      dimP[0]=cDim[0];
    
    for (int i=1;i<NDIM;i++) 
      dimP[i]=cDim[i];
    
    for (int i=1;i<=NDIM;++i)
      {
	strideP[i] = strideP[i-1]*dimP[i-1];
	deltaP[i] = dimP[i-1]-cDim[i-1];
      }

    return strideP[NDIM];
  }

  template <typename T>
  static long actualSize(const T dim[NDIM], int FFTWPadding, int periodic)
  {
    T strideP[NDIM+1]; 
    T dimP[NDIM+1]; 
    T deltaP[NDIM+1]; 
    return actualSize(dim,dimP,deltaP,strideP,FFTWPadding,periodic);   
  }

  template <typename T>
  static Data *rearrange(Data *src, Data *dest,const T dim[NDIM], 
			 int FFTWPadding, int periodic, bool toFFTW=true)
  {   
    Data *newData=dest;
    T strideP[NDIM+1];
    T deltaP[NDIM+1];
    T dimP[NDIM]; 
   
    long nval=1;
    long i,j;

    for (i=0;i<NDIM;i++) nval*=dim[i];   
    
    //long size = 
      actualSize(dim,dimP,deltaP,strideP,FFTWPadding,periodic);

    // if (!FFTWPadding)
    //   {
    // 	if (src==dest) return dest;
    //   }

    //newData=dest;
    long lineSize=sizeof(Data)*dim[0];
    if (toFFTW)
      {	
	if ((deltaP[1]>0)||(newData!=src))
	  {
	    long dj=dim[0];
	    long di=dim[0] + deltaP[1];
	    j=nval;
	    i=strideP[NDIM]-deltaP[1];
	    for (;i>0;)
	      {	
		memmove(&newData[i-dim[0]],&src[j-dim[0]],lineSize);
		i-=di;j-=dj;		
	      }
	  }
      }
    else
      {
	if ((deltaP[1]>0)||(newData!=src))
	  {	    
	    long dj=dim[0];
	    long di=dim[0] + deltaP[1];
	    j=0;
	    i=0;
	    for (;i<strideP[NDIM];)
	      {
		memmove(&newData[j],&src[i],lineSize);
		j+=dj;i+=di;	   	
	      }
	  }
      }

    return newData;
  }

  template <typename T>
  static Data *rearrangeToFFTW(Data *src,Data *dst, const T dim[NDIM], 
			       int FFTWPadding, int periodic)
  {
    return rearrange(src,dst,dim,FFTWPadding,periodic,true);
  }

  template <typename T>
  static Data *rearrangeFromFFTW(Data *src,Data *dst, const T dim[NDIM],
				 int FFTWPadding, int periodic)
  {
    return rearrange(src,dst,dim,FFTWPadding,periodic,false);
  }
  
  template <typename T>
  static T indGen(T k,T dim)
  {
    return (k<=(dim/2))?k:(k-dim);  
  }

  template <typename T>
  static T indGen(T k[NDIM],T dim[NDIM], T out[NDIM])
  {
    for (int i=0;i<NDIM;i++)
      out[i]=(k[i]<=(dim[i]/2))?k[i]:(k[i]-dim[i]);  
  }

  template <typename T>
  static T indGen(T *k,T *dim, T *out, int ndims)
  {
    for (int i=0;i<ndims;i++)
      out[i]=(k[i]<=(dim[i]/2))?k[i]:(k[i]-dim[i]);  
  }

  template <typename T,typename T2,typename T3>
  static double norm(const T k, const T2 delta, const T3 size)
  {
    double a=(k<=(delta/2))?k:(k-delta);
    return (a*size)/delta;
  }

  template <typename T,typename T2>
  static double norm(const T k, const T2 delta)
  {
    double a=(k<=(delta/2))?k:(k-delta);
    return a/delta;
  }
 
  template <typename T,typename T2>
  static double normV(const T k[NDIM], const T2 delta[NDIM])
  {
    double a[NDIM];
    double result=0;

    for (long i=0;i<NDIM;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2>
  static double norm(const T *k, const T2 *delta, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2,typename T3>
  static double normV(const T k[NDIM], const T2 delta[NDIM], const T3 size[NDIM])
  {
    double a[NDIM];
    double result=0;

    for (long i=0;i<NDIM;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2,typename T3>
  static double norm(const T *k, const T2 *delta, const T3 *size, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2>
  static double norm2V(const T k[NDIM], const T2 delta[NDIM])
  {
    double a[NDIM];
    double result=0;

    for (long i=0;i<NDIM;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2>
  static double norm2(const T *k, const T2 *delta, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2,typename T3>
  static double norm2V(const T k[NDIM], const T2 delta[NDIM], const T3 size[NDIM])
  {
    double a[NDIM];
    double result=0;

    for (long i=0;i<NDIM;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2,typename T3>
  static double norm2(const T *k, const T2 *delta, const T3 *size, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]/2))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  /** \brief implements a hanning filter in fourier space
   *  \param k1Norm The norm of the normalized wave vector (i.e. the norm of k1 where each
   *  component of k1 is -0.5<k1_i<0.5
   */
  static double hanning(double k1Norm)
  {
    static const double twopi=atan(1.0L)*8.0L;    
    return (k1Norm>0.5)?0:(0.5L*(1.0L+cos(twopi*k1Norm)));
  }

};

/** \}*/
#include "../internal/namespace.footer"
#endif

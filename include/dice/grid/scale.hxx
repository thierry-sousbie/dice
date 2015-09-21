#ifndef __SCALE_HXX__
#define __SCALE_HXX__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#include <algorithm>
#include <string>

#include "valLocationType.hxx"

#include "../tools/types/typeSelect.hxx"

#include "../internal/namespace.header"

template <typename T>
struct ScaleT {

  struct ScaleTypeV {
    enum Type {LINEAR=0,LOGARITHMIC=1,QUADRATIC=2,CUBIC=3,QUARTIC=4,UNDEFINED=-1};
  };
  typedef typename ScaleTypeV::Type ScaleType;

  struct ScaleTypeSelect : public TypeSelectT<ScaleTypeV> {
    ScaleTypeSelect()
    {
      this->insert("linear",ScaleTypeV::LINEAR);
      this->insert("logarithmic",ScaleTypeV::LOGARITHMIC);
      this->insert("quadratic",ScaleTypeV::QUADRATIC);
      this->insert("cubic",ScaleTypeV::CUBIC);
      this->insert("quartic",ScaleTypeV::QUARTIC);
    }
    std::string name() {return "scale_type";}
  };

  // typedef ValLocationTypeV ValLocationTypeV;
  // typedef ValLocationType  ValLocationType;  
  typedef T Data;

  /*
  static scaleTypeT str2Type(std::string str)
  {
    
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);

    if (str=="LINEAR") return scaleTypeV::LINEAR;
    if (str=="LOGARITHMIC") return scaleTypeV::LOGARITHMIC;
    if (str=="QUADRATIC") return scaleTypeV::QUADRATIC;
    if (str=="CUBIC") return scaleTypeV::CUBIC;
    if (str=="QUARTIC") return scaleTypeV::QUARTIC;
    
    fprintf(stderr,"ERROR in scale::str2Type : unknow type '%s'.\n",str.c_str());
    exit(-1);
  }

  static std::string type2Str(scaleTypeT type)
  {
    if (type==scaleTypeV::LINEAR) return "LINEAR";
    if (type==scaleTypeV::LOGARITHMIC) return "LOGARITHMIC";
    if (type==scaleTypeV::QUADRATIC) return "QUADRATIC";
    if (type==scaleTypeV::CUBIC) return "CUBIC";
    if (type==scaleTypeV::QUARTIC) return "QUARTIC";
    

    fprintf(stderr,"ERROR in scale::str2Type : invalid type .\n");
    exit(-1);
  }
  
  static std::string valLocation2Str(valLocationT vl)
  {
    if (vl==valLocationV::CELL) return "C";
    if (vl==valLocationV::VERTEX) return "V";
    
    fprintf(stderr,"ERROR in valLocation2Str : invalid type .\n");
    exit(-1);
  }
  */

  static T valueAt(double start, double stop, long N, long index, 
		   ScaleType sct=ScaleTypeV::LINEAR, 
		   ValLocationType ltype=ValLocationTypeV::VERTEX,
		   bool periodic = false)
  {
    double val=0;

    if (sct==ScaleTypeV::LINEAR) { 
      double dx=(stop-start)/(N);
      val =start+dx*index;
    }
    else if (sct==ScaleTypeV::LOGARITHMIC) {
      if (start*stop<=0) {
	fprintf(stderr,"error: cannot generate logscale from a 0-crossing range.\n");
	exit(0);
      }
      double s=(start<0)?-1:1;
      val=valueAt(log(fabs(start)),log(fabs(stop)),N,index,
		  ScaleTypeV::LINEAR,ValLocationTypeV::VERTEX);
      val = s*exp(val);
    }
    else if ((int)sct>=2) {
      int pp=(int)sct;
      val=valueAt(pow(start,1./pp),pow(stop,1./pp),N,index,
		  ScaleTypeV::LINEAR,ValLocationTypeV::VERTEX);
      val=pow(val,double(pp));
    }

    if (ltype==ValLocationTypeV::CELL) {
      double val2 = valueAt(start,stop,N,index+1,sct,ValLocationTypeV::VERTEX);
      val = (val+val2)*0.5;
    }
    else
      {
	if (index==0) return start;
	if (index==N) return stop;
      }

    return val;
  }
  
  static std::vector<T> genScale(double start, double stop, long N, 
				 ScaleType sct, ValLocationType ltype)
  {
    std::vector<double> scale;
    long i;
     
    if (sct==ScaleTypeV::LINEAR) { 
      double dx=(stop-start)/(N);
      scale.resize(N+1);
      for (i=0;i<N+1;i++) scale[i]=start+dx*i;
    } 
    else if (sct==ScaleTypeV::LOGARITHMIC) {
      if (start*stop<=0) {
	fprintf(stderr,"error: cannot generate logscale from a 0-crossing range.\n");
	exit(0);
      }
      double s=(start<0)?-1:1;
      scale=genScale(log(fabs(start)),log(fabs(stop)),N,
		     ScaleTypeV::LINEAR,ValLocationTypeV::VERTEX);
      for (i=0;i<N+1;i++) scale[i]=s*exp(scale[i]);
    }
    else if ((int)sct>=2) {
      int pp=(int)sct;
      scale=genScale(pow(start,1./pp),pow(stop,1./pp),N,
		     ScaleTypeV::LINEAR,ValLocationTypeV::VERTEX);
      for (i=0;i<N+1;i++) scale[i]=pow(scale[i],double(pp));
    }
  
    if (ltype==ValLocationTypeV::CELL) {
      int j;
      for (j=0;j<scale.size()-1;j++)  scale[j]=(scale[j]+scale[j+1])*0.5;
      scale.resize(scale.size()-1);
    }
    else
      {
	scale.front()=start;
	scale.back()=stop;	
      }

    return scale;
  }

  static std::vector<T> genScaleDelta(double start, double stop, long N, 
				      ScaleType sct, ValLocationType ltype)
  {
    std::vector<double> scale=genScale(start,stop,N,sct,ValLocationTypeV::VERTEX);
    long j;
        
    if (ltype==ValLocationTypeV::CELL) {
      for (j=0;j<scale.size()-1;j++)  scale[j]=(scale[j+1]-scale[j]);
      scale.resize(scale.size()-1);
    }
    else
      {
	std::vector<T> tmp=genScale(start,stop,N,sct,ValLocationTypeV::CELL);
	scale.front()-=start;
	for (j=1;j<scale.size()-1;j++)  
	  scale[j]=(tmp[j]-tmp[j-1]);
	scale.back()-=tmp.back();
      }

    return scale;
  }


  template <typename dT>
  static bool rescale(dT &x,ScaleType sct)
  { 
    switch (sct)
      {
      case ScaleTypeV::LOGARITHMIC:
	x=log((x));
	return true;
      case ScaleTypeV::LINEAR:
	return false;
      }

    x=pow(x,1./sct);
    return true;
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,ScaleType xsct,ScaleType ysct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,dT &z,ScaleType xsct,ScaleType ysct,ScaleType zsct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
    rescale(z,zsct);
  }

  template <typename dT>
  static void unrescale(dT &x, ScaleType sct)
  {
    if (sct==ScaleTypeV::LOGARITHMIC) x=exp(x);
    else if (sct>=2) x=pow(x,sct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,ScaleType xsct,ScaleType ysct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,dT &z,ScaleType xsct,ScaleType ysct,ScaleType zsct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
    unrescale(z,zsct);
  }
};

#include "../internal/namespace.footer"
#endif

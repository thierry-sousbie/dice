#ifndef __BASIC_SIMPLEX_DATA_POLICIES_HXX__
#define __BASIC_SIMPLEX_DATA_POLICIES_HXX__


#include "../../internal/namespace.header"

namespace simplexInitDataPolicy {

  template <typename T, int N, class M, class S, class DT>
  struct Zero
  {
    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {
      for (int i=0;i<N;++i) 
	result[i] = T();
    }
  };

  template <typename T, int N, class M, class S, class DT>
  struct Copy
  {
    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {
      for (int i=0;i<N;++i) 
	result[i] = (initVal==NULL)?T():initVal[i];
    }  
  };

  template <typename T, int N, class M, class S, class DT>
  struct VolumeWeighted
  {
    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {
      double v=mesh->computeVolume(simplex);
      for (int i=0;i<N;++i) 
	result[i]=(initVal==NULL)?(v):(v*initVal[i]);
    }
  };

  template <typename T, int N, class M, class S, class DT>
  struct ProjectedVolumeWeighted
  {
    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {
      double v=mesh->computeProjectedVolume(simplex);
      for (int i=0;i<N;++i) 
	result[i]=(initVal==NULL)?(v):(v*initVal[i]);
    }
  };

  template <typename T, int N, class M, class S, class DT>
  struct Barycenter
  {
    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {    
      typedef typename S::Vertex Vertex;
      typedef typename Vertex::Coord Coord;
      const double factor = 1.0L/S::NVERT;
      const Coord *refCoords=simplex->getVertex(0)->getCoordsConstPtr();    

      std::copy(refCoords,refCoords+N,result);
      
      for (int i=1;i<S::NVERT;++i) 
	{
	  Vertex *vertex = simplex->getVertex(i);
	  const Coord *coords=vertex->getCoordsConstPtr();
	  
	  for (int j=0;j<N;++j)
	    result[j] += mesh->getGeometry()->checkCoordConsistency(coords[j],refCoords[j],j);
	}
      for (int i=0;i<N;++i) result[i]*=factor;
      mesh->getGeometry()->checkBoundary(result);
    }  
  };

  template <typename T, int N, class M, class S, class DT>
  struct SegTracers
  {
    typedef typename S::Coord Coord;
    typedef typename S::SegmentHandle SegmentHandle;

    static void init(T *result, const M *mesh, const S *simplex, const DT *initVal)
    {     
      for (int i=0;i<S::NSEG;++i)
	{
	  const SegmentHandle sh=const_cast<S*>(simplex)->getSegmentHandle(i);
	  T *res = &result[S::NDIM_W*i];

	  const Coord *c0=sh->getVertex(0)->getCoordsConstPtr();
	  const Coord *c1=sh->getVertex(1)->getCoordsConstPtr();
	  
	  std::copy_n(c0,S::NDIM_W,res);
	  for (int j=0;j<S::NDIM_W;++j)
	    {
	      res[j] += mesh->getGeometry()->checkCoordConsistency(c1[j],c0[j],j);
	      res[j] *= 0.5;
	    }
	  mesh->getGeometry()->checkBoundary(res);
	}
    }
  };


}

namespace simplexRefineDataPolicy {
  
  template <int PASS, typename T, int N, class M, class V, class S>
  struct Dummy
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex,
		       const T* refValue, void *buffer)
    {}
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct Copy
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, 
		       const T* refValue, void *buffer)
    {
      for (int i=0;i<N;++i)
	result[i]=refValue[i];
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct HalfHalf
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, 
		       const T* refValue, void *buffer)
    {
      for (int i=0;i<N;++i)
	result[i]=0.5*refValue[i];
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct VolumeWeighted
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, 
		       T* refValue, void *buffer)
    {
      double v1 = mesh->computeVolume(simplex);
      double v2 = mesh->computeVolume(otherSimplex);
      double fac = (v1==0)?0:(v1/(v1+v2));

      for (int i=0;i<N;++i)
	{
	  result[i]=refValue[i]*fac;
	  refValue[i]-=result[i];
	}
    }
  };

  template <typename T, int N, class M, class V, class S>
  struct VolumeWeighted<1,T,N,M,V,S>
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex,T* refValue,
		       void *buffer)
    {
      for (int i=0;i<N;++i)
	result[i]=refValue[i];
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct ProjectedVolumeWeighted
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, 
		       T* refValue, void *buffer)
    {
      double v1 = mesh->computeProjectedVolume(simplex,false);
      double v2 = mesh->computeProjectedVolume(otherSimplex,false);
      double fac = (v1==0)?0:(v1/(v1+v2));

      for (int i=0;i<N;++i)
	{
	  result[i]=refValue[i]*fac;
	  refValue[i]-=result[i];
	}
    }
  };

  template <typename T, int N, class M, class V, class S>
  struct ProjectedVolumeWeighted<1,T,N,M,V,S>
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex,T* refValue,
		       void *buffer)
    {
      for (int i=0;i<N;++i)
	result[i]=refValue[i];
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct ProjectedDensityWeighted
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, 
		       T* refValue, void *buffer)
    {
      double v1 = mesh->computeProjectedVolume(simplex,false);
      double v2 = mesh->computeProjectedVolume(otherSimplex,false);
      double d1=0;
      double d2=0;
      
      for (int i=0;i<S::NVERT;++i)
	{
	  d1+=simplex->getVertex(i)->projectedDensity.getValue();
	  d2+=otherSimplex->getVertex(i)->projectedDensity.getValue();
	}
      d1/=S::NVERT;d2/=S::NVERT;
      double m1=v1*d1;
      double m2=v2*d2;
      
      double fac = (*refValue)/(m1+m2);
      
      (*result)=m1*fac;
      (*refValue)=m2*fac;    
    }
  };

  template <typename T, int N, class M, class V, class S>
  struct ProjectedDensityWeighted<1,T,N,M,V,S>
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex,T* refValue,
		       void *buffer)
    {
      for (int i=0;i<N;++i)
	result[i]=refValue[i];
    }
  };
  

  template <int PASS, typename T, int N, class M, class V, class S>
  struct Barycenter
  {
    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, T* refValue,
		       void *buffer)
    {    
      typedef typename S::Vertex Vertex;
      typedef typename Vertex::Coord Coord;
      const double factor = 1.0L/S::NVERT;
      const Coord *refCoords=simplex->getVertex(0)->getCoordsConstPtr();    

      std::copy(refCoords,refCoords+N,result);
      
      for (int i=1;i<S::NVERT;++i) 
	{
	  Vertex *vertex = simplex->getVertex(i);
	  const Coord *coords=vertex->getCoordsConstPtr();
	  
	  for (int j=0;j<N;++j)
	    result[j] += mesh->getGeometry()->checkCoordConsistency(coords[j],refCoords[j],j);
	}
      for (int i=0;i<N;++i) result[i]*=factor;
      mesh->getGeometry()->checkBoundary(result);
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct SegmentsBarycenter
  {
    static const int NDIM=S::NDIM;
    static const int NDIM_W=S::NDIM_W;
    typedef typename V::Coord Coord;
    typedef typename S::SegmentHandle SegmentHandle;

    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, T* refValue,
		       void *buffer)
    {
      for (int i=0;i<S::NSEG;++i)
	{
	  const SegmentHandle sh=const_cast<S*>(simplex)->getSegmentHandle(i);
	  T *res = &result[S::NDIM_W*i];

	  const Coord *c0=sh->getVertex(0)->getCoordsConstPtr();
	  const Coord *c1=sh->getVertex(1)->getCoordsConstPtr();
	  
	  std::copy_n(c0,S::NDIM_W,res);
	  for (int j=0;j<S::NDIM_W;++j)
	    {
	      res[j] += mesh->getGeometry()->checkCoordConsistency(c1[j],c0[j],j);
	      res[j] *= 0.5;
	    }
	  mesh->getGeometry()->checkBoundary(res);
	}
    }
  };

  template <int PASS, typename T, int N, class M, class V, class S>
  struct RegressionWithTracerAndSegTracers_segTracers
  {
    static const int NDIM=S::NDIM;
    static const int NDIM_W=S::NDIM_W;
    typedef typename V::Coord Coord;
    typedef typename S::SegmentHandle SegmentHandle;

    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, T* refValue,
		       void *buffer)
    {    
      T *tBuffer=static_cast<T*>(buffer) + (2*NDIM_W);   
      V **v=static_cast<V**>(static_cast<void*>(&tBuffer[S::NVERT*NDIM_W]));
      
      // Slower but more generic ...
      for (int i=0;i<S::NVERT;++i)
	{
	  int sid=simplex->findSegmentIndex(v[i],newVertex);
	  if (sid>=0)
	    std::copy_n(&tBuffer[i*NDIM_W],NDIM_W,&result[sid*NDIM_W]);	      
	}   
    }
  };  
  
  template <int PASS, typename T, int N, class M, class V, class S>
  struct RegressionWithTracerAndSegTracers_tracer
  {
    static const int NDIM=S::NDIM;
    static const int NDIM_W=S::NDIM_W;
    typedef typename V::Coord Coord;

    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, T* refValue,
		       void *buffer)
    {        
      T *tBuffer=static_cast<T*>(buffer);      
      V **v=static_cast<V**>(static_cast<void*>(&tBuffer[2*NDIM_W + S::NVERT*NDIM_W]));
      
      if (simplex->getVertexIndex(v[S::NVERT])<0)
	std::copy_n(&tBuffer[NDIM_W],NDIM_W,result);
      else
	std::copy_n(tBuffer,NDIM_W,result);      
    }  
  };
  
  template <int PASS, typename T, int N, class M, class V, class S>
  struct RegressionWithTracer
  {
    static const int NDIM=S::NDIM;
    static const int NDIM_W=S::NDIM_W;
    typedef typename V::Coord Coord;

    static void refine(T *result, const M *mesh, const V* newVertex, 
		       const S* simplex, const S* otherSimplex, T* refValue,
		       void *buffer)
    {    
      // return Barycenter<PASS,T,N,M,V,S>::refine(result,mesh,newVertex,
      // 						simplex,otherSimplex,
      // 						refValue,buffer);
      T *tBuffer=static_cast<T*>(buffer);
      V *v0=*static_cast<V**>(static_cast<void*>(&tBuffer[2*NDIM_W]));
      
      // The first tracer in buffer belongs to the simplex containing vertex v0
      if (simplex->getVertexIndex(v0)<0)
	std::copy_n(&tBuffer[NDIM_W],NDIM_W,result);
      else 
	std::copy_n(tBuffer,NDIM_W,result);      

      /*
      Coord base[NDIM][NDIM_W];
      Coord origin[NDIM_W];
      mesh->getBaseVectors(simplex,base);
      std::fill_n(origin,NDIM_W,0);      
      for (int j=0;j<S::NVERT;++j) 
	{
	  const Coord *coords=simplex->getVertex(j)->getCoordsConstPtr();
	  
	  for (int k=0;k<NDIM_W;++k)
	    origin[k] += coords[k];
	}
      for (int j=0;j<NDIM_W;++j) origin[j]/=S::NVERT;

      Coord tmp[NDIM_W];
      std::copy_n(result,NDIM_W,tmp);
      for (int i=0;i<NDIM_W;++i) tmp[i]-=origin[i];
      mesh->getGeometry()->normalize(tmp);
      for (int i=0;i<NDIM;++i)
	{
	  mesh->getGeometry()->normalize(&base[i][0]);
	  double dot = mesh->getGeometry()->dot_noCheck(tmp,base[i]);
	  double angle=acos(dot)*180.0/3.14159;
	  if (angle<85)
	    printf("Angle = %e\n",angle);
	}
      */
    }  
  };
}

namespace simplexCoarsenDataPolicy {
  
  template <typename T, int N>
  struct Dummy
  {
    static void coarsen(T *keepValue, T* rmValue)
    {}
  };

  template <typename T, int N>
  struct Keep
  {
    static void coarsen(T *keepValue, T* rmValue)
    {}
  };

  template <typename T, int N>
  struct Copy
  {
    static void coarsen(T *keepValue, T* rmValue)
    {
      for (int i=0;i<N;++i) 
	keepValue[i]=rmValue[i];
    }
  };
  
  template <typename T, int N>
  struct Add
  {
    static void coarsen(T *keepValue, T* rmValue)
    {
      for (int i=0;i<N;++i) 
	keepValue[i]+=rmValue[i];
    }
  };

  template <typename T, int N>
  struct Average
  {
    static void coarsen(T *keepValue, T* rmValue)
    {
      for (int i=0;i<N;++i) 
	keepValue[i]=(keepValue[i]+rmValue[i])*0.5;
    }
  };

}

#include "../../internal/namespace.footer"
#endif

#ifndef __INTERNAL_PROJECTOR_WEIGHT_FUNCTOR_HXX__
#define __INTERNAL_PROJECTOR_WEIGHT_FUNCTOR_HXX__

#include <cstdlib>

#include "../../tools/wrappers/standardDrand48ReentrantWrapper.hxx"
#include "../../geometry/simplexUniformSampler.hxx"

#include "../../internal/namespace.header"

namespace internal {

  struct Tag_ProjectAll{};
  struct Tag_ExcludeIfTagged{};

  template <int order, class M, class WF, class WDF, class T, class HT, 
	    class SimplexCache=Tag_ProjectAll> 
  class ProjectorWeightFunctorT;
  
  // Implementation without possiblity to exclude tagged simplices
  template <class M, class WF, class T, class HT>
  class ProjectorWeightFunctorT<0,M,WF,WF,T,HT,Tag_ProjectAll>
  {
  public:
    static const int ORDER = 0;
    typedef typename M::Simplex Simplex;
    typedef typename M::Vertex Vertex;

    typedef ProjectorWeightFunctorT<0,M,WF,WF,T,HT,Tag_ProjectAll> MyType;

    ProjectorWeightFunctorT(M *mesh, const WF &wf_, int nThreads, bool init=true):
      m(mesh),wf(wf_)
    {
      // Precompute the weight of all local simplices and assign it to the cache
      if (init)
	mesh->template visitSimplices<MyType>(*this,true,false,false,nThreads);	
      seed.resize(nThreads);
      long initSeed=1321;
      for (int i=0;i<seed.size();++i) 
	DRand48_rWapper::srand48_r(initSeed*i,&seed[i]);
    }

    template <typename C, typename OUT = T>
    OUT compute(Simplex *simplex,const C* coords) const 
    {return get<OUT>(simplex);}
    
    void operator()(Simplex *simplex, int th) const {simplex->cache.d=wf(simplex);}    
    template <typename OUT=T>
    OUT get(Simplex *simplex) const {return hlp::numericStaticCast<OUT>(simplex->cache.d);}
    template <typename OUT>
    int getGradient(const Simplex *simplex, OUT *result) const {return 0;}
    
    static char getTag(Simplex *simplex){return 0;}
    static void setTag(Simplex *simplex, char tag){}

    template <typename CT, typename IT1, typename IT2>
    long sample(Simplex *simplex, CT N, IT1 coordsOut, IT2 valueOut) const
    {
      int th=omp_get_thread_num();
      long count = 
	m->template generateRandomSample<Simplex,CT,IT1,Simplex::NDIM>
	(simplex,N,coordsOut,&seed[th]);
      double m=getMass(simplex)/count;
      std::fill_n(valueOut,count,m);
      return count;
      /*
      typedef typename Simplex::Coord Coord;
      Coord base[Simplex::NVERT-1][Simplex::NDIM];
      
      m->template getBaseVectors<Simplex,Coord,Simplex::NDIM>(simplex,base);  
      const Coord *p0=simplex->getVertex(0)->getCoordsConstPtr();
      
      long count=SimplexUniformSamplerT<Simplex::NVERT-1,Simplex::NDIM>::
	generate(p0,base,N,coordsOut);   

      // Normalize the values so that the sum is the total mass !
      double m=getMass(simplex)/count;
      std::fill_n(valueOut,count,m);
      return count;
      */
    }

    double getMass(const Simplex *simplex) const
    {
      return m->computeProjectedVolume(simplex)*wf(simplex);
    }

  protected:
    const M *m;
    const WF &wf;
    mutable std::vector<typename DRand48_rWapper::RandData> seed;
  };

  template <class M, class WF, class WDF,class T, class HT>
  class ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ProjectAll>
  {
  public:
    static const int ORDER = 1;
    typedef typename M::Simplex Simplex;
    typedef typename M::Vertex Vertex;

    static const bool useWeightVector=(sizeof(double)<sizeof(T));

    typedef ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ProjectAll> MyType;

    ProjectorWeightFunctorT(M *mesh, const WF &wf_, const WDF &wdf_, int nThreads,
			    bool init=true):
      m(mesh),
      wf(wf_),
      wdf(wdf_)
    {
      if (useWeightVector) {weightVector.assign(mesh->getNSimplices(),0);}
      // Precompute the weight of all local simplices and assign it to the cache
      if (init)	  
	mesh->template visitSimplices<MyType>(*this,true,false,false,nThreads);
	
      seed.resize(nThreads);
      long initSeed=1321;
      for (int i=0;i<seed.size();++i) 
	DRand48_rWapper::srand48_r(initSeed*i,&seed[i]);

      
    }
    
    template <typename C, typename OUT = T>
    OUT compute(Simplex *simplex,const C* coords) const 
    {
      double grad[M::NDIM];      
      //const typename M::Coord* c0 = simplex->getVertex(0)->getCoordsConstPtr();
      wdf(simplex,grad);
      HT wei=wf(simplex->getVertex(0));
      for (int i=0;i<M::NDIM;++i)
	{
	  HT g=grad[i];
	  HT c=coords[i];
	  wei-=g*c;
	}
      return hlp::numericStaticCast<OUT>(wei);
    }
    
    void operator()(Simplex *simplex, int th) const 
    {
      typedef typename M::Coord Coord;
      if (useWeightVector)
	{
	  weightVector[simplex->getLocalIndex()]=
	    compute<Coord,T>(simplex,simplex->getVertex(0)->getCoordsConstPtr());
	}
      else
	{
	  simplex->cache.d=
	    compute<Coord,double>(simplex,simplex->getVertex(0)->getCoordsConstPtr());
	}

      /*
      double grad[M::NDIM];      
      const typename M::Coord* c0 = simplex->getVertex(0)->getCoordsConstPtr();
      wdf(simplex,grad);

      if (useWeightVector)
	{
	  weightVector[simplex->getLocalIndex()]=wf(simplex->getVertex(0));
	  for (int i=0;i<M::NDIM;++i)
	    {
	      T g=grad[i];
	      T c=c0[i];
	      weightVector[simplex->getLocalIndex()]-=g*c;
	    }
	}
      else
	{
	  simplex->cache.d=wf(simplex->getVertex(0));
	  for (int i=0;i<M::NDIM;++i)
	    simplex->cache.d -= grad[i]*c0[i];
	}
      */    
    }

    template <typename OUT = T>
    OUT get(Simplex *simplex) const 
    {
      if (useWeightVector)
	return hlp::numericStaticCast<OUT>(weightVector[simplex->getLocalIndex()]);
      else
	return hlp::numericStaticCast<OUT>(simplex->cache.d);
    }

    template <typename OUT>
    int getGradient(Simplex *simplex, OUT *result) const
    {
      return getGradientT<OUT>
	(simplex,result,typename hlp::SameType<OUT,double>::Result());      
    }

    // Generate a sampling of the weight inside the simplex
    template <typename CT, typename IT1, typename IT2>
    long sample(Simplex *simplex, CT N, IT1 *coordsOut, IT2 *valueOut) const
    {
      int th=omp_get_thread_num();
      double value[Simplex::NVERT];
      for (int i=0;i<Simplex::NVERT;++i)
	value[i]=wf(simplex->getVertex(i));
      long count = 
	m->template generateRandomSample<Simplex,double,CT,IT1*,IT2*,Simplex::NDIM>
	(simplex,value,N,coordsOut,valueOut,&seed[th]);
      double factor=0;
      for (int i=0;i<count;++i)
	factor+=valueOut[i];
      factor = getMass(simplex)/factor;
      for (int i=0;i<count;++i)
	valueOut[i]*=factor;
      return count;

      /*
      typedef typename Simplex::Coord Coord;
      Coord base[Simplex::NVERT-1][Simplex::NDIM];
      
      m->template getBaseVectors<Simplex,Coord,Simplex::NDIM>(simplex,base);  
      const Coord *p0=simplex->getVertex(0)->getCoordsConstPtr();
      double value[Simplex::NVERT];
      // double grad[Simplex::NDIM];
      // getGradient(simplex,grad);
      for (int i=0;i<Simplex::NVERT;++i)
	value[i]=wf(simplex->getVertex(i));
      
      long count=SimplexUniformSamplerT<Simplex::NVERT-1,Simplex::NDIM>::
	generate(p0,base,value,N,coordsOut,valueOut);
      
      // Normalize the values so that the sum is the total mass !
      double factor=0;
      for (int i=0;i<count;++i)
	factor+=valueOut[i];
      factor = getMass(simplex)/factor;
      for (int i=0;i<count;++i)
	valueOut[i]*=factor;
      return count;
*/
    }

    static char getTag(Simplex *simplex){return 0;}
    static void setTag(Simplex *simplex, char tag) {}

    double getMass(const Simplex *simplex) const
    {
      double mass=wf(simplex->getVertex(0));
      for (int i=1;i<Simplex::NVERT;++i) mass+=wf(simplex->getVertex(i));
      return m->computeProjectedVolume(simplex)*(mass/Simplex::NVERT);
    }

  protected:
    const M *m;
    const WF &wf;
    const WDF &wdf;
    mutable std::vector<typename DRand48_rWapper::RandData> seed;
    mutable std::vector<T> weightVector;

    template <class TT>
    int getGradientT(Simplex *simplex, TT *result, const hlp::IsTrue) const
    {
      return wdf(simplex,result);
    }

    template <class TT>
    int getGradientT(Simplex *simplex, TT *result, const hlp::IsFalse) const
    {
      double tmp[M::NDIM];
      int n=wdf(simplex,tmp);
      std::copy_n(tmp,M::NDIM,result);
      return n;
    }
  };
  
  // Implementation with possiblity to exclude tagged simplices
  template <class M, class WF, class WDF, class T, class HT>
  class ProjectorWeightFunctorT<0,M,WF,WDF,T,HT,Tag_ExcludeIfTagged> :
    public ProjectorWeightFunctorT<0,M,WF,WDF,T,HT,Tag_ProjectAll>
  {
    typedef ProjectorWeightFunctorT<0,M,WF,WDF,T,HT,Tag_ProjectAll> Base;
  public:
    typedef typename M::Simplex Simplex;
    typedef typename M::Vertex Vertex;
    typedef ProjectorWeightFunctorT<0,M,WF,WDF,T,HT,Tag_ExcludeIfTagged> MyType;

    ProjectorWeightFunctorT(M *mesh, const WF &wf_, int nThreads, bool dontTag=false):
      Base(mesh,wf_,nThreads,false)
    {
      // Precompute the weight of all local simplices and assign it to the cache
      cache.assign(mesh->getNSimplices(),0);
      if (!dontTag)
	mesh->template visitSimplices<MyType>(*this,true,false,false,nThreads);
    }

    template <typename C, typename OUT = T>
    OUT compute(Simplex *simplex,const C* coords) const 
    {
      if (simplex->isLocal()&&(cache[simplex->getLocalIndex()]==0))
	return Base::template compute<C,OUT>(simplex,coords);
      else return hlp::numericStaticCast<OUT>(0);      
    }

    void operator()(Simplex *simplex, int th) const
    {
      Base::operator()(simplex,th);
      if (simplex->isLocal()&&(simplex->cache.c[0]!=0))
	cache[simplex->getLocalIndex()]=simplex->cache.c[0];
    }

    template <typename OUT=T>
    OUT get(Simplex *simplex) const 
    {
      if (simplex->isLocal()&&(cache[simplex->getLocalIndex()]==0))
	return Base::template get<OUT>(simplex);
      else
	return hlp::numericStaticCast<OUT>(0);
    }

    template <typename OUT>
    int getGradient(Simplex *simplex, OUT *result) const 
    {return Base::getGradient(simplex,result);}
  
    char getTag(Simplex *simplex) const
    {
      if (!simplex->isLocal()) 
	return -1;
      else
	return cache[simplex->getLocalIndex()];
    }

    void setTag(Simplex *simplex, char tag) const
    {
      if (simplex->isLocal()) 
	cache[simplex->getLocalIndex()]|=tag;
    }

  private:
    mutable std::vector<char> cache;
  };

  template <class M, class WF, class WDF, class T, class HT>
  class ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ExcludeIfTagged> :
    public ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ProjectAll>
  {
    typedef ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ProjectAll> Base;
  public:
    typedef typename M::Simplex Simplex;
    typedef typename M::Vertex Vertex;
    typedef ProjectorWeightFunctorT<1,M,WF,WDF,T,HT,Tag_ExcludeIfTagged> MyType;

    ProjectorWeightFunctorT(M *mesh, const WF &wf_, const WDF &wdf_, int nThreads, 
			    bool dontTag=false):
      Base(mesh,wf_,wdf_,nThreads,false)
    {
      // Precompute the weight of all local simplices and assign it to the cache
      cache.assign(mesh->getNSimplices(),0);
      if (!dontTag)
	mesh->template visitSimplices<MyType>(*this,true,false,false,nThreads);
    }

    template <typename C, typename OUT = T>
    OUT compute(Simplex *simplex,const C* coords) const 
    {
      if (simplex->isLocal()&&(cache[simplex->getLocalIndex()]==0))
	return Base::template compute<C,OUT>(simplex,coords);
      else return hlp::numericStaticCast<OUT>(0);
    }

    void operator()(Simplex *simplex, int th) const 
    {      
      Base::operator()(simplex,th);
      if (simplex->isLocal()&&(simplex->cache.c[0]!=0))
	cache[simplex->getLocalIndex()]=simplex->cache.c[0];
      /*
      if ((!simplex->isLocal())||(simplex->cache.c[0]==0))
	Base::operator()(simplex,th);
      else
	{
	  cache[simplex->getLocalIndex()]=simplex->cache.c[0];
	  simplex->cache.d=0;
	}
      */
    }

    template <typename OUT=T>
    OUT get(Simplex *simplex) const 
    {
      if (simplex->isLocal()&&(cache[simplex->getLocalIndex()]==0))
	return Base::template get<OUT>(simplex);
      else return hlp::numericStaticCast<OUT>(0);
    }

    template <typename OUT>
    int getGradient(Simplex *simplex, OUT *result) const 
    {
      if (simplex->isLocal()&&(cache[simplex->getLocalIndex()]==0))
	Base::getGradient(simplex,result);
      else
	std::fill_n(result,M::NDIM,0);
      return M::NDIM;
    }

    char getTag(Simplex *simplex) const
    {
      if (!simplex->isLocal()) 
	return -1;
      else
	return cache[simplex->getLocalIndex()];
    }

    void setTag(Simplex *simplex, char tag) const
    {
      if (simplex->isLocal()) 
	cache[simplex->getLocalIndex()]|=tag;
    }

  private:
    mutable std::vector<char> cache;
  };
}

#include "../../internal/namespace.footer"
#endif

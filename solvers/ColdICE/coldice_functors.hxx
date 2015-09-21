#ifndef __COLDICE_FUNCTORS_HXX__
#define __COLDICE_FUNCTORS_HXX__

#include <functional>
#include <iostream>
#include <map>

#include <string.h>
#include <dice/mesh/cellData/cellDataFunctors_interface.hxx>

// Use to identify negative valued voxels in AMR grids
template <class T=double>
class LessOrNanT:
  public std::unary_function<T,bool>
{
public:
  LessOrNanT(T ref=T()):refVal(ref){}
  bool operator() (const T &a) const {return (a<refVal)||(a!=a);}
private:
  const T refVal;
};

// Computes projected density of a simplex with constant density distribution
template <class M>
class ProjectedDensityFunctorT: 
  public dice::cellDataFunctors::InterfaceT<M,typename M::Simplex>
{
public:
  typedef dice::cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
  typedef typename Base::Cell Cell;        
  ProjectedDensityFunctorT(const M* m, int flags=dice::cellDataFunctors::F_NO_FLAG):
    Base(m,flags){}    
  double get(const Cell *c) const
  {
    double result=Base::mesh->template computeProjectedVolume<Cell,long double>(c);
    if (result!=0) result=c->mass.getValue()/result; // 1.0L
    return result;
  }
  double get(const Cell *c, int n) const {return get(c);}
  std::string getName() const {return std::string("projectedDensity");}
  int getSize() const {return 1;}
};

// Computes projected density of a simplex with linear density distribution (i.e. defined at
// vertices)
template <class M>
class SmoothedProjectedDensityFunctorT: 
  public dice::cellDataFunctors::InterfaceT<M,typename M::Simplex>
{
public:
  typedef dice::cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
  typedef typename Base::Cell Cell;  
  typedef typename M::Vertex Vertex;
  typedef typename M::Simplex Simplex;  
  SmoothedProjectedDensityFunctorT(const M* m,int flags=dice::cellDataFunctors::F_NO_FLAG):
    Base(m,flags){}    
  double get(const Cell *c) const
  {    
    static const double factor = 1.0L/(M::NDIM+1);
    double result=0;
    for (int i=0;i<Simplex::NVERT;++i)
      result+=c->getVertex(i)->projectedDensity.getValue();
    return result*factor;
  }
  double get(const Cell *c, int n) const {return get(c);}
  std::string getName() const {return std::string("smoothedProjectedDensity");}
  int getSize() const {return 1;}
};

// Associate the index of the non refined original simplex to any simplex that belongs to
// its refinement
/*
template <class M>
class DomainIndexFunctorT: 
  public dice::cellDataFunctors::InterfaceT<M,typename M::Simplex>
{
public:
  typedef dice::cellDataFunctors::InterfaceT<M,typename M::Simplex> Base;
  typedef typename Base::Cell Cell;  
  typedef typename M::Vertex Vertex;
  typedef typename M::Simplex Simplex;  
  DomainIndexFunctorT(const M* m,int flags=dice::cellDataFunctors::F_NO_FLAG):
    Base(m,flags){}    
  double get(const Cell *c) const
  {    
    return static_cast<double>(c->getConstRoot()->getLocalIndex()) + 
      static_cast<double>(Base::mesh->getMpiCom()->rank())/
      Base::mesh->getMpiCom()->size();
  }
  double get(const Cell *c, int n) const {return get(c);}
  std::string getName() const {return std::string("domainIndex");}
  int getSize() const {return 1;}
};
*/

template <class M, long N=100>
class FilterBadSimplicesT
{
  typedef typename M::Simplex Simplex;
public:
  FilterBadSimplicesT(M *mesh, int nThreads=dice::glb::num_omp_threads)
  {
    this->mesh=mesh;

    anisotropy.resize(nThreads);
    minHeight.resize(nThreads);
    volume.resize(nThreads);

    anisotropyF = mesh->getSimplexFunctorPtr("projectedAnisotropy");
    minHeightF = mesh->getSimplexFunctorPtr("projectedMinHeight");
    volumeF = mesh->getSimplexFunctorPtr("signedProjectedVolume");
  }

  void operator()(Simplex *s, int th) const
  {
    std::map<double,Simplex *>& anisotropy=this->anisotropy[th];
    std::map<double,Simplex *>& minHeight=this->minHeight[th];
    std::map<double,Simplex *>& volume=this->volume[th];

    double val;
    val=(*anisotropyF)(s);
    if ((anisotropy.size()<N)||(anisotropy.begin()->first <  val))
      {
	anisotropy.insert(std::make_pair(val,s));
	if (anisotropy.size()>N) anisotropy.erase(anisotropy.begin());
      }
    
    val=(*minHeightF)(s);
    if ((minHeight.size()<N)||(minHeight.rbegin()->first >  val))
      {
	minHeight.insert(std::make_pair(val,s));
	if (minHeight.size()>N)
	  minHeight.erase(std::prev(minHeight.end()));
      }
    
    val=std::abs((*volumeF)(s));
    if ((volume.size()<N)||(volume.rbegin()->first >  val))
      {
	volume.insert(std::make_pair(val,s));	
	if (volume.size()>N) 
	  volume.erase(std::prev(volume.end()));
      }
  }

  void print()
  {
    for (int i=1;i<anisotropy.size();++i)
      {
	anisotropy[0].insert(anisotropy[i].begin(),anisotropy[i].end());
	minHeight[0].insert(minHeight[i].begin(),minHeight[i].end());
	volume[0].insert(volume[i].begin(),volume[i].end());
      }
    anisotropy.resize(1);
    minHeight.resize(1);
    volume.resize(1);
    while (anisotropy[0].size()>N) {anisotropy[0].erase(anisotropy[0].begin());}
    /*
    anisotropy[0].resize((anisotropy.size()<N)?anisotropy.size():N);
    minHeight[0].resize((minHeight.size()<N)?minHeight.size():N);
    volume[0].resize((volume.size()<N)?volume.size():N);
    */
    std::map<double,Simplex *>& anisotropy=this->anisotropy[0];
    std::map<double,Simplex *>& minHeight=this->minHeight[0];
    std::map<double,Simplex *>& volume=this->volume[0];

    std::cout << std::endl << "Anisotropy (minHeight/Volume):" << std::endl;
    int count=0;
    for (auto it=anisotropy.rbegin();(it!=anisotropy.rend())&&(count<N);++it,++count)
      std::cout << it->first <<"("<<(*minHeightF)(it->second) <<","<<(*volumeF)(it->second)<<")"<<std::endl;


    std::cout << std::endl << "MinHeight (anisotropy/volume):" << std::endl;
    count=0;
    for (auto it=minHeight.begin();(it!=minHeight.end())&&(count<N);++it,++count)
      std::cout << it->first <<"("<<(*anisotropyF)(it->second) <<","<<(*volumeF)(it->second)<<")"<< std::endl;

    std::cout << std::endl << "Volume (anisotropy/minHeight):" << std::endl;
    count=0;
    for (auto it=volume.begin();(it!=volume.end())&&(count<N);++it,++count)
      std::cout << it->first <<"("<<(*anisotropyF)(it->second) <<","<<(*minHeightF)(it->second)<<")"<< std::endl;
  }
  
private:
  M *mesh;
  mutable std::vector< std::map<double,Simplex *> > anisotropy;
  mutable std::vector< std::map<double,Simplex *> > minHeight;
  mutable std::vector< std::map<double,Simplex *> > volume;

  const typename M::SimplexFunctor *anisotropyF;
  const typename M::SimplexFunctor *minHeightF;
  const typename M::SimplexFunctor *volumeF;
};

#endif

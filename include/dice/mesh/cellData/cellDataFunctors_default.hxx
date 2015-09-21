#ifndef __CELL_DATA_FUNCTORS_DEFAULT_HXX__
#define __CELL_DATA_FUNCTORS_DEFAULT_HXX__

#include <limits>
#include "../../mesh/cellData/cellDataFunctors_interface.hxx"

#include "../../internal/namespace.header"

namespace cellDataFunctors {  

  template <class M, int W>
  class VertexDataT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;     
    typedef typename Cell::DataElementInfo DataInfo;
    VertexDataT(const M* m, int flags=F_NO_FLAG):
      Base(m,flags),
      dataInfo(Cell::template DeclaredCellData<W>::getDataInfo())
    {}
    double get(const Cell *c) const
    {return (double)c->template getDataElementPtr<W>()->getValue();}     
    double get(const Cell *c, int n) const
    {return (double)c->template getDataElementPtr<W>()->getValueAt(n);}     
    std::string getName() const {return dataInfo.name;}
    int getSize() const {return Cell::template DeclaredCellData<W>::Result::SIZE;}
  private:
    const DataInfo dataInfo;
  };

  template <class M, int W>
  class SimplexDataT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;     
    typedef typename Cell::DataElementInfo DataInfo;
    SimplexDataT(const M* m, int flags=F_NO_FLAG):
      Base(m,flags),
      dataInfo(Cell::template DeclaredCellData<W>::getDataInfo())
    {}
    double get(const Cell *c) const
    {return (double)c->template getDataElementPtr<W>()->getValue();}
    double get(const Cell *c, int n) const
    {return (double)c->template getDataElementPtr<W>()->getValueAt(n);}
    std::string getName() const {return dataInfo.name;}
    int getSize() const {return Cell::template DeclaredCellData<W>::Result::SIZE;}
  private:
    const DataInfo dataInfo;
  };

  // data definitions
  template <class M>
  class VertexFlagsT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;        
    VertexFlagsT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getFlags();}
    std::string getName() const {return std::string("flags");}
    int getSize() const {return 1;}
  };

  template <class M>
  class SimplexFlagsT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    SimplexFlagsT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getFlags();}
    std::string getName() const {return std::string("flags");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class VertexGenerationT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;  
    typedef typename M::Vertex::GlobalIdentity GlobalIdentity;
    VertexGenerationT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const
    {
      GlobalIdentity id=c->getGeneration();
      if (id.rank()==0) 
	return static_cast<double>(id.id());
      else
	return static_cast<double>(-(int)id.rank());
    }
    std::string getName() const {return std::string("generation");} 
    int getSize() const {return 1;}
  };

  template <class M>
  class VertexRankT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;        
    VertexRankT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getGlobalIdentity().rank();}
    std::string getName() const {return std::string("rank");} 
    int getSize() const {return 1;}
  };

  template <class M>
  class SimplexRankT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    SimplexRankT(const M* m, int flags=F_NO_FLAG):Base(m,flags)
    {
      myRank=m->getMpiCom()->rank();
    } 
    double get(const Cell *c, int n) const 
    {return (double)c->getGlobalIdentity(myRank).rank();}
    std::string getName() const {return std::string("rank");}
    int getSize() const {return 1;}
  private:
    int myRank;
  };

  template <class M>
  class VertexGlobalIndexT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;        
    VertexGlobalIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getGlobalIdentity().id();}
    std::string getName() const {return std::string("globalIndex");}
    int getSize() const {return 1;}
  };

  template <class M>
  class SimplexGlobalIndexT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    SimplexGlobalIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags)
    {
      myRank=m->getMpiCom()->rank();
    }    
    double get(const Cell *c, int n) const 
    {return (double)c->getGlobalIdentity(myRank).id();}
    std::string getName() const {return std::string("globalIndex");}
    int getSize() const {return 1;}
  private:
    int myRank;
  }; 

  template <class M>
  class VertexLocalIndexT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;        
    VertexLocalIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getLocalIndex();}
    std::string getName() const {return std::string("localIndex");}
    int getSize() const {return 1;}
  };

  template <class M>
  class SimplexLocalIndexT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    SimplexLocalIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const {return (double)c->getLocalIndex();}
    std::string getName() const {return std::string("localIndex");}
    int getSize() const {return 1;}
  }; 
  /*
  template <class M>
  class DomainIndexT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;      
    DomainIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c) const
    {    
      if (c->isShadowOrGhost())
	return -1;
      else
	return static_cast<double>(c->getConstRoot()->getLocalIndex()) + 
	  static_cast<double>(Base::mesh->getMpiCom()->rank())/
	  Base::mesh->getMpiCom()->size();
    }
    double get(const Cell *c, int n) const {return get(c);}
    std::string getName() const {return std::string("domainIndex");}
    int getSize() const {return 1;}
  };
  */

  template <class M>
  class DomainIndexT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;      
    DomainIndexT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c) const
    {    
      return (c->isShadow())?0:(double)c->getGeneration().id();
    }
    double get(const Cell *c, int n) const {return get(c);}
    std::string getName() const {return std::string("domainIndex");}
    int getSize() const {return 1;}
  };

  template <class M>
  class LevelT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;      
    LevelT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c) const
    {    
      return (c->isShadow())?0:(double)c->getGeneration().rank();
    }
    double get(const Cell *c, int n) const {return get(c);}
    std::string getName() const {return std::string("level");}
    int getSize() const {return 1;}
  };

  template <class M>
  class WeightT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    WeightT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const 
    {return (c->isShadowOrGhost())?0:(double)c->getWeight();}
    std::string getName() const {return std::string("weigth");}
    int getSize() const {return 1;}
  };

  template <class M>
  class DepthT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    DepthT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const
    {return (c->isShadowOrGhost())?0:(double)c->getDepth();}
    std::string getName() const {return std::string("depth");}
    int getSize() const {return 1;}
  };

  template <class M>
  class ProjectedVolumeT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    ProjectedVolumeT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const
    {return Base::mesh->computeProjectedVolume(c);}
    std::string getName() const {return std::string("projectedVolume");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class VolumeT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;
    VolumeT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}
    double get(const Cell *c, int n) const {return Base::mesh->computeVolume(c);}
    std::string getName() const {return std::string("volume");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class SignedProjectedVolumeT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;        
    SignedProjectedVolumeT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}    
    double get(const Cell *c, int n) const 
    {return Base::mesh->computeProjectedVolume(c,true);}
    std::string getName() const {return std::string("signedProjectedVolume");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class SignedVolumeT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;
    SignedVolumeT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}
    double get(const Cell *c, int n) const
    {return Base::mesh->computeVolume(c,true);}
    std::string getName() const {return std::string("signedVolume");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class ExtraDimsT: 
    public InterfaceT<M,typename M::Vertex>
  {
  public:
    typedef InterfaceT<M,typename M::Vertex> Base;
    typedef typename Base::Cell Cell;
    ExtraDimsT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}
    double get(const Cell *c, int n) const {return c->getCoord(M::NDIM+n);}
    std::string getName() const {return std::string("U");}
    int getSize() const {return M::NDIM_W - M::NDIM;}
  }; 

  template <class M>
  class ProjectedAnisotropyT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;
    ProjectedAnisotropyT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}
    double get(const Cell *c, int n) const
    {      
      typedef typename M::Simplex Simplex;
      
      // compute the simplex extent
      auto fh=c->getConstFacetHandle(0);
      double sMax=Base::mesh->computeProjectedVolume(fh.asPointer());
      double sMin=sMax;

      for (int i=1;i<Simplex::NFACET;++i)
	{//const_cast<Cell*>(c)
	  fh=c->getConstFacetHandle(i);
	  double tmp=Base::mesh->computeProjectedVolume(fh.asPointer());
	  if (tmp>sMax) sMax=tmp;
	  if (tmp<sMin) sMin=tmp;
	}
      if (sMin==0) return std::numeric_limits<double>::max();
      else return sMax/sMin;
    }
    std::string getName() const {return std::string("projectedAnisotropy");}
    int getSize() const {return 1;}
  }; 

  template <class M>
  class ProjectedMinHeightT: 
    public InterfaceT<M,typename M::Simplex>
  {
  public:
    typedef InterfaceT<M,typename M::Simplex> Base;
    typedef typename Base::Cell Cell;
    ProjectedMinHeightT(const M* m, int flags=F_NO_FLAG):Base(m,flags){}
    double get(const Cell *c, int n) const
    {      
      typedef typename M::Simplex Simplex;
      
      // compute the simplex extent
      auto fh=c->getConstFacetHandle(0);
      double sMax=Base::mesh->computeProjectedVolume(fh.asPointer());

      for (int i=1;i<Simplex::NFACET;++i)
	{//const_cast<Cell*>(c)
	  fh=c->getConstFacetHandle(i);
	  double tmp=Base::mesh->computeProjectedVolume(fh.asPointer());	  
	  if (tmp>sMax) sMax=tmp;
	}

      double v=Base::mesh->computeProjectedVolume(c);
      return (v*M::NDIM)/sMax;
    }
    std::string getName() const {return std::string("projectedMinHeight");}
    int getSize() const {return 1;}
  }; 


} //namespace

#include "../../internal/namespace.footer"
#endif

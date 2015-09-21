#ifndef __LOCAL_AMR_GRID_VISITORS_HXX__
#define __LOCAL_AMR_GRID_VISITORS_HXX__

#include "../internal/namespace.header"

namespace localAmrGridVisitors
{

  template <class AMR, class F>
  class CountIfT 
  {
    typedef typename AMR::Voxel Voxel;
    typedef typename AMR::ICoord ICoord;
    typedef typename AMR::Data Data;
    static const long NDIM=AMR::NDIM;
  public:

    CountIfT():
      count(0)      
    {}

    CountIfT(const F &func):
      functor(func),
      count(0)
    {}

    static void initialize(Voxel *rootVoxel) {}

    bool visit(Voxel *voxel)
    { 
      if (voxel->isLeaf())
	{
	  if (functor(voxel->data))
	    {
#pragma omp atomic
	      count++;
	    }
	  return false;
	}
      return true;
    }

    static bool visit(Voxel *voxel, int i)
    {return true;}

    static void visited(Voxel *voxel, int i=0) 
    {}
      
    long getCount() {return count;}
  private:
    F functor;
    //const Data val;
    long count;
  };

  template <class AMR, class F>
  class GetVoxelIfT 
  {
    typedef typename AMR::Voxel Voxel;
    typedef typename AMR::ICoord ICoord;
    typedef typename AMR::Data Data;
    static const long NDIM=AMR::NDIM;
  public:

    GetVoxelIfT()
    {}

    GetVoxelIfT(const F &func):
      functor(func)      
    {}

    static void initialize(Voxel *rootVoxel) {}

    bool visit(Voxel *voxel)
    { 
      if (voxel->isLeaf())
	{
	  if (functor(voxel->data))
	    {
#pragma omp critical
	      result.push_back(voxel);
	    }
	  return false;
	}
      return true;
    }

    static bool visit(Voxel *voxel, int i)
    {return true;}

    static void visited(Voxel *voxel, int i=0) 
    {}
      
    const std::vector<Voxel*> &getResult() {return result;}
  private:
    F functor;    
    std::vector<Voxel*> result;
  };

  template <class AMR>
  class SetValueT 
  {
    typedef typename AMR::Voxel Voxel;    
    typedef typename AMR::Data Data;
  public:

    SetValueT(const Data &val=Data()):value(val)
    {}

    static void initialize(Voxel *rootVoxel) {}

    bool visit(Voxel *voxel) const
    { 
      if (voxel->isLeaf())
	{
	  voxel->data = value;
	  return false;
	}
      return true;
    }

    static bool visit(Voxel *voxel, int i)
    {return true;}

    static void visited(Voxel *voxel, int i=0)
    {}
   
  private:
    const Data value;
  };

  template <class AMR>
  class GetMaxLevelT 
  {
    typedef typename AMR::Voxel Voxel;
  public:

    GetMaxLevelT():
      maxLevel(0)      
    {}

    static void initialize(Voxel *rootVoxel) {}
  
    bool visit(Voxel *voxel)
    { 
      if (voxel->isLeaf())
	{
	  if (voxel->getLevel()>maxLevel)
	    maxLevel=voxel->getLevel();
	  return false;
	}
      return true;
    }

    static bool visit(Voxel *voxel, int i)
    {return true;}

    static void visited(Voxel *voxel, int i=0) 
    {}
      
    long getMaxLevel() const {return maxLevel;}
  private:
    int maxLevel;
  };

  template <class AMR, class Visitor>
  class LeavesVisitor
  {
    typedef typename AMR::Voxel Voxel;    
  public:

    LeavesVisitor(Visitor &visitor):
      v(visitor)
    {}

    static void initialize(Voxel *root){}

    bool visit(Voxel *voxel)
    { 
      if (voxel->isLeaf()) 
	{
	  v.visit(voxel);
	  return false;
	}
	 
      return true;
    }

    static bool visit(Voxel *voxel,int i) 
    {
      return true;
    }

    static void visited(Voxel *voxel) {}
    static void visited(Voxel *voxel, int i) {}    
  private:
    Visitor &v;
  };

}

#include "../internal/namespace.footer"
#endif

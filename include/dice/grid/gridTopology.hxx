#ifndef __GRID_TOPOLOGY_HXX__
#define __GRID_TOPOLOGY_HXX__

#include "../dice_globals.hxx"

#include "./regularGridNavigation.hxx"

#include "../internal/namespace.header"

template <int ND>
class GridTopologyT
{
public:
  typedef GridTopologyT<ND> MyType;
  static const int NDIM = ND;

  typedef RegularGridNavigation GridNav;
  typedef typename GridNav::Direction Direction;  

  static std::string classHeader() {return "grid_topology";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  GridTopologyT()
  {}

  GridTopologyT(int N)
  {
    reset(N);
  }

  ~GridTopologyT()
  {}

  void clear()
  {
    topology.clear();
  }

  void reset(int N)
  {
    clear();
    topology.resize(N);
  }

  unsigned long size() const
  {
    return topology.size();
  }

  template <class Slicer>
  void setNeighbors(int index, const Slicer &slicer)
  {
    int dim,dir,id;

    for (int j=0;j<slicer.neighborsCount();++j)
      {
	slicer.getNeighborInfo(j,dim,dir,id);
	setNeighbor(index,dim,dir,id);
      }
  }
  
  template <class Slicer>
  void setFromSlicer(const typename Slicer::Params &params,
		     long size,int nThreads=glb::num_omp_threads)
  {
    Slicer slicer;
    reset(size);
    for (long i=0;i<size;++i)
      {
	slicer.slice(params, i, size, nThreads);
	setNeighbors(i,slicer);
      }
  }

  void setNeighbor(int index, int dim, int dir, int neiIndex)
  {
    if (dir>=0)
      topology[index].nei[dim][1] = neiIndex;
    else
      topology[index].nei[dim][0] = neiIndex;
  }

  int getNeighbor(int index, int dim, int dir) const
  {
    int result;
    if (dir>=0)
      result=topology[index].nei[dim][1];
    else
      result=topology[index].nei[dim][0];
    return result;    
  }

  int getNeighbor(int index, Direction dir) const
  {
    for (int i=0;i<NDIM;++i)
      {
	if (index>=0)
	  {
	    int d = GridNav::getDir(dir,i);
	    if (d>0) index = topology[index].nei[i][1];
	    else if (d<0) index = topology[index].nei[i][0];
	  }
      }
    return index;
  }

  template <class R>
  void read(R* reader)
  {
    float version;
    R::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);

    clear();
    unsigned long sz;
    
    reader->read(&sz);
    topology.resize(sz);
    for (unsigned long i=0;i<sz;++i)
      reader->read(&(topology[i].nei[0][0]),NDIM*2);
  }

  template <class W>
  void write(W* writer) const
  {
    writer->writeHeader(classHeader(),classVersion());

    unsigned long sz = topology.size();
    writer->write(&sz);
    for (unsigned long i=0;i<sz;++i)
      writer->write(&(topology[i].nei[0][0]),NDIM*2);
  }

private:

  struct Neighborhood
  {
    static const int NDIM = MyType::NDIM;

    Neighborhood()
    {
      std::fill(&nei[0][0],&nei[0][0]+2*NDIM,-1);
    }

    int nei[NDIM][2];
  };

  std::vector<Neighborhood> topology;
};

#include "../internal/namespace.footer"
#endif

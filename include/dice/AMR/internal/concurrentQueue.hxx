#ifndef __INTERNAL_CONCCURENT_QUEUE_HXX__
#define __INTERNAL_CONCCURENT_QUEUE_HXX__

#include <queue>
#include "../../tools/OMP/openMP_interface.hxx"

#include "../../internal/namespace.header"

namespace internal {

  class ConccurentQueue {
  public:
    ConccurentQueue():
      nextAvail(-1)
    {}

    void clear()
    {
#pragma omp critical(CQUEUE)
      {	
	enroled.clear();
	std::queue<int> empty;
	q.swap(empty);
      }
    }

    void resize(long sz)
    {
#pragma omp critical(CQUEUE)
      {
	enroled.resize(sz,0);
      }
    }

    long nWaiting() const
    {
      return q.size();
    }

    bool isEnroled(int index) const
    {
      return enroled[index];
    }

    bool enrol(int index)
    {
      if (enroled[index]) return false;
#pragma omp critical(CQUEUE)
      {	
	q.push(index);
	if (q.size()==1) nextAvail=index;
	enroled[index]=1;
      }
      return true;
    }

    int availableFor() const
    {            
      return nextAvail;      
    }
    
    bool done(int index)
    {
      if (index==nextAvail)
	{
#pragma omp critical(CQUEUE)
	  {	    
	    q.pop();
	    if (q.size()==0) nextAvail=-1;
	    else nextAvail=q.front();
	    enroled[index]=0;
	  }
	  return true;
	}
      return false;
    }

  private:
    volatile int nextAvail;
    std::queue<int> q;
    std::vector<int> enroled;
  };

}

#include "../../internal/namespace.footer"
#endif

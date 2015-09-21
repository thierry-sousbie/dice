#ifndef __GET_BALL_HELPER_HXX__
#define __GET_BALL_HELPER_HXX__

#include "../../internal/namespace.header"

namespace internal {
  namespace mesh {
  
    // Gets simplices adjacent to a vertex
    // front and newFront are vectors of simplices
    // Ball is a set of simplices
    template <class S, class V, class SV, class B>
    void getBall(V *vertex, SV &front, SV &newFront, B &ball)
    {     
      newFront.clear();
      for (auto it=front.begin();it!=front.end();++it)
	{
	  for (int i=0;i<S::NNEI;++i)
	    {
	      auto nei=(*it)->getNeighbor(i);
					  
	      if (((*it)->getVertex(i)!=vertex)&&
		  (ball.find(nei)==ball.end()))
		{
		  newFront.push_back(nei);
		  ball.insert(nei);
		}					
	    }
	}
      if (newFront.size()>0)
	{
	  front.swap(newFront);
	  getBall<S>(vertex,front,newFront,ball);
	}
    };

  }
}

#include "../../internal/namespace.footer"
#endif

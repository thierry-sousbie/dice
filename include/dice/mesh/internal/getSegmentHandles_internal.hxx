#ifndef __GET_SEGMENT_HANDLES_INTERNAL_HXX__
#define __GET_SEGMENT_HANDLES_INTERNAL_HXX__

#include "../../internal/namespace.header"

namespace internal {

  template <int ND, class SEGH, class SEG>
  struct GetSegmentHandle
  {
    template <class V, class S>
    static inline SEGH get(int i, S* me, V *vertices)
    {
      int j=0;    
      while (i>=(ND-j))
	{
	  i-=(ND-j);
	  j++;
	};
      return SEGH(SEG(vertices[j],vertices[j+i+1],me));
    }

    static inline std::pair<int,int> getIndices(int i)
    {
      int j=0;    
      while (i>=(ND-j))
	{
	  i-=(ND-j);
	  j++;
	};
     
      return std::make_pair(j,j+i+1);
    }
  };

  template <class SEGH, class SEG>
  struct GetSegmentHandle<1,SEGH,SEG>
  {
    template <class V, class S>
    static inline SEGH get(int i, S* me, V *vertices)
    {
      return SEGH(SEG(vertices[0],vertices[1],me));
    }

    static inline std::pair<int,int> getIndices(int i)
    {
      return std::make_pair(0,1);
    }
  };

  template <class SEGH, class SEG>
  struct GetSegmentHandle<2,SEGH,SEG>
  {
    template <class V, class S>
    static inline SEGH get(int i, S* me, V *vertices)
    {
      switch (i)
	{
	case 0: return SEGH(SEG(vertices[0],vertices[1],me));
	case 1: return SEGH(SEG(vertices[0],vertices[2],me));
	case 2: return SEGH(SEG(vertices[1],vertices[2],me));
	}
      return SEGH(SEG());
    }

    static inline std::pair<int,int> getIndices(int i)
    {
      switch (i)
	{
	case 0: return std::make_pair(0,1);
	case 1: return std::make_pair(0,2);
	case 2: return std::make_pair(1,2);
	}
      return std::make_pair(0,0);      
    }
  };

  template <class SEGH, class SEG>
  struct GetSegmentHandle<3,SEGH,SEG>
  {
    template <class V, class S>
    static inline SEGH get(int i, S* me, V *vertices)
    {
      switch (i)
	{
	case 0: return SEGH(SEG(vertices[0],vertices[1],me));
	case 1: return SEGH(SEG(vertices[0],vertices[2],me));
	case 2: return SEGH(SEG(vertices[0],vertices[3],me));
	case 3: return SEGH(SEG(vertices[1],vertices[2],me));
	case 4: return SEGH(SEG(vertices[1],vertices[3],me));
	case 5: return SEGH(SEG(vertices[2],vertices[3],me)); 
	}
      return SEGH(SEG());     
    }

    static inline std::pair<int,int> getIndices(int i)
    {
      switch (i)
	{
	case 0: return std::make_pair(0,1);
	case 1: return std::make_pair(0,2);
	case 2: return std::make_pair(0,3);
	case 3: return std::make_pair(1,2);
	case 4: return std::make_pair(1,3);
	case 5: return std::make_pair(2,3);
	}
      return std::make_pair(0,0);      
    }
  };

} // internal

#include "../../internal/namespace.footer"
#endif

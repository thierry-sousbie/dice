#ifndef __FACET_IMPLEMENTATION_HXX__
#define __FACET_IMPLEMENTATION_HXX__

#include "../../internal/namespace.header"

namespace internal {

  template <int ND>
  struct FacetImplementationT;

  template <>
  struct FacetImplementationT<2>
  {
    template <class S>
    static typename S::Vertex *getVertex(S *simplex, int index, int which)
    {
      static int id[3][2]={{1,2},{2,0},{0,1}};
      return simplex->getVertex(id[index][which]);
    }   
    
    template <class S, class OutputIterator>
    static void getVertices(S *simplex, int index, OutputIterator out)
    {
      static int id[3][2]={{1,2},{2,0},{0,1}};

      (*out)=simplex->getVertex(id[index][0]);++out;
      (*out)=simplex->getVertex(id[index][1]);
    }

    template <class C, class OutputIterator>
    static void computeProjectedNormal(const C vCoord[2][2], OutputIterator out)
    {
      (*out)=vCoord[1][1]-vCoord[0][1];++out;
      (*out)=vCoord[0][0]-vCoord[1][0];
    }
  };

  template <>
  struct FacetImplementationT<3>
  {
    template <class S>
    static typename S::Vertex *getVertex(S *simplex, int index, int which)
    {
      static int id[4][3]={{1,2,3},{3,2,0},{0,1,3},{2,1,0}};
      return simplex->getVertex(id[index][which]);
    }   

    template <class S, class OutputIterator>
    static void getVertices(S *simplex, int index, OutputIterator out)
    {
      static int id[4][3]={{1,2,3},{3,2,0},{0,1,3},{2,1,0}};

      (*out)=simplex->getVertex(id[index][0]);++out;
      (*out)=simplex->getVertex(id[index][1]);++out;
      (*out)=simplex->getVertex(id[index][2]);      
    }

    template <class C, class OutputIterator>
    static void computeProjectedNormal(const C vCoord[3][3], OutputIterator out)
    {           
      (*out)=
	(vCoord[1][1]-vCoord[0][1])*(vCoord[2][2]-vCoord[0][2])-
	(vCoord[1][2]-vCoord[0][2])*(vCoord[2][1]-vCoord[0][1]);
      ++out;
      (*out)=
	(vCoord[1][2]-vCoord[0][2])*(vCoord[2][0]-vCoord[0][0])-
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][2]-vCoord[0][2]); 
      ++out;
      (*out)=
	(vCoord[1][0]-vCoord[0][0])*(vCoord[2][1]-vCoord[0][1])-
	(vCoord[1][1]-vCoord[0][1])*(vCoord[2][0]-vCoord[0][0]);
      
    }
  };

}

#include "../../internal/namespace.footer"
#endif

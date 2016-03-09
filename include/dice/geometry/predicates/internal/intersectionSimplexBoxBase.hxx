#ifndef INTERSECTION_SIMPLEX_BASE_HXX__
#define INTERSECTION_SIMPLEX_BASE_HXX__

#include "../../geometricProperties.hxx"

#include "../../../tools/wrappers/boostMultiprecisionFloat128.hxx"

#include "../pointInSimplex.hxx"
#include "../pointInBox.hxx"
#include "../intersectionSegmentFacet.hxx"


#ifdef HAVE_BOOST
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#endif
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#endif


#include "../../../internal/namespace.header"

namespace internal {
  template <int ND, int filterType>
  class IntersectionSimplexBoxBaseT
  {
    typedef predicate::PointInBoxT<ND,filterType> PointInBox;
    typedef predicate::PointInSimplexT<ND,filterType> PointInSimplex;
    typedef predicate::IntersectionSegmentFacetT<ND,filterType> IntersectionSegmentFacet;

    static const int VERT_PER_FACE = (1<<(ND-1));
  public:
    
    template <class T>
    static int test(const T vCoord[ND+1][ND], const T box[2][ND])
    {
      // check simplex vertices in box
      for (int i=0;i<ND+1;++i)
	if (PointInBox::test(box,vCoord[i])) return 1;

      // check box vertices in simplex, we store points on a signle facet of the cube
      T point[VERT_PER_FACE][ND];
      for (int j=0;j<VERT_PER_FACE;++j)
	{
	  for (int k=0;k<(ND-1);++k)
	    {
	      if (j&(1<<k))
		point[j][k+1]=box[1][k+1];
	      else
		point[j][k+1]=box[0][k+1];
	    }
		
	  point[j][0]=box[1][0];
	  if (PointInSimplex::test(vCoord,point[j])) return 1;
	  
	  point[j][0]=box[0][0];
	  if (PointInSimplex::test(vCoord,point[j])) return 1;
	}

      // Get all the facets
      T fCoord[ND+1][ND][ND];
      for (int k=0;k<ND+1;++k)
	for (int j=0,i=0;j<ND+1;++j)
	  if (j!=k) std::copy_n(vCoord[k],ND,fCoord[k][i++]);     

      // Check facet / voxel edge intersections
      for (int f=0;f<VERT_PER_FACE;f++)
	{
	  for (int i=0;i<ND;++i)
	    if (box[1][i] != point[f][i])
	      {
		// test segments departing from points on the stored box facet
		for (int j=0;j<ND+1;++j)
		  if (IntersectionSegmentFacet::test(fCoord[j],point[f],i,box[1][i]))
		    return 1;
		
		// and also test segments departing from the other facet
		if (i>0) 
		{
		  point[f][0]=box[1][1];
		  
		  for (int j=0;j<ND+1;++j)
		    if (IntersectionSegmentFacet::test(fCoord[j],point[f],i,box[1][i]))
		      return 1;

		  point[f][0]=box[1][0];
		}
	      }
	}
      
      return 0;
    }

    template <class T, class G>
    static int test(const T vCoord[ND+1][ND], const T box[2][ND], const G* geometry)  
    {
      
      // check simplex vertices in box
      for (int i=0;i<ND+1;++i)
	if (PointInBox::test(box,vCoord[i],geometry)) return 1;
      
      // check box vertices in simplex, we store points on a signle facet of the cube
      T point[VERT_PER_FACE][ND];
      for (int j=0;j<VERT_PER_FACE;++j)
	{
	  for (int k=0;k<(ND-1);++k)
	    {
	      if (j&(1<<k))
		point[j][k+1]=box[1][k+1];
	      else
		point[j][k+1]=box[0][k+1];
	    }
		
	  point[j][0]=box[1][0];
	  if (PointInSimplex::test(vCoord,point[j],geometry)) return 1;
	  
	  point[j][0]=box[0][0];
	  if (PointInSimplex::test(vCoord,point[j],geometry)) return 1;
	}
      
      
      // Get all the facets
      T fCoord[ND+1][ND][ND];
      for (int k=0;k<ND+1;++k)
	for (int j=0,i=0;j<ND+1;++j)
	  if (j!=k) std::copy_n(vCoord[k],ND,fCoord[k][i++]);     

      // Check facet / voxel edge intersections
      for (int f=0;f<VERT_PER_FACE;f++)
	{
	  for (int i=0;i<ND;++i)
	    if (box[1][i] != point[f][i])
	      {
		// test segments departing from points on the stored box facet
		for (int j=0;j<ND+1;++j)
		  if (IntersectionSegmentFacet::test(fCoord[j],point[f],i,box[1][i],geometry))
		    return 1;
		
		// and also test segments departing from the other facet
		if (i>0) 
		{
		  point[f][0]=box[1][0];
		  
		  for (int j=0;j<ND+1;++j)
		    if (IntersectionSegmentFacet::test(fCoord[j],point[f],i,box[1][i],geometry))
		      return 1;

		  point[f][0]=box[1][0];
		}
	      }
	}
      
      return 0;
    }
    
  };

} // namespace internal
 
#include "../../../internal/namespace.footer"
#endif

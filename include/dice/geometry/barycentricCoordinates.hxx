#ifndef __BARYCENTRIC_COORDINATES_HXX__
#define __BARYCENTRIC_COORDINATES_HXX__

#include "../geometry/geometricProperties.hxx"
#include "../tools/helpers/helpers.hxx"

#include "./internal/barycentricCoordinates_implementation.hxx"

/**
 * @file 
 * @brief  A class to compute barycentric coordinates of a simplex and test if a point 
 * is inside/outside/on boundary
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Geometry
 *   \{
 */

/**
 * \class BarycentricCoordinatesT
 * \brief A class to compute barycentric coordinates of a simplex and test if a point 
 * is inside/outside/on boundary.
 *
 * LIMITATIONS: for now, the simplex cannot be embeded in a higher dimensional space.
 *
 * \tparam ND the number of dimensions.
 * \tparam GP a class to handle boundary conditions (see GeometricPropertiesT)
 * \todo Embed in NDW dimensions space with NDW != ND
 */

template < int ND , typename CT=long double, class GP = DummyGeometricPropertiesT<ND> >
class BarycentricCoordinatesT:
  protected internal::BarycentricCoordinatesBaseT<ND,CT>
{
public:
  typedef GP GeometricProperties;
  typedef double Coord;
  
  static const int NDIM = ND;
  
  /** \brief An enum representing different possible locations of a point with respect 
   *  to a simplex 
   */
  enum Location {
    Inside=0,             /**< The point is inside the simplex */
    Outside=(1<<0),       /**< The point is outside the simplex */
    OnVertex=(1<<1),      /**< The point falls on a vertex of the simplex */
    OnSegment=(1<<2),     /**< The point falls on a segment of the (2+)-simplex */
    OnTriangle=(1<<3),    /**< The point falls on a triangle of the (3+)-simplex*/
    OnTetrahedron=(1<<4), /**< The point falls on a tetrahedron of the (4+)-simplex */
    On4Simplex=(1<<5),    /**< The point falls on a 4-simplex face of the (5+)-simplex */
    On5Simplex=(1<<6),    /**< The point falls on a 5-simpelx face of the (6+)-simplex */
    OnBoundary=(1<<1)|(1<<2)|(1<<3)|(1<<4)|(1<<5)|(1<<6) 
    /**< The point falls on a boundary of any sort of the simplex (if 'location' is 
     *   of type Location, a point is on boundary if location&OnBoundary is true)
     */
  };

  /** 
   *  \param coords The ND+1 coordinates of the simplex vertices
   *  \param geometry_ boundary conditions geometry handling class 
   *         (see GeometricPropertiesT)
   */
  template <class T>
  BarycentricCoordinatesT(const T* const coords[ND+1],const GP &geometry_):
    geometry(geometry_)
  {
    //geometry=geometry_;
    std::copy(coords[0], coords[0]+ND, origin);

    T base[ND][ND];
    geometry.template getBaseVectors<T,T,ND,ND>(coords,base);   
    degenerate = Base::computeCoeffs(base);
  }

  /** \param coords The ND+1 coordinates of the simplex vertices
   *  \param geometry_ boundary conditions geometry handling class 
   *         (see GeometricPropertiesT)
   */
  template <class T>
  BarycentricCoordinatesT(const T coords[ND+1][ND],const GP &geometry_):
    geometry(geometry_)
  {
    //geometry=geometry_;
    std::copy(coords[0], coords[0]+ND, origin);

    T base[ND][ND];
    geometry.template getBaseVectors<T,T,ND,ND>(coords,base);   
    degenerate = Base::computeCoeffs(base);
  }
  

  /** \param coords The ND+1 coordinates of the simplex vertices  
   */
  template <class T>
  BarycentricCoordinatesT(const T* const coords[ND+1]):
    geometry(GP())
  {
    //geometry=geometry_;
    std::copy(coords[0], coords[0]+ND, origin);

    T base[ND][ND];
    geometry.template getBaseVectors<T,T,ND,ND>(coords,base);   
    degenerate = Base::computeCoeffs(base);
  }

  /** \param coords The ND+1 coordinates of the simplex vertices  
   */
  template <class T>
  BarycentricCoordinatesT(const T coords[ND+1][ND]):
    geometry(GP())
  {
    //geometry=geometry_;
    std::copy(coords[0], coords[0]+ND, origin);

    T base[ND][ND];
    geometry.template getBaseVectors<T,T,ND,ND>(coords,base);   
    degenerate = Base::computeCoeffs(base);
  }

  /* \param simplex A simplex (see SimplexT)  
   *  \tparam S A simplex class with method getVerticesCoordsConstPtr()
   */
  /*
  template <class S>
  BarycentricCoordinatesT(const S* simplex):
    geometry(GP())
  {     
    typedef typename S::Coord SCoord;
    SCoord base[ND][ND];
    const SCoord *p[S::NVERT];
    simplex->getVerticesCoordsConstPtr(p);
    geometry.template getBaseVectors<SCoord,ND,ND>(p,base);
    std::copy(p[0],p[0]+ND,origin);  
    degenerate = Base::computeCoeffs(base);
  }
  */

  /** \brief Compute the (NDIM+1) barycentric coordinates
   *  \param coords The cartesian coordinates of the point
   *  \param[out] lambda The barycentric coordinates of the point
   */
  template <class T>
  void computeLambda(const T coords[], T lambda[]) const
  {
    T vec[NDIM];
    geometry.template getVector<T,ND>(origin,coords,vec);    
    Base::computeLambda(vec,lambda);
  }

  /** \brief Compute the gradient of a function over the simplex
   *  \param[in] fVal The value of the function at each vertex
   *  \param[out] gradient The gradient of the function over the simplex
   *  \return true if the simplex is degenerate, in which case the gradient 
   *  is set to 0
   */
  template <class T, class T2>
  bool computeGradient(const T fVal[], T2 gradient[]) const //,bool debug=false) const
  {
    if (degenerate)
      {
	for (int i=0;i<NDIM;++i)
	  gradient[i]=0;
      }
    else
      {
	/*
	if (debug)
	  {
	    for (int i=0;i<NDIM;++i) gradient[i]=0;
	for (int i=0;i<NDIM;++i)
	  for (int j=0;j<NDIM;++j)
	    {
	      printf("[%d %d]: gradient[%d]= %g += %g * (%g-%g)\n",i,j,j,gradient[j],
		     Base::matrix[i][j],fVal[i+1],fVal[0]);
		     
	      gradient[j]+=Base::matrix[i][j]*(fVal[i+1]-fVal[0]);
	    }
	  }
	*/
	/*
	for (int i=0;i<NDIM;++i)
	  {
	    gradient[i]=0;
	    for (int j=0;j<NDIM;++j)
	      gradient[i]+=Base::matrix[i][j]*(fVal[j+1]-fVal[0]);
	  }
	*/

	CT tmp[NDIM];
	for (int i=0;i<NDIM;++i) tmp[i]=0;
	for (int i=0;i<NDIM;++i)
	  for (int j=0;j<NDIM;++j)
	    tmp[j]+=Base::matrix[i][j]*(fVal[i+1]-fVal[0]);

	for (int i=0;i<NDIM;++i)
	  gradient[i]=hlp::numericStaticCast<T2>(tmp[i]);
	
      }
    return degenerate;
  }


  /** \brief Get the location of a point with respect to the simplex. Use getLocation 
   *  instead if the barycentric coordinates are already known.
   *  \param coords Coordinates of the point
   *  \param tolerance how close do we have to be from a boundary to be on it ?
   *  \return a Location type value, indicating if the point is inside the simplex 
   *   (Location::Inside), outside (Location::Outside) or on one of its 
   *   vertex/segment/face ... (Location::Vertex, Location::Segment, ....)
   *   The point is on the boundary (of any type) if (return&OnBoundary) is true.
   * \todo We miss a consistant way to define the tolerance ...
   * \warning This function maynot be used as a geometric predicate as its result are not
   * guaranteed to be consistent !
   */
  template <class T>
  Location getLocationFromCoords(const T *coords, double tolerance=0) const
				 
  //std::numeric_limits<double>::epsilon()*10000.0F) const
  {
    T lambda[NDIM+1];
    computeLambda(coords,lambda);
    return getLocation(lambda,tolerance);
  }

  /** \brief Get the location of a point with respect to the simplex. This is quite 
   *  inprecise and should not be used as a geometric predicate
   *  \param lambda Barycentric coordinates of the point
   *  \param tolerance how close do we have to be from a boundary to be on it ?
   *  \return a Location type value, indicating if the point is inside the simplex 
   *   (Location::Inside), outside (Location::Outside) or on one of its 
   *   vertex/segment/face ... (Location::Vertex, Location::Segment, ....)
   *   The point is on the boundary (of any type) if (return&Location::Boundary) is true.
   * \todo We miss a consistent way to define the tolerance ...
   * \warning This function maynot be used as a geometric predicate as its result are not
   * guaranteed to be consistent !
   */
  template <class T>
  Location getLocation(const T* lambda, double tolerance=0) const
		       // double tolerance=
		       // std::numeric_limits<double>::epsilon()*10000.0F) const
  {
    Location result=Inside;

    int nZeroes=0;
    int nOnes=0;
    int nOut=0;
    int nIn=0;

    for (int i=0;i<=NDIM;++i)
      {
	if (fabs(lambda[i])<=tolerance) nZeroes++;
	else if (fabs(lambda[i]-1)<=tolerance) nOnes++;
	else if ((lambda[i]<0)||(lambda[i]>1)) nOut++;
	else nIn++;
      }
    
    if (nIn<=NDIM)
      {
	if (nOut>0) result=Outside;
	else result=Location(1<<(NDIM+1-nZeroes));
      }

    return result;      
  }

  /** \brief check if the simplex is degenerate
   *  \return true if the simplex is degenerate, false otherwise
   */
  bool isDegenerate()
  {
    return degenerate;
  }

private:
  typedef internal::BarycentricCoordinatesBaseT<ND,CT> Base;
    
  const GeometricProperties &geometry;
  Coord origin[NDIM];
  bool degenerate;
};

/** \}*/
#include "../internal/namespace.footer"
#endif

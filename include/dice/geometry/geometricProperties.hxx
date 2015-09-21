#ifndef __GEOMETRIC_PROPERTIES_HXX__
#define __GEOMETRIC_PROPERTIES_HXX__

#include <cmath>        // for std::abs
#include <limits>
#include "../tools/helpers/helpers.hxx"
#include "../geometry/boundaryType.hxx"

/**
 * @file 
 * @brief  Definiton of a class for computing generic geometric properties with different 
 * boundary conditions
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup Geometry
 *   \{
 */

/**
 * \class GeometricPropertiesT
 * \brief A class for computing generic geometric properties under different 
 * boundary conditions.
 *
 * \tparam ND  Number of spatial dimensions
 * \tparam NDW Total number of dimensions in the embedding space
 * \tparam BD  Boundary conditions along the first ND dimensions (of type 
 * BoundaryType::Type)
 * \tparam BW  Boundary conditions in the NDW-ND last embeding dimensions (of type 
 * BoundaryType::Type)
 */

template <typename CT, int ND, int NDW, int BD,int BW>
class GeometricPropertiesT
{

public:
  typedef CT Coord;
  static const int BOUNDARY_TYPE = BD;
  static const int WORLD_BOUNDARY_TYPE = BW;
  static const int NDIM=ND;
  static const int NDIM_W=NDW;

  /** 
   * \param x0_ An iterator to the vector of leftmost coordinates of the bounding box 
   * along each of the NDW dimension.
   * \param delta_ An iterator to the size of the bounding box along each NDW dimension
   * \tparam InputIterator Input iterator type
   */
  template <class InputIterator>
  GeometricPropertiesT(InputIterator x0_,InputIterator delta_)
  {
    /*
    x0.assign(x0_,x0_+NDIM_W);
    xmax.assign(x0_,x0_+NDIM_W);
    //xmid.assign(x0_,x0_+NDIM_W);
    delta.assign(delta_,delta_+NDIM_W);
    halfDelta.assign(delta_,delta_+NDIM_W);

    for (int i=0;i<NDIM_W;i++) 
      {
	halfDelta[i]=0.5*delta[i];
	//xmid[i]+=halfDelta[i];
	xmax[i]+=delta[i];
      }
*/
    for (int i=0;i<NDIM_W;i++) 
      {
	x0[i]=*x0_;++x0_;
	delta[i]=*delta_;++delta_;	
	halfDelta[i]=0.5*delta[i];
	xmax[i]=x0[i]+delta[i];
      }
  }
  friend class GeometricPropertiesT<CT,ND-1,NDW-1,BD,BW>;
  GeometricPropertiesT(const GeometricPropertiesT<CT,ND+1,NDW+1,BD,BW> &geom, int dim)
  { 
    // x0.resize(NDIM_W);
    // xmax.resize(NDIM_W);
    // delta.resize(NDIM_W);
    // halfDelta.resize(NDIM_W);
    for (int i=0;i<dim;++i)
      {
	x0[i]=geom.x0[i];
	xmax[i]=geom.xmax[i];
	delta[i]=geom.delta[i];
	halfDelta[i]=geom.halfDelta[i];
      }
    for (int i=dim;i<NDIM_W;++i)
      {
	x0[i]=geom.x0[i+1];
	xmax[i]=geom.xmax[i+1];
	delta[i]=geom.delta[i+1];
	halfDelta[i]=geom.halfDelta[i+1];
      }
  }

  GeometricPropertiesT()
  {}

  CT getEpsilon(int i=0)
  {
    return 0;
  }
  
  
  /** \brief Test if a point with coordinates x is on the boundary of the bounding box.   
   * \param x the coordinates of the point to test
   * \return true if it is on the bounary, false otherwise.
   * \note there is no boundary for peridoic boundary conditions, so this function 
   *  return false in that case.
   */
  template <class T>
  int onBoundary(const T * x) const
  {
    int result=0;
    int out=(1<<NDIM_W)-1;
    for (int i=0; i<NDIM; ++i)
      {	
	if ((x[i]>=xmax[i])||(x[i]<=x0[i]))
	  {
	    if ((x[i]==xmax[i])||(x[i]==x0[i])) 
	      result|=(1<<i);
	    else out=0;
	  }	
      }
    return result&out;
  }

  /** \brief compute the square of the distance between two points |xb-xa|**2
   * \param xa A pointer to the D coordinates of the first point
   * \param xb A pointer to the D coordinates of the second point
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the square of the distance between xa and xb
   */
  template <class T, int D = NDIM_W>
  static T distance2(const T * __restrict xa,const T * __restrict xb)
  {
    return hlp::Distance2<D,T>::compute(xa,xb);   
  }
  
  /** \brief compute the distance between two points |xb-xa|**2
   * \param xa A pointer to the D coordinates of the first point
   * \param xb A pointer to the D coordinates of the second point
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the distance between xa and xb
   */
  template <class T, int D = NDIM_W>
  static T distance(const T * __restrict xa,const T * __restrict xb) 
  {
    return sqrt(distance2<T,D>(xa,xb));
  }

  /** \brief compute the squared norm of a vector
   * \param xa the coordinates of the vector   
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the squared norm of the vector
   */
  template <class T, int D = NDIM_W>
  static T norm2(const T * xa)
  {
    return hlp::Norm2<D,T>::compute(xa);    
    /*
    double nrm=0;
    for (int i=0;i<D;++i)
      nrm+=xa[i]*xa[i];
    return nrm;
    */
  }

  /** \brief compute the norm of a vector
   * \param xa the coordinates of the vector   
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the norm of the vector
   */
  template <class T, int D = NDIM_W>
  static T norm(const T * xa)
  {
    return sqrt(norm2<T,D>(xa));
  }

  /** \brief compute the norm of a vector without checking anything.
   * \param xa the coordinates of the vector   
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the norm of the vector
   */
  template <class T, int D = NDIM_W>
  static T norm_noCheck(const T * xa)
  {
    return sqrt(norm2<T,D>(xa));
  }

  /** \brief Normalize a vector
   * \param[in,out] xa the coordinates of the vector   
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the norm of the vector
   */
  template <class T, int D = NDIM_W>
  static T normalize(T * xa)
  {
    T nrm = norm<T,D>(xa);
    if (nrm!=0)
      {
	T nrm_inv = T(1.0)/nrm;
	for (int i=0;i<D;++i) xa[i]*=nrm_inv;
      }
    return nrm;
  }
  
  /** \brief Normalize a vector without checking anything.
   * \param[in,out] xa the coordinates of the vector   
   * \tparam T data type
   * \tparam D the number of dimensions
   * \return the norm of the vector
   */
  template <class T, int D = NDIM_W>
  static T normalize_noCheck(T * xa)
  {
    return normalize<T,D>(xa);
  }

  /** \brief check if a coordinate is within the bounding box. Always true for periodic 
   * boundaries)   
   */
  template <class T>
  bool inBound(T x, int dim) const
  {
    return ((x>=x0[dim])&&(x<xmax[dim]));
  }

  /** \brief Check that the coordinates of a point lie within the box when boundary 
   * conditions are periodic and update them if necessary. Moreover, the coordinates
   * of the point are changed if they lie at a distance of the boundary less than the size
   * of the box times the numerical precision (being exaclty on the left boundary is OK)
   * The coordinates remain unchanged for non-periodic boundary conditions.
   * \param[in,out] xa the coordinates of the point  
   * \tparam T data type   
   */
  template <class T>
  static void sanitizeBoundary(T *xa)
  {
  }

  /** \brief Check that the coordinates of a point lie within the box when boundary 
   * conditions are periodic and update them if necessary. 
   * The coordinates remain unchanged for non-periodic boundary conditions.
   * \param[in,out] xa the coordinates of the point  
   * \tparam T data type   
   */
  template <class T>
  static void checkBoundary(T *xa)
  {
  }

  /** 
   * \brief Check that the 'which' coordinate \a c of a point lies within the boundaries
   * when the boundary conditions are periodic. The returned coordinate is always \a c if 
   * boundary conditions are non-periodic. 
   * \param c the coordinates of the point along dimensions 'which'
   * \param which the index of the dimensions to check
   * \tparam T data type   
   * \return the corrected coordinate
   */
  template <class T>
  static T checkBoundary(T c, int which)
  {
    return c;
  }

  /** \brief Check that a set of points coordinates (usually a simplex's vertices) are as 
   * close as possible from the coordinates of a reference point. The coordinates 
   * \a vCoords are never changed for non-periodic boundary conditions.
   * \param[in,out] vCoords The D coords of the NVERT points
   * \param refPoint the reference point
   * \return true if a coordinates was modified, false otherwise
   * \tparam T data type 
   * \tparam NVERT the number of points to check
   * \tparam D the number of dimensions for each point   
   */
  template <class T,int  NVERT=NDIM+1, int D=NDIM_W>
  static bool checkCoordsConsistency(T (&vCoords)[NVERT][D],const T * refPoint)
  {    
    return false;
  }

  /** \brief Check that a point coordinate is as close as possible from another 
   * coordinates of a reference point when boundary conditions are periodic. The returned
   * value is always equal to \a vCoord for non-periodic boundary conditions.
   * \param vCoord the coordinate along axis \a dim of the vertex to check
   * \param refPoint the coordinate along axis \a dim of the reference point
   * \param dim the dimension we are checking
   * \return the new value of vCoord
   * \tparam T data type 
   */
  template <class T, class T2>
  static T checkCoordConsistency(T vCoord, const T2 refPoint, int dim)
  { 
    return vCoord;
  } 

  /** \brief Get the sign of the difference of two coords (b-a) along dimension dim
   *  or 0 if equal
   */
  template <class T>
  T coordDiffSign(T a, T b, int dim) const
  { 
    if (b==a) return T(0);
    return (b>a)?T(1):T(-1);
  } 

  /** \brief returns true if coords are consistent (see checkCoordConsistency)
   * \a vCoords are never changed for non-periodic boundary conditions.
   * \param[in,out] vCoords The D coords of the NVERT points
   * \param refPoint the reference point
   * \tparam NVERT the number of points to check
   * \tparam D the number of dimensions for each point
   */
  template <class T, int  NVERT=NDIM+1, int D=NDIM_W>
  bool coordsAreConsistent(const T vCoords[NVERT][D],const T *refPoint) const
  {
    return true;
  }
 
  /** \brief Correct a coordinate difference \a len along axis \a which so that it lies 
   * between -boxSize/2 and boxSize/2 when boundary conditions are periodic. The returned
   * value is always equal to \a len for non-periodic boundary conditions.
   * \param len the coordinate difference along dimensions 'which'
   * \param which the index of the dimensions to check
   * \tparam T data type   
   * \return the corrected difference
   */
  template <class T>
  static T correctCoordsDiff(T len, int which)
  {
    return len;
  }

  /** \brief Displace a point with coordinates xa by vector dx.
   * \param[in,out] xa the coordinates of the point
   * \param dx the coordinates of the displacement vector
   * \tparam T data type     
   */
  template <class T>
  static void displace(T * __restrict xa,const T *__restrict dx)
  {
    for (int i=0;i<NDIM_W;i++) xa[i]+=dx[i];    
  }

  /** \brief Compute the coordinate of the barycenter of two points
   * \param xa the coordinates of the first point
   * \param xb the coordinates of the second point
   * \param[out] out the coordinates of barycenter (mid-point)
   * \tparam T data type     
   */
  template <class T>
  static void midPointCoords(const T * __restrict xa, const T * __restrict xb, 
			     T * __restrict out)
  {  
    for (int i=0;i<NDIM_W;i++)
      {
	out[i] = 0.5*(xa[i]+xb[i]);
      }   
  }

  /** \brief Compute the vector (xb-xa) defined by the difference of two points 
   *   coordinates. The difference and eventual poriodization are achieve using
   *   the type with highest precision.
   * \param xa the coordinates of the first point
   * \param xb the coordinates of the second point
   * \param[out] vec the resulting vector (xb-xa)
   * \tparam T data type  
   * \tparam D the number of dimensions of the vector 
   */
  template <class T, class T2, int D=NDIM_W>
  static void getVector(const T *xa, const T *xb, T2 *vec) 
  {
    typedef typename hlp::MostPreciseFP<T,T2>::Result MPT;
    for (long i=0;i<D;i++) 
      vec[i]= (hlp::numericStaticCast<MPT>(xb[i])-
	       hlp::numericStaticCast<MPT>(xa[i]));
  }

  /** \brief Transform the coordinates of a point into a vector
   * \param xa the coordinates of the point   
   * \param[out] vec the vector representing point xa
   * \tparam T data type  
   * \tparam D the number of dimensions of the vector 
   */
  template <class T, class T2, int D=NDIM_W>
  static void getVector(const T *xa, T2 *vec) 
  {
    for (long i=0;i<D;i++) 
      vec[i]=hlp::numericStaticCast<T2>(xa[i]);
  }

  /** \brief Compute the vector base formed by (D+1) points of dimensions DW. This can 
   * be used to compute the base vectors defined by the (D+1) points of a D-simplex 
   * embeded in a DW-dimensional space.
   * \param p The coordinates of the (D+1) DW-dimensional points
   * \param[out] base the (D) DW-dimensional vectors forming the base.
   * \tparam T data type  
   * \tparam D the number of spatial dimensions (=number of base vector = number of point-1)
   * \tparam DW the number of dimensions in the embedding space (=number of coordinates 
   * for each point)
   */
  template <class T, class T2, int D=NDIM, int DW=NDIM_W>
  static void getBaseVectors(const T * const p[D+1], T2 (&base)[D][DW])
  {    
    const T* p0=p[0]; // __restrict ??
    for (long i=0;i<D;i++) 
      {
	const T* p1=p[i+1]; // __restrict ??
	for (long j=0;j<DW;j++) 
	  {
	    base[i][j]=p1[j]-p0[j];
	    //vec[i][j]=correctCoordsDiff(p1[j]-p0[j],j);
	  }
      }  
  }

  /** \brief Compute the vector base formed by (D+1) points of dimensions DW. This can 
   * be used to compute the base vectors defined by the (D+1) points of a D-simplex 
   * embeded in a DW-dimensional space.
   * \param p The coordinates of the (D+1) DW-dimensional points
   * \param[out] base the (D) DW-dimensional vectors forming the base.
   * \tparam T data type  
   * \tparam D the number of spatial dimensions (=number of base vector = number of point-1)
   * \tparam DW the number of dimensions in the embedding space (=number of coordinates 
   * for each point)
   */
  template <class T, class T2, int D=NDIM, int DW=NDIM_W>
  static void getBaseVectors(const T p[D+1][D], T2 (&base)[D][DW])
  {        
    for (long i=0;i<D;i++) 
      for (long j=0;j<DW;j++) 
	base[i][j]=p[i+1][j]-p[0][j];
  }

  /** \brief Compute the dot product of two vector
   * \param a the coordinates of the first vector
   * \param b the coordinates of the second vector   
   * \tparam T data type  
   * \tparam D the number of dimensions of the vectors
   * \return the result of the dot product
   */
  template <class T, class T2, int D=NDIM_W>
  static T dot(const T *a, const T2 *b)
  {
    T dt=0;
    for (long i=0;i<D;i++) 
      dt+=a[i]*b[i];
    return dt;
  }

  /** \brief Compute the dot product of two vector without checking the validity of 
   *  the vectors
   * \param a the coordinates of the first vector
   * \param b the coordinates of the second vector   
   * \tparam T data type  
   * \tparam D the number of dimensions of the vectors
   * \return the result of the dot product
   */
  template <class T, class T2, int D=NDIM_W>
  static T dot_noCheck(const T *a, const T2 *b)
  {
    return dot<T,T2,D>(a,b);
  }

  /** \brief returns the size of the bounding box along dimension \a i
   */
  CT getBBoxSize(int i=0) const
  {
    return delta[i];
  }

 
  

protected:
  CT x0[NDIM_W];
  CT xmax[NDIM_W];
  CT delta[NDIM_W];
  CT halfDelta[NDIM_W];
  // std::vector<CT> x0;
  // std::vector<CT> xmax;
  // std::vector<CT> delta;
  // std::vector<CT> halfDelta;
};

// FIXME: need to optimize distance for instance ...
template <typename CT,int ND, int NDW, int BW>
class GeometricPropertiesT<CT,ND,NDW,BoundaryType::PERIODIC,BW>
{
private:
  typedef GeometricPropertiesT<CT,ND,NDW,BoundaryType::NONE,BW> 
  GeometricPropertiesNoBoundary;

  //NB: a must be positive !
  template <class T, class U>
  static inline T square_P(const T /*__restrict*/ a, const U /*__restrict*/ h, const U /*__restrict*/ d)
  {return (a>=h)?(hlp::square(a-d)):(hlp::square(a));}

  // Distance2 with periodic boundary conditions
  template <int D, class T, class U>
  struct Distance2_P {                           
    static inline T compute(const T* /*__restrict*/ a, const T* /*__restrict*/ b, 
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return Distance2_P<D-1,T,U>::compute(a+1,b+1,h+1,d+1) + 
	square_P(std::abs(*b-*a),*h,*d); }
  }; 
  template <class T, class U>
  struct Distance2_P<1,T,U> {                        
    static inline T compute(const T* /*__restrict*/ a, const T* /*__restrict*/ b,
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return square_P(std::abs(*b-*a),*h,*d); }                        
  }; 
  template <class T, class U>
  struct Distance2_P<0,T,U> {                        
    static inline T compute(const T* /*__restrict*/ a, const T* /*__restrict*/ b, 
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return 0; }                        
  }; 

  // Norm2 with periodic boundary conditions
  template <int D, class T, class U>
  struct Norm2_P {                           
    static inline T compute(const T* a,
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return Norm2_P<D-1,T,U>::compute(a+1,h+1,d+1) + 
	square_P(std::abs(*a),*h,*d); }
  }; 
  template <class T, class U>
  struct Norm2_P<1,T,U> {                        
    static inline T compute(const T* a,
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return square_P(std::abs(*a),*h,*d); }                        
  }; 
  template <class T, class U>
  struct Norm2_P<0,T,U> {                        
    static inline T compute(const T* a,
			    const U* /*__restrict*/ h, const U* /*__restrict*/ d)
    { return 0; }                        
  }; 

public:
  typedef CT Coord;
  static const int BOUNDARY_TYPE = BoundaryType::PERIODIC;
  static const int WORLD_BOUNDARY_TYPE = BW;
  static const int NDIM=ND;
  static const int NDIM_W=NDW;

  template <class InputIterator>
  GeometricPropertiesT(InputIterator x0_,InputIterator delta_)
  {
    /*
    x0.assign(x0_,x0_+NDIM_W);
    xmax.assign(x0_,x0_+NDIM_W);
    //xmid.assign(x0_,x0_+NDIM_W);
    delta.assign(delta_,delta_+NDIM_W);
    halfDelta.assign(delta_,delta_+NDIM_W);
    for (int i=0;i<NDIM_W;i++) 
      {
	halfDelta[i]=0.5*delta[i];
	//xmid[i]+=halfDelta[i];
	xmax[i]+=delta[i];
      }
*/
    for (int i=0;i<NDIM_W;i++) 
      {
	x0[i]=*x0_;++x0_;
	delta[i]=*delta_;++delta_;
	halfDelta[i]=0.5*delta[i];
	xmax[i]=x0[i]+delta[i];
	epsilon[i]=delta[i]*std::numeric_limits<CT>::epsilon();
      }
  }

  friend class GeometricPropertiesT<CT,ND-1,NDW-1,BoundaryType::PERIODIC,BW>;
  GeometricPropertiesT(const GeometricPropertiesT<CT,ND+1,NDW+1,BoundaryType::PERIODIC,BW> &geom, int dim)
  { 
    // x0.resize(NDIM_W);
    // xmax.resize(NDIM_W);
    // delta.resize(NDIM_W);
    // halfDelta.resize(NDIM_W);
    for (int i=0;i<dim;++i)
      {
	x0[i]=geom.x0[i];
	xmax[i]=geom.xmax[i];
	delta[i]=geom.delta[i];
	halfDelta[i]=geom.halfDelta[i];	
	epsilon[i]=geom.epsilon[i];
      }
    for (int i=dim;i<NDIM_W;++i)
      {
	x0[i]=geom.x0[i+1];
	xmax[i]=geom.xmax[i+1];
	delta[i]=geom.delta[i+1];
	halfDelta[i]=geom.halfDelta[i+1];
	epsilon[i]=geom.epsilon[i+1];
      }
  }

  CT getEpsilon(int i=0)
  {
    return epsilon[i];
  }

  template <class T>
  static int onBoundary(const T * x)
  {
    return 0;
  }

  template <class T, int D = NDIM_W>
  T distance2(const T * /*__restrict*/ xa,const T * /*__restrict*/ xb) const
  {
    typedef typename hlp::IF_< (D<NDIM),	
      hlp::ConstantValue<D>,		
      hlp::ConstantValue<NDIM> >::Result DIM1;  

    return 
      Distance2_P<DIM1::value,T,CT>::compute(xa,xb,&halfDelta[0],&delta[0]) +
      hlp::Distance2<D-DIM1::value,T>::compute(xa+DIM1::value,xb+DIM1::value);
  }

  /*
    // FIXME: need to add template param int DIM = NDIM_W without slowing down !
  template <class T, int DIM = NDIM_W>
  double distance2(const T * __restrict xa,const T * __restrict xb) const
  {  
    double result=0;
    for (int i=0;i<NDIM;i++)
      {
	double dx=xb[i]-xa[i];

	if (dx>=halfDelta[i]) dx-=delta[i];
	if (dx<-halfDelta[i]) dx+=delta[i];
	result += dx*dx;
      }    
    
    for (int i=NDIM;i<NDIM_W;i++)
      {
	double dx=xb[i]-xa[i];
	result += dx*dx;
      }
    
    return result;    
  }
  */

  template <class T, int D = NDIM_W>
  T distance(const T * /*__restrict*/ xa,const T * /*__restrict*/ xb) const
  {
    return sqrt(distance2<T,D>(xa,xb));
  }

  template <class T, int D = NDIM_W>
  T norm2(const T * xa)
  {
    typedef typename hlp::IF_< (D<NDIM),	
      hlp::ConstantValue<D>,		
      hlp::ConstantValue<NDIM> >::Result DIM1;  

    return 
      Norm2_P<DIM1::value,T,CT>::compute(xa,&halfDelta[0],&delta[0]) +
      hlp::Norm2<D-DIM1::value,T>::compute(xa+DIM1::value);    
  }

  template <class T, int D = NDIM_W>
  T norm(const T * xa)
  {
    return sqrt(norm2<T,D>(xa));
  }

  template <class T, int D = NDIM_W>
  T norm_noCheck(const T * xa)
  {
    return GeometricPropertiesNoBoundary::template norm<T,D>(xa);    
  }

  template <class T, int D = NDIM_W>
  T normalize(T * xa)
  {
    T nrm = norm<T,D>(xa);
    if (nrm!=0)
      {
	T nrm_inv = T(1.0)/nrm;
	for (int i=0;i<D;++i) 
	  xa[i] = correctCoordsDiff(*xa,i) * nrm_inv;
      }
    return nrm;
  }

  template <class T, int D = NDIM_W>
  static T normalize_noCheck(T * xa)
  {
    return GeometricPropertiesNoBoundary::template normalize_noCheck<T,D>(xa);
  }

  template <class T>
  static bool inBound(T x, int dim)
  {
    return true;
  }
  /*
  template <class T>
  void checkBoundary(T * __restrict xa) const
  {   
    for (int i=0;i<NDIM;i++)
      {
	if (xa[i]>=xmax[i]) xa[i]-=delta[i];
	else if (xa[i]<x0[i]) xa[i]+=delta[i];
      }
  }

  template <class T>
  T checkBoundary(T c, int which) const
  {   
    if (c>=xmax[which]) c-=delta[which];
    else if (c<x0[which]) c+=delta[which];
    return c;
  }
*/

  template <class T>
  void sanitizeBoundary(T * __restrict xa) const
  {   
    for (int i=0;i<NDIM;i++)
      {
	if (xa[i]>=xmax[i]-epsilon[i]) {xa[i]-=delta[i];if (xa[i]<x0[i]) xa[i]+=delta[i];}
	else if (xa[i]<x0[i]+epsilon[i]) {xa[i]+=delta[i];if (xa[i]>=xmax[i]) xa[i]-=delta[i];}
      }
  }
  
  // The implementation is that way to prevent bugs due to limited precision on one side when  one ofthe boundaries is 0. Keep it that way !.
  template <class T>
  void checkBoundary(T * __restrict xa) const
  {   
    for (int i=0;i<NDIM;i++)
      {
	if (xa[i]>=xmax[i]) {xa[i]-=delta[i];if (xa[i]<x0[i]) xa[i]+=delta[i];}
	else if (xa[i]<x0[i]) {xa[i]+=delta[i];if (xa[i]>=xmax[i]) xa[i]-=delta[i];}
      }
  }

  // The implementation is that way to prevent bugs due to limited precision on one side when  one ofthe boundaries is 0. Keep it that way !.
  template <class T>
  T checkBoundary(T c, int which) const
  {   
    if (c>=xmax[which]) {c-=delta[which];if (c<x0[which]) c+=delta[which];}
    else if (c<x0[which]) {c+=delta[which];if (c>=xmax[which]) c-=delta[which];}
    return c;
  }

  template <class T,int  NVERT=NDIM+1, int D=NDIM_W>
  bool checkCoordsConsistency(T vCoords[NVERT][D],const T * refPoint) const
  {   
    bool modified=false;
    /*
    //T tmp[D];
    // First choose a side (e.g. the side of the first vertex)
    for (int i=0;i<D;++i)
      vCoords[0][i] = checkBoundary<T>(vCoords[0][i],i);
    //tmp[i]=checkBoundary<T>(vCoords[0][i],i);
    
    for (int i=1;i<NVERT;++i)
      for (int j=0;j<D;++j)
	vCoords[i][j] = vCoords[0][j]+correctCoordsDiff<T>(vCoords[i][j]-vCoords[0][j],j);
    */
    for (int i=0;i<NVERT;++i)
      for (int j=0;j<hlp::MinT<D,NDIM>::value;++j)
	vCoords[i][j] = checkCoordConsistency(vCoords[i][j],refPoint[j],j,modified);
    //vCoords[i][j] = refPoint[j]+correctCoordsDiff<T>(vCoords[i][j]-refPoint[j],j);
    return modified;
  }

  /*
  template <class T>
  T checkCoordConsistency(T vCoord, const T refPoint, int dim) const
  {       
    return refPoint+correctCoordsDiff<T>(vCoord-refPoint,dim);  
  }
  */
  template <class T, class T2>
  T checkCoordConsistency(T vCoord, const T2 &refPoint, int which) const
  {
    // typedef typename hlp::MostPreciseFP<T,T2>::Result MPT_;
    // typedef typename hlp::MostPreciseFP<MPT_,CT>::Result MPT;
    
    T result=vCoord;    
    if (which<NDIM)
      {
	//T len=vCoord-refPoint;
	if (vCoord-halfDelta[which]>=refPoint)
	  {
	    // This must NOT be changed !
	    T d=xmax[which];
	    d-=x0[which];
	    result= vCoord-d;//-delta[which];
	  }
	else if (vCoord+halfDelta[which]<refPoint)
	  {
	    // This must NOT be changed !
	    T d=xmax[which];
	    d-=x0[which];
	    result= vCoord+d;//+delta[which];
	  }
      }
    return result;    
  }

  template <class T, class T2>
  T checkCoordConsistency(T vCoord, const T2 &refPoint, int which, bool &modified) const
  {
    // typedef typename hlp::MostPreciseFP<T,T2>::Result MPT_;
    // typedef typename hlp::MostPreciseFP<MPT_,CT>::Result MPT;
    
    T result=vCoord;    
    if (which<NDIM)
      {
	//T len=vCoord-refPoint;
	if (vCoord-halfDelta[which]>=refPoint)
	  {
	    // This must NOT be changed !
	    T d=xmax[which];
	    d-=x0[which];
	    result= vCoord-d;//-delta[which];
	    modified=true;
	  }
	else if (vCoord+halfDelta[which]<refPoint)
	  {
	    // This must NOT be changed !
	    T d=xmax[which];
	    d-=x0[which];
	    result= vCoord+d;//+delta[which];
	    modified=true;
	  }
      }
    return result;    
  }

  template <class T, int  NVERT=NDIM+1, int D=NDIM_W>
  bool coordsAreConsistent(const T vCoords[NVERT][D],const T *refPoint) const
  {
    //return false;
    bool result=true;
    for (int i=0;i<NVERT;++i)
      for (int j=0;j<hlp::MinT<D,NDIM>::value;++j)
	{
	  if ((vCoords[i][j]-halfDelta[j]>=refPoint[j])||
	      (vCoords[i][j]+halfDelta[j]<refPoint[j]))
	    result=false;
	}
    
    return result;
  }

  template <class T>
  T correctCoordsDiff(T len, int which) const
  {   
    T result=len;
    if (which<NDIM)
      {
	if (len>=halfDelta[which])
	  {
	    T d=hlp::numericStaticCast<T>(xmax[which]);
	    d-=x0[which];
	    result=(len-d);//(len-delta[which]);
	  }
	else if (len<-halfDelta[which])
	  {
	    T d=hlp::numericStaticCast<T>(xmax[which]);
	    d-=x0[which];
	    result=(len+d);//delta[which]+len;
	  }
      }
    return result;    
  }  

  /** \brief Get the sign of the difference of two coords (b-a) along dimension dim
   *  or 0 if they are equal. 
   */
  template <class T>
  T coordDiffSign(T a, T b, int dim) const
  {   
    if (b>a)
      {
	if ((dim<NDIM)&&(b>halfDelta[dim]+a)) return -1;
	else return 1;
      }
    else if (b<a)
      {
	if ((dim<NDIM)&&(a>halfDelta[dim]+b)) return 1;
	else return -1;
      }

    return 0;
  } 

  template <class T>
  void displace(T * __restrict xa,const T * __restrict dx) const
  {
    for (int i=0;i<NDIM_W;i++) xa[i]+=dx[i];
    checkBoundary(xa);
  }

  // FIXME: optimize this ?
  template <class T>
  void midPointCoords(const T * __restrict xa, const T * __restrict xb,
		      T * __restrict out) const
  {    
    for (int i=0;i<NDIM;i++)
      {
	T d=fabs(xb[i]-xa[i]);
	if (d>halfDelta[i]) 
	  {	    
	    if (xa[i]<xb[i]) 
	      d=xb[i]+0.5*(delta[i]-d);
	    else
	      d=xa[i]+0.5*(delta[i]-d);

	    if (d<x0[i]) 
	      d+=delta[i];
	    else if (d>=xmax[i]) 
	      d-=delta[i];

	    (*out)=d;
	  }
	else (*out)=(xa[i]+xb[i])*0.5;
	
	++out;
      }
    for (int i=NDIM;i<NDIM_W;i++)
      {
	(*out)=(xa[i]+xb[i])*0.5;
	++out;
      }
  }

  template <class T, class T2, int D=NDIM_W>
  void getVector(const T *a, const T *b, T2* vec) const
  {
    typedef typename hlp::MostPreciseFP<T,T2>::Result MPT;
    static const int midD = hlp::MinT<D,NDIM>::value;
    
    for (int i=0;i<midD;i++) 
      vec[i]=correctCoordsDiff<MPT>(hlp::numericStaticCast<MPT>(b[i])-
				    hlp::numericStaticCast<MPT>(a[i]),
				    i);

    for (int i=midD;i<D;i++) 
      vec[i]=hlp::numericStaticCast<MPT>(b[i])-hlp::numericStaticCast<MPT>(a[i]);

    // hlp::numericStaticCast<T2>(b[i])-
    //   hlp::numericStaticCast<T2>(a[i]),i);
  }

  template <class T, class T2, int D=NDIM_W>
  void getVector(const T *a, T2* vec) const
  {
    typedef typename hlp::MostPreciseFP<T,T2>::Result MPT;
    static const int midD = hlp::MinT<D,NDIM>::value;

    for (int i=0;i<midD;i++) 
      vec[i]=correctCoordsDiff(hlp::numericStaticCast<MPT>(a[i]),i);

    for (int i=midD;i<D;i++) 
      vec[i]=hlp::numericStaticCast<MPT>(a[i]);
  }

  template <class T, class T2, int D, int DW>
  void getBaseVectors(const T * const p[D+1], T2 (&vec)[D][DW]) const
  {
    const T* p0=p[0]; // __restrict ??
    for (int i=0;i<D;i++) 
      {
	const T* p1=p[i+1]; // __restrict ??
	//getVector<T,T2,DW>(p0,p1,vec[i]);
	for (int j=0;j<DW;j++) 
	  {
	    vec[i][j]=correctCoordsDiff(p1[j]-p0[j],j);
	  }
      }  
  }

  template <class T, class T2, int D, int DW>
  void getBaseVectors(const T p[D+1][D], T2 (&vec)[D][DW]) const
  {  
    for (int i=0;i<D;i++) 
      for (int j=0;j<DW;j++) 
	//getVector<T,T2,DW>(p0,p1,vec[i]);
	vec[i][j]=correctCoordsDiff(p[i+1][j]-p[0][j],j);
  }

  template <class T, class T2, int D=NDIM_W>
  T dot(const T *a, const T2 *b)
  {
    T dt=0;
    for (long i=0;i<D;i++) 
      dt+=correctCoordsDiff(a[i],i)*correctCoordsDiff(b[i],i);
    return dt;
  }

  template <class T, class T2, int D=NDIM_W>
  static T dot_noCheck(const T *a, const T2 *b)
  {
    return GeometricPropertiesNoBoundary::template dot_noCheck<T,T2,D>(a,b);
  }

  CT getBBoxSize(int i=0) const
  {
    return delta[i];
  }

protected:
  CT x0[NDIM_W];
  CT xmax[NDIM_W];
  CT delta[NDIM_W];
  CT halfDelta[NDIM_W];
  CT epsilon[NDIM_W];
  // std::vector<CT> x0;
  // std::vector<CT> xmax;
  // std::vector<CT> delta;
  // std::vector<CT> halfDelta;
};


template <typename CT,int ND, int NDW,int BW>
class GeometricPropertiesT<CT,ND,NDW,BoundaryType::NONE,BW>
{
public:
  typedef CT Coord;
  static const int BOUNDARY_TYPE = BoundaryType::NONE;
  static const int WORLD_BOUNDARY_TYPE = BW;
  static const int NDIM=ND;
  static const int NDIM_W=NDW;

  GeometricPropertiesT()
  {}
  
  GeometricPropertiesT(const GeometricPropertiesT<CT,ND-1,NDW,BoundaryType::NONE,BW> &geom, int dim)
  {}

  CT getEpsilon(int i=0)
  {
    return 0;
  }
  
  template <class T>
  static int onBoundary(const T * x)
  {
    return false;
  }

  template <class T, int D = NDIM_W>
  static T distance2(const T * __restrict xa,const T * __restrict xb)
  {
    return hlp::Distance2<D,T>::compute(xa,xb);   
  }

  template <class T, int D = NDIM_W>
  static T distance(const T * __restrict xa,const T * __restrict xb)
  {
    return sqrt(distance2<T,D>(xa,xb));
  }

  template <class T, int D = NDIM_W>
  static T norm2(const T * xa)
  {
    return hlp::Norm2<D,T>::compute(xa);     
  }
  
  template <class T, int D = NDIM_W>
  static T norm(const T * xa)
  {
    return sqrt(norm2<T,D>(xa));
  }

  template <class T, int D = NDIM_W>
  static T norm_noCheck(const T * xa)
  {
    return sqrt(norm2<T,D>(xa));
  }

  template <class T, int D = NDIM_W>
  static T normalize(T * xa)
  {
    T nrm = norm<T,D>(xa);
    if (nrm!=0)
      {
	T nrm_inv = T(1.0)/nrm;
	for (int i=0;i<D;++i) xa[i]*=nrm_inv;
      }
    return nrm;
  }

  template <class T, int D = NDIM_W>
  static T normalize_noCheck(T * xa)
  {
    return normalize<T,D>(xa);
  }

  template <class T>
  static bool inBound(T x, int dim)
  {
    return true;
  }

  template <class T>
  static void sanitizeBoundary(T *xa)
  {
  }

  template <class T>
  static void checkBoundary(T *xa)
  {
  }

  template <class T>
  static T checkBoundary(T c, int which)
  {
    return c;
  }
  
  template <class T,int  NVERT=NDIM+1, int D=NDIM_W>
  static bool checkCoordsConsistency(T vCoords[NVERT][D],const T * refPoint)
  {    
    return false;
  }

  template <class T,class T2>
  static T checkCoordConsistency(T vCoord, const T2 &refPoint, int dim)
  { 
    return vCoord;
  } 

  template <class T, int  NVERT=NDIM+1, int D=NDIM_W>
  bool coordsAreConsistent(const T vCoords[NVERT][D],const T *refPoint) const
  {
    return true;
  }

  template <class T>
  static T correctCoordsDiff(T len, int which)
  {
    return len;
  }
 
  template <class T>
  T coordDiffSign(T a, T b, int dim) const
  { 
    if (b==a) return T(0);
    return (b>a)?T(1):T(-1);
  } 

  template <class T>
  static void displace(T * __restrict xa,const T *__restrict dx)
  {
    for (int i=0;i<NDIM_W;i++) xa[i]+=dx[i];    
  }

  template <class T>
  static void midPointCoords(const T * __restrict xa, const T * __restrict xb, 
			     T * __restrict out)
  {  
    for (int i=0;i<NDIM_W;i++)
      {
	out[i] = 0.5*(xa[i]+xb[i]);
      }   
  }

  template <class T, class T2, int D=NDIM_W>
  static void getVector(const T *xa, const T *xb, T2 *vec) 
  {
    typedef typename hlp::MostPreciseFP<T,T2>::Result MPT;
    for (long i=0;i<D;i++) 
      vec[i]= (hlp::numericStaticCast<MPT>(xb[i])-
	       hlp::numericStaticCast<MPT>(xa[i]));

    // for (long i=0;i<D;i++) 
    //   vec[i]=hlp::numericStaticCast<T2>(xb[i])-hlp::numericStaticCast<T2>(xa[i]);
  }

  template <class T, class T2, int D=NDIM, int DW=NDIM_W>
  static void getBaseVectors(const T * const p[D+1], T2 (&base)[D][DW])
  {    
    const T* p0=p[0]; // __restrict ??
    for (long i=0;i<D;i++) 
      {
	const T* p1=p[i+1]; // __restrict ??
	for (long j=0;j<DW;j++) 
	  {
	    base[i][j]=p1[j]-p0[j];
	    //vec[i][j]=correctCoordsDiff(p1[j]-p0[j],j);
	  }
      }  
  }

  template <class T, class T2, int D=NDIM, int DW=NDIM_W>
  static void getBaseVectors(const T p[D+1][D], T2 (&base)[D][DW])
  {    
    for (long i=0;i<D;i++) 
      for (long j=0;j<DW;j++) 
	base[i][j]=p[i+1][j]-p[0][j];
  }

  template <class T, class T2, int D=NDIM_W>
  static T dot(const T *a, const T2 *b)
  {
    T dt=0;
    for (long i=0;i<D;i++) 
      dt+=a[i]*b[i];
    return dt;
  }

  template <class T, class T2, int D=NDIM_W>
  static T dot_noCheck(const T *a, const T2 *b)
  {
    return dot<T,T2,D>(a,b);
  }

  static CT getBBoxSize(int i=0)
  {
    return 0;
  }

};

template <int ND>
class DummyGeometricPropertiesT : 
  public GeometricPropertiesT<double,ND,ND,BoundaryType::NONE,BoundaryType::NONE>
{};

/** \}*/
#include "../internal/namespace.footer"
#endif

#ifndef __QUADRATIC_SIMPLEX_HXX__
#define __QUADRATIC_SIMPLEX_HXX__

#ifdef NDEBUG
#define EIGEN_NO_DEBUG
#endif

#include "finiteElementTypeEnum.hxx"

#include <Eigen/Dense>
//#include "../algebra/internal/surfaceFitter_implementation.hxx"
#include "internal/quadraticSimplex_implementation.hxx"

/**
 * @file 
 * @brief Defines a class for computing different properties of a quadratic simplex finite 
 * element.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup FiniteElements
 *   \{
 */

/**
 * \class QuadraticSimplexT
 * A class used for computing different properties of a quadratic simplex finite 
 * element defined by its coordinates at vertices and center segments.
 * \tparam ND  Number of dimensions of the simplex
 * \tparam NDW Number of embedding space dimensions
 * \tparam CT  The number type to use for computations
 */
template <int ND, int NDW, typename CT = double>
class QuadraticSimplexT
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef QuadraticSimplexT<ND,NDW,CT> MyType;

  static const int NDIM = ND;
  static const int NDIM_W = NDW;
  static const int NVERT = NDIM+1;
  static const int NSEG = (NDIM*(NDIM+1))/2;
  static const int NCOEFS = NVERT+NSEG;

  static const bool HIERARCHICAL = (NDIM!=NDIM_W);

  static const fETypeE::Type FE_TYPE = fETypeE::simplex;

  typedef CT CompType;

public:

  QuadraticSimplexT()
  {
    reset();
  }

  /** \brief initialize the quadratic simplex from a simplex and tracers for each segment 
   *  midpoint
   */
  template <class Simplex, typename T, typename G>
  QuadraticSimplexT(Simplex *s, const T *segTracers, const G *geometry)
  {
    setPointsCount=0;
    set(s,segTracers,geometry);
  }

  /** \brief initialize the quadratic simplex from a simplex and tracers for each segment 
   *  midpoint
   */
  template <class Simplex, typename T>
  QuadraticSimplexT(Simplex *s, const T *segTracers)
  {
    setPointsCount=0;
    set(s,segTracers);
  }
  
  /** \brief resets the stored points coordinates.
   */
  void reset()
  {
    std::fill_n(pointIsSet,NCOEFS,0);
    setPointsCount=0;
  }

  /** \brief initialize the quadratic simplex from a simplex and tracers for each segment 
   *  midpoint. Segtracers should be given in the order the segments are indexed in Simplex:
   *  If [a,b] represents segment [vertex(a),vertex(b)],
   *  then the order is {[0,1];[0,2];[1,2]} in 2D and
   *  {[0,1];[0,2];[0,3];[1,2];[1,3];[2,3]} in 3D.
   */
  template <class Simplex, typename T, typename G>
  void set(Simplex *simplex, const T *segTracers, const G *geometry)
  {
    if (setPointsCount!=NCOEFS) 
      {
	std::fill_n(pointIsSet,NCOEFS,1);
	setPointsCount=NCOEFS;
      }

    typedef typename Simplex::Coord Coord;
    const Coord *cRef=simplex->getVertex(0)->getCoordsConstPtr();
    // vertices
    for (int j=0;j<Simplex::NVERT;++j)
      {
	const Coord *c=simplex->getVertex(j)->getCoordsConstPtr();
	for (int k=0;k<NDIM_W;++k)
	  coordsMatrix.col(j)[k]=
	    geometry->checkCoordConsistency(c[k],cRef[k],k);
      }

    // and segment tracers    
    const T* c=segTracers;//simplex->segTracers.getConstPointer();
    for (int j=0;j<NSEG;++j,c+=NDIM_W)
      {
	for (int k=0;k<NDIM_W;++k)
	  coordsMatrix.col(j+Simplex::NVERT)[k]=
	    geometry->checkCoordConsistency(c[k],cRef[k],k);
      }

    // In that case, we should adopt the hierarchical representation !
    if (HIERARCHICAL)
      internal::quadSimplex::standardToHierarchical<NDIM>(coordsMatrix);
    
  }

  /** \brief initialize the quadratic simplex from a simplex and tracers for each segment 
   *  midpoint
   */
  template <class Simplex, typename T>
  void set(Simplex *simplex, const T *segTracers)
  {
    if (setPointsCount!=NCOEFS) 
      {
	std::fill_n(pointIsSet,NCOEFS,1);
	setPointsCount=NCOEFS;
      }

    typedef typename Simplex::Coord Coord;
    const Coord *cRef=simplex->getVertex(0)->getCoordsConstPtr();
    // vertices
    for (int j=0;j<Simplex::NVERT;++j)
      {
	const Coord *c=simplex->getVertex(j)->getCoordsConstPtr();
	for (int k=0;k<NDIM_W;++k)
	  coordsMatrix.col(j)[k]=c[k];
      }

    // and segment tracers    
    const T* c=segTracers;//simplex->segTracers.getConstPointer();
    for (int j=0;j<NSEG;++j,c+=NDIM_W)
      {
	for (int k=0;k<NDIM_W;++k)
	  coordsMatrix.col(j+Simplex::NVERT)[k]=c[k];
      }

    // In that case, we should adopt the hierarchical representation !
    if (HIERARCHICAL)
      internal::quadSimplex::standardToHierarchical<NDIM>(coordsMatrix);
  }

  /** \brief Specifies the coordinates of a vertex or segment midpoint. All points
   *  (i.e. 6 in 2D, 10 in 3D) must be specified before any other non static functions 
   *  is used. Note that vertices should be set before segments !
   *  \param start an iterator to the coordinates of the point
   *  \param vi the index of the vertex corresponding to this point if si<0, or the index
   *   of one of the segment extremity if si>=0
   *  \param si the index of the vertex at the other extremity of the segment. If si>=0, 
   *  the tracer corresponds to segment (vi,si), while it corresponds to vertex vi 
   * otherwise.
   *  \tparam Hierarchical if false (default), segments points coordinates are interpreted
   *  as regular coordinates. If Hierarchical is true, then the segments tracers 
   *  coordinates are in the system with origin the middle point of the segment.
   * \warning Vertices (i.e. vi={0,..NDIM} and si<0) must be set before the segments are !
   */
  template <class InputIterator, bool Hierarchical=false>
  void setPointCoord(InputIterator start, int vi, int si=-1)
  {    
    int index = internal::quadSimplex::GetIndexHelper<NDIM,2>::get(vi,si);

    if (!pointIsSet[index]) 
      {
	++setPointsCount;
	pointIsSet[index]=true;
      }
   
    auto it=start;
    for (int i=0;i<NDIM_W;++i,++it)
      coordsMatrix.col(index)[i]=*it;

    if (index>NVERT)
      {
	if (Hierarchical)
	  {
	    if (!HIERARCHICAL)
	      internal::quadSimplex::hierarchicalToStandard<NDIM>(coordsMatrix,index);
	  }
	else
	  {
	    if (HIERARCHICAL)
	      internal::quadSimplex::standardToHierarchical<NDIM>(coordsMatrix,index);
	  }
      }
  }

  /** \brief Converts barycentric coordinates within the simplex to spatial coordinates
   */
  template <class T>
  void barycentricToPosition(T barycentricCoords[NDIM+1], T result[NDIM_W]) const
  {
   
    if (setPointsCount != NCOEFS)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	printf("Some points where not set: %d valid among %d\n",
	       setPointsCount,NCOEFS);
	exit(-1);
      }    

    Eigen::Matrix<CompType,NCOEFS,1> shape;    
    internal::quadSimplex::ShapeFunctionsHelper<NDIM,HIERARCHICAL,2>::
	  evalShapeFunctions(barycentricCoords,shape);    

    Eigen::Matrix<CompType,NDIM_W,1> out = coordsMatrix * shape;
    
    for (int i=0;i<NDIM_W;++i)
      result[i]=out[i];
  }

   /** \brief Converts barycentric coordinates to deviation vector (i.e. difference between
    *   the spatial coordinate on the quadratic simplex and that on the flat linear one).
   */
  template <class T>
  void barycentricToDeviation(T barycentricCoords[NDIM+1], T result[NDIM_W]) const
  {   
    if (setPointsCount != NCOEFS)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	printf("Some points where not set: %d valid among %d\n",
	       setPointsCount,NCOEFS);
	exit(-1);
      }    

    
    T vals[NDIM_W][NCOEFS];
    Eigen::Map< Eigen::Matrix<CompType,NDIM_W,NCOEFS,Eigen::RowMajor> > valsM(&vals[0][0]);

    if (HIERARCHICAL)
      {
	valsM.template rightCols<NSEG>()=coordsMatrix.template rightCols<NSEG>();
      }
    else
      {
	valsM=coordsMatrix;
	internal::quadSimplex::standardToHierarchical(valsM);	
      }

    valsM.template leftCols<NVERT>().array()=0;
    /*
    if (HIERARCHICAL)
      {
	for (int i=0;i<NDIM_W;++i)
	  for (int j=NVERT;j<NVERT+NSEG;++j)
	    vals[i][j] = coordsMatrix.col(j)[i];	
      }
    else
      {
	Mat tmpCM = coordsMatrix;
	internal::quadSimplex::standardToHierarchical(tmpCM);

	for (int i=0;i<NDIM_W;++i)
	  for (int j=NVERT;j<NVERT+NSEG;++j)
	    vals[i][j] = tmpCM.col(j)[i];
      }
    */
    interpolateVec(barycentricCoords,vals,result,NDIM_W);    
  }

   /** \brief Interpolate a function with values val[] given at each tracer to a point
    *  with barycentric coordinates barycentricCoords
    *  \param barycentricCoords barycentric coordinates of the point to interpolate
    *  \param val value of the function at tracer points. The first NDIM+1 points 
    *  correspond to the first NDIM+1 vertices, while the rest correspond to segment 
    *  tracers given in the right order. If [a,b] represent segment [vertex(a),vertex(b)],
    *  then the order is {[0,1];[0,2];[1,2]} in 2D and
    *  {[0,1];[0,2];[0,3];[1,2];[1,3];[2,3]} in 3D.
   */
  template <class TC, class TV> 
  static TV interpolate(TC barycentricCoords[NDIM+1], TV val[NCOEFS])
  {   
    Eigen::Matrix<CompType,NCOEFS,1> shape;
    internal::quadSimplex::ShapeFunctionsHelper<NDIM,false,2>::
      evalShapeFunctions(barycentricCoords,shape);
    
    double result=shape(0)*val[0];
    for (int i=1;i<NCOEFS;++i) result+=shape(i)*val[i];
    
    return result;
  }

  /** \brief Interpolate a vector function with values val[] given at each tracer to a
   *   point with barycentric coordinates barycentricCoords.
    *  \param barycentricCoords barycentric coordinates of the point to interpolate
    *  \param val value of the function at tracer points. The first NDIM+1 points 
    *  correspond to the first NDIM+1 vertices, while the rest correspond to segment 
    *  tracers given in the right order. If [a,b] represent segment [vertex(a),vertex(b)],
    *  then the order is {[0,1];[0,2];[1,2]} in 2D and
    *  {[0,1];[0,2];[0,3];[1,2];[1,3];[2,3]} in 3D.
    *  \param result used to store the interpolated vector.
    *  \param nDims The dimension of the vector to interpolate. 
   */
  template <class TC, class TV, class TR> 
  static void interpolateVec(TC barycentricCoords[NDIM+1], 
			     TV val[][NCOEFS], TR result[], int nDims)
  {
    Eigen::Matrix<CompType,NCOEFS,1> shape;
    internal::quadSimplex::ShapeFunctionsHelper<NDIM,false,2>::
      evalShapeFunctions(barycentricCoords,shape);
    
    for (int i=0;i<nDims;++i)
      {
	result[i]=shape(0)*val[i][0];
	for (int j=1;j<NCOEFS;++j) 
	  result[i]+=shape(j)*val[i][j];
      }   
  }

  /** \brief Evaluates the determinant of the jacobian matrix at each NDIM+1 vertices. The 
   *  returned value already includes the 1/(NDIM!) factor necessary for changing 
   *  variables from barycentric coordinates.
   *  \return the jacobian determinant at vertices times 1/(NDIM!)
   */
  void evalJacobianDetAtVertices(CompType result[NDIM+1]) const
  {
    for (int i=0;i<NVERT;++i)
      {
	CompType bc[NDIM+1]={0};
	bc[i]=1;
	result[i]=evalJacobianDet(bc);
      }
      /*
    static const double factor = 
      static_cast<CT>(hlp::IntPower<4,NDIM+1>::value)/
      static_cast<CT>(hlp::FactorialT<NDIM>::value);

    internal::quadSimplex::TransformHelper<NDIM,NDIM_W,HIERARCHICAL,2>::
      evalJacobianDetAtVertices(coordsMatrix,result)*factor;
      */
    // DEBUG
    /*
    for (int i=0;i< NDIM+1;++i)
      {
	CompType c[NDIM+1]={0};
	c[i]=1;
	printf("result(%d) = (%e == %e) ?\n",i,result[i],evalJacobianDet(c));
      }
    */
    /* */
  }
  
  /** \brief Evaluates the determinant of the jacobian matrix at a point with given
   *  barycentric coordinates. The returned value already includes the 1/(NDIM!) factor 
   *  necessary for changing variables from barycentric coordinates.
   *  \return the jacobian determinant times 1/(NDIM!)
   */
  template <typename T>
  CompType evalJacobianDet(T barycentricCoords[NDIM+1]) const
  {
    Eigen::Matrix<CompType,NDIM_W+1,NDIM+1> J;
    internal::quadSimplex::TransformHelper<NDIM,NDIM_W,HIERARCHICAL,2>::
      evalQuarterJacobian(barycentricCoords,coordsMatrix,J);
    
    CompType result;
    if (NDIM==NDIM_W)
      {
	static const CompType factor = 
	  static_cast<CompType>(hlp::IntPower<4,NDIM+1>::value)/
	  static_cast<CompType>(hlp::FactorialT<NDIM>::value);

	// This is faster ...
	result= J.determinant()*factor;
	/*
	Eigen::Matrix<CompType,NDIM_W,NDIM> K;
	for (int i=0;i<NDIM;++i)
	  K.col(i) = (J.col(i+1).template segment<NDIM_W>(1)-
		      J.col(0).template segment<NDIM_W>(1));
	CompType det2= (K.transpose()*K).determinant();
	*/
	//std::cout << " Result : " << result <<" == "<< sqrt(std::abs(det2))*factor/4;	
      }
    else
      {
	static const CompType factor = 
	  static_cast<CompType>(hlp::IntPower<4,NDIM>::value)/
	  static_cast<CompType>(hlp::FactorialT<NDIM>::value);
	/*
	#pragma omp critical
	{
	  std::cout << "\nJ=" << J <<"\n";
	  std::cout << "C=" << coordsMatrix <<"\n";

	  Eigen::Matrix<CompType,NDIM_W,NDIM> K;
	  for (int i=0;i<NDIM;++i)
	    K.col(i) = (J.col(i+1).segment(1,NDIM_W)-
			J.col(0).segment(1,NDIM_W))*4;

	  std::cout << "K=" << K <<"\n";
	  std::cout << "K2=" << (K.transpose()*K) <<"\n";
	  std::cout << "det =" << sqrt(std::abs((K.transpose()*K).determinant()))/2 <<"\n";
	}
	*/
	Eigen::Matrix<CompType,NDIM_W,NDIM> K;
	for (int i=0;i<NDIM;++i)
	  K.col(i) = (J.col(i+1).template segment<NDIM_W>(1)-
		      J.col(0).template segment<NDIM_W>(1));
	CompType det2= (K.transpose()*K).determinant();
	
	/*
	CompType det2=
	  (J.transpose()*J).determinant();	   
	*/

	if (det2<0)
	  result=-sqrt(-det2)*factor;
	else
	  result=sqrt(det2)*factor;	   
      }
    
    return result;
  }
   
  /*
  void printMat() const
  {
    std::cout << "mat = \n" << coordsMatrix <<"\n";   
  }
  */
private:
  typedef Eigen::Matrix<CompType,NDIM_W,NCOEFS> Mat;

  int pointIsSet[NCOEFS];
  Mat coordsMatrix;
  
  int setPointsCount;
};

/** \}*/
#include "../internal/namespace.footer"
#endif

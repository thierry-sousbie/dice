#ifndef __SURFACE_FITTER_HXX__
#define __SURFACE_FITTER_HXX__

#ifdef NDEBUG
#define EIGEN_NO_DEBUG
#endif

#ifdef HAVE_BOOST
#include "../tools/wrappers/boostMultiprecisionForEigen3.hxx"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/SVD>
#include <Eigen/Dense>
#pragma GCC diagnostic pop

#include "../tools/helpers/helpers.hxx"

#include "internal/surfaceFitter_implementation.hxx"

#include "../internal/namespace.header"

template <int ND, int ND_W, int DEG=2, typename CT = double>
class SurfaceFitterT
{
  template <int _ND, int _ND_W, int _DEG, typename _CT> class FitFunctorT;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  enum Method {fit_method_qr=0, 
	       fit_method_svd=1, 	       
	       fit_method_fast = fit_method_qr,
	       fit_method_precise = fit_method_svd,
	       fit_method_default = fit_method_fast};

  typedef SurfaceFitterT<ND,ND_W,DEG,CT> MyType;
  typedef FitFunctorT<ND,ND_W,DEG,CT> FitFunctor;  
  typedef CT CompType;

  static const int NDIM = ND;
  static const int NDIM_W = ND_W;
  static const int DEGREE = DEG;
  static const int NCOEFS = internal::SurfaceFitterHelperT<ND,DEG>::NCOEFS;

  SurfaceFitterT():
    base(BaseMat::Identity()),
    baseIsCanonical(true),
    allowOverDetermined(true),
    allowUnderDetermined(false)
  {}

  template <class T>
  void addPoint(T point)
  {
   BaseVec tmp;
    for (int i=0;i<NDIM_W;++i)
      tmp[i]=point[i];
    //tmp.t().print("Adding:");
    points.push_back(tmp);
  } 

  template <class InputIterator>
  void addPoints(InputIterator start, InputIterator stop)
  {
    BaseVec tmp;
    int nPoints = std::distance(start,stop);
    
    points.reserve(points.size()+nPoints);   
 
    for (InputIterator it=start; it!=stop; ++it)
      {
	for (int i=0;i<NDIM_W;++i)
	  tmp[i]=(*it)[i];
	points.push_back(tmp);
      }    
  } 

  template <class T>
  void getBase(T baseOut[NDIM_W][NDIM_W]) const
  {
    Eigen::Map< Eigen::Matrix<T,NDIM_W,NDIM_W> > out(baseOut);
    out=base;
  }

  void resetBase()
  {
    if (!baseIsCanonical)
      {
	base=BaseMat::Identity();
	baseIsCanonical=true;
      }
  }

  template <class T>
  void setTangentSpaceBase(const T base[NDIM][NDIM_W])
  {
    T *ptr[NDIM];
    for (int i=0;i<NDIM;++i)
      ptr(i)=&base[i][0];
    setTangentSpaceBase(ptr);
  }

  template <class InputIterator>
  void setTangentSpaceBase(InputIterator start, bool debug=false)
  {
    BaseMat rawBase;
    InputIterator it=start;
    int removed[NDIM_W]={0};
    
    // Set the tangent space first
    for (int i=0;i<NDIM;++i,++it)
      {
	for (int j=0;j<NDIM_W;++j)
	  rawBase(j,i)=(*it)[j];

	rawBase.col(i).normalize();
	// Now we want to find the unit vector most aligned with this vector 
	// so that we do not use it in the normal space base first guess
	
	CompType sp=0;
	int rm_id=0;
	for (int j=0;j<NDIM_W;++j)
	  {
	    if ((!removed[j]) && (fabs(rawBase(j,i))>sp))
	      {
		// Highest scalar product -> remove it
		sp=fabs(rawBase(j,i));//*norm_inv;
		rm_id=j;
	      }
	  }
	removed[rm_id]=1;
      }

    // Now complete our base with the unit vectors that are the least aligned to the 
    // tangent space ones
    int cur=-1;
    for (int i=NDIM;i<NDIM_W;++i)
      {
	while (removed[++cur]);
	rawBase.col(i).fill(0);	
	rawBase(cur,i)=1;
      }

    //rawBase.print("rawBase:");
     
    // Now do some Grahm-schmidt.
    // rawBase = base.R with R upper triangular and 'base' the orthonormalized base
    base.fill(0);
    BaseMat R; 
    R.fill(0);          
    for (int j=0;j<rawBase.cols();j++) 
      {
	BaseVec v =rawBase.col(j);
	if (j>0) 
	  {
	    for(int i=0;i<j;i++) 
	      {
		R(i,j) = (base.col(i).transpose() *  rawBase.col(j))(0);
		v = v - R(i,j) * base.col(i);
	      }
	  }
	R(j,j) = v.norm();
	base.col(j) = v / R(j,j);
      }

    // We do a second Grahm-schmidt as the first one may give poor quality vectors if
    // the input vectors were far from being orthogonal ...
    /* */
    rawBase=base;
    base.fill(0);  
    R.fill(0);    
    for (int j=0;j<rawBase.cols();j++) 
      {
	BaseVec v =rawBase.col(j);
	if (j>0) 
	  {
	    for(int i=0;i<j;i++) 
	      {
		R(i,j) = (base.col(i).transpose() *  rawBase.col(j))(0);
		v = v - R(i,j) * base.col(i);
	      }
	  }
	R(j,j) = v.norm();
	base.col(j) = v / R(j,j);
      }
    /* */

    //if (debug) base.print("Tangent base:");
    baseIsCanonical=false;
  }

  template <class InputIterator>  
  FitFunctor setBaseAndFit(InputIterator start,
			   Method method=fit_method_default, 
			   bool debug=false)
  {
    setTangentSpaceBase(start,debug);
    return fit(method);//debug);
  }

  template <class T>
  FitFunctor setBaseAndFit(const T base[NDIM][NDIM_W],
			   Method method=fit_method_default,  
			   bool debug=false)
  {
    setTangentSpaceBase(base,debug);
    return fit(method);//debug);
  }  
  
  FitFunctor fit(Method method=fit_method_default) const//bool debug=false) const
  {
    static const int NC = NCOEFS;
    int nUsedPoints = points.size(); // NC;   

    if ((!allowUnderDetermined)&&(nUsedPoints<NC)) 
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("System is under-determined !\n");
	exit(-1);
      }
    if ((!allowOverDetermined)&&(nUsedPoints>NC)) 
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("System is over-determined !\n");
	exit(-1);
      }
    
    DMatrix A(nUsedPoints,NC);
    DMatrix f(nUsedPoints,NDIM_W-NDIM);
  
    CompType taylor[NC];
    Eigen::Map< Eigen::Matrix<CompType,NC,1> > taylorVec(taylor);
     
    for (int i=0;i<nUsedPoints;++i)
      {
	// Projected coordinates
	BaseVec pCoords;

	if (!baseIsCanonical)
	  pCoords=base.transpose()*points[i];
	else
	  pCoords=points[i];

	// pCoords.t().print("p:");
	// Taylor coefs for that point
	internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(pCoords,&taylor[0]);
	
	A.row(i) = taylorVec.transpose();
	f.row(i) = pCoords.segment(NDIM,NDIM_W-NDIM).transpose();	
      }
   
    // Diagonal scaling matrix used to improve the conditon number of A
    Eigen::Matrix<CompType,NC,1> S;
    for (int i=0;i<NC;++i)
      {
	CompType tmp=A.col(i).norm();

	if (tmp!=0) 
	  S[i]=1.0/tmp;
	else 
	  S[i]=1;
	
	A.col(i) *= S[i];	  
      }
    
    DVec coefs[NDIM_W-NDIM];   
    doSolve(A,f,coefs,method);

    // Don't forget to unscale the solution
    for (int i=0;i<NDIM_W-NDIM;++i)
      for (int j=0;j<NC;++j)
	coefs[i][j]*=S[j];
    
    /*
    DVec coefs[NDIM_W-NDIM];        
    CompType error=0;
    int m=static_cast<int>(method);
    int mMax=static_cast<int>(fit_method_precise);
    do {
      doSolve(A,f,coefs,static_cast<Method>(m++)); 

      // Don't forget to unscale the solution
      for (int i=0;i<NDIM_W-NDIM;++i)
	for (int j=0;j<NC;++j)
	  coefs[i][j]*=S[j];

      FitFunctor fit(base,coefs);     	
      CompType center[NDIM_W]={0};
      CompType dist=0;
      for (int i=0;i<points.size();++i)
	for (int j=0;j<NDIM_W;++j)
	  center[j]+=points[i][j];

      for (int j=0;j<NDIM_W;++j)
	center[j]/=points.size();

      for (int i=0;i<points.size();++i)
	{
	  CompType tmpDist=0;
	  for (int j=0;j<NDIM_W;++j)
	    tmpDist += (points[i][j] - center[j])*(points[i][j] - center[j]);
	  dist += sqrt(tmpDist);
	}
      dist /= points.size();

      error=0;
      for (int i=0;i<points.size();++i)
	{
	  CompType tmpIn[NDIM_W];
	  CompType tmpOut[NDIM_W];
	  for (int j=0;j<NDIM_W;++j)
	    tmpIn[j]=points[i][j];

	  fit.projectToSurface(tmpIn,tmpOut);

	  CompType err=0;
	  for (int j=0;j<NDIM_W;++j)
	    {
	      if (dist!=0)
		{
		  CompType tmp = fabs((tmpOut[j]-tmpIn[j])/dist);// / tmpIn[i]);
		  if (tmp>err) err=tmp;
		} else printf("NULL DIST\n");
	    }
	  if (err>error) error=err;

	  if (err>4.E-1) 
	    {
	      printf("%e @ %d\n",error,m-1);
	      printf(" Point (%.12e %.12e %.12e %.12e %.12e %.12e)\n",
		     tmpIn[0],tmpIn[1],tmpIn[2],tmpIn[3],tmpIn[4],tmpIn[5]);
	      printf("    -> (%.12e %.12e %.12e %.12e %.12e %.12e)\n",
		     tmpOut[0],tmpOut[1],tmpOut[2],tmpOut[3],tmpOut[4],tmpOut[5]);
	    }
	}
    } while ( (m<=mMax) && (error>4.E-1) );
    */

    /*
    if (debug) 
      {
	A.print("A:");
	f.print("f:");	


	for (int i=0;i<NDIM_W-NDIM;++i)
	  coefs[i].t().print("C:");
      }
    */
      return FitFunctor(base,coefs,baseIsCanonical);
  }

  // Returns the cosine of the minimum angle between any old tangent space base vector and
  // any new normal space base vector
  // => 0 means they are orthogonal so refitting has converged.
  template <class T>
  double refit(FitFunctor &fitF, const T point[][NDIM_W], int nPoints,
	       Method method=fit_method_default,
	       bool evalAngle=false)
  {
    CompType projectedTangentSpaceBase[NDIM][NDIM_W];  
    CompType pt[NDIM_W];

    // Average the normal space basis at each point (each vector should point toward 
    // similar directions)
    std::copy_n(point[0],NDIM_W,pt);
    fitF.computeTangentSpaceBase(pt,projectedTangentSpaceBase);
    for (int i=1;i<nPoints;++i)
      {
	CompType tmpBase[NDIM][NDIM_W];
	std::copy_n(point[i],NDIM_W,pt);
	fitF.computeTangentSpaceBase(pt,tmpBase);
	for (int j=0;j<NDIM;++j)
	  for (int k=0;k<NDIM_W;++k)
	    projectedTangentSpaceBase[j][k]+=tmpBase[j][k];
      }

    if (nPoints>1)
      {
	for (int j=0;j<NDIM;++j)
	  for (int k=0;k<NDIM_W;++k)
	    projectedTangentSpaceBase[j][k]/=nPoints;
      }

    BaseMat oldBase;
    if (evalAngle) {oldBase=base;/*oldBase.print("-> oldBase:");*/}
    
    fitF=setBaseAndFit(projectedTangentSpaceBase,method);
    //base.print("-> newBase:");
    if (evalAngle) 
      {
	//(oldBase.cols(0,NDIM-1).t()*base.cols(NDIM,NDIM_W-1)).print("-> Q:");
	CompType retVal = (oldBase.template block<NDIM_W,NDIM>(0,0).transpose()*
			   base.template block<NDIM_W,NDIM_W-NDIM>(0,NDIM).norm()).norm();
	return hlp::numericStaticCast<double>(retVal);	  
      }

    return 0;
  }
  
  void reset(bool keepBase=false)
  {
    points.clear();
    //vals.clear();
    if (!keepBase) 
      resetBase();
  }

  long getPointsCount() const
  {
    return points.size();
  }
  /*
  static long getFitCoeffsCount(int degree) const
  {
    return ((degree+1)*(degree+2))/2;
  }
  */
  void allowUnderDeterminedSystems(bool val=true)
  {
    allowUnderDetermined=val;
  }

  void allowOverDeterminedSystems(bool val=true)
  {
    allowOverDetermined=val;
  }

private:
  typedef Eigen::Matrix<CompType,Eigen::Dynamic,Eigen::Dynamic> DMatrix;
  typedef Eigen::Matrix<CompType,Eigen::Dynamic,1> DVec;
  typedef Eigen::Matrix<CompType,NDIM_W,NDIM_W> BaseMat;
  typedef Eigen::Matrix<CompType,NDIM_W,1> BaseVec;

  // we need that to have a container of vectors, could using boost align instead, 
  // but we need to require v1.56+ then ...
  // Using unaligned here should not be noticeably slow.
  typedef Eigen::Matrix<CompType,NDIM_W,1,Eigen::DontAlign> UnalignedBaseVec;

  void doSolve(const DMatrix &A, const DMatrix &f, 
	       DVec coefs[NDIM_W-NDIM],
	       Method method) const
  {
    if ( method == fit_method_svd )//qr.rank() < A.cols() )
      {
	// Use the more expensive SVD decomposition only when the matrix does not have 
	// full rank ?
	// glb::console->print<LOG_DEBUG>
	//   ("Using SVD: effective rank is %d < %d (should be >= %d)\n",
	//    rank,A.cols(),(int)A.rows());
	// Eigen::JacobiSVD<DMatrix> svd(qr.matrixQR(),
        // 	Eigen::ComputeThinU|Eigen::ComputeThinV);
	Eigen::JacobiSVD<DMatrix> svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);
	for (int i=0;i<NDIM_W-NDIM;++i)
	  coefs[i]=svd.solve(f.col(i));
      }
    else
      {
	Eigen::ColPivHouseholderQR<DMatrix> qr(A);
	for (int i=0;i<NDIM_W-NDIM;++i)
	  coefs[i]=qr.solve(f.col(i));
      }    
    /*
    CompType res[NDIM_W-NDIM];
    CompType sigma[NDIM_W-NDIM];
    
    for (int i=0;i<NDIM_W-NDIM;++i)
      {
	//CompType ref=f.col(i).lpNorm<Eigen::Infinity>();
	//CompType ref=f.col(i).maxCoeff() - f.col(i).minCoeff()
	DVec tmp = (f.col(i).array() - f.col(i).mean());
	sigma[i] = tmp.template lpNorm<2>();
	if (sigma[i]!=0)
	  {
	    tmp = (f.col(i) - A*coefs[i]);	    
	    res[i]=tmp.template lpNorm<2>() ;
	    //(A*coefs[i]-f.col(i)).lpNorm<Eigen::Infinity>() / ref;
	  }
	else res[i]=0;
      }

    CompType error=res[0]/sigma[0];
    for (int i=1;i<NDIM_W-NDIM;++i)
      {
	CompType tmp=res[i]/sigma[i];
	if (tmp>error) error = tmp;
      }
    if (error>=0.5)
      {
#pragma omp critical 
	printf ("ERRORS(%d): %.12e %.12e %.12e (%.12e %.12e %.12e)\n",(int)method,
		res[0],res[1],res[2],sigma[0],sigma[1],sigma[2]);
      }
    return error;
    */    
  }

  template <int MYDEG>
  void fitRec(hlp::IsTrue,
	      const DMatrix &A,
	      const DMatrix &f, 
	      DVec coefs[NDIM_W-NDIM],
	      Method method,
	      double condLimit=1.E3) const
  {  
    doSolve(A,f,coefs,method);
    /*
    Eigen::JacobiSVD<DMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    for (int i=0;i<NDIM_W-NDIM;++i)    
      coefs(i)=svd(f.col(i));
    */  
  }

  template <int MYDEG>
  void fitRec(hlp::IsFalse,
	      const DMatrix &A,
	      const DMatrix &f, 
	      DVec coefs[NDIM_W-NDIM],
	      Method method,
	      double condLimit=1.E3) const
  {
    static int NC_CUR=internal::SurfaceFitterHelperT<NDIM,MYDEG>::NCOEFS;
    static int NC_LOW=internal::SurfaceFitterHelperT<NDIM,MYDEG-1>::NCOEFS;

    if (false)
      {
	DMatrix A2(points.size(),NC_LOW);	
	DVec coefs2[NDIM_W-NDIM];

	long deg[NC_CUR];
	internal::SurfaceFitterHelperT<NDIM,MYDEG>::computeTaylorDegree(deg);
	// printf("DEG%d: %ld",MYDEG,deg[0]);
	// for (int i=1;i<NC_CUR;++i) printf(" %ld",deg[i]);
	// printf("\n");

	internal::removeHighOrderCols<NDIM,MYDEG>(A,A2,deg);	  

	typedef typename hlp::IsTrueT<(MYDEG<=1)>::Result Status;
	fitRec<MYDEG-1>(Status(),A2,f,coefs2,method);

	for (int i=0;i<NDIM_W-NDIM;++i)	  
	  internal::addHighOrderElements<NDIM,MYDEG>(coefs2[i],coefs[i],deg);
      }
    else 
      {
	doSolve(A,f,coefs,method);
	/*
	Eigen::JacobiSVD<DMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	for (int i=0;i<NDIM_W-NDIM;++i)      
	  coefs(i)=svd(f.col(i));
	*/
      }
  }

  void recursiveFit(const DMatrix&A,
		    const DMatrix&f, 
		    DVec coefs[NDIM_W-NDIM],
		    Method method,
		    double condLimit=1.E3) const
  {
    typedef typename hlp::IsTrueT<(DEGREE<=0)>::Result Status;
    fitRec<DEGREE>(Status(),A,f,coefs,method,condLimit);
  }

  template <int _ND, int _ND_W, int _DEG, typename _CT>
  class FitFunctorT
  {        
  public:
    static const int NDIM = _ND;
    static const int NDIM_W = _ND_W;
    static const int DEGREE = _DEG;    

  private:
    typedef SurfaceFitterT<NDIM,NDIM_W,DEGREE,_CT> SF;
    friend class SurfaceFitterT<NDIM,NDIM_W,DEGREE,_CT>;
    
  public:
    static const int NCOEFS = SF::NCOEFS;
    typedef _CT CompType;
    FitFunctorT() {}

    template <typename T>
    void getValue(const T coords[NDIM], T point[NDIM_W-NDIM]) const
    {
      Eigen::Matrix<CompType,NDIM_W,1> vec;	  
      for (int i=0;i<NDIM;++i) 
	vec[i]=hlp::numericStaticCast<CompType>(coords[i]);
      for (int i=NDIM;i<NDIM_W;++i)
	vec[i]=0;

      // Compute taylor coefficients
      CompType taylor[NCOEFS];
      internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(vec,taylor);
      Eigen::Map< Eigen::Matrix<CompType,NCOEFS,1> > taylorVec(taylor);
      
      // Compute normal space coordinates
      vec.segment(NDIM,NDIM_W-NDIM) = coefs*taylorVec;   
    
      if (!baseIsCanonical)
	vec = base*vec;
	  
      for (int i=0;i<NDIM_W-NDIM;++i) 
	point[i]=hlp::numericStaticCast<T>(vec[i+NDIM]);
    }

    template <typename T>
    void projectToSurface(const T point[NDIM_W], T projectedPoint[NDIM_W], 
			  bool debug=false) const
    {         
      // vecOut stores data at projectedPoint
      //Eigen::Matrix<T,NDIM_W,1> vecOut;//(projectedPoint);

      //for (int i=0;i<NDIM_W;++i) vecOut[i]=point[i];

      // Set vecOut to the coordinates of the point
      //std::copy_n(point,NDIM_W,projectedPoint);
      
      Eigen::Matrix<CompType,NDIM_W,1> vec;
      // This is slower but allows casting to Boost::multiprecision types ...
      for (int i=0;i<NDIM_W;++i) 
	vec[i]=hlp::numericStaticCast<CompType>(point[i]);//vecOut[i]);
      
      //if (debug) vec.t().print("-> oldCoord:");

      // Project point onto the tangent space
      Eigen::Matrix<CompType,NDIM_W,1> pVec;
      if (baseIsCanonical)
	pVec.segment(0,NDIM) = vec.segment(0,NDIM);
      else
	pVec.segment(0,NDIM) = base.template block<NDIM_W,NDIM>(0,0).transpose()*vec;      

      // Compute taylor coefficients
      CompType taylor[NCOEFS];
      internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(pVec,taylor);
      Eigen::Map< Eigen::Matrix<CompType,NCOEFS,1> > taylorVec(taylor);

      // Compute normal space coordinates
      pVec.segment(NDIM,NDIM_W-NDIM) = coefs*taylorVec;

      // And project back
      if (baseIsCanonical)
	vec.noalias() = pVec;
      else
	vec.noalias() = base*pVec;

      for (int i=0;i<NDIM_W;++i) 
	projectedPoint[i]=hlp::numericStaticCast<T>(vec[i]);

      //if (debug) (base.t()*base).print("base_T * base:");
      /*
      if (debug) 
	{
	  //coefs.t().print("coefs:");
	  pVec.t().print("-> pVec:");
	  vec.t().print("-> newCoord:");
	}
      */
    }

    template <typename T>
    void computeTangentSpaceBase(const T point[NDIM_W], 
				 T projectedTangentSpace[NDIM][NDIM_W]) const
    {     
      // Point coordinates as a vector      
      Eigen::Matrix<CompType,NDIM_W,1> vec;
      for (int i=0;i<NDIM_W;++i) 
	vec(i)=hlp::numericStaticCast<CompType>(point[i]);
      
      // Project point onto the tangent space
      Eigen::Matrix<CompType,NDIM,1> pVec;
      if (baseIsCanonical)
	pVec=vec;
      else
	pVec = base.template block<NDIM_W,NDIM>(0,0).transpose()*vec;
      /*
      {
	double t[NCOEFS];
	internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(pVec,&t[0]);
	printf("normal Taylor : %e",t[0]);
	for (int i=1;i<NCOEFS;++i)
	  printf(", %e",t[i]);
	printf("\n");
      }
      */
      // Compute taylor coefficients for the gradient
      CompType taylor[NDIM][NCOEFS];
      internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylorGradient(pVec,taylor);
      /*
      for (int j=0;j<NDIM;++j)
	{
	  printf("grad Taylor[%d]: %e",j,taylor[j][0]);
	  for (int i=1;i<NCOEFS;++i)
	    printf(", %e",taylor[j][i]);
	  printf("\n");
	}
      */
      for (int i=0;i<NDIM;++i)
	{
	  Eigen::Map< Eigen::Matrix<CompType,NCOEFS,1> > taylorVec(&taylor[i][0]);	  
	  Eigen::Matrix<CompType,NDIM_W,1> tmpT;
	  
	  // Compute gradient along the ith base vector
	  tmpT.fill(0);
	  tmpT[i]=1;
	  tmpT.segment(NDIM,NDIM_W-NDIM) = coefs*taylorVec;
	  
	  //Project back and store the result in projectedTangentSpace[i]	 
	  Eigen::Matrix<CompType,NDIM_W,1> tVec;

	  if (baseIsCanonical)
	    tVec=tmpT;
	  else
	    tVec=base*tmpT;

	  for (int j=0;j<NDIM_W;++j) 
	    projectedTangentSpace[i][j]=hlp::numericStaticCast<T>(tVec(j));
	}
    }

    FitFunctorT(const FitFunctorT& other):
      baseIsCanonical(other.baseIsCanonical)
    {
      coefs=other.coefs;
      if (!baseIsCanonical)
	base=other.base;
    }

    FitFunctor &operator=(const FitFunctorT& other)
    {
      if (this==&other)
	return *this;
      
      coefs=other.coefs;
      baseIsCanonical=other.baseIsCanonical;
      if (!baseIsCanonical)
	base=other.base;

      return *this;
    }

  private:
    typedef Eigen::Matrix<CompType,NDIM_W,NDIM_W> BaseMat;        
    typedef Eigen::Matrix<CompType,NDIM_W-NDIM,NCOEFS> CoefMat;    
        
    FitFunctorT(const BaseMat &base_,
		const Eigen::Matrix<CompType,Eigen::Dynamic,1> coefs_[NDIM_W-NDIM],
		bool canonicalBase=false):
      baseIsCanonical(canonicalBase)
    {     
      if (!baseIsCanonical) base=base_;

      for (int i=0;i<NDIM_W-NDIM;++i)
	coefs.row(i)=coefs_[i].transpose();
    }
        
    BaseMat base; // tangent base
    CoefMat coefs;
    bool baseIsCanonical;
  };  

  //std::vector<BaseVec> points;
  // because eigen may segfault when using AVX
  // FIXME: use Boost::align instead (need version 1.56+)
  std::vector<UnalignedBaseVec> points; 
  BaseMat base; 
  bool baseIsCanonical;
  bool allowOverDetermined;
  bool allowUnderDetermined;
};

#include "../internal/namespace.footer"
#endif

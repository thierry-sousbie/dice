#ifndef __SURFACE_FITTER_ARMA_HXX__
#define __SURFACE_FITTER_ARMA_HXX__

//#ifdef HAVE_ARMADILLO
#define ARMA_DONT_USE_CXX11
#ifdef NDEBUG
//#define ARMA_NO_DEBUG
#endif
#include <armadillo>
//#endif

#include "internal/surfaceFitter_implementation.hxx"

#include "../internal/namespace.header"
/*
#ifndef HAVE_ARMADILLO

// Produce an explicit compile time error when trying to instanciate without ARMADILLO
class ARMADILLO_REQUIRED_BUT_NOT_FOUND;
template <int ND, int ND_W>
class SurfaceFitter: ARMADILLO_REQUIRED_BUT_NOT_FOUND {};

#else
*/
template <int ND, int ND_W, int DEG=2, typename CT = double>
class SurfaceFitterT
{
  template <int _ND, int _ND_W, int _DEG, typename _CT> class FitFunctorT;

public:
  typedef SurfaceFitterT<ND,ND_W> MyType;
  typedef FitFunctorT<ND,ND_W,DEG,CT> FitFunctor;  
  typedef CT CompType;

  static const int NDIM = ND;
  static const int NDIM_W = ND_W;
  static const int DEGREE = DEG;
  static const int NCOEFS = internal::SurfaceFitterHelperT<ND,DEG>::NCOEFS;

  SurfaceFitterT():
    baseIsSet(false),
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
    arma::Mat<T> out(baseOut,NDIM_W,NDIM_W,false);
    out=base;
  }

  template <class T>
  void setTangentSpaceBase(const T base[NDIM][NDIM_W])
  {
    T *ptr[NDIM];
    for (int i=0;i<NDIM;++i)
      ptr[i]=&base[i][0];
    setTangentSpaceBase(ptr);
  }

  template <class InputIterator>
  void setTangentSpaceBase(InputIterator start, bool debug=false)
  {
    typename arma::Mat<CompType>::template fixed<NDIM_W,NDIM_W> rawBase;
    InputIterator it=start;
    int removed[NDIM_W]={0};
    
    // Set the tangent space first
    for (int i=0;i<NDIM;++i,++it)
      {
	for (int j=0;j<NDIM_W;++j)
	  rawBase.col(i)[j]=(*it)[j];

	rawBase.col(i) /= arma::norm(rawBase.col(i));
	// Now we want to find the unit vector most aligned with this vector 
	// so that we do not use it in the normal space base
	
	//double norm_inv = 1.0L/arma::norm(rawBase.col(i));
	CompType sp=0;
	int rm_id=0;
	for (int j=0;j<NDIM_W;++j)
	  {
	    if ((!removed[j]) && (fabs(rawBase.col(i)[j])>sp))
	      {
		// Highest scalar product -> remove it
		sp=fabs(rawBase.col(i)[j]);//*norm_inv;
		rm_id=j;
	      }
	  }
	removed[rm_id]=1;
      }

    // And complete our base with the unit vectors that are the least aligned to the 
    // tangent space ones
    int cur=-1;
    for (int i=NDIM;i<NDIM_W;++i)
      {
	while (removed[++cur]);
	rawBase.col(i).fill(0);
	rawBase.col(i)[cur]=1;
      }

    //rawBase.print("rawBase:");
     
    // Now do some Grahm-schmidt.
    // rawBase = base.R with R upper triangular and 'base' the orthonormalized base
    base.fill(0);
    BaseMat R; 
    R.fill(0);          
    for (int j=0;j<rawBase.n_cols;j++) 
      {
	BaseVec v =rawBase.col(j);
	if (j>0) 
	  {
	    for(int i=0;i<j;i++) 
	      {
		R(i,j) = arma::as_scalar(base.col(i).t() *  rawBase.col(j));
		v = v - R(i,j) * base.col(i);
	      }
	  }
	R(j,j) = arma::norm(v,2);
	base.col(j) = v / R(j,j);
      }

    // We do a second Grahm-schmidt as the first one may give poor quality vectors if
    // the input vectors were far from being orthogonal ...
    
    rawBase=base;
    base.fill(0);
    //BaseMat R; 
    R.fill(0);          
    for (int j=0;j<rawBase.n_cols;++j) 
      {
	BaseVec v=rawBase.col(j);
	if (j>0) 
	  {
	    for(int i=0;i<j;++i) 
	      {
		R(i,j) = arma::as_scalar(base.col(i).t() * rawBase.col(j));
		v = v - R(i,j) * base.col(i);
	      }
	  }
	R(j,j) = arma::norm(v,2);
	base.col(j) = v / R(j,j);
      }
    

    if (debug) base.print("Tangent base:");

    baseIsSet=true;
  }

  template <class InputIterator>  
  FitFunctor setBaseAndFit(InputIterator start, double &cond, double limitCondNumber=-1.0, bool debug=false)
  {
    setTangentSpaceBase(start,debug);
    return fit(cond,limitCondNumber);//debug);
  }

  template <class T>
  FitFunctor setBaseAndFit(const T base[NDIM][NDIM_W], double &cond, double limitCondNumber=-1.0, bool debug=false)
  {
    setTangentSpaceBase(base,debug);
    return fit(cond,limitCondNumber);//debug);
  }  
  
  FitFunctor fit(double &cond, double limitCondNumber=1.E-3) const//bool debug=false) const
  {
    if (!baseIsSet)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("Tangent space base needs to be set before calling fit() !\n");
	exit(-1);
      }
    if ((!allowUnderDetermined)&&(points.size()<NCOEFS)) 
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>
	  ("System is under-determined !\n");
	exit(-1);
      }
    if ((!allowOverDetermined)&&(points.size()>NCOEFS)) 
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("System is over-determined !\n");
	exit(-1);
      }

    arma::Mat<CompType> A(points.size(),NCOEFS);
    arma::Mat<CompType> f(points.size(),NDIM_W-NDIM);        

    CompType taylor[NCOEFS];
    arma::Col<CompType> taylorVec(&taylor[0],NCOEFS,false);
    
    for (int i=0;i<points.size();++i)
      {
	// Projected coordinates
	BaseVec pCoords = base.t()*points[i];
	// pCoords.t().print("p:");
	// Taylor coefs for that point
	internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(pCoords,&taylor[0]);
	
	A.row(i) = taylorVec.t();
	f.row(i) = pCoords.subvec(NDIM,NDIM_W-1).t();	
      }
   
    // Diagonal scaling matrix used to improve the conditon number of A
    typename arma::Col<CompType>::template fixed<NCOEFS> S;    
    for (int i=0;i<NCOEFS;++i)
      {
	CompType tmp=arma::norm(A.col(i),2);

	if (tmp!=0) 
	  S[i]=1.0/tmp;
	else 
	  S[i]=1;
	
	A.col(i) *= S[i];	  
      }
     
    arma::Col<CompType> coefs[NDIM_W-NDIM];    
    /*
    for (int i=0;i<NDIM_W-NDIM;++i)      
      arma::solve(coefs[i],A,f.col(i));
    */
    if (limitCondNumber>0)
      recursiveFit(A,f,coefs,limitCondNumber);
    else
      {
	cond = arma::cond(A);
	for (int i=0;i<NDIM_W-NDIM;++i)      
	  arma::solve(coefs[i],A,f.col(i));
      }
    
    
    // Don't forget to unscale the solution
    for (int i=0;i<NDIM_W-NDIM;++i)
      for (int j=0;j<NCOEFS;++j)
	coefs[i][j]*=S[j];
    /*
    if (arma::cond(A)>1.E3)
      {
	printf(" Condition number: %e\n",arma::cond(A));
      }
    */
    /*
    if (debug) 
      {
	A.print("A:");
	f.print("f:");	

	printf(" Condition number: %e\n",arma::cond(A));

	for (int i=0;i<NDIM_W-NDIM;++i)
	  coefs[i].t().print("C:");
      }
    */
    return FitFunctor(base,coefs);
  }

  // Returns the cosine of the minimum angle between any old tangent space base vector and
  // any new normal space base vector
  // => 0 means they are orthogonal so refitting has converged.
  template <class T>
  double refit(FitFunctor &fitF, const T point[][NDIM_W], int nPoints, 
	       double &cond, double limitCondNumber=-1.0, bool evalAngle=false)
  {
    T projectedTangentSpaceBase[NDIM][NDIM_W];  
     
    // Average the normal space basis at each point
    fitF.computeTangentSpaceBase(&point[0][0],projectedTangentSpaceBase);
    for (int i=1;i<nPoints;++i)
      {
	T tmpBase[NDIM][NDIM_W];
	fitF.computeTangentSpaceBase(&point[i][0],tmpBase);
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
    
    fitF=setBaseAndFit(projectedTangentSpaceBase,cond,limitCondNumber);
    //base.print("-> newBase:");
    if (evalAngle) 
      {
	//(oldBase.cols(0,NDIM-1).t()*base.cols(NDIM,NDIM_W-1)).print("-> Q:");
	return arma::norm(oldBase.cols(0,NDIM-1).t()*base.cols(NDIM,NDIM_W-1));
      }

    return 0;
  }
  
  void reset(bool keepBase=false)
  {
    points.clear();
    //vals.clear();
    if (!keepBase) baseIsSet=false;
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
  /*
  template <class T1, class T2>
  void toTangentSpace(T1 *coords, T2 *newCoords) const
  {
       arma::vec::fixed<NDIM_W> vec;
    for (int j=0;j<NDIM_W;++j)
      vec[i]=coords[i];
      
    // Project point
    arma::vec::fixed<NDIM> pc=base.t()*vec;
    for (int j=0;j<NDIM_W;++j) 
      newCoords[j]=pc[j];    
  }
  */

  template <int MYDEG>
  void fitRec(hlp::IsTrue,
	      const arma::Mat<CompType> &A,
	      const arma::Mat<CompType> &f, 
	      arma::Col<CompType> coefs[NDIM_W-NDIM],
	      double condLimit=1.E3) const
  {  
    for (int i=0;i<NDIM_W-NDIM;++i)      
      arma::solve(coefs[i],A,f.col(i));
  }

  template <int MYDEG>
  void fitRec(hlp::IsFalse,
	      const arma::Mat<CompType> &A,
	      const arma::Mat<CompType> &f, 
	      arma::Col<CompType> coefs[NDIM_W-NDIM],
	      double condLimit=1.E3) const
  {
    static int NC_CUR=internal::SurfaceFitterHelperT<NDIM,MYDEG>::NCOEFS;
    static int NC_LOW=internal::SurfaceFitterHelperT<NDIM,MYDEG-1>::NCOEFS;

    if (arma::cond(A)>condLimit)
      {
	arma::Mat<CompType> A2(points.size(),NC_LOW);	
	arma::Col<CompType> coefs2[NDIM_W-NDIM];

	long deg[NC_CUR];
	internal::SurfaceFitterHelperT<NDIM,MYDEG>::computeTaylorDegree(deg);
	// printf("DEG%d: %ld",MYDEG,deg[0]);
	// for (int i=1;i<NC_CUR;++i) printf(" %ld",deg[i]);
	// printf("\n");

	internal::removeHighOrderCols<NDIM,MYDEG>(A,A2,deg);	  

	typedef typename hlp::IsTrueT<(MYDEG<=1)>::Result Status;
	fitRec<MYDEG-1>(Status(),A2,f,coefs2);

	for (int i=0;i<NDIM_W-NDIM;++i)	  
	  internal::addHighOrderElements<NDIM,MYDEG>(coefs2[i],coefs[i],deg);
      }
    else 
      {
	for (int i=0;i<NDIM_W-NDIM;++i)      
	  arma::solve(coefs[i],A,f.col(i));
      }
  }

  void recursiveFit(const arma::Mat<CompType> &A,
		    const arma::Mat<CompType> &f, 
		    arma::Col<CompType> coefs[NDIM_W-NDIM],
		    double condLimit=1.E3) const
  {
    typedef typename hlp::IsTrueT<(DEGREE<=0)>::Result Status;
    fitRec<DEGREE>(Status(),A,f,coefs,condLimit);
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
    void projectToSurface(const T point[NDIM_W], T projectedPoint[NDIM_W], 
			  bool debug=false) const
    {         
      // vec stores data at projectedPoint
      arma::Col<T> vec(projectedPoint,NDIM_W,false);

      // Set vec to the coordinats of the point
      std::copy_n(point,NDIM_W,projectedPoint);
      if (debug) vec.t().print("-> oldCoord:");

      // Project point onto the tangent space
      typename arma::Col<CompType>::template fixed<NDIM_W> pVec;
      pVec.subvec(0,NDIM-1) = base.cols(0,NDIM-1).t()*vec;      

      // Compute taylor coefficients
      CompType taylor[NCOEFS];
      internal::SurfaceFitterHelperT<NDIM,DEGREE>::computeTaylor(pVec,taylor);
      arma::Col<CompType> taylorVec(taylor,NCOEFS,false);

      // Compute normal space coordinates
      pVec.subvec(NDIM,NDIM_W-1) = coefs*taylorVec;

      // And project back
      vec = base*pVec;

      //if (debug) (base.t()*base).print("base_T * base:");
      if (debug) 
	{
	  //coefs.t().print("coefs:");
	  pVec.t().print("-> pVec:");
	  vec.t().print("-> newCoord:");
	}
    }

    template <typename T>
    void computeTangentSpaceBase(const T point[NDIM_W], 
				 T projectedTangentSpace[NDIM][NDIM_W]) const
    {     
      // Point coordinates as a vector
      const arma::Col<T> vec(const_cast<T*>(&point[0]),NDIM_W,false);
      
      // Project point onto the tangent space
      typename arma::Col<CompType>::template fixed<NDIM> pVec;
      pVec = base.cols(0,NDIM-1).t()*vec;
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
	  arma::Col<CompType> taylorVec(&taylor[i][0],NCOEFS,false);	  
	  typename arma::Col<CompType>::template fixed<NDIM_W> tmpT;
	  
	  // Compute gradient along the ith base vector
	  tmpT.fill(0);
	  tmpT[i]=1;
	  tmpT.subvec(NDIM,NDIM_W-1) = coefs*taylorVec;
	  
	  //Project back and store the result in projectedTangentSpace[i]
	  arma::Col<T> tVec(&projectedTangentSpace[i][0],NDIM_W,false);
	  tVec=base*tmpT;
	}
    }

  private:
    
    FitFunctorT(const typename arma::Mat<CompType>::template fixed<NDIM_W,NDIM_W> &base_,
		const typename arma::Col<CompType> coefs_[NDIM_W-NDIM])
    {     
      base=base_;
      for (int i=0;i<NDIM_W-NDIM;++i)
	coefs.row(i)=coefs_[i].t();
    }
        
    typename arma::Mat<CompType>::template fixed<NDIM_W,NDIM_W> base; // tangent base
    typename arma::Mat<CompType>::template fixed<NDIM_W-NDIM,NCOEFS> coefs;
  };

  typedef typename arma::Mat<CompType>::template fixed<NDIM_W,NDIM_W> BaseMat;
  typedef typename arma::Col<CompType>::template fixed<NDIM_W> BaseVec;

  std::vector<BaseVec> points;  
  BaseMat base; 
  bool baseIsSet;
  bool allowOverDetermined;
  bool allowUnderDetermined;
};
//#endif // HAVE_ARMADILLO

#include "../internal/namespace.footer"
#endif

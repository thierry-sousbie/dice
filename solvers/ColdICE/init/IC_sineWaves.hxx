#ifndef __IC_SINE_WAVES_HXX__
#define __IC_SINE_WAVES_HXX__

#include <math.h>

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include <dice/algebra/inverseMatrix.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_sineWaves : 
  public InitialConditionsInterfaceT<M>
{
private:
  static const int NDIM=M::NDIM;
  static const int NDIM_W=M::NDIM_W;

public:
  typedef InitialConditionsInterfaceT<M> Interface;
  typedef typename Interface::ImplicitTesselation ImplicitTesselation;
  typedef typename Interface::Params Params;
  typedef typename Interface::Vertex Vertex;

  static std::string name() {return "sineWaves";}
  std::string getName() const {return name();}

  template <class PM, class R>
  IC_sineWaves(PM &paramsManager, R *reader, 
	       dice::MpiCommunication *com=0) //nullptr
  {
    parse(paramsManager,reader);
  }

  void initialize(const double *x0, const double *delta, double t, 
		  const Units &units, int pass)
  {
    std::copy_n(x0,NDIM,x0_);
    std::copy_n(delta,NDIM,delta_);

    vFact[0]=0;
    vFact[1]=0;

    xFact[0]=1.0;
    xFact[1]=0.0;

    if (cosmo)
      {
	double aStart=units.cosmology.a_of_tau(t,0);	

	vFact[0]=
	  units.cosmology.h_over_h0(aStart)*
	  units.cosmology.f1(aStart)*
	  units.H*units.length*aStart/units.velocity; 	

	vFact[1]=
	  units.cosmology.h_over_h0(aStart)*
	  units.cosmology.f2(aStart)*
	  units.H*units.length*aStart/units.velocity; 

	xFact[0]=units.cosmology.d1_plus(aStart);
	xFact[1]=units.cosmology.d2_plus(aStart);
      }    
  }

  void displace(double *coords, double &density, Vertex *v)
  {
    double x[NDIM];
    double u[NDIM];
    double uco[NDIM];
    double dx[NDIM];
    double dv[NDIM];
    
    for (int i=0;i<NDIM;++i)
      x[i]=(coords[i]-x0_[i])/delta_[i];

    // Change basis
    for (int i=0;i<NDIM;++i)
      {
	uco[i]=0;
	for (int j=0;j<NDIM;++j)
	  uco[i]+=kInv[j][i]*x[j];
	u[i]=uco[i];
      }
    /*
    // Multiply by metric tensor to get contravariant coords (not working)
    for (int i=0;i<NDIM;++i)
      {
	u[i]=0;
	for (int j=0;j<NDIM;++j)
	  u[i]+=g[i][j]*uco[j];
	//u[i]=uco[i];
      }
    */
    // Scale and peridoize
    for (int i=0;i<NDIM;++i)
      {
	u[i]*=scaleFactor[i]; //if (u[i]>1) printf("u[i]=%lg\n",u[i]);
	if (u[i]<0)  u[i]+=int(-u[i])+1;
	if (u[i]>=1) u[i]-=int(u[i]);
      }
    /*
    for (int i=0;i<NDIM;++i) 
      {	     
	double dx = xFact[0]*dx1(x,i) + xFact[1]*dx2(x,i);
	double dv = xFact[0]*vFact[0]*dx1(x,i) + xFact[1]*vFact[1]*dx2(x,i);

	coords[i]+=dx*delta_[i];
	coords[i+NDIM]+=dv*delta_[i];
      }
    */
    for (int i=0;i<NDIM;++i) 
      {
	dx[i]= xFact[0]*dx1(u,i) + xFact[1]*dx2(u,i);
	dv[i]= xFact[0]*vFact[0]*dx1(u,i) + xFact[1]*vFact[1]*dx2(u,i);
      }

    for (int j=0;j<NDIM;++j)
      coords[j]+=dx[j]*delta_[j];
    for (int j=0;j<NDIM;++j)
      coords[j+NDIM]+=dv[j]*delta_[j];

    // Staggered grid !
    //for (int j=0;j<NDIM;++j)
    //  coords[j]+=1.0/256.0*0.123;
    //coords[j]+=0.5*delta_[j]/sgParams.resolution[j];
    
  }

  Params getParams()
  {
    Params p;
    std::copy_n(sgParams.x0,NDIM,&p.x0.front());
    std::copy_n(sgParams.delta,NDIM,&p.delta.front());

    p.useCosmo=cosmo;
    p.cosmoParams.set(1.0,0,0,1.0);
    p.mass=1;
    p.aStart=1.e-2;

    p.unitLength=1.0;
    p.unitVelocity=1.0;
    p.unitMass=1.0;

    p.defaultUnitLength="-";
    p.defaultUnitVelocity="-";
    p.defaultUnitMass="-";

    p.G=1;
    p.H=1;
    
    return p;
  }

  ImplicitTesselation *createImplicitTesselation()
  { 
    SimplicialGrid *sg = new SimplicialGrid(sgParams);
    return reinterpret_cast<ImplicitTesselation*>(sg);
  }

  bool useLagrangianMass() const {return true;}
private:
  typedef typename ImplicitTesselation::Cell Cell;
  typedef dice::SimplicialGridT<NDIM,NDIM_W,Cell,M::IS_PERIODIC> SimplicialGrid;
  typedef typename SimplicialGrid::Params SimplicialGridParams;
  
  SimplicialGridParams sgParams;  

  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {
    cosmo = 1;
    cosmo = paramsManager.
      get("cosmo",Interface::parserCategory(),cosmo,
	  reader,PM::FILE_FIRST,
	  "Whether to use cosmology");

    if (cosmo)
      {
	std::fill_n(sgParams.x0,NDIM,-0.5);
	std::fill_n(sgParams.delta,NDIM,1.0);
      }
    else
      {
	std::fill_n(sgParams.x0,NDIM,-0.5);
	std::fill_n(sgParams.delta,NDIM,1.0);
      }
    sgParams.parse(reader,paramsManager,Interface::parserCategory());

    amplitude[0] = (cosmo)?2:0.1;
    amplitude[0] = paramsManager.
      get("amplitude",Interface::parserCategory(),amplitude[0],0,
	  reader,PM::FILE_FIRST,
	  "Amplitude of the sine waves");

    for (int i=1;i<NDIM;++i)
      {
	amplitude[i] = amplitude[i-1];
	amplitude[i] = paramsManager.
	  get("amplitude",Interface::parserCategory(),amplitude[i],i,
	      reader,PM::FILE_FIRST,
	      "Amplitude of the sine waves");
      }

    for (int k=0;k<NDIM;++k)
      {
	int orMin=10000;
	buildNextIVec(orientation,k);
	for (int i=0;i<NDIM;++i)
	  {
	    orientation[k][i] = paramsManager.
	      get("orientation",Interface::parserCategory(),orientation[k][i],k*NDIM+i,
		  reader,PM::FILE_FIRST,
		  "A vector of -INTEGER- components indicating the direction orthogonal to the plane wave W.R.T each axis (integer components ensure correct reconnection of the waves with periodic boundary conditions). If you give a single vector, an arbitrary orthonormal basis will be built.");
	    if ((orMin>orientation[k][i])&&(orientation[k][i]>0)) orMin=orientation[k][i];
	  }
	
	// Make sure the components have no common factors
	// Could use sqrt(orMin)+1 and test orMin separately, but who cares ...
	for (int i=2;i<=orMin;++i)
	  {
	    bool isFactor=true;
	    for (int j=0;j<NDIM;++j)
	      if (int(orientation[k][j]/i)*i != orientation[k][j])
		isFactor=false;
	    
	    if (isFactor)
	      {
		for (int j=0;j<NDIM;++j) orientation[k][j]/=i;
		orMin/=i;
		if (i>1) i--; //account for factors multiplicity 
	      }
	  }
      }
    
    for (int j=0;j<NDIM;++j)
      {
	scaleFactor[j]=0;
	for (int i=0;i<NDIM;++i)
	  {
	    kVec[j][i]=orientation[j][i];
	    scaleFactor[j] += kVec[j][i]*kVec[j][i];
	  }
	scaleFactor[j]=sqrt(scaleFactor[j]);

	for (int i=0;i<NDIM;++i)
	  kVec[j][i]/=scaleFactor[j];
      }
    /*
    if (NDIM==2)
      {
	double tmp=kVec[1][0];
	kVec[1][0]=-kVec[1][1];
	kVec[1][1]=tmp;
      }
    else 
      {
	for (int j=1;j<NDIM;++j)
	  for (int i=0;i<NDIM;++i)
	    kVec[j][i]=(i==j);

	// Graam schmidt
	for (int k=1;k<NDIM;k++)
	  {
	    for (int i=0;i<k;++i)
	      for (int j=0;j<NDIM;j++)
		kVec[k][j]-=kVec[i][j]*kVec[k][j] * kVec[i][j];

	    double dot=0;
	    for (int j=0;j<NDIM;++j) dot+=kVec[k][j];
	    dot=1.0/sqrt(dot);
	    for (int j=0;j<NDIM;++j) kVec[k][j]*=dot;
	  }
      }
    */
    dice::InverseMatrixT<NDIM,double,double>::compute(kVec,kInv);

    for (int i=0;i<NDIM;++i)
      for (int j=0;j<NDIM;++j)
	{
	  g[i][j]=0;
	  for (int k=0;k<NDIM;++k)
	    g[i][j]+=kVec[i][k]*kVec[j][k];
	}
    //for (int i=0;i<NDIM;++i)
    //  std::copy_n(kVec[i],NDIM,kVec[i]+NDIM);

    // for (int k=0;k<NDIM;++k)
    //   {
    // 	scaleFactor[k]=kNorm[k];
    // 	 for (int i=0;i<NDIM;++i)
    // 	   scaleFactor[k]+=orientation[k][i]*orientation[k][i];
    // 	 scaleFactor[k] = sqrt(scaleFactor[k]);
    //   }
    
    dice::glb::console->printFlush<dice::LOG_INFO>
      ("Initial conditions scale factor is f = (%lg, %lg).\n",scaleFactor[0],scaleFactor[1]);
    dice::glb::console->printFlush<dice::LOG_INFO>
      ("Initial conditions scale factor is k = (%lg, %lg; %lg %lg).\n",kVec[0][0],kVec[0][1],kVec[1][0],kVec[1][1]);
    /*
    rescaleBox = 0;
    rescaleBox = paramsManager.
      get("rescaleBox",Interface::parserCategory(),rescaleBox,
	  reader,PM::FILE_FIRST,
	  "Whether to rescale the box in order to conserve distances between waves when 'orientation' vector is not unitary.");

    if (rescaleBox)
      {
	for (int i=0;i<NDIM;++i)
	  {
	    sgParams.x0[i]*=scaleFactor;
	    sgParams.delta[i]*=scaleFactor;
	  }
      }

    std::cout << "Orientation: " ;
    for (int i=0;i<NDIM;++i) std::cout << orientation[i] << " ";
    std::cout << std::endl;
    std::cout << "Orientation (k): " ;
    for (int i=0;i<NDIM_W;++i) std::cout << kVec[0][i] << " ";
    std::cout << std::endl;
    std::cout << "Scale factor : " << scaleFactor << std::endl;;
    //exit(0);
    */
  }

  double dx1(double x[NDIM], int index)
  {
    static const double twoPi = 2.0*M_PI;    
    return -sin(2.0*M_PI*(x[index]-0.5))*amplitude[index]/twoPi;
  }

  double dx2(double x[NDIM], int index)
  {
    // Removed second order !
    return 0;

    static const double twoPi = 2.0*M_PI;
    
    if (!cosmo) return 0;

    double result = 1.0/(4.0*twoPi)*sin(2.0L*M_PI*(x[index]-0.5))*amplitude[index];

    double tmp=0;
    for (int i=0;i<NDIM;++i)
      if (i!=index) tmp+=cos(2.0L*M_PI*(x[i]-0.5))*amplitude[i];
    
    return result*tmp;
  }

  template <typename I>
  void buildNextIVec(I iVec[NDIM][NDIM], int which)
  {
    if (which>=NDIM) return;
    
    if (which==0)
      {
	for (int i=0;i<NDIM;++i)
	  iVec[0][i]=(i==0);
      }
    else if (which==1)
      {
	if (NDIM==2)
	  {
	    iVec[1][0]=-iVec[0][1];
	    iVec[1][1]= iVec[0][0];
	  }
	else
	  {
	    std::vector<int> id;
	    for (int i=0;i<NDIM;++i)
	      {
		iVec[1][i]=0;
		if (iVec[0][i]!=0) id.push_back(i);
	      }
	    if (id.size()==1) iVec[1][(id[0]+1)%NDIM]=iVec[0][id[0]];
	    else
	      {
		iVec[1][id[0]]=-iVec[0][id[1]];
		iVec[1][id[1]]= iVec[0][id[0]];
	      }
	  }
      }
    else if (which==2)
      {
	iVec[2][0]=iVec[0][1]*iVec[1][2]-iVec[0][2]*iVec[1][1];
	iVec[2][1]=iVec[0][2]*iVec[1][0]-iVec[0][0]*iVec[1][2];
	iVec[2][2]=iVec[0][0]*iVec[1][1]-iVec[0][1]*iVec[1][0];
      }
  }

  int cosmo;
  double amplitude[NDIM];
  int orientation[NDIM][NDIM]; // integer valued vector corresponding to X vector
  //int rescaleBox;
  double scaleFactor[NDIM];
  double x0_[NDIM];
  double delta_[NDIM]; 
  double vFact[2];
  double xFact[2];
  double kVec[NDIM][NDIM]; // The rotated unit vectors basis
  double kInv[NDIM][NDIM];
  double g[NDIM][NDIM]; // metric tensor
};

#endif

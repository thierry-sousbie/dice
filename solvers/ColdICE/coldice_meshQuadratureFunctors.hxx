#ifndef __COLDICE_MESH_QUADRATURE_FUNCTORS_HXX__
#define __COLDICE_MESH_QUADRATURE_FUNCTORS_HXX__

#include <stdlib.h>
#include <functional>
#include <dice/finiteElements/finiteElementTypeEnum.hxx>
#include <dice/finiteElements/quadraticSimplex.hxx>
#include <dice/tools/wrappers/standardDrand48ReentrantWrapper.hxx>

// Functors for quadrature a piecewise linear and quadratic mesh

// Kinetic energy
template <class M, int ELEMENT_ORDER=1>
class KineticEnergyMQT
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = ELEMENT_ORDER;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = 2;
  
  KineticEnergyMQT()
  {}
 
  void set(const Mesh *mesh, const Simplex *s)
  {
    double halfVolume = mesh->computeProjectedVolume(s)*0.5;
    
    for (int i=0;i<NVERT;++i)
      {
	Vertex *v=s->getVertex(i);
	for (int j=0;j<NDIM_W-NDIM;++j)
	  vc[j][i]=v->getCoord(NDIM+j);
	m[i]=halfVolume*v->projectedDensity.getValue();
      }
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    double vx=a*vc[0][0]+b*vc[0][1]+c*vc[0][2];
    double vy=a*vc[1][0]+b*vc[1][1]+c*vc[1][2];
    double rho=a*m[0]+b*m[1]+c*m[2];
    return rho*(vx*vx+vy*vy);
  } 

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    double vx=a*vc[0][0]+b*vc[0][1]+c*vc[0][2]+d*vc[0][3];
    double vy=a*vc[1][0]+b*vc[1][1]+c*vc[1][2]+d*vc[1][3];
    double vz=a*vc[2][0]+b*vc[2][1]+c*vc[2][2]+d*vc[2][3];
    double rho=a*m[0]+b*m[1]+c*m[2]+d*m[3];
    return rho*(vx*vx+vy*vy+vz*vz);
  }

  std::string getName() const {return std::string("kineticEnergy");}
private:
  Coord vc[NDIM_W-NDIM][NVERT];
  double m[NVERT];
};

// Kinetic energy over a quadratic mesh
template <class M>
class KineticEnergyMQT<M,2>
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;
  static const int NSEG = Simplex::NSEG;

  static const int FE_ORDER = 2;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = 4;

  typedef dice::QuadraticSimplexT<NDIM,NDIM> QS;

  KineticEnergyMQT()
  {}
  
  void set(const Mesh *mesh, const Simplex *s)
  {    
    for (int i=0;i<NVERT;++i)
      {
	Vertex *vertex=s->getVertex(i);
	for (int j=0;j<NDIM_W-NDIM;++j)
	  v[j][i]=vertex->getCoord(NDIM+j);	
      }

    const Coord *coord = s->segTracers.getConstPointer();
    for (int i=0;i<NSEG;++i)
      {
	const Coord *c=&coord[i*NDIM_W+NDIM];
	for (int j=0;j<NDIM_W-NDIM;++j)
	  v[j][i+NVERT]=c[j];
      }
    halfMass=s->mass.getValue()*0.5;
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    double bc[3]={a,b,c};
    double vInterp[NDIM_W-NDIM];
    QS::template interpolateVec(bc,v,vInterp,NDIM_W-NDIM);
    
    return halfMass*(vInterp[0]*vInterp[0]+vInterp[1]*vInterp[1]);
  } 

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    double bc[4]={a,b,c,d};
    double vInterp[NDIM_W-NDIM];
    QS::interpolateVec(bc,v,vInterp,NDIM_W-NDIM);
    
    return halfMass*(vInterp[0]*vInterp[0]+vInterp[1]*vInterp[1]+vInterp[2]*vInterp[2]);
  }

  std::string getName() const {return std::string("kineticEnergy_Q");}
private:
  Coord v[NDIM_W-NDIM][NVERT+NSEG];
  double halfMass;
};

// Volume of the mesh
template <class M, int ELEMENT_ORDER=1>
class VolumeMQT
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = ELEMENT_ORDER;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = 1;
  
  typedef dice::QuadraticSimplexT<NDIM,NDIM_W> QS;
  
  
  VolumeMQT()
  {}
  
  void set(const Mesh *mesh, Simplex *s)
  {  
    volume = mesh->computeVolume(s);        
  }

  // 2D
  double operator()(double a, double b, double c) const
  {   
    return volume;
  }

  // 3D
  double operator()(double a, double b, double c, double d) const
  {  
    return volume;
  }

  std::string getName() const {return std::string("volume");}
private:
  double volume;
};

//  Volume of the quadratic mesh
template <class M>
class VolumeMQT<M,2>
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = 2;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = NDIM_W;

  // Here the embedding space has NDIM_W dims, we don't want the projected volume !
  typedef dice::QuadraticSimplexT<NDIM,NDIM_W> QS;
  
  
  VolumeMQT()
  {}

  void set (const Mesh *mesh, Simplex *s)
  {    
    qSimplex.set(s,s->segTracers.getConstPointer(),mesh->getGeometry());    
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    double bc[3]={a,b,c};
    return std::abs(qSimplex.evalJacobianDet(bc));
  }

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    double bc[4]={a,b,c,d};
    return std::abs(qSimplex.evalJacobianDet(bc));
  }

  std::string getName() const {return std::string("volume_Q");}
private:
  QS qSimplex; 
};

// Projected volume of the mesh
template <class M, int ELEMENT_ORDER=1>
class ProjectedVolumeMQT
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = ELEMENT_ORDER;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = 0;
  
  ProjectedVolumeMQT()
  {}

  void set(const Mesh *mesh, const Simplex *s)
  {
    volume = mesh->computeProjectedVolume(s);
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    return volume;
  }

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    return volume;
  }

  std::string getName() const {return std::string("projectedVolume");}
private:
  double volume;
};

// Projected Volume of the quadratic mesh
template <class M>
class ProjectedVolumeMQT<M,2>
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = 2;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = NDIM;

  typedef dice::QuadraticSimplexT<NDIM,NDIM> QS;
  
  
  ProjectedVolumeMQT()
  {}
  
  void set(const Mesh *mesh, const Simplex *s)
  {
    qSimplex.set(s,s->segTracers.getConstPointer(),mesh->getGeometry());
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    double bc[3]={a,b,c};
    return std::abs(qSimplex.evalJacobianDet(bc));
  }

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    double bc[4]={a,b,c,d};
    return std::abs(qSimplex.evalJacobianDet(bc));
  }

  std::string getName() const {return std::string("projectedVolume_Q");}
private:
  QS qSimplex; 
};

// Not correctly implemented yet
template <class M>
class InterVolumeMQT
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;
  static const int NSEG = Simplex::NSEG;

  static const int FE_ORDER = 2;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = NDIM+1;

  typedef dice::QuadraticSimplexT<NDIM,NDIM_W> QS;

  InterVolumeMQT()
  {}
  
  void set(const Mesh *mesh, const Simplex *s)
  {
    volume = mesh->computeVolume(s);
    qSimplex.set(s,s->segTracers.getConstPointer(),mesh->getGeometry());   
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    double bc[3]={a,b,c};
    Coord delta[NDIM_W];
    double d=0;

    qSimplex.barycentricToDeviation(bc,delta);
    for (int i=0;i<NDIM_W;++i)
      d+=delta[i]*delta[i];
    d=sqrt(d);

    double avgV = (volume + std::abs(qSimplex.evalJacobianDet(bc)))/2;
    // Need to multiple by sine of the angle here ...
    return d*avgV;
  }   

private:
  QS qSimplex; 
  double volume;
};

// Mass from density values at vertices (to check interpolation works correctly) 
template <class M>
class MassFromRhoMQT
{
public:
  typedef M Mesh;
  typedef typename Mesh::Simplex Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Coord Coord;
  
  static const int NDIM = Simplex::NDIM;
  static const int NDIM_W = Simplex::NDIM_W;
  static const int NVERT = Simplex::NVERT;

  static const int FE_ORDER = 1;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::simplex;
  static const int QUADRATURE_DEGREE = 2;

  MassFromRhoMQT()
  {}
  
  void set(const Mesh *mesh, const Simplex *s)
  {    
    double volume = mesh->computeProjectedVolume(s);
    for (int i=0;i<NVERT;++i)
      m[i]=volume*s->getVertex(i)->projectedDensity.getValue();      
  }

  // 2D
  double operator()(double a, double b, double c) const
  {
    return a*m[0]+b*m[1]+c*m[2];
  } 

  // 3D
  double operator()(double a, double b, double c, double d) const
  {
    return a*m[0]+b*m[1]+c*m[2]+d*m[3];
  }

  std::string getName() const {return std::string("massFromRho");}
private:
  double m[NVERT];  
};

// Compute radial profile F(r)
// Potential energy
template <class M>
class RadialProfileMT
{
  static const int NDIM=M::NDIM;
  static const int NDIM_W=M::NDIM_W;
  typedef typename M::Simplex Simplex;
  typedef typename M::Vertex Vertex;
public:
  RadialProfileMT(M* mesh, int samplesPerOverlap=10, 
		  double r0=0, bool useLog=false,
		  int nThreads=dice::glb::num_omp_threads)
  {    
    double x0[NDIM];
    double delta[NDIM];
    this->mesh=mesh;
    this->samplesPerOverlap=samplesPerOverlap;

    mesh->getBoundingBox(x0,delta,false);
    
    for (int i=0;i<M::NDIM;++i)
      center[i]=x0[i]+delta[i]*0.5;
    
    y.resize(nThreads);
    count.resize(nThreads);
    nIntervals=pow(mesh->getGlobalNCells(NDIM),1.0/NDIM);
    if (nIntervals>100) nIntervals=100;

    xMax=0;
    for (int i=0;i<M::NDIM;++i)
      {
	if (xMax<delta[i]) xMax=delta[i];
	//if (xMax<abs(x0[i])) xMax=abs(x0[i]);
	//if (xMax<abs(x0[i]+delta[i])) xMax=abs(x0[i]+delta[i]);
      }
    xMax/=2;

    for (int i=0;i<nThreads;++i)
      {
	y[i].resize(nIntervals,0);
	count[i].resize(nIntervals,0);
      }

    x.resize(nIntervals+1);        
    if (useLog)
      {	
	double logRMax=log10(xMax);
	double logR0=(r0<=0)?(-3):log10(r0);//logRMax-nIntervals/10;
	//printf("logR0 = %g %g\n",logR0,logRMax);
	for (int i=0;i<nIntervals;++i)
	  x[i]=pow(10.0,logR0+static_cast<double>(i)*(logRMax-logR0)/nIntervals);
      }
    else
      {	
	for (int i=0;i<nIntervals;++i)
	  x[i]=r0+static_cast<double>(i)*(xMax-r0)/nIntervals;
      }

    x.back()=xMax;
    xMaxInv=1.0/xMax;   

    seed.resize(nThreads);
    long initSeed=1321;
    for (int i=0;i<seed.size();++i) 
      dice::DRand48_rWapper::srand48_r(initSeed*i,&seed[i]);
     
    volume=delta[0];
    for (int i=1;i<NDIM;++i) volume*=delta[i];


    for (long i=0;i<x.size();++i)
      x[i]*=x[i];
    xMax*=xMax;
    xMaxInv*=xMaxInv;

  }

  void setCenterCoordinates(const double *c)
  {
    for (int i=0;i<M::NDIM;++i)
      center[i]=c[i];
  }

  void operator()(Simplex *s, int th) const
  {
    static const int N=samplesPerOverlap;
    static const int NMax=10*N;
    double c[NMax*NDIM];
    double v[NMax];
    double vertexValue[Simplex::NVERT];

    // Find the interval the simplex belongs to
    double rmin=x.back();
    double rmax=x.front();
    for (int i=0;i<Simplex::NVERT;++i)
      {
	auto coords=s->getVertex(i)->getCoordsConstPtr();
	double cCoords[NDIM];
	for (int j=0;j<NDIM;++j)
	  cCoords[j]=coords[j]-center[j];

	double r=cCoords[0]*cCoords[0];
	for (int j=1;j<NDIM;++j)
	  r+=cCoords[j]*cCoords[j];
	if (r<rmin) rmin=r;
	if (r>rmax) rmax=r;	
      }
    
    // The simplex is outside the range
    if (rmin>=x.back()) return;
    
    auto xEnd=std::upper_bound(x.begin(),x.end(),rmin);
    if (xEnd==x.begin()) ++xEnd;
    auto xBegin=xEnd-1;
   
    if (*xEnd > rmax)
      {
	// the simplex is fully contained in a single bin
	long index=std::distance(x.begin(),xBegin);
	count[th][index]+=N;
	y[th][index]+=s->mass.getValue();
	return;
      }
    else
      {
	do {++xEnd;} 
	while ((xEnd!=x.end())&&
	       ((*xEnd)<rmax));
      }    

    // Generate the sampling points
    for (int i=0;i<Simplex::NVERT;++i)
      vertexValue[i]=s->getVertex(i)->projectedDensity.getValue();

    int delta=std::distance(xBegin,xEnd)-1;
    long nSamples=N*delta;
    if (nSamples>NMax) nSamples=NMax;

    mesh->template generateRandomSample<Simplex,double,int,double*,double*,Simplex::NDIM>
      (s,vertexValue,nSamples,c,v,&seed[th]);
    
    // The factor to get sampling points mass from their interpolated density
    double factor=0;
    for (long i=0;i<nSamples;++i)
      factor+=v[i];
    factor = s->mass.getValue()/factor;

    // And make an histogram from the samples ...
    double *value=v;
    long index=0;
    for (long i=0;i<NDIM*nSamples;i+=NDIM)
      {
	double r=(c[i]-center[0])*(c[i]-center[0]);
	for (int j=1;j<M::NDIM;++j)
	  r+=(c[i+j]-center[j])*(c[i+j]-center[j]);
	
	if ((r<x.back())&&(r>x.front()))
	  {
	    auto it=std::upper_bound(xBegin,xEnd,r);
	    index=std::distance(x.begin(),it)-1;	  
	    ++count[th][index];
	    y[th][index]+=(*value)*factor;
	    ++value;
	  }
      }   
  }

  template <class Com>
  void toFile(const std::string &fName, 
	      Com *mpiCom, 
	      const std::string &comment="") const
  {
    std::vector<long> ct_sum(nIntervals,0);
    std::vector<double> val_sum(nIntervals,0);

    for (long i=0;i<nIntervals;++i)
      {
	ct_sum[i]=count[0][i];
	val_sum[i]=y[0][i];
	for (int j=1;j<y.size();++j)
	  {
	    val_sum[i]+=y[j][i];
	    ct_sum[i]+=count[j][i];
	  }
      }

    if (mpiCom->size()>1)
      {
	mpiCom->Reduce_inplace(ct_sum,0,MPI_SUM);
	mpiCom->Reduce_inplace(val_sum,0,MPI_SUM);
      }

    if (mpiCom->rank()==0) 
      {
	FILE *f=fopen(fName.c_str(),"w");

	if (comment.size()>0)
	  fprintf(f,"# radius value count (%s)\n",comment.c_str());
	else
	  fprintf(f,"# radius value count\n");

	for (long i=0;i<nIntervals;++i)
	  {
	    if (ct_sum[i]!=0)
	      {
		val_sum[i]/=shellVolume(sqrt(x[i]),sqrt(x[i+1]));
		fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x[i]),val_sum[i],ct_sum[i]);
	      }
	    else fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x[i]),0.0,ct_sum[i]);
	  }
	fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x.back()),0.0,0L);
	fclose(f);
      }
  }

  void toFile(const std::string &fName, const std::string &comment="") const
  {
    FILE *f=fopen(fName.c_str(),"w");

    if (comment.size()>0)
      fprintf(f,"# radius value count (%s)\n",comment.c_str());
    else
      fprintf(f,"# radius value count\n");

    for (long i=0;i<nIntervals;++i)
      {
	long ct=count[0][i];
	double val=y[0][i];
	for (int j=1;j<y.size();++j)
	  {
	    val+=y[j][i];
	    ct+=count[j][i];
	  }

	if (ct!=0)
	  {
	    /*
	    if (i==0) 
	      {
		
		printf("x = %g %g, vol=%g\n val=%g ct=%ld",x[i],x[i+1],volume,val,ct);
	      }
*/
	    val/=shellVolume(sqrt(x[i]),sqrt(x[i+1]));
	    fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x[i]),val,ct);
	  }
	else fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x[i]),0.0,ct);
      }
    fprintf(f,"%15.15lg %15.15lg %ld\n",sqrt(x.back()),0.0,0L);
    fclose(f);
  }

private:
  double shellVolume(double r0, double r1) const
  {
    if (M::NDIM==3)
      return 4.0/3.0*M_PI*(r1*r1*r1-r0*r0*r0);
    else
      return M_PI*(r1*r1-r0*r0);
  }

  M *mesh;
  std::vector<double> coords[M::NDIM];
  mutable std::vector< std::vector<double> > y;
  mutable std::vector< std::vector<long> > count;
  mutable std::vector<typename dice::DRand48_rWapper::RandData> seed;
  std::vector<double> x;
  double xMax;
  double xMaxInv;
  long nIntervals;
  double volume;
  int samplesPerOverlap;
  double center[M::NDIM];
};

#endif

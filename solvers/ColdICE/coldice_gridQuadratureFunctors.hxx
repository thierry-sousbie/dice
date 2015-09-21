#ifndef __COLDICE_GRID_QUADRATURE_FUNCTORS_HXX__
#define __COLDICE_GRID_QUADRATURE_FUNCTORS_HXX__

#include <stdio.h>
#include <functional>
#include <dice/finiteElements/finiteElementTypeEnum.hxx>

// Mass
template <class G, class K, int DEGREE=1>
class MassGQT
{
public:
  static const int FE_ORDER = 1;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::cuboid;
  static const int QUADRATURE_DEGREE = DEGREE;
  
  MassGQT(G *grid):
    g(grid)
  {
    grid->initializeInterpolationKernel(k);
    volume=grid->getConstantVolumeElement();
  }

  double operator()(int field, double x, double y)
  {
    double c[2]={x,y};
    double result;
    g->applyKernel(k,c,&result,field);
    return result*volume;   
  }

  double operator()(int field, double x, double y, double z)
  {
    double c[3]={x,y,z};
    double result;
    g->applyKernel(k,c,&result,field);
    return result*volume;
  }

private:
  double volume;
  K k;
  G *g;
};

// Potential energy
template <class G, class K, int DEGREE=2>
class PotentialEnergyGQT
{
public:
  static const int FE_ORDER = 1;
  static const dice::fETypeE::Type FE_TYPE = dice::fETypeE::cuboid;
  static const int QUADRATURE_DEGREE = DEGREE;
  
  PotentialEnergyGQT(G *density, G *potential):
    d(density),
    p(potential)
  {
    potential->initializeInterpolationKernel(pKernel);
    density->initializeInterpolationKernel(dKernel);

    alpha=-0.5*potential->getConstantVolumeElement();
  }

  double operator()(int field, double x, double y)
  {
    double c[2]={x,y};
    double density;
    double potential;
    
    d->applyKernel(dKernel,c,&density,field);
    p->applyKernel(pKernel,c,&potential,field);
    
    return alpha*potential*density;
  }

  double operator()(int field, double x, double y, double z)
  {
    double c[3]={x,y,z};
    double density;
    double potential;
    
    d->applyKernel(dKernel,c,&density,field);
    p->applyKernel(pKernel,c,&potential,field);
    
    return alpha*potential*density;    
  }

private:  
  K pKernel;
  K dKernel;
  G *d;
  G *p;
  double alpha;
};

// Compute radial profile F(r)
// Potential energy
template <class G>
class RadialProfileGT
{
public:
  RadialProfileGT(G* grid, int nThreads=dice::glb::num_omp_threads)
  {
    this->grid=grid;
    for (int i=0;i<G::NDIM;++i)
      coords[i]=grid->getCellCoord(i);

    y.resize(nThreads);
    count.resize(nThreads);

    nIntervals=0;
    xMax=0;
    //double xMin=0;
    for (int i=0;i<G::NDIM;++i)
      {	
	if (nIntervals<grid->getCellCoord(i).size())
	  nIntervals=grid->getCellCoord(i).size();
	
	double a=grid->getVertexCoord(i).front();
	double b=grid->getVertexCoord(i).back();

	if (xMax<(b-a)) xMax=b-a;
	center[i]=(a+b)/2;
      }
    xMax/=2;
    nIntervals/=2;

    for (int i=0;i<nThreads;++i)
      {
	y[i].resize(nIntervals,0);
	count[i].resize(nIntervals,0);
      }    

    x.resize(nIntervals+1);    
    for (int i=0;i<nIntervals;++i)
      x[i]=static_cast<double>(i)*xMax/nIntervals;
    x.back()=xMax;
    xMaxInv=1.0/xMax;
  }

  void setCenterCoordinates(const double *c)
  {
    for (int i=0;i<G::NDIM;++i)
      center[i]=c[i];
  }

  template <class IT>
  void operator()(IT &it, int th=0) const
  {
    double c[G::NDIM];
    getCoords(it,c);
    double r=(c[0]-center[0])*(c[0]-center[0]);
    for (int i=1;i<G::NDIM;++i)
      r+=(c[i]-center[i])*(c[i]-center[i]);
    r=sqrt(r);
    //printf("R=%g / %g += %g\n",r,xMax,(*it));
    if (r<xMax)
      {
	int index=static_cast<int>(r*xMaxInv*nIntervals);
	++count[th][index];
	y[th][index]+=(*it);
      }
  }
  
  template <class Com>
  void toFile(const std::string &fName,
	      Com *mpiCom,
	      const std::string &comment="")
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
		val_sum[i]/=ct_sum[i];//shellVolume(x[i],x[i+1])
		fprintf(f,"%15.15lg %15.15lg %ld\n",x[i],val_sum[i],ct_sum[i]);
	      }
	    else fprintf(f,"%15.15lg %15.15lg %ld\n",x[i],0.0,ct_sum[i]);
	  }
	fprintf(f,"%15.15lg %15.15lg %ld\n",x.back(),0.0,0L);
	fclose(f);
      }
    
  }
  
  void toFile(const std::string &fName, 
	      const std::string &comment="") const
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
	    val/=ct;//shellVolume(x[i],x[i+1])
	    fprintf(f,"%15.15lg %15.15lg %ld\n",x[i],val,ct);
	  }
	else fprintf(f,"%15.15lg %15.15lg %ld\n",x[i],0.0,ct);
      }
    fprintf(f,"%15.15lg %15.15lg %ld\n",x.back(),0.0,0L);
    fclose(f);
  }

private:  

  double shellVolume(double r0, double r1) const
  {
    if (G::NDIM==3)
      return 4.0/3.0*M_PI*(r1*r1*r1-r0*r0*r0);
    else
      return 4.0*M_PI*(r1*r1-r0*r0);
  }

  template <class IT>
  void getCoords(IT &it, double *c) const
  {
    const int *w=it.get_w();
    for (int i=0;i<G::NDIM;++i)
      c[i]=coords[i][w[i]];
  }
  
  G *grid;
  std::vector<double> coords[G::NDIM];
  mutable std::vector< std::vector<double> > y;
  mutable std::vector< std::vector<long> > count;
  std::vector<double> x;
  double xMax;
  double xMaxInv;
  long nIntervals;
  double center[G::NDIM];

};

#endif

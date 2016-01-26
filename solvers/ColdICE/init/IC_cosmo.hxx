#ifndef __IC_COSMO_HXX__
#define __IC_COSMO_HXX__

#include <fstream>

#include <math.h>
#include <dice/tools/IO/IOHelpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>
#include <dice/mesh/tesselation/simplicialGrid.hxx>
#include "initialConditionsInterface.hxx"

template <class M>
class IC_cosmo : 
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

  static std::string name() {return "cosmo";}
  std::string getName() const {return name();}

  template <class PM, class R>
  IC_cosmo(PM &paramsManager, R *reader, 
	   dice::MpiCommunication *com=0): //nullptr
    com_(com),
    velocityFactor_(0)
  { 
    parse(paramsManager,reader);
  }

  long getNPasses() const
  {
    long n=static_cast<long>(cum.back()/chunkSize_);
    
    if (chunkSize_*n < cum.back())
      return n+1;
    else
      return n;
  }

  void initialize(const double *x0, const double *delta, double t, 
		  const Units &units, int pass) 
  {
    if (pass==0)
      {
	if (readHeaders(fname_)<1)
	  {
	    dice::PRINT_SRC_INFO(dice::LOG_ERROR);
	    dice::glb::console->print<dice::LOG_ERROR>
	      ("Invalid filename: '%s'\n",fname_.c_str());
	    throw std::runtime_error("Invalid filename");
	  }
	
	std::copy_n(x0,NDIM,x0_);
	std::copy_n(delta,NDIM,delta_);     
	for (int i=0;i<NDIM;++i)
	  {
	    epsilon_[i]=0.01*(delta[i]/iLen_);
	    deltaInv_[i]=double(iLen_)/delta[i];
	    //x0_[i]-=epsilon_[i];
	  }

	velocityFactor_=sqrt(headers[0].time);
      }
    
    startIndex_=pass*chunkSize_;
    deltaIndex_=readChunks(startIndex_,chunkSize_);
    stopIndex_=startIndex_+deltaIndex_;
  }
  
  ImplicitTesselation *createImplicitTesselation()
  { 
    typedef typename ImplicitTesselation::Cell Cell;
    typedef dice::SimplicialGridT<NDIM,NDIM_W,Cell,M::IS_PERIODIC> SimplicialGrid;
    typedef typename SimplicialGrid::Params SimplicialGridParams;

    if (readHeaders(fname_)<1)
      {
	dice::PRINT_SRC_INFO(dice::LOG_ERROR);
	dice::glb::console->print<dice::LOG_ERROR>
	  ("Invalid filename: '%s'\n",fname_.c_str());
	throw std::runtime_error("Invalid filename");
      }

    GadgetHeader &h=headers[0];

    SimplicialGridParams sgParams;
    std::copy_n(x0_,NDIM,sgParams.x0);
    //std::fill_n(sgParams.x0,NDIM,0);
    std::fill_n(sgParams.delta,NDIM,h.BoxSize);
    std::fill_n(sgParams.resolution,NDIM,iLen_/2);
    sgParams.t=initTesselationType;

    SimplicialGrid *sg = new SimplicialGrid(sgParams);
    return static_cast<ImplicitTesselation*>(sg);
  }

  void release()
  {
    std::vector< GadgetHeader > eHeaders;
    std::vector< std::string > eFullName;
    std::vector< unsigned long > eCum;
    
    std::vector<double> ePos;    
    std::vector<double> eVel;
    
    eHeaders.swap(headers);
    eFullName.swap(fullName);
    eCum.swap(cum);
    
    ePos.swap(pos);
    eVel.swap(vel);
  }

  Params getParams()
  {
    // It is OK to query parameters when IC files do not exist as we may be 
    // restarting a run, in which case they are not needed.
    if (readHeaders(fname_)<1)
      return Params();
    
    GadgetHeader &h=headers[0];
    Params p;

    p.useCosmo=1;
    p.aStart=h.time;
    p.mass=h.mass[1]*cum.back();
 
    p.unitLength=3.085678e22;    // 1.0 Mpc
    p.unitVelocity=1.e3;         // 1km/s
    p.unitMass=1.989e40;         // 1e10 solar masses

    p.defaultUnitLength="1Mpc";
    p.defaultUnitVelocity="1km/s";
    p.defaultUnitMass="1E10 M_sol";

    p.G=6.67384e-11; // in m3 kg-1 s-2
    p.H=3.2409e-18;  // value in mks for h0=1

    p.cosmoParams.set(h.Omega0,h.OmegaLambda,0.0,h.HubbleParam);    
    
    std::fill_n(&p.x0.front(),NDIM,0);
    std::fill_n(&p.delta.front(),NDIM,h.BoxSize);
    for (int i=0;i<NDIM;++i) p.x0[i]=x0_[i];
    // iLen is divided by 2 because we also need to displace tracers.
    // std::fill_n(&p.resolution.front(),NDIM,iLen_/2);
    return p;
  }
  
  void displace(double *coords, double &density, Vertex *v) 
  {
    // Sorry, that's ugly ;(
    // but we need a way to know if the particle was already displaced during a previous pass!
    if ((coords[NDIM]==0)&&(coords[NDIM+1]==0))
      {
	long index = coords2index(coords)-startIndex_;
	//printf("(%g %g %g) => I=%ld\n",coords[0],coords[1],coords[2],index);
	
	if ((index>=0)&&(index<deltaIndex_))
	  {
	    //std::copy_n(&pos[index*3],NDIM,coords);
	    //std::copy_n(&vel[index*3],NDIM,coords+NDIM);
	    for (int i=0;i<NDIM;++i) 
	      coords[i]=pos[index*3+i]+x0_[i];
	    for (int i=0;i<NDIM;++i) 
	      coords[i+NDIM]=vel[index*3+i]*velocityFactor_;			     
	  }
	/*
	printf("->(%g %g %g %g %g %g)\n",
	       coords[0],coords[1],coords[2],
	       coords[3],coords[4],coords[5]);
	*/
      }
  } 

  bool useLagrangianMass() const {return true;}
  
protected:
  dice::TesselationType initTesselationType;
  std::string fname_;
  double maxMemPerNode_;
  long chunkSize_;
  bool rowMajor_;
  long startIndex_;
  long stopIndex_;
  long deltaIndex_;

  template <class PM, class R>
  void parse(PM &paramsManager, R *reader)
  {
    fname_ = "ics_gadget.dat";    
    fname_ = paramsManager.
      get("fileName",Interface::parserCategory(),fname_,
	  reader,PM::FILE_FIRST,
	  "Gagdet file name (without .N extension for multiple files).");

    maxMemPerNode_ = 1;
    maxMemPerNode_ = paramsManager.
      get("maxMemPerNode",Interface::parserCategory(),maxMemPerNode_,
	  reader,PM::FILE_FIRST,
	  "The maximum amount of memory to use on each MPI node, expressed in GigaBytes.");

    rowMajor_ = 1;
    rowMajor_ = paramsManager.
      get("rowMajor",Interface::parserCategory(),rowMajor_,
	  reader,PM::FILE_FIRST,
	  "Use row major particle ordering if true (i.e. index of perticle changes faster with X coordinate). Default is rowMinor.");

    initTesselationType = dice::TesselationType::ANY;
    std::string initTesselationTypeStr = paramsManager.template 
      get<std::string>("initTesselationType",Interface::parserCategory(),
		       dice::TesselationTypeSelect().getString(initTesselationType,true),
		       reader,PM::FILE_FIRST,
		       dice::TesselationTypeSelect().getAllString
		       ("Topology of the initial grid simplicial tesselation (%s)"));
    initTesselationType = 
      dice::TesselationTypeSelect().getVal(initTesselationTypeStr,true); 
    
    bool haveHeaders=false;
    if (readHeaders(fname_)>=1)
      haveHeaders=1;

    for (int i=0;i<NDIM;++i)
      {
	x0_[i] = (haveHeaders)?(-headers[0].BoxSize/2):0;
	x0_[i] = paramsManager.
	  get("x0",Interface::parserCategory(),x0_[i],i,
	      reader,PM::FILE_FIRST,
	      "Position of the lower left corner of the box");
      }
    
  }

  long readChunks(long start, long sz)
  {
    pos.resize(sz*3);
    vel.resize(sz*3);
    
    long nRead=0;
  
    if ((com_==0)||(com_->rank()==0)) //nullptr
      {
	readChunk(start,sz,&pos[0],0);
	nRead=readChunk(start,sz,&vel[0],1);
      }
    
    if ((com_!=0)&&(com_->size()>1)) //nullptr
      {
	com_->Bcast(&nRead,0,1);
	com_->Bcast(pos,0,nRead*3);
	com_->Bcast(vel,0,nRead*3);     
      }

    return nRead;
  }

  unsigned long coords2index(double *coords) const
  {
    unsigned long index=stride_[0]*static_cast<unsigned long>
      ((coords[0]-x0_[0]+epsilon_[0])*deltaInv_[0]);
    for (int i=1;i<NDIM;++i)
      index+=stride_[i]*static_cast<unsigned long>
	((coords[i]-x0_[i]+epsilon_[i])*deltaInv_[i]);
    
    return index;
  }
  
private:
  dice::MpiCommunication *com_;
  double x0_[3];
  double delta_[3];  
  double deltaInv_[3];  
  long iLen_;
  long stride_[3];
  double epsilon_[3];
  double velocityFactor_;

  static const int fillBytes = 60;
  struct GadgetHeader
  {
    int npart[6];                        
    double mass[6];                      
    double time;                         
    double redshift;                     
    int flag_sfr;                        
    int flag_feedback;                   
    unsigned int npartTotal[6];          
    int flag_cooling;                    
    int num_files;                       
    double BoxSize;                      
    double Omega0;                       
    double OmegaLambda;                  
    double HubbleParam;                  
    int flag_stellarage;                 
    int flag_metals;                     
    unsigned int npartTotalHighWord[6];  
    int  flag_entropy_instead_u;         
    char fill[fillBytes];                       
  };                       
  
  bool useFloats;
  long floatSize;

  std::vector< GadgetHeader > headers;
  std::vector< std::string > fullName;
  std::vector< unsigned long > cum;

  std::vector<double> pos;
  std::vector<double> vel;

  long getFileNames(const std::string &name, int index=-1)
  {
    std::stringstream ss;
    ss<<name;
    if (index>=0) ss << "." << index;
    else fullName.clear();

    if (!dice::myIO::exists(ss.str()))
      {
	if (index<0) return getFileNames(name,index+1);
	else return index;
      }
    else 
      {
	fullName.push_back(ss.str());	      
	if (index<0) return 1;
      }
 
    return getFileNames(name,index+1);    
  }

  long readHeaders(const std::string &name)
  {
    long N = getFileNames(name);
    headers.resize(N);
    
    if (N<1) return N;

    for (long i=0;i<fullName.size();++i)
      {
	std::ifstream is(fullName[i].c_str(),std::ifstream::binary);
	
	int dummy;
	is.read((char*)&dummy,sizeof(int));
	if (dummy!=256)
	  {
	    dice::glb::console->print<dice::LOG_WARNING>
	      ("File '%s' is an invalid gadget file !\n",
	       fullName[i].c_str());
	    if (dummy==(256L<<24)) 
	      dice::glb::console->print<dice::LOG_WARNING>
	      ("File and machine endianness differ.\n");
	    return 0;
	  }
	is.read((char*)&headers[i],sizeof(GadgetHeader));
	is.read((char*)&dummy,sizeof(int));
      } 

    unsigned long npartTotal=
      (static_cast<unsigned long>(headers[0].npartTotalHighWord[1])<<32)+ 
      headers[0].npartTotal[1];

    cum.resize(N+1);cum[0]=0;
    if (N==1) 
      cum[1]=npartTotal;
    else
      {
	for (long i=0;i<headers.size();++i)
	  cum[i+1]=cum[i]+headers[i].npart[1];
      }
    
    if (cum.back() != npartTotal)
      {
	dice::PRINT_SRC_INFO(dice::LOG_ERROR);
	dice::glb::console->print<dice::LOG_ERROR>
	  ("Inconsistent header: total npart = %ld, counted = %ld\n",
	   npartTotal,cum.back());
	throw std::runtime_error("Invalid filename");
      }

    std::ifstream is(fullName[0].c_str(),std::ifstream::binary);
    is.seekg(0, is.end);
    long fileSize=is.tellg();
    
    long nPart1=cum[1]-cum[0];
    long dLenI=sizeof(int)*8+nPart1*(6*sizeof(double)+sizeof(int))+sizeof(GadgetHeader);
    long dLenL=sizeof(int)*8+nPart1*(6*sizeof(double)+sizeof(long))+sizeof(GadgetHeader);
    
    if ((fileSize==dLenI)||(fileSize==dLenL))
      {
	// dice::glb::console->print<dice::LOG_INFO>
	//   ("Using doubles: fsize=%ld (%ld/%ld)\n",fileSize,dLenI,dLenL);
	useFloats=false;
	floatSize=sizeof(double);
	chunkSize_ = ((1ul<<30)*maxMemPerNode_)/(sizeof(double)*6);
      }
    else
      {
	// dice::glb::console->print<dice::LOG_INFO>
	//   ("Using floats: fsize=%ld (%ld/%ld)\n",fileSize,dLenI,dLenL);
	useFloats=true;
	floatSize=sizeof(float);
	chunkSize_ = ((1ul<<30)*maxMemPerNode_)/((sizeof(double)*6+sizeof(float)*3));
      }
    
    if (cum.back()<chunkSize_)
      chunkSize_=cum.back();
    
    iLen_=pow(cum.back(),1.0/NDIM);
    while (iLen_*iLen_*iLen_<cum.back()) ++iLen_;
    if (iLen_*iLen_*iLen_<cum.back())
      {
	dice::PRINT_SRC_INFO(dice::LOG_ERROR);
	dice::glb::console->print<dice::LOG_ERROR>
	  ("The number of particle should be a power of 3! (N=%ld)",cum.back());
	throw std::runtime_error("Invalid particles count filename");
      }

    stride_[0]=1;
    for (int i=1;i<NDIM;++i) stride_[i]=stride_[i-1]*iLen_;
    
    if (rowMajor_)
      {
	for (int i=0;i<NDIM/2;++i)
	  std::swap(stride_[i],stride_[NDIM-i-1]);
      }

    return N;
  }

  long readChunk(long start, long sz, double *ptr, long which)
  {
    if ((start>=cum.back())||(start<0)) return 0;

    auto it=std::upper_bound(cum.begin(),cum.end(),start);
    long fIndex=std::distance(cum.begin(),it)-1;
    long nSkip=(start-cum[fIndex]);
    long nLeft=sz;

    while ((nLeft>0)&&(fIndex<headers.size()))
      {
	long nCurLeft=(cum[fIndex+1]-cum[fIndex]);
	long delta=sizeof(int)*(3+2*which)+sizeof(GadgetHeader)+
	  floatSize*(cum[fIndex+1]-cum[fIndex])*3*which;
	std::ifstream is(fullName[fIndex].c_str(),std::ifstream::binary);
	is.seekg(delta+nSkip*floatSize*3, is.beg);
	nCurLeft -= nSkip;
	nSkip=0;
	
	long nRead=std::min(nCurLeft,nLeft);
	if (useFloats)
	  {
	    float *tmp=new float[nRead*3];
	    is.read((char*)tmp,nRead*floatSize*3);
	    std::copy(tmp,tmp+nRead*3,ptr);
	    delete[] tmp;
	  }
	else is.read((char*)ptr,nRead*floatSize*3);
	ptr+=nRead*3;
	nLeft-=nRead;
	++fIndex;
      }

    return sz-nLeft;
  }

  
};

#endif

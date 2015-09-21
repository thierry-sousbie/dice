#ifndef __MESH_PARAMS_HXX__
#define __MESH_PARAMS_HXX__

#include <string>
#include "../dice_globals.hxx"

//#include "../tools/IO/paramsManager.hxx"
//#include "../tools/IO/myIO.hxx"

#include "../partition/partitionType.hxx"

#include "../mesh/tesselation/tesselationType.hxx"

/**
 * @file 
 * @brief  A class used to define and manage parameters used to setup meshT class.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/** 
 * \class MeshParamsT
 * \brief  A class used to define and manage parameters used to setup meshT class.
 */

template <int D, int DW, typename C> 
struct MeshParamsT
{
  typedef MeshParamsT<D,DW,C> MyType;

  static const int NDIM = D;
  static const int NDIM_W = DW;

  static std::string parserCategory() {return "mesh";}
  static std::string classHeader() {return "mesh_params";}
  static float classVersion() {return 0.10;}
  static float compatibleSinceClassVersion() {return 0.10;}

  typedef C Coord;
  //typedef ParamsManagerT<ParamsParser,Console,LOG_INFO> ParamsManager;

  double allocFactor; //!< fraction of new object to allocate when container is full
  Coord x0[NDIM_W];    //!< box origin
  Coord delta[NDIM_W]; //!< box size
  //long resolution[NDIM];  //!< initial box resolution in pixels
  PartitionType initPartitionType; //!< initial partition type
  RefinePartitionType refinePartitionType; //!< load balancing method
  //TesselationType initTesselationType; //!< initial tesselation topology
  double initPartitionTolerance; //!< tolerance on work load for initial partition
  double repartTolerance; //!< tolerance on work load for partitions (after repartitionning)
  double repartThreshold; //!< Repartition is triggered when max(nSimplex)/min(nSimplex) > repartThreshold
  
  MeshParamsT()
  {
    setDefault();
    serializedVersion=classVersion();
  }

  void setSerializedVersion(float ver)
    {
      serializedVersion=ver;
    }
  
  /*
  MeshParamsT(const MyType &defaultParams, 
	      ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {
    setDefault(defaultParams);
    parse(myIO::BinaryReaderT<>::nullReader()); 
  }
  */
  /*
  template <class R>
  MeshParamsT(const MyType &defaultParams, R *reader, 
	      ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {
    setDefault(defaultParams);
    parse(reader);
  }
    
  MeshParamsT(ParamsParser *parser):
    manager(parser, glb::console)
  {
    setDefault();
    parse(myIO::BinaryReaderT<>::nullReader());
  }
 
  template <class R>
  MeshParamsT(R *reader, ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {    
    setDefault();
    parse(reader);
  }
  */
  /*
  MeshParamsT()
  {
    setDefault();      
  }
 
  template <class R>
  MeshParamsT(const MyType &defaultParams)
  {
    setDefault(defaultParams);   
  }
  */
  ~MeshParamsT()
  {
  }
 
  /*
  template <class LT>
  void report() const
  {
    manager.report<LT>();
  }
  */

  void setDefault()
  {
    for (int i=0;i<NDIM_W;i++)
      {
	x0[i]=(i<NDIM)?-1:0; 
	delta[i]=(i<NDIM)?2:0;
	//if (i<NDIM) resolution[i]=64;
      }
    initPartitionTolerance=1.01; // Allow 1% imbalance on initial partiton
    repartTolerance=initPartitionTolerance;
    repartThreshold=1.15; // Allow 15% imbalance max
    allocFactor=1.0; // Alloc 100% new objects when a container is full (->double the size)

    //initTesselationType = TesselationType::ANY;
    initPartitionType = PartitionType::KWAY;
    refinePartitionType = RefinePartitionType::ADAPTIVE;  
  }
  /*
  void setDefault(const MyType &def)
  {
    for (int i=0;i<NDIM_W;i++)
      {
	x0[i]=def.x0[i];
	delta[i]=def.delta[i];
	if (i<NDIM) resolution[i]=def.resolution[i];
      }
    initPartitionTolerance=def.initPartitionTolerance;
    repartTolerance=def.repartTolerance;
    repartThreshold=def.repartThreshold;
    allocFactor=def.allocFactor;
    initTesselationType=def.initTesselationType;
    initPartitionType=def.initPartitionType;
    refinePartitionType=def.refinePartitionType;
  }
  */
  template <class BW>
  void write(BW *writer) const
  {
    writer->writeHeader(classHeader(),classVersion());
    writer->write(x0,NDIM_W);
    writer->write(delta,NDIM_W);
    //writer->write(resolution,NDIM);

    int tmp;

    tmp=initPartitionType;writer->write(&tmp);
    //tmp=refinePartitionType;writer->write(&tmp);
    //tmp=initTesselationType;writer->write(&tmp);

    writer->write(&initPartitionTolerance);
    //writer->write(&repartTolerance);
    //writer->write(&repartThreshold);
  }

  template <class BR>
  void read(BR *reader)
  {
    float version;
    BR::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
      (glb::console,reader,version,true);
    
    reader->read(x0,NDIM_W);
    reader->read(delta,NDIM_W);
    //reader->read(resolution,NDIM);

    int tmp;

    reader->read(&tmp);initPartitionType=static_cast<PartitionType>(tmp);
    //reader->read(&tmp);refinePartitionType=tmp;
    //reader->read(&tmp);initTesselationType=static_cast<TesselationType>(tmp);

    reader->read(&initPartitionTolerance);
    //reader->read(&repartTolerance);
    //reader->read(&repartThreshold);
  }

  template <class BR, class PM>
  void parse(BR *reader, PM &manager)
  {    
    // float version;
    // BR::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
    //   (glb::console,reader,version,true);

    for (int i=0;i<NDIM_W;i++)
      {
	x0[i]=manager.
	  get("x0",parserCategory(),x0[i],i,
	      reader,PM::FILE_FIRST,
	      "Initial coordinates of the bounding box lower left corner");
	
	delta[i]=manager.
	  get("delta",parserCategory(),delta[i],i,
	      reader,PM::FILE_FIRST,
	      "Size of the bounding box");			 
	/*
	if (i<NDIM) 
	  resolution[i]=manager.
	    get("resolution",parserCategory(),resolution[i],i,
		reader,PM::FILE_FIRST,
		"Resolution of the initial mesh (in pixels before tesselation)");
	*/
      }   

    allocFactor=manager.
      get("allocFactor",parserCategory(),allocFactor,
	  reader,PM::PARSER_FIRST,
	  "What fraction of the current number of element to reallocate when a pool is full.");
    
    initPartitionTolerance=manager.
      get("initPartitionTolerance",parserCategory(),initPartitionTolerance,
	  reader,PM::FILE_FIRST,
	  "Tolerance on the load balance of the initial partition");    
    
    repartTolerance=manager.
      get("repartTolerance",parserCategory(),repartTolerance,
	  reader,PM::PARSER_FIRST,
	  "Tolerance on the load balance after repartitionning");
		     
    repartThreshold=manager.
      get("repartThreshold",parserCategory(),repartThreshold,
	  reader,PM::PARSER_FIRST,
	  "Inbalance factor that triggers repartitionning");
    /*	     
    std::string initTesselationTypeStr = manager.template 
      get<std::string>("initTesselationType",parserCategory(),
		       TesselationTypeSelect().getString(initTesselationType,true),
		       reader,PM::FILE_FIRST,
		       TesselationTypeSelect().getAllString
		       ("Type of initial grid tesselation into simplices (%s)"));
    initTesselationType = TesselationTypeSelect().getVal(initTesselationTypeStr,true); 
    */
    std::string initPartitionTypeStr = manager.template 
      get<std::string>("initPartitionType",parserCategory(),
		       PartitionTypeSelect().getString(initPartitionType,true), 
		       reader,PM::FILE_FIRST,
		       PartitionTypeSelect().getAllString
		       ("The method to use for initial partitionning (%s)"));
    initPartitionType = PartitionTypeSelect().getVal(initPartitionTypeStr,true); 

    std::string refinePartitionTypeStr = manager.template 
      get<std::string>("refinePartitionType",parserCategory(),
		       RefinePartitionTypeSelect().getString(refinePartitionType,true),
		       reader,PM::PARSER_FIRST,
		       RefinePartitionTypeSelect().getAllString
		       ("The method to use for repartitionning (%s)"));
    refinePartitionType = RefinePartitionTypeSelect().getVal(refinePartitionTypeStr,true); 
  }

  template <class PP>
  void parse(PP &parser)
  {    
    for (int i=0;i<NDIM_W;i++)
      {
	x0[i]=parser.
	  get("x0",parserCategory(),x0[i],i,
	      "Initial coordinates of the bounding box lower left corner");
	
	delta[i]=parser.
	  get("delta",parserCategory(),delta[i],i,
	      "Size of the bounding box");			 
	/*
	if (i<NDIM) 
	  resolution[i]=parser.
	    get("resolution",parserCategory(),resolution[i],i,
		"Resolution of the initial mesh (in pixels before tesselation)");
	*/
      }   

    allocFactor=parser.
      get("allocFactor",parserCategory(),allocFactor,
	  "What fraction of the current number of element to reallocate when a pool is full.");
    
    initPartitionTolerance=parser.
      get("initPartitionTolerance",parserCategory(),initPartitionTolerance,
	  "Tolerance on the load balance of the initial partition");    
    
    repartTolerance=parser.
      get("repartTolerance",parserCategory(),repartTolerance,
	  "Tolerance on the load balance after repartitionning");
		     
    repartThreshold=parser.
      get("repartThreshold",parserCategory(),repartThreshold,
	  "Inbalance factor that triggers repartitionning");
    /*		     
    std::string initTesselationTypeStr = parser.template 
      get<std::string>("initTesselationType",parserCategory(),
		       TesselationTypeSelect().getString(initTesselationType,true),
		       TesselationTypeSelect().getAllString
		       ("Topology of the initial simplicial tesselation of the grid (%s)"));
    initTesselationType = TesselationTypeSelect().getVal(initTesselationTypeStr,true); 
    */
    std::string initPartitionTypeStr = parser.template 
      get<std::string>("initPartitionType",parserCategory(),
		       PartitionTypeSelect().getString(initPartitionType,true),
		       PartitionTypeSelect().getAllString
		       ("The method to use for initial partitionning (%s)"));
    initPartitionType = PartitionTypeSelect().getVal(initPartitionTypeStr,true); 

    std::string refinePartitionTypeStr = parser.template 
      get<std::string>("refinePartitionType",parserCategory(),
		       RefinePartitionTypeSelect().getString(refinePartitionType,true),
		       RefinePartitionTypeSelect().getAllString
		       ("The method to use for repartitionning (%s)"));
    refinePartitionType = RefinePartitionTypeSelect().getVal(refinePartitionTypeStr,true); 
  }
private:
  float serializedVersion;
};

/** \}*/
#include "../internal/namespace.footer"
#endif

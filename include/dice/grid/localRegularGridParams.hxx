#ifndef __LOCAL_REGULAR_GRID_PARAMS_HXX__
#define __LOCAL_REGULAR_GRID_PARAMS_HXX__

#include <string>
#include "../dice_globals.hxx"

#include "../grid/scale.hxx"
#include "../grid/valLocationType.hxx"

#include "../tools/IO/paramsManager.hxx"
#include "../tools/IO/myIO.hxx"

#include "../internal/namespace.header"

template <int ND>
struct LocalRegularGridParamsT
{
  typedef LocalRegularGridParamsT<ND> MyType;

  static const int NDIM = ND;

  static std::string parserCategory() { return "localGrid"; }
  static std::string classHeader() { return "local_regular_grid__params"; }
  static float classVersion() { return 0.10; }
  static float compatibleSinceClassVersion() { return 0.10; }

  // typedef ParamsManagerT<ParamsParser,Console,LOG_INFO> ParamsManager;

  typedef ScaleT<double> Scale;
  typedef typename Scale::ScaleType ScaleType;
  typedef typename Scale::ScaleTypeV ScaleTypeV;

  double x0[NDIM];                   //!< coordinate of the leftmost vertex
  double delta[NDIM];                //!< size of the box
  int resolution[NDIM];              //!< resolution along each axis
  ScaleType scale[NDIM];             //!< scale type along each axis
  ValLocationType valLocation[NDIM]; //!< location of values (cells or vertex)

  int lowMargin[NDIM];  //!< Number of voxels in the low margin
  int highMargin[NDIM]; //!< Number of voxels in the high margin
  int nFields;          //!< Number of fields per voxel

  int haveParentGrid;         //!< True if the grid is a subgrid of a larger one
  int parentNDim;             //!< Number of dimensions in the parent grid
  int position[NDIM];         //!< position of the lower left corner within the parent grid
  int parentResolution[NDIM]; //!< resolution of the parent grid
  long minElementsCount;      //!< The minimum size that should be allocated for the grid

  double parentX0[NDIM];    //!< coordinate of the leftmost vertex of parent grid (if haveParentGrid)
  double parentDelta[NDIM]; //!< size fo the parent grid box (if haveParentGrid)

  template <class LOG>
  void print() const
  {
    if (NDIM == 3)
    {
      glb::console->print<LOG>("Resolution : [%g %g %g] [%g %g %g] [%d,%d,%d]\n",
                               x0[0], x0[1], x0[2], delta[0], delta[1], delta[2],
                               resolution[0], resolution[1], resolution[2]);
      glb::console->print<LOG>("Parents : [%d %d %d] [%d,%d,%d]\n",
                               position[0], position[1], position[2],
                               resolution[0], resolution[1], resolution[2]);
    }
    else
    {
      glb::console->print<LOG>("Resolution : [%g %g] [%g %g] [%d %d]\n",
                               x0[0], x0[1], delta[0], delta[1],
                               resolution[0], resolution[1]);
      glb::console->print<LOG>("Parents : [%d %d] [%d %d]\n",
                               position[0], position[1],
                               resolution[0], resolution[1]);
    }
  }

  LocalRegularGridParamsT()
  {
    setDefault();
  }
  /*
    // We don't want to specify a copy constructor ...
  LocalRegularGridParamsT(const MyType &defaultParams,
        const ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {
    setDefault(defaultParams);
    parse(myIO::BinaryReaderT<>::nullReader());
  }
  */
  /*
  template <class R>
  LocalRegularGridParamsT(const MyType &defaultParams, R *reader,
        const ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {
    setDefault(defaultParams);
    parse(reader);
  }

  LocalRegularGridParamsT(const ParamsParser *parser):
    manager(parser, glb::console)
  {
    setDefault();
    parse(myIO::BinaryReaderT<>::nullReader());
  }

  template <class R>
  LocalRegularGridParamsT(R *reader, const ParamsParser *parser=glb::dummyPParser):
    manager(parser, glb::console)
  {
    setDefault();
    parse(reader);
    }
*/

  ~LocalRegularGridParamsT()
  {
  }

  void setDefault()
  {
    for (int i = 0; i < NDIM; i++)
    {
      x0[i] = 0;
      delta[i] = 1.0;
      resolution[i] = 64;
      scale[i] = ScaleTypeV::LINEAR;
      valLocation[i] = ValLocationTypeV::CELL;
      lowMargin[i] = 0;
      highMargin[i] = 0;
      parentResolution[i] = resolution[i];

      position[i] = 0;
      parentX0[i] = x0[i];
      parentDelta[i] = delta[i];
    }
    nFields = 1;
    haveParentGrid = 0;
    parentNDim = NDIM;
    minElementsCount = 0;
  }

  template <class W>
  void write(W *writer) const
  {
    writer->writeHeader(classHeader(), classVersion());
    writer->write(x0, NDIM);
    writer->write(delta, NDIM);
    writer->write(resolution, NDIM);
    for (int i = 0; i < NDIM; ++i)
    {
      int k;
      k = static_cast<int>(scale[i]);
      writer->write(&k);
      k = static_cast<int>(valLocation[i]);
      writer->write(&k);
    }
    writer->write(lowMargin, NDIM);
    writer->write(highMargin, NDIM);
    writer->write(&nFields);
    writer->write(&haveParentGrid);
    writer->write(&parentNDim);
    writer->write(parentResolution, NDIM);
    writer->write(&minElementsCount);

    writer->write(position, NDIM);
    writer->write(parentX0, NDIM);
    writer->write(parentDelta, NDIM);
  }

  template <class R>
  void read(R *reader)
  {
    float version;
    R::template checkHeaderAndReport<LOG_ERROR, LOG_WARNING, MyType>(glb::console, reader, version, true);

    reader->read(x0, NDIM);
    reader->read(delta, NDIM);
    reader->read(resolution, NDIM);
    for (int i = 0; i < NDIM; ++i)
    {
      int k;
      reader->read(&k);
      scale[i] = static_cast<ScaleType>(k);
      reader->read(&k);
      valLocation[i] = static_cast<ValLocationType>(k);
    }
    reader->read(lowMargin, NDIM);
    reader->read(highMargin, NDIM);
    reader->read(&nFields);
    reader->read(&haveParentGrid);
    reader->read(&parentNDim);
    reader->read(parentResolution, NDIM);
    reader->read(&minElementsCount);

    reader->read(position, NDIM);
    reader->read(parentX0, NDIM);
    reader->read(parentDelta, NDIM);
  }
  /*
  template <class W>
  void writeManaged(W *writer) const
  {
    writer->writeHeader(classHeader(),classVersion());
    manager.write(writer);
  }

  template <class LT>
  void report() const
  {
    manager.report<LT>();
  }

  void setDefault(const MyType &def)
  {
    for (int i=0;i<NDIM;i++)
      {
  x0[i]=def.x0[i];
  delta[i]=def.delta[i];
  resolution[i]=def.resolution[i];
  scale[i]=def.scale[i];
  valLocation[i]=def.valLocation[i];
  lowMargin[i]=def.lowMargin[i];
  highMargin[i]=def.highMargin[i];
  parentResolution[i]=def.parentResolution[i];

  position[i]=def.position[i];
  parentX0[i]=def.parentX0[i];
  parentDelta[i]=def.parentDelta[i];
      }
    nFields=def.nFields;
    haveParentGrid=def.haveParentGrid;
    parentNDim=def.parentNDim;
    minElementsCount=def.minElementsCount;
  }
  */
  /*
  template <class BR>
  void read(BR *reader)
  {
    return parse(reader,false);
  }
  */
  template <class BR, class PM>
  void parse(BR *reader, PM &manager)
  {
    // float version;
    // BR::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
    //   (glb::console,reader,version,true);

    for (int i = 0; i < NDIM; i++)
    {
      x0[i] = manager.get("x0", parserCategory(), x0[i], i,
                          reader, PM::FILE_FIRST);

      delta[i] = manager.get("delta", parserCategory(), delta[i], i,
                             reader, PM::FILE_FIRST);

      resolution[i] = manager.get("resolution", parserCategory(), resolution[i], i,
                                  reader, PM::FILE_FIRST);

      lowMargin[i] = manager.get("lowMargin", parserCategory(), lowMargin[i], i,
                                 reader, PM::FILE_FIRST);

      highMargin[i] = manager.get("highMargin", parserCategory(), highMargin[i], i,
                                  reader, PM::FILE_FIRST);

      parentResolution[i] = manager.get("parentResolution", parserCategory(), parentResolution[i], i,
                                        reader, PM::FILE_FIRST);

      parentX0[i] = manager.get("parentX0", parserCategory(), parentX0[i], i,
                                reader, PM::FILE_FIRST);

      parentDelta[i] = manager.get("parentDelta", parserCategory(), parentDelta[i], i,
                                   reader, PM::FILE_FIRST);

      position[i] = manager.get("position", parserCategory(), position[i], i,
                                reader, PM::FILE_FIRST);

      std::string scaleStr = manager.template get<std::string>("scale", parserCategory(),
                                                               typename Scale::ScaleTypeSelect().getString(scale[i], true), i,
                                                               reader, PM::FILE_FIRST);
      scale[i] = typename Scale::ScaleTypeSelect().getVal(scaleStr, true);

      std::string valLocationStr = manager.template get<std::string>("valLocation", parserCategory(),
                                                                     ValLocationTypeSelect().getString(valLocation[i], true), i,
                                                                     reader, PM::FILE_FIRST);
      valLocation[i] = ValLocationTypeSelect().getVal(valLocationStr, true);
    }

    haveParentGrid = manager.get("haveParentGrid", parserCategory(), haveParentGrid,
                                 reader, PM::FILE_FIRST);

    parentNDim = manager.get("parentNDim", parserCategory(), parentNDim,
                             reader, PM::FILE_FIRST);

    minElementsCount = manager.get("minElementsCount", parserCategory(), minElementsCount,
                                   reader, PM::FILE_FIRST);
  }

  template <class PP>
  void parse(PP &parser)
  {
    // float version;
    // BR::template checkHeaderAndReport<LOG_ERROR,LOG_WARNING,MyType>
    //   (glb::console,reader,version,true);

    for (int i = 0; i < NDIM; i++)
    {
      x0[i] = parser.get("x0", parserCategory(), x0[i], i);

      delta[i] = parser.get("delta", parserCategory(), delta[i], i);

      resolution[i] = parser.get("resolution", parserCategory(), resolution[i], i);

      lowMargin[i] = parser.get("lowMargin", parserCategory(), lowMargin[i], i);

      highMargin[i] = parser.get("highMargin", parserCategory(), highMargin[i], i);

      parentResolution[i] = parser.get("parentResolution", parserCategory(), parentResolution[i], i);

      parentX0[i] = parser.get("parentX0", parserCategory(), parentX0[i], i);

      parentDelta[i] = parser.get("parentDelta", parserCategory(), parentDelta[i], i);

      position[i] = parser.get("position", parserCategory(), position[i], i);

      std::string scaleStr = parser.template get<std::string>("scale", parserCategory(),
                                                              typename Scale::ScaleTypeSelect().getString(scale[i], true), i);
      scale[i] = typename Scale::ScaleTypeSelect().getVal(scaleStr, true);

      std::string valLocationStr = parser.template get<std::string>("valLocation", parserCategory(),
                                                                    ValLocationTypeSelect().getString(valLocation[i], true), i);
      valLocation[i] = ValLocationTypeSelect().getVal(valLocationStr, true);
    }

    haveParentGrid = parser.get("haveParentGrid", parserCategory(), haveParentGrid);

    parentNDim = parser.get("parentNDim", parserCategory(), parentNDim);

    minElementsCount = parser.get("minElementsCount", parserCategory(), minElementsCount);
  }

  /*
  template <class subParamT>
  void getSubGridParams(subParamT &sp, int which[subParamT::NDIM], long minElCount=0) const
  {
    int i;
    sp.nFields=nFields;
    for (i=0;i<subParamT::NDIM;i++)
      {
  int id = which[i];
  sp.x0[i]          = x0[id];
  sp.delta[i]       = delta[id];
  sp.resolution[i]  = resolution[id];
  sp.scale[i]       = scale[id];
  sp.valLocation[i] = valLocation[id];
  sp.lowMargin[i]   = lowMargin[id];
  sp.highMargin[i]  = highMargin[id];
      }

    sp.nFields=nFields;
    sp.haveParentGrid=1;
    sp.parentNDim=NDIM;
    sp.minElementsCount=minElCount;
    for (i=0;i<subParamT::NDIM;i++)
      {
  sp.position[i]=0;
  sp.parentResolution[i]=sp.resolution[i];
  sp.parentX0[i]=sp.x0[i];
  sp.parentDelta[i]=sp.delta[i];
      }
  }

  ParamsManager manager;
  */
};

#include "../internal/namespace.footer"
#endif

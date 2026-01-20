#ifndef __COLDICE_FILE_DUMPS_HXX__
#define __COLDICE_FILE_DUMPS_HXX__

#include <vector>
#include <map>

#include <dice/dice_globals.hxx>
#include <dice/tools/helpers/eventManager.hxx>

class FileDumps
{
public:
  static std::string parserCategory() { return "dump"; }

  enum Type
  {
    Mesh,
    FatMesh,
    Caustics,
    Lines,
    Subsets,
    Density,
    Amr,
    Potential,
    RadialGridDensity,
    RadialMeshDensity
  };

  FileDumps()
  {
    // Version is the version of VlasovPoisson !
    types.push_back(Info("mesh", Mesh, "mesh", 0.15));
    types.push_back(Info("fatMesh", FatMesh, "mesh with ghost/shadows and tracers info", 0.23));
    types.push_back(Info("caustics", Caustics, "caustics surface sub-mesh", 0.15));
    types.push_back(Info("lines", Lines, "lagrangian 1D lines sub-mesh", 0.15));
    types.push_back(Info("subsets", Subsets, "mesh subsets", 0.15));
    types.push_back(Info("density", Density, "density grid", 0.15));
    types.push_back(Info("amr", Amr, "AMR grid", 0.15));
    types.push_back(Info("potential", Potential, "potential grid", 0.15));

    types.push_back(Info("radialGridDensity", RadialGridDensity,
                         "radial density profile from grid", 0.16));
    types.push_back(Info("radialMeshDensity", RadialMeshDensity,
                         "radial density profile from mesh", 0.16));
  }

  template <class PM, class R>
  void parseFromManager(PM &paramsManager, R *reader,
                        int classVersion, float serializedVersion)
  {
    typename PM::Parser *parser = paramsManager.getParser();
    int reset = parser->get("reset", parserCategory(), 0,
                            "Set this parameters to reset any default file dumps (i.e. reset dumps programmed from restart file)");

    if (reset)
      doParse<PM, PM::IGNORE_FILE>(paramsManager, reader, serializedVersion);
    else
      doParse<PM, PM::PARSER_FIRST>(paramsManager, reader, serializedVersion);

    setTypesMap();
  }

  bool checkEvent(Type e, bool update = false)
  {
    // some versions of icpc do not support operator[] correctly for maps ...
    const Info &info = *typesMap.at(static_cast<int>(e));

    if ((info.dumpEvery >= 1) && (step_ - info.lastDumpStep >= info.dumpEvery))
    {
      if (update)
        info.lastDumpStep = step_;
      return true;
    }

    if (info.eventManager.checkEventTimeFrame(time_, time_ + dt_))
    {
      return true;
    }

    return false;
  }

  bool recheckEvent(Type e)
  {
    // some versions of icpc do not support operator[] correctly for maps ...
    const Info &info = *typesMap.at(static_cast<int>(e));

    if ((info.dumpEvery >= 1) && ((step_ - info.lastDumpStep >= info.dumpEvery) ||
                                  (step_ == info.lastDumpStep)))
    {
      return true;
    }

    if (info.eventManager.checkEventTimeFrame(time_, time_ + dt_))
    {
      return true;
    }

    return false;
  }

  template <class W>
  void write(W *writer)
  {
    typename std::vector<Info>::size_type sz = types.size();
    writer->write(&sz);
    for (long i = 0; i < types.size(); ++i)
    {
      writer->write(&types[i].name);
      writer->write(&types[i].id);
      writer->write(&types[i].lastDumpStep);
    }
  }

  template <class R>
  void read(R *reader)
  {
    typename std::vector<Info>::size_type sz;
    reader->read(&sz);
    for (long i = 0; i < sz; ++i)
    {
      std::string name;
      int id;
      reader->read(&name);
      reader->read(&id);
      auto it = typesMap.find(id);
      if ((it != typesMap.end()) && (name == (*it).second->name))
        reader->read(&((*it).second->lastDumpStep));
      else
      {
        int dummy;
        reader->read(&dummy);
      }
    }
  }

  void setOutputDir(const std::string &oDir)
  {
    fileNameMaker.setOutputDir(oDir);
  }

  std::string getOutputDir()
  {
    return fileNameMaker.getOutputDir();
  }

  void updateState(int step, double time, double dt, dice::MpiCommunication *mpiCom = 0)
  {
    step_ = step;
    time_ = time;
    dt_ = dt;
    fileNameMaker.update(step, mpiCom);
  }

  std::string getGlobalFName(Type e)
  {
    const Info &info = *typesMap[static_cast<int>(e)];
    return fileNameMaker.getGlobal(info.name.c_str());
  }

  std::string getLocalFName(Type e)
  {
    const Info &info = *typesMap[static_cast<int>(e)];
    return fileNameMaker.getLocal(info.name.c_str());
  }

  std::string getLocalFNameFormat(Type e)
  {
    const Info &info = *typesMap[static_cast<int>(e)];
    return fileNameMaker.getLocalFormat(info.name.c_str());
  }

private:
  struct Info
  {
    Info(const std::string &n, int i, const std::string &d, float introducedAt) : name(n), id(i), description(d),
                                                                                  lastDumpStep(-1), dumpEvery(0),
                                                                                  introducedAtVersion(introducedAt)
    {
    }

    std::string name;
    int id;
    std::string description;

    dice::EventManager eventManager;
    mutable int lastDumpStep;
    int dumpEvery;
    float introducedAtVersion;
  };

  std::vector<Info> types;
  std::map<int, Info *> typesMap;

  int step_;
  double time_;
  double dt_;

  void setTypesMap()
  {
    typesMap.clear();
    for (unsigned long i = 0; i < types.size(); ++i)
      typesMap.insert(std::make_pair(types[i].id, &types[i]));
  }

  template <class PM, typename PM::Priority PM_PRIORITY, class R>
  void doParse(PM &paramsManager, R *reader, float serializedVersion)
  {
    for (long i = 0; i < types.size(); ++i)
    {
      Info &info = types[i];
      char description[1024];
      char name[255];

      sprintf(description,
              "Sets %s file dumps time interval (fmt: 'start:stop:every' or a filename containg discrete dump times).",
              info.description.c_str());
      sprintf(name, "%sAt", info.name.c_str());
      // printf("Parsing %s (%f>=%f)\n",name,serializedVersion,info.introducedAtVersion);
      std::string at;
      at = paramsManager.get(name, parserCategory(), at, reader, PM_PRIORITY, description,
                             serializedVersion >= info.introducedAtVersion);

      info.eventManager.addEvents(at);

      sprintf(description,
              "Sets %s file dumps time-steps interval (dump every N timesteps, never if N<=0).",
              info.description.c_str());
      sprintf(name, "%sEvery", info.name.c_str());
      int every = 0;
      every = paramsManager.get(name, parserCategory(), every, reader, PM_PRIORITY, description,
                                serializedVersion >= info.introducedAtVersion);

      info.dumpEvery = every;
    }
  }

  class FileNameMaker
  {
  public:
    FileNameMaker()
    {
      strcpy(outputDir, "./");
    }
    /*
    FileNameMaker(const char *oDir, long curStep, dice::MpiCommunication *mpiCom=0)
    {
      update(oDir,curStep,mpiCom);
    }
    */
    FileNameMaker(long curStep, dice::MpiCommunication *mpiCom = 0)
    {
      strcpy(outputDir, "./");
      update(curStep, mpiCom);
    }
    /*
    void update(const char *oDir, long curStep, dice::MpiCommunication *mpiCom=0)
    {
      strcpy(outputDir,oDir);
      update(curStep,mpiCom);
    }
    */
    void update(long curStep, dice::MpiCommunication *mpiCom = 0)
    {
      sprintf(localFNameFormat, "%s_S%6.6ld", "R%5.5d", curStep);

      if ((mpiCom != 0) && (mpiCom->size() > 1))
        sprintf(localFName, "_R%5.5d_S%6.6ld", mpiCom->rank(), curStep);
      else
        sprintf(localFName, "_S%6.6ld", curStep);

      sprintf(globalFName, "_S%6.6ld", curStep);
    }

    const std::string getLocal(const char *name) const
    {
      char buffer[1024];
      sprintf(buffer, "%s%s%s", outputDir, name, localFName);
      return std::string(buffer);
    }

    const std::string getLocalFormat(const char *name) const
    {
      char buffer[1024];
      sprintf(buffer, "%s%s%s", outputDir, name, localFNameFormat);
      return std::string(buffer);
    }

    const std::string getGlobal(const char *name) const
    {
      char buffer[1024];
      sprintf(buffer, "%s%s%s", outputDir, name, globalFName);
      // printf("%s+%s+%s = %s\n",outputDir,name,globalFName,buffer);
      return std::string(buffer);
    }

    void setOutputDir(const std::string &oDir)
    {
      if (oDir.size() == 0)
        strcpy(outputDir, "./");
      else if (oDir[oDir.size() - 1] != '/')
        sprintf(outputDir, "%s/", oDir.c_str());
      else
        strcpy(outputDir, oDir.c_str());
    }

    std::string getOutputDir()
    {
      return std::string(outputDir);
    }

  private:
    char localFName[256];
    char globalFName[256];
    char localFNameFormat[256];
    char outputDir[256];
  } fileNameMaker;
};

#endif

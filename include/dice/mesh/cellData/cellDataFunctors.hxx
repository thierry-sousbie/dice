#ifndef __CELL_DATA_FUNCTORS_HXX__
#define __CELL_DATA_FUNCTORS_HXX__

#include "../../mesh/cellData/cellDataFunctors_interface.hxx"

#include "../../internal/namespace.header"

// Generic definitions
template <class M>
class CellDataFunctorsT
{
public:
  typedef M Mesh;
  typedef CellDataFunctorsT<M> MyType;
  typedef typename M::Simplex Simplex;
  typedef typename M::Vertex Vertex;
  typedef cellDataFunctors::InterfaceT<M, Vertex> VertexFunctor;
  typedef cellDataFunctors::InterfaceT<M, Simplex> SimplexFunctor;

  CellDataFunctorsT(const Mesh *m) : mesh(m)
  {
    setDataElements();
  }
  ~CellDataFunctorsT()
  {
    for (unsigned long i = 0; i < simplexFunctor.size(); ++i)
      delete simplexFunctor[i];
    for (unsigned long i = 0; i < vertexFunctor.size(); ++i)
      delete vertexFunctor[i];
    /*
      for (typename SimplexMap::iterator it=simplexFunctor.begin();
      it!=simplexFunctor.end();++it)
      delete it->second;
      for (typename VertexMap::iterator it=vertexFunctor.begin();
      it!=vertexFunctor.end();++it)
      delete it->second;
    */
  }

private:
  // FIXME: copy not implemented yet (ever)...
  CellDataFunctorsT() {}
  CellDataFunctorsT(const MyType &init) {}
  MyType &operator=(const MyType &other) {}

  void setDataElements()
  {
    // Create a functor for data associated to simplices ...
    Simplex::createAllDataElementFunctorPtr(mesh, std::back_inserter(simplexFunctor));
    for (unsigned long i = 0; i < simplexFunctor.size(); ++i)
      simplexFunctorMap.insert(std::make_pair(simplexFunctor[i]->getName(), simplexFunctor[i]));

    // and vertices
    Vertex::createAllDataElementFunctorPtr(mesh, std::back_inserter(vertexFunctor));
    for (unsigned long i = 0; i < vertexFunctor.size(); ++i)
      vertexFunctorMap.insert(std::make_pair(vertexFunctor[i]->getName(), vertexFunctor[i]));
  }

public:
  template <template <class TM> class T>
  bool insert(int flags = cellDataFunctors::F_NO_FLAG, bool replaceIfExists = false)
  {
    return add__<T>(flags, replaceIfExists, cell_traits<typename T<M>::Cell>());
  }

  VertexFunctor *getVertexFunctor(unsigned int i)
  {
    if (i < vertexFunctor.size())
      return vertexFunctor[i];
    return NULL;
  }

  SimplexFunctor *getSimplexFunctor(unsigned int i)
  {
    if (i < simplexFunctor.size())
      return simplexFunctor[i];
    return NULL;
  }

  VertexFunctor *getVertexFunctor(const std::string &str)
  {
    typename std::map<std::string, VertexFunctor *>::iterator it =
        vertexFunctorMap.find(str);
    if (it == vertexFunctorMap.end())
      return NULL;
    else
      return it->second;
  }

  SimplexFunctor *getSimplexFunctor(const std::string &str)
  {
    typename SimplexMap::iterator it =
        simplexFunctorMap.find(str);
    if (it == simplexFunctorMap.end())
      return NULL;
    else
      return it->second;
  }

  int getVertexFunctorCount() const
  {
    return vertexFunctor.size();
  }

  int getSimplexFunctorCount() const
  {
    return simplexFunctor.size();
  }

  int getTotalCount() const
  {
    return simplexFunctor.size() + vertexFunctor.size();
  }

private:
  template <class CT>
  struct cell_traits
  {
    typedef CT Cell;
  };

  const Mesh *mesh;

  typedef std::map<std::string, SimplexFunctor *> SimplexMap;
  typedef std::map<std::string, VertexFunctor *> VertexMap;
  typedef std::vector<SimplexFunctor *> SimplexVec;
  typedef std::vector<VertexFunctor *> VertexVec;

  SimplexMap simplexFunctorMap;
  VertexMap vertexFunctorMap;

  SimplexVec simplexFunctor;
  VertexVec vertexFunctor;

  template <template <class TM> class T>
  bool add__(int flags, bool replaceIfExists, cell_traits<Simplex>)
  {
    SimplexFunctor *data = new T<M>(mesh, flags);
    auto it = simplexFunctorMap.find(data->getName());
    if (it != simplexFunctorMap.end())
    {
      if (replaceIfExists)
        std::swap(it->second, data);

      delete data;
      return replaceIfExists;
    }
    simplexFunctorMap.insert(std::make_pair(data->getName(), data));
    simplexFunctor.push_back(data);
    return true;
  }
  template <template <class TM> class T>
  bool add__(int flags, bool replaceIfExists, cell_traits<Vertex>)
  {
    VertexFunctor *data = new T<M>(mesh, flags);
    auto it = vertexFunctorMap.find(data->getName());
    if (it != vertexFunctorMap.end())
    {
      if (replaceIfExists)
        std::swap(it->second, data);

      delete data;
      return replaceIfExists;
    }
    vertexFunctorMap.insert(std::make_pair(data->getName(), data));
    vertexFunctor.push_back(data);
    return true;
  }
};

#include "../../internal/namespace.footer"
#endif

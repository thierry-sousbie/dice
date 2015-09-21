#ifndef __IMPLICIT_TESSELATION_HXX__
#define __IMPLICIT_TESSELATION_HXX__

#include "../../tools/types/cell.hxx"

/**
 * @file 
 * @brief  Defines an interface for implicit tesselations used to initialize a mesh
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class ImplicitTesselationT
 * \brief Defines an interface for implicit tesselations used to initialize a mesh. An 
 * implicit tesselation should as much as possible not store explicitely the topology 
 * and geometry of the mesh so that its memory usage is constant (i.e. it does not 
 * depend on the size of the tesselation). Implicit tesselations can be used to 
 * initialize explicit meshes that spans numerous MPI nodes without having to store 
 * the whole tesselation on each MPI node at any time.
 * \tparam ND the number of dimensions for the tesselation
 * \tparam NDW the number of the embedding space
 * \tparam CT a struct that stores cells identity. See e.g. CellT
 */
template <int ND, int NDW, class CT = CellT<> >
class ImplicitTesselationT
{
public:
  typedef CT Cell;  
  static const int NDIM = ND;
  static const int NDIM_W = NDW;
  
  virtual ~ImplicitTesselationT() {}
  
  virtual int getNeighbors(Cell cell, std::vector<Cell> &out, int type=-1) const=0;
  virtual int getNeighbors(Cell cell, Cell *out, int type=-1) const=0;
  virtual int getVertices(Cell cell, Cell *out) const=0;
  virtual void getPosition(Cell cell,double *out) const=0;
  virtual int getComponentIndex(Cell cell) {return 0;}
  virtual void getX0(double *out) const=0;
  virtual void getDelta(double *out) const=0;
  virtual void getResolution(int *out) const=0;
  virtual std::vector<unsigned long> getNCells() const=0;
  virtual std::vector<unsigned long> getNCellsPerVoxel() const=0;
  virtual int getNDims() const=0;
  virtual int getNDimsW() const=0;  
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

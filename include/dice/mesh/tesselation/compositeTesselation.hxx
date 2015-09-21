#ifndef __COMPOSITE_TESSELATION_HXX__
#define __COMPOSITE_TESSELATION_HXX__

/**
 * @file 
 * @brief  A class to define the implicit tesselation resulting from the composition of 
 * other tesselations
 * @author Thierry Sousbie
 */
#include <vector>
#include "./implicitTesselation.hxx"

#include "../../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class CompositeTesselationT
 * \tparam ND the number of dimensions for the tesselation
 * \tparam NDW the number of the embedding space
 * \tparam CT a struct that stores cells identity. See e.g. CellT
 */
template <int ND, int NDW ,class CT>
class CompositeTesselationT : public ImplicitTesselationT<ND,NDW,CT>
{
public:
  typedef CT Cell;  
  static const int NDIM = ND;
  static const int NDIM_W = NDW;
  
  typedef ImplicitTesselationT<ND,NDW,CT> Interface;
  
  template <typename InputIterator>
  CompositeTesselationT(InputIterator x0_, InputIterator delta_)
  {
    for (int i=0; i<NDIM+1; ++i)
      nCellsCum[i].push_back(0);
    setBBox(x0_,delta_);
  }

  CompositeTesselationT()
  {
    for (int i=0; i<NDIM+1; ++i)
      nCellsCum[i].push_back(0);

    std::fill(x0,x0+NDIM,-1);
    std::fill(delta,delta+NDIM,2);
  }

  template <typename InputIterator>
  void setBBox(InputIterator x0_, InputIterator delta_)
  {
    std::copy(x0_,x0_+NDIM,x0);
    std::copy(delta_,delta_+NDIM,delta);
  }

  void add(Interface *it)
  {
    tesselations.push_back(it);
    auto newNCells=it->getNCells();
    for (int i=0;i<NDIM+1;++i)
      nCellsCum[i].push_back(newNCells[i]+nCellsCum[i].back());
  }

  long getComponentId(Cell cell) const
  {
    return getTesselationId(cell);
  }
  
  /* interface implementation */
  int getNeighbors(Cell cell, std::vector<Cell> &out, int type=-1) const
  {
    long id=getTesselationId(cell);
    Interface *t=tesselations[id];
    int result;
    if (id>0)
      {
	cell.setId(cell.id()-nCellsCum[cell.type()][id]);
	result=t->getNeighbors(cell,out,type);
	for (int i=0;i<result;++i) 
	  out[i].setId(out[i].id()+nCellsCum[out[i].type()][id]);
      }
    else result=t->getNeighbors(cell,out,type);
    return result;
  }
  
  int getNeighbors(Cell cell, Cell *out, int type=-1) const
  {
    long id=getTesselationId(cell);
    Interface *t=tesselations[id];
    int result;
    if (id>0)
      {
	cell.setId(cell.id()-nCellsCum[cell.type()][id]);
	result=t->getNeighbors(cell,out,type);
	for (int i=0;i<result;++i) 
	  out[i].setId(out[i].id()+nCellsCum[out[i].type()][id]);
      }
    else result=t->getNeighbors(cell,out,type);
    return result;
  }
  
  int getVertices(Cell cell, Cell *out) const
  {
    long id=getTesselationId(cell);
    Interface *t=tesselations[id];
    int result;
    if (id>0)
      {
	cell.setId(cell.id()-nCellsCum[cell.type()][id]);
	result=t->getVertices(cell,out);
	for (int i=0;i<result;++i) 
	  out[i].setId(out[i].id()+nCellsCum[out[i].type()][id]);
      }
    else result=t->getVertices(cell,out);
    return result;
  }
  
  void getPosition(Cell cell,double *out) const
  {
    long id=getTesselationId(cell);
    Interface *t=tesselations[id];
    if (id>0)
      {
	cell.setId(cell.id()-nCellsCum[cell.type()][id]);
	t->getPosition(cell,out);
      }
    else t->getPosition(cell,out);
  }
  
  void getX0(double *out) const
  {
    std::copy(x0,x0+NDIM,out);
  }
  
  void getDelta(double *out) const
  {
    std::copy(delta,delta+NDIM,out);
  }
  
  void getResolution(int *out) const
  {
    if (tesselations.size())
      {
	for (int i=0;i<NDIM;++i)
	  {
	    (*out)=0;
	    ++out;
	  }
      }
    else tesselations.back()->getResolution(out);
  }
  
  std::vector<unsigned long> getNCells() const
  {
    std::vector<unsigned long> result(NDIM+1);
    for (int i=0;i<NDIM+1;++i)
      result[i]=nCellsCum[i].back();
    return result;
  }
  
  std::vector<unsigned long> getNCellsPerVoxel() const
  {
    return getNCells();
  }

  int getNDims() const
  {
    return NDIM;
  }

  int getNDimsW() const
  {
    return NDIM_W;
  }
  /* end of interface */

private:
  double x0[NDIM];
  double delta[NDIM];
  std::vector<Interface *> tesselations;
  std::vector< unsigned long > nCellsCum[NDIM+1];     

  long getTesselationId(Cell cell) const
  {
    typename Cell::Type type=cell.type();
    typename Cell::Id id=cell.id();
    
    if (id>=nCellsCum[type].back()) return -1;
    auto it=std::upper_bound(nCellsCum[type].begin(),nCellsCum[type].end(),id);
    return std::distance(nCellsCum[type].begin(),it)-1;
  }
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

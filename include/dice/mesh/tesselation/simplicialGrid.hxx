#ifndef __SIMPLICIAL_GRID_HXX__
#define __SIMPLICIAL_GRID_HXX__

/**
 * @file 
 * @brief  A class to define the implicit tesselation of a grid into simplices
 * @author Thierry Sousbie
 */


#include <stdio.h>
#include <assert.h>

#include <set>
#include <map>
#include <vector>

#include "../../dice_globals.hxx"
#include "./tesselationType.hxx"
#include "./implicitTesselation.hxx"
#include "./simplicialGridParams.hxx"

#include "../../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class SimplicialGridT
 * A class to define the implicit tesselation of a grid into simplices.
 * \tparam ND the number of dimensions for the tesselation
 * \tparam NDW the number of the embedding space
 * \tparam CT a struct that stores cells identity. See e.g. CellT
 * \tparam P whether the grid is periodic (true) or not (false)
 * \warning no access to intermediate cells (i.e 0<dim<NDIM) for alternate
 * tesselations.
 */
template <int ND, int NDW ,class CT ,bool P>
class SimplicialGridT : public ImplicitTesselationT<ND,NDW,CT>
{
public:
  typedef CT Cell;  
  typedef typename CT::Type Type;
  typedef typename CT::Id Id;

  typedef ImplicitTesselationT<ND,NDW,CT> Interface;
  typedef SimplicialGridParamsT<ND,NDW,CT> Params;
 
  static const int NDIM_MAX = 6;
  static const int NDIM = ND;
  static const int NDIM_W = NDW;
  static const bool PERIODIC=P;
  
  SimplicialGridT()
    {
      initialized=false;
      createVertices();
    }
  /*
  template <typename inputIteratorF, typename inputIteratorI>
  SimplicialGridT(inputIteratorF _x0,inputIteratorF _delta, inputIteratorI _size, 
		  TesselationType _t = TesselationTypeV::ANY, bool forceAlternate_=false)
  {
    initialized=false;
    createVertices();
    build(_x0,_delta,_size,_t,forceAlternate_);
  }
*/

  SimplicialGridT(const Params &p)
  {
    initialized=false;
    createVertices();
    build(p);
  }
  
  ~SimplicialGridT()
  {
    
  } 
  
  /*
  template <typename inputIteratorF, typename inputIteratorI>
  void build(inputIteratorF _x0,inputIteratorF _delta, inputIteratorI resolution, 
	     TesselationType t = TesselationTypeV::ANY, bool forceAlternate_=false)
  */
  void build(const Params &p, bool forceAlternate_=false)
  {
    std::copy(p.x0,p.x0+NDIM,x0);
    std::copy(p.delta,p.delta+NDIM,delta);
    std::copy(p.resolution,p.resolution+NDIM,dims);
    for (int i=0;i<NDIM;i++) xmax[i]=x0[i]+delta[i];
    for (int i=0;i<NDIM;i++) deltaHalf[i]=delta[i]/2;
    for (int i=0;i<NDIM;i++) stepSize[i]=delta[i]/dims[i];
    forceAlternate=forceAlternate_;
    build(p.t);
    initialized=true;   
  }  
  
  template <class OutputIterator>
  int getNeighbors(Cell cell, OutputIterator out, int type=-1) const
  {
    ReferenceCell ref;
    cellToReferenceCell(cell,ref);    
    long config=Coords2BoundaryConfig(ref.coords);
    if (type<0) type=cell.type();
    if (config<0)
      {
	cellToReferenceCell(cell,ref,true);    
	std::vector<long> result;
	long count=getNeighborsIndex(ref,type,std::back_inserter(result));
	for (long i=0;i<count;i++)
	  {
	    *out = Cell(type,result[i]);
	    ++out;
	  }
	return count;
      }
    else
      {
	const std::vector<long> &shift = shiftTable[config][cell.type()][type][ref.simplex->index];
	ReferenceCell ref2(const_cast<Simplex*>(&allCells[type][0]),ref.coords);
	long refIndex=referenceCellToIndex(ref2);
	for (unsigned long i=0;i<shift.size();i++)
	  {	  
	    *out = Cell(type,refIndex + shift[i]);
	    ++out;
	  }
	return shift.size();
      }
    return 0;
  }

  template <class OutputIterator>
  int getVertices(Cell cell, OutputIterator out) const 
  { 
    return getNeighbors(cell,out,0);
  }
  
  template <class OutputIterator>
  void getPosition(Cell cell,OutputIterator out) const
  {
    if (cell.type()==0)
      {
	ReferenceCell ref;
	cellToReferenceCell(cell,ref);
	for (int i=0;i<NDIM;i++) 
	  {
	    (*out)=ref.coords[i]*stepSize[i]+x0[i];
	    ++out;
	  }
      }
    else
      {
	std::vector<double> tmp(NDIM,0);
	Cell vertices[NDIM+1];
	int n = getVertices(cell,vertices);
	double nInv=1.0/n;
	for (int i=0;i<n;i++)
	  {
	    ReferenceCell ref;
	    cellToReferenceCell(vertices[i],ref);
	    for (int j=0;j<NDIM;j++)
	      tmp[j]+=ref.coords[j]*stepSize[j]+x0[j];
	  }
	for (int j=0;j<NDIM;j++)
	  {
	    *out = tmp[j]*nInv;
	    ++out;
	  }	
      }
  }

  template <class OutputIterator>
  void getX0(OutputIterator out) const
  {
    for (int i=0;i<NDIM;++i,++out)
      (*out)=x0[i];
  }

  template <class OutputIterator>
  void getDelta(OutputIterator out) const
  {
    for (int i=0;i<NDIM;++i,++out)
      (*out)=delta[i];
  }

  template <class OutputIterator>
  void getResolution(OutputIterator out) const
  {
    for (int i=0;i<NDIM;++i,++out)
      (*out)=dims[i];
  }

  /* */
  // this is to fulfill the implicitTesselation interface requirements
  int getNeighbors(Cell cell, Cell *out, int type=-1) const
  {return getNeighbors<Cell*>(cell,out,type);}
  int getNeighbors(Cell cell, std::vector<Cell> &out, int type=-1) const
  {
    out.clear();
    return getNeighbors(cell,std::back_inserter(out),type);
  }
  int getVertices(Cell cell, Cell *out) const
  {return getVertices<Cell*>(cell,out);}
  void getPosition(Cell cell,double *out) const
  {return getPosition<double*>(cell,out);}
  // void getPosition(Cell cell,float *out) const
  // {return getPosition<float*>(cell,out);}
  void getX0(double *out) const
  {getX0<double*>(out);}
  void getDelta(double *out) const
  {getDelta<double*>(out);}
  void getResolution(int *out) const
  {getResolution<int*>(out);}  
  /* */

  std::vector<unsigned long> getNCells() const
  {
    return std::vector<unsigned long>(nCells,nCells+NDIM+1);
    //if (perVox) return nCellsPerVox;
    //else return nCells;
  }

  std::vector<unsigned long> getNCellsPerVoxel() const
  {
    return std::vector<unsigned long>(nCellsPerVoxel,nCellsPerVoxel+NDIM+1);
    //return nCellsPerVox;
  }

  int getNDims() const
  {
    return NDIM;
  }

  int getNDimsW() const
  {
    return NDIM_W;
  }

  void addSimplex(int v0=0, int v1=0, int v2=0, int v3=0, int v4=0, int v5=0, int v6=0)
  {
    int v[7]={v0,v1,v2,v3,v4,v5,v6};
    addSimplex(v);
  }

  template <class InputIterator>
  void addSimplex(InputIterator verticesId)
  {
    int v[NDIM+1];
    for (int i=0;i<NDIM+1;++i)
      {
	v[i]=(*verticesId);
	++verticesId;
      }    
    addSimplex(v);
  }


protected:
  double x0[NDIM];
  double delta[NDIM];
  double deltaHalf[NDIM];
  double stepSize[NDIM];
  double xmax[NDIM];
  long dims[NDIM];
  bool initialized;
  bool isAlternate;
  bool forceAlternate;

  template <int VND, int VNDW>
  struct VertexT
  { 
    static const int NDIM = VND;
    static const int NDIM_W = VNDW;
 
    VertexT(long index_)
    {
      for (int j=0;j<NDIM_W;j++)
	{
	  if (index_&(1<<j)) 
	    coords[j]=true;
	  else 
	    coords[j]=false;
	}
      index=index_;
    }   

    bool coords[NDW];
    int index;

    template <class T>
    long getSwappedIndex(T swap[NDIM]) const
    {
      long result=index;
      for (int i=0;i<NDIM;i++)
	{
	  if (swap[i])
	    {
	      if (index&(1<<i)) 
		result &= ~(1ul<<i);
	      else
		result |= (1ul<<i);
	    }
	}
      return result;
    }
  };

  template <class S>
  class SimplexConfigT
  {
  public:
    typedef SimplexConfigT<S> MyType;
    typedef S Simplex;
    static const int NDIM = Simplex::NDIM;
    static const int NDIM_W = Simplex::NDIM_W;

    SimplexConfigT(Simplex *s)
    {
      simplex=s;
      for (int i=0;i<NDIM;i++)
	shift[i]=0;
    }

    template <typename T>
    SimplexConfigT(Simplex *s,const T shift_[NDIM])
    {
      simplex=s;
      for (int i=0;i<NDIM;i++)
	shift[i]=shift_[i];      
    }

    long getConfig() const
    {
      long config=0;
      
      for (int i=0;i<NDIM;i++)
	if (shift[i]!=0) 
	  config|=(1l<<i);
	
      return config;
    }

    std::string toString() const
    {
      std::string str;
      str = simplex->toString();
      char tmp[1000];
      sprintf(tmp,"CONFIG(%s [%d,%d,%d])",str.c_str(),shift[0],shift[1],shift[2]);
      return std::string(tmp);
    }

    Simplex *simplex;
    int shift[NDIM];  
  };

  template <int SND, int SNDW>
  class SimplexT
  {
  public:
    typedef SimplexT<ND,NDW> MyType;
    typedef VertexT<ND,NDW> Vertex;
    typedef SimplexConfigT<MyType> Config;
    
    static const int NDIM = SND;
    static const int NDIM_W = SNDW;
    
    template <class InputIterator>
    SimplexT(InputIterator vBegin, InputIterator vEnd, long index_=-1)
    {
      std::set<Vertex*> vSet(vBegin,vEnd);
      vertices.assign(vSet.begin(),vSet.end());
      nVertices=vertices.size();
      if (std::distance(vBegin,vEnd) != nVertices)
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Ill defined simplex (some vertices are identical).\n");
	  exit(-1);
	}
      type=nVertices-1;
      index=index_;
    }
    
    void addNeighbor(MyType *c)
    {
      neighbors[c->type].push_back(c);
    }

    long getRightBoundaryConfig() const
    {
      long result=~(0l);
      for (long i=0;i<nVertices;i++)
	result &= vertices[i]->index;
      return result;
    }
    
    bool contains(const MyType &other) const
    {
      if (other.type>type) return false;
      
      int nFound=0;
      int i=0;int j=0;
      while ((j<other.nVertices)&&(i<nVertices))
	{
	  if (vertices[i]<other.vertices[j]) i++;
	  else if (other.vertices[j]<vertices[i]) j++;
	  else {nFound++;i++;j++;}
	};
      return (nFound==other.nVertices);
    }

    bool isNeighbor(const MyType &other) const
    {
      if (other.type!=type) return false;
      if (type==0) return false;

      int nFound=0;
      int i=0;int j=0;
      while ((j<nVertices)&&(i<nVertices))
	{
	  if (vertices[i]<other.vertices[j]) i++;
	  else if (other.vertices[j]<vertices[i]) j++;
	  else {nFound++;i++;j++;}
	};
      return (nFound==(nVertices-1));
    }

    bool operator<(const MyType &other) const 
    {
      if (other.nVertices!=nVertices) return nVertices<other.nVertices;
      for (int i=0;i<nVertices;i++)
	{
	  if (vertices[i] == other.vertices[i]) continue;
	  return (vertices[i] < other.vertices[i]);
	}
      return false;
    }

    bool operator==(const MyType &other) const
    {
      if (other.nVertices!=nVertices) return false;
      for (int i=0;i<nVertices;i++)
	if (vertices[i] != other.vertices[i]) return false;
      return true;
    }

    std::string neighborsToString() const
    {
      char tmp[10000];
      strcpy(tmp,"");

      for (int i=0;i<=NDIM;i++)
	{
	  for (int j=0;j<neighbors[i].size();j++)
	    {
	      sprintf(tmp,"%s  ->%s.\n",tmp,neighbors[i][j]->toString().c_str());
	    }
	}
      return std::string(tmp);
    }

    std::string toString() const
    {
      char tmp[255];
      if (type==0)
	sprintf(tmp,"%d-cell(%d):{%d}",type,index,
		vertices[0]->index);
      if (type==1)
	sprintf(tmp,"%d-cell(%d):{%d,%d}",type,index,
		vertices[0]->index,
		vertices[1]->index);
      else if (type==2)
	sprintf(tmp,"%d-cell(%d):{%d,%d,%d}",type,index,
		vertices[0]->index,
		vertices[1]->index,
		vertices[2]->index);
      else if (type==3)
	sprintf(tmp,"%d-cell(%d):{%d,%d,%d,%d}",type,index,
		vertices[0]->index,
		vertices[1]->index,
		vertices[2]->index,
		vertices[3]->index);
      return std::string(tmp);
    }

    int nVertices;
    int type;
    int index;
    std::vector<MyType *> neighbors[NDIM+1];
    std::vector<Vertex*> vertices;
  };

  template <class S>
  struct ReferenceCellT
  {
    typedef SimplexConfigT<S> MyType;
    typedef S Simplex;
    static const int NDIM = Simplex::NDIM;
    static const int NDIM_W = Simplex::NDIM_W;    

    ReferenceCellT()
    {

    }

    template <typename T>
    ReferenceCellT(Simplex *s,const T coords_[NDIM])
    {
      simplex=s;
      for (int i=0;i<NDIM;i++) coords[i]=coords_[i];
    }

    Simplex *simplex;
    int coords[NDIM];

    std::string toString() const
    {
      std::string str;
      str = simplex->toString();
      char tmp[1000];
      sprintf(tmp,"REF(%s [%d,%d,%d])",str.c_str(),coords[0],coords[1],coords[2]);
      return std::string(tmp);
    }

    template <class OutputIterator> 
    long getEquivalentConfigs(OutputIterator out, bool isAlternate=false)
    {  
      long nSym=0;
      int symDir[NDIM];
      int symConf[NDIM];

      //glb::console->print<LOG_DEBUG>("looking for equivalent configs for %s:\n",toString().c_str());

      for (int i=0;i<NDIM;i++)
	{
	  int j=1;
	  for (;j<simplex->nVertices;j++)
	    {	    
	      if ((simplex->vertices[j]->index&(1<<i)) != (simplex->vertices[j-1]->index&(1<<i)))
		break;
	    }
	  if (j==simplex->nVertices) 	    
	    {	      	     
	      if (simplex->vertices[0]->index&(1<<i))
		symConf[i]=1;
	      else
		symConf[i]=-1;
	
	      symDir[nSym++]=i;
	    }	 
	  else symConf[i]=0;
	  if ((isAlternate)&(coords[i]&1)) symConf[i]=-symConf[i];	
	}

      //glb::console->print<LOG_DEBUG>("Found %ld symmetries : [%d %d %d]\n",nSym,symConf[0],symConf[1],symConf[2]);
      int disp[NDIM];
      for (int i=0;i<(1<<nSym);i++)
	{
	  memset(disp,0,sizeof(int)*NDIM);
	  for (int j=0;j<nSym;j++)
	    {
	      if (i&(1<<j))
		disp[symDir[j]]=symConf[symDir[j]];
	      else 
		disp[symDir[j]]=0;		
	    }
	  *out = SimplexConfig(simplex,disp);
	  ++out;
	}

      return (1<<nSym);
    }  
    
  };

  typedef SimplexT<NDIM,NDIM_W> Simplex;
  typedef typename Simplex::Vertex Vertex;
  typedef typename Simplex::Config SimplexConfig;

  typedef ReferenceCellT<Simplex> ReferenceCell;
  
  unsigned long nCellsPerVoxel[NDIM+1];
  unsigned long nCells[NDIM+1];
  unsigned long nVoxels;
  
  unsigned long cellStride[NDIM+1][NDIM+1];
  unsigned long int addStride[NDIM+1][NDIM+1];
  unsigned long facesPerVoxel[NDIM+1][1<<NDIM];//for a type, in each face of the cube
  unsigned long facesCellStride[NDIM+1][1<<NDIM][NDIM+1];
  unsigned long voxelStride[NDIM+1];
  //std::vector<Simplex> simplices;

  std::vector<Vertex> vertices;
  std::vector<Simplex> allCells[NDIM+1];
  std::set<Simplex> allCellsSet[NDIM+1];
  /*
  std::vector<Vertex> allVerticesConfigs[1<<NDIM];
  std::vector<Simplex> allCellsConfigs[1<<NDIM][NDIM+1];
  std::vector<Set> allCellsSetConfigs[1<<NDIM][NDIM+1];
  */
  typedef typename std::set<Simplex>::iterator CellSetIt;

  
  void createVertices()
  {
    for (long i=0;i<(1l<<NDIM);i++)
      vertices.push_back(Vertex(i));

    for (unsigned long i=0;i<vertices.size();i++)
      {
	std::vector<Vertex *> vec(1,&vertices[i]);
	Simplex s(vec.begin(),vec.end(),i);
	allCellsSet[0].insert(s);
	allCells[0].push_back(s);
      }
  }

  void addSimplex(int a[NDIM+1])
  {
    std::vector< Vertex *> v;
    for (int i=0;i<NDIM+1;i++)
      v.push_back(&vertices[a[i]]);

    Simplex s(v.begin(),v.end(),allCells[NDIM].size());
    if (addCells(s)==false)
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("There is a duplicated simplex.\n");
	exit(-1);
      }    
  }

  bool addCells(Simplex &s)
  { 
    typedef std::pair<CellSetIt,bool> InsertResult;    
        
    InsertResult result=allCellsSet[s.type].insert(s);
    
    if (result.second == false) 
      {
	//printf("Skipped %s.\n",s.toString().c_str());
	return false;
      }
    //glb::console->print<LOG_DEBUG>("Inserted %s.\n",s.toString().c_str());
    allCells[s.type].push_back(s);
    
    int nvert = s.nVertices;
    int type = s.type;
    if (type<=1) return true;

    Vertex *newCellV[NDIM+1];
    for (long i=0;i<(1l<<nvert);i++)
      {
	int count=0;
	for (int j=0;j<nvert;j++) 
	  {
	    if (i&(1<<j)) 
	      newCellV[count++]=s.vertices[j];
	  }
	if (count!=nvert-1) continue;
	Simplex cell(newCellV,newCellV+count,allCells[type-1].size());
	addCells(cell);
      }

    return true;
  }

  void setNeighbors()
  {
    for (int type1=NDIM;type1>=0;type1--)
      {
	std::vector<Simplex> &cells =  allCells[type1];
	for (unsigned long i=0;i<cells.size();i++)
	  {
	    for (int type2=0;type2<=type1;type2++)
	      {
		std::vector<Simplex> &cells2 =  allCells[type2];
		for (unsigned long j=0;j<cells2.size();j++)
		  {
		    /*
		    if (type1==type2)
		      {
			if (j<=i) {j=i;continue;}
			if (cells[i].isNeighbor(cells2[j]))
			  {
			    cells[i].addNeighbor(&cells2[j]);
			    cells2[j].addNeighbor(&cells[i]);
			  }
		      }
		    */
		    if (type1==type2) continue;
		    else if (cells[i].contains(cells2[j]))
		      {
			cells[i].addNeighbor(&cells2[j]);
			cells2[j].addNeighbor(&cells[i]);
		      }		    
		  }
	      }
	  }
      }

    // add vertices neighbors (two vertices are neighbor if they share and edge)
    std::vector<Simplex> &cells=allCells[1];
    for (unsigned long i=0;i<cells.size();i++)
      {
	cells[i].neighbors[0][0]->addNeighbor(cells[i].neighbors[0][1]);
	cells[i].neighbors[0][1]->addNeighbor(cells[i].neighbors[0][0]);
      }
  }

  void setSimplicesRegular()
  {
    if (NDIM==1)
      {
	addSimplex(0,1);
      }
    else if (NDIM==2)
      {
	addSimplex(0,1,2);
	addSimplex(1,2,3);
      }
    else if (NDIM==3)
      {
	addSimplex(0,7,1,3);
	addSimplex(0,7,3,2);
	addSimplex(0,7,2,6);
	addSimplex(0,7,6,4);
	addSimplex(0,7,4,5);
	addSimplex(0,7,5,1);
      }   
  }

  void setSimplicesMinimal()
  {
    if (NDIM<3) return setSimplicesRegular();

    isAlternate=true;
    if (NDIM==3)
      {
	addSimplex(0,1,2,4);
	addSimplex(1,2,3,7);
	addSimplex(2,4,6,7);
	addSimplex(1,2,4,7);
	addSimplex(1,4,5,7);
      }
    else if (NDIM==4)
      {
	addSimplex(0,1,2,4,8);
	addSimplex(4,8,12,13,14);
	addSimplex(2,8,10,11,14);
	addSimplex(2,4,6,7,14);
	addSimplex(1,8,9,11,13);
	addSimplex(1,4,5,7,13);
	addSimplex(1,2,3,7,11);
	addSimplex(7,11,13,14,15);
	addSimplex(1,2,4,8,14);
	addSimplex(1,4,8,13,14);
	addSimplex(1,2,8,11,14);
	addSimplex(1,2,4,7,14);
	addSimplex(1,8,11,13,14);
	addSimplex(1,4,7,13,14);
	addSimplex(1,2,7,11,14);
	addSimplex(1,7,11,13,14);
      }
    
  }

  bool dimsAreEven()
  {
    for (int i=0;i<NDIM;i++)
      if (dims[i]&1) return false;
    return true;
  }

  TesselationType setSimplices(TesselationType t = TesselationTypeV::ANY)
  {
    if (t==TesselationTypeV::ANY) 
      {
	if (NDIM<4) 
	  {
	    if ((dimsAreEven())||(!PERIODIC)) forceAlternate=true;
	    setSimplicesRegular();
	    return TesselationTypeV::REGULAR;
	  }
	else 
	  {
	    setSimplicesMinimal();
	    return TesselationTypeV::MINIMAL;
	  }
      }
    else if (t==TesselationTypeV::REGULAR) 
      {
	setSimplicesRegular();
	return TesselationTypeV::REGULAR;
      }
    else if (t==TesselationTypeV::ALTERNATE) 
      {
	forceAlternate=true;
	setSimplicesRegular();
	return TesselationTypeV::ALTERNATE; //FIXME: CHECK THIS !
      }
    else if (t==TesselationTypeV::MINIMAL) 
      {
	setSimplicesMinimal();
	return TesselationTypeV::MINIMAL;
      }

    return TesselationTypeV::UNDEFINED;
  }
 
  // if config&(1<<i) is true => mirror along the ith dimension
  Vertex *getNeighborVertexPtr(Vertex *v, long config, bool alternate=false) const
  {
    long index = v->index;
    if (config==0) return const_cast<Vertex *>(&vertices[index]);
    if (alternate) return const_cast<Vertex *>(&vertices[index]);

    long newIndex=index;
    for (int i=0;i<NDIM;i++)
      {
	if (config&(1l<<i))
	  {
	    long val=(1l<<i);
	    if (newIndex&val)
	      newIndex &= (~val);
	    else
	      newIndex |= val;
	  }
      }
    return const_cast<Vertex *>(&vertices[newIndex]);
  }

  // if config&(1<<i) is true => mirror along the ith dimension
  Simplex *getNeighborCellPtr(Simplex *s, long config, bool alternate=false) const
  {
    if (config==0) return s;
    if (alternate) return s;

    Vertex* vert[s->nVertices];
    for (int i=0;i<s->nVertices;i++)
      vert[i]=getNeighborVertexPtr(s->vertices[i],config,alternate);
    Simplex newSimplex(vert,vert+s->nVertices);
    CellSetIt it=allCellsSet[s->type].find(newSimplex);
    
    if (it==allCellsSet[s->type].end())
      {
	//printf("changed %s -> NULL (%s).\n",s->toString().c_str(),newSimplex.toString().c_str());
	return NULL;
      }
    else 
      {
	//printf("changed %s -> %s.\n",s->toString().c_str(),allCells[s->type][it->index].toString().c_str());
	return const_cast<Simplex *>(&allCells[s->type][it->index]);
      }
  }

  bool checkAlternate()
  {
    if (NDIM<=2) return false;  
    std::vector<Simplex> &cells = allCells[NDIM-1];
    for (unsigned long j=0;j<cells.size();j++)
      {
	Simplex *s=&cells[j];
	long config = s->getRightBoundaryConfig();
	if (config==0) continue;
	if (getNeighborCellPtr(s,config)==NULL) return true;
      }       
  
    return false;
  }

  template <typename OutputIterator>
  long getLocalNeighbors(Simplex *reference, int type, OutputIterator out) const
  {    
    //glb::console->print<LOG_DEBUG>("Local neighbors for %s:\n",reference->toString().c_str());
    long count=0;
    for (unsigned long i=0;i<reference->neighbors[type].size();i++)
      {
	//glb::console->print<LOG_DEBUG>("    ->%s\n",reference->neighbors[type][i]->toString().c_str());
	*out = reference->neighbors[type][i];
	++out;
	++count;
      }
    return count;
  }  
  
  // WARNING if reference->type==type, the reference is incorectly included
  // in the neighbors
  template <typename OutputIterator>
  long getAllNeighborsConfigs(ReferenceCell &reference, int type, OutputIterator outConfig) const
  {
    //glb::console->print<LOG_DEBUG>("looking for %d-neighbors for  %s\n",type,reference.toString().c_str());
    if ((reference.simplex->type==type)&&(type>0))
      {	
	std::vector<Simplex *> nei;
	getLocalNeighbors(reference.simplex,type-1,std::back_inserter(nei));
	//glb::console->print<LOG_DEBUG>("Found %ld %d-neighbors:\n",nei.size(),type-1);
	//for (int i=0;i<nei.size();i++) glb::console->print<LOG_DEBUG>("    ->%s\n",nei[i]->toString().c_str());
	long count=0;
	for (unsigned long i=0;i<nei.size();i++)
	  {
	    ReferenceCell neiRef(nei[i],reference.coords);
	    count+=getAllNeighborsConfigs(neiRef,type,outConfig);
	  }
	return count;
      }

    long count=0;
    //SimplexConfig refConf[1<<NDIM];
    std::vector<SimplexConfig> refConf;   
    long nRefConf = reference.getEquivalentConfigs(std::back_inserter(refConf),isAlternate);
    //glb::console->print<LOG_DEBUG>("Found %ld equivalent configs:\n",refConf.size());
   
      
    for (int i=0;i<nRefConf;i++)
      {
	////glb::console->print<LOG_DEBUG>("    ->%s\n",refConf[i].toString().c_str());
	//glb::console->print<LOG_DEBUG>("    eq conf %d::%s\n",i,refConf[i].toString().c_str());
	Simplex *s=getNeighborCellPtr(refConf[i].simplex,
				      refConf[i].getConfig(),
				      isAlternate);
	//glb::console->print<LOG_DEBUG>("-->eq conf %d::%s\n",i,s->toString().c_str());
	std::vector<Simplex *> nei;
	getLocalNeighbors(s,type,std::back_inserter(nei));
	//glb::console->print<LOG_DEBUG>("        ->%ld local neighbors.\n",nei.size());
	for (unsigned long j=0;j<nei.size();j++)
	  {
	    (*outConfig)=SimplexConfig(nei[j],refConf[i].shift);
	    ++outConfig;
	    count++;
	  }
      }
    return count;
  } 

  template <class T>
  long swapIndex(long index,const T coords[NDIM]) const 
  {
    long result=index;
    if (!isAlternate) return index;
    for (int i=0;i<NDIM;i++)
      {
	long c=clipCoord(coords[i],i);
	if (c&1)
	  {
	    if (index&(1<<i)) 
	      result &= ~(1ul<<i);
	    else
	      result |= (1ul<<i);
	  }
      }
    return result;
  }
 
  template <class T>
  SimplexConfig getMainConfig(SimplexConfig &conf, const T coords[NDIM]) const
  {
    SimplexConfig result=conf;
    if (isAlternate)
      {
	// FIXME: NOT working  for 0<type<NDIM ? Or is it ?
	Simplex *s=result.simplex;	
	unsigned long isOdd[NDIM];
	for (int j=0;j<NDIM;j++) 
	  {
	    long c=clipCoord(coords[j]+conf.shift[j],j);
	    if (c&1)
	      isOdd[j]=1<<j;
	    else
	      isOdd[j]=0;
	  }

	unsigned long sym=0;
	for (int i=0;i<s->nVertices;i++)
	  for (int j=0;j<NDIM;j++)
	    {
	      // true if side is main side
	      if ( (isOdd[j]) == (s->vertices[i]->index&(1ul<<j)) )
		sym |= (1ul<<j);
	    }
	//printf("sym=~%lud (%ld %ld)\n",sym,isOdd[0],isOdd[1]);
	sym=~sym;
	result.simplex = getNeighborCellPtr(s,sym,isAlternate);
	for (int i=0;i<NDIM;i++)
	  if (sym&(1ul<<i)) result.shift[i]+=1;
      }
    else
      {
	Simplex *s=result.simplex;
	long sym=~(0l);
	for (int i=0;i<s->nVertices;i++)
	  sym &= s->vertices[i]->index;
	  
	result.simplex = getNeighborCellPtr(s,sym,isAlternate);
	for (int i=0;i<NDIM;i++)
	  if (sym&(1l<<i)) result.shift[i]+=1;
      }
    return result;
  }

  template <class T>
  long simplexToIndex(Simplex *simplex,const T refCoords[NDIM]) const
  {
    SimplexConfig c(simplex);
    return simplexConfigToIndex(c,refCoords);
  }

  template <typename T>
  T clipCoord(T val, int dir) const
  {
    T result=val;
    if (result<0) 
      {
	if (PERIODIC)
	  return result+dims[dir];
	else 
	  return -1;
      }
    if (result>=dims[dir]) 
      {
	if (PERIODIC) 
	  return result-dims[dir];
	else if (result==dims[dir])
	  return result;
	else return -1;
      }
    return result;    
  }
  
  long referenceCellToIndex(ReferenceCell &ref) const
  {
    long index=0;    
    int type=ref.simplex->type;
    for (int i=0;i<NDIM;i++) 
      {
	long c=clipCoord(ref.coords[i],i);
	index+=c*cellStride[type][i];
      }
    return index + ref.simplex->index;
  }

  template <class T>
  long simplexConfigToIndex(SimplexConfig &config,const T refCoords[NDIM]) const
  {
    long index=0;
    long coord[NDIM];
    int type=config.simplex->type;
    //glb::console->print<LOG_DEBUG>("computing index for %s:\n",config.toString().c_str());
    SimplexConfig mainConfig = getMainConfig(config,refCoords);    
    //glb::console->print<LOG_DEBUG>("  main config is %s\n",mainConfig.toString().c_str());
    //glb::console->print<LOG_DEBUG>("   Coords = [%ld,%ld,%ld] -> ",(long)refCoords[0],(long)refCoords[1],(long)refCoords[2]);
    for (int i=0;i<NDIM;i++)
      {
	coord[i]=refCoords[i]+mainConfig.shift[i];
	if (coord[i]<0) 
	  {
	    if (PERIODIC)
	      coord[i]+=dims[i];
	    else 
	      return -1;
	  }

	if (coord[i]>=dims[i]) 
	  {
	    if (PERIODIC) 
	      coord[i]-=dims[i];
	    else if (coord[i]==dims[i])
	      {
		// FIXME: CHECK IMPLEMENTATION FOR ALTERNATE !!!
		Simplex *tmp=mainConfig.simplex;		
		if ((isAlternate)&&(coord[i]&1))
		  {
		    for (int j=0;j<tmp->nVertices;j++)
		      {
			int vIndex = tmp->vertices[j]->index;
			if ((~vIndex) & (1<<i)) 
			  return -1;
		      }
		  }
		else
		  {
		    for (int j=0;j<tmp->nVertices;j++)
		      {
			int vIndex = tmp->vertices[j]->index;
			if (vIndex & (1<<i)) 
			  return -1;
		      }
		  }
	      }
	    else return -1;
	  }
      }    
    //glb::console->print<LOG_DEBUG>("[%ld,%ld,%ld]\n",(long)coord[0],(long)coord[1],(long)coord[2]);
    for (int i=0;i<NDIM;i++) index+=coord[i]*cellStride[type][i];
    
    long simplexIndex;
    if (isAlternate)
      {
	//int isOdd[NDIM];
	for (int j=0;j<NDIM;j++) 
	  {
	    clipCoord(refCoords[j]+mainConfig.shift[j],j);
	    //long c=clipCoord(refCoords[j]+mainConfig.shift[j],j);
	    //if (c&1) isOdd[j]=1;
	    //else isOdd[j]=0;
	  }

	if (mainConfig.simplex->type==0)
	  simplexIndex=0;//mainConfig.simplex->vertices[0]->getSwappedIndex(isOdd);
	else if (mainConfig.simplex->type==NDIM)
	  simplexIndex = mainConfig.simplex->index;	
	else 
	  {
	    simplexIndex = mainConfig.simplex->index;	
	    // FIXME: WONT WORK FOR INTERMEDIATE CELLS !!!!
	    // how to retrieve the local simplex index ????? 
	    
	    // Vertex *v[mainConfig.simplex->nVertices];
	    // for (int i=0;i<mainConfig.simplex->nVertices;i++)
	    //   v[i]=&vertices[mainConfig.simplex->vertices[i]->getSwappedIndex(isOdd)];
	    // Simplex simplex(v,v+mainConfig.simplex->nVertices);
	    // printf("Alternate simplex: %s\n",simplex.toString().c_str());
	    // simplexIndex=allCellsSet[mainConfig.simplex->type].find(simplex)->index;
	    
	  }


	//simplexIndex = mainConfig.simplex->index;	
	
      }
    else simplexIndex = mainConfig.simplex->index;
    //glb::console->print<LOG_DEBUG>("  -> Reference Index=%ld (+%d)\n",index,simplexIndex);

    return index + simplexIndex;
  }

  template <class OutputIterator>
  long getNeighborsIndex(ReferenceCell &reference,int type, OutputIterator out) const			
  {
    if ((reference.simplex->type==0)&&(type==NDIM))
      {
	
	/*
	std::vector<SimplexConfig> refConf;   
	long nRefConf = reference.getEquivalentConfigs(std::back_inserter(refConf),isAlternate);
	for (int i=0;i<nRefConf;i++)
	  {
	    glb::console->print<LOG_DEBUG>("    eq conf %d::%s\n",i,refConf[i].toString().c_str());
	    Simplex *s=getNeighborCellPtr(refConf[i].simplex,
				      refConf[i].getConfig(),
				      isAlternate);
	    glb::console->print<LOG_DEBUG>("-->eq conf %d::%s\n",i,s->toString().c_str());
	  }
	exit(0);
	*/
	/*
	ReferenceCell ref(s,refCoords);
	getNeighborsIndex(ReferenceCell &reference,int type, OutputIterator out);
	getNeighborsIndex(ReferenceCell &reference,int type, OutputIterator out);
	getNeighborsIndex(ReferenceCell &reference,int type, OutputIterator out);
	getNeighborsIndex(ReferenceCell &reference,int type, OutputIterator out)
	*/
      }

    long refIndex=-1;
    if (reference.simplex->type == type)
      refIndex = simplexToIndex(reference.simplex,reference.coords);
      
    std::set<long> result;
    std::vector<SimplexConfig> all;
    getAllNeighborsConfigs(reference,type,std::back_inserter(all));
    //glb::console->print<LOG_DEBUG>("found %ld configs.\n",all.size());
    for (unsigned long i=0;i<all.size();i++)
      {
	//glb::console->print<LOG_DEBUG>("conf %d: %s\n",i,all[i].toString().c_str());
	long index = simplexConfigToIndex(all[i],reference.coords);
	if (index>=0) result.insert(index);
      }

    long count=0;
    for (std::set<long>::iterator it=result.begin();it!=result.end();++it)
      {
	if (*it!=refIndex)
	  {
	    *out=*it;
	    ++out;
	    ++count;
	  }
      }
    return count;
  }
  /*
  template <class OutputIterator>
  long getNeighborsIndex(ReferenceCell cell,int type,OutputIterator out)
  {
    return getNeighborsIndex(cell.simplex,cell.coords,type,out);
  }
  */
  void cellToReferenceCell(const Cell cell,ReferenceCell &ref, bool correct=false) const
  {
    long id=cell.id();
    int type=cell.type();
    long cum=0;
    //FIXME
    // Here check if we are on the right boundary, 
    // and get a correct cell index (seems difficult ...)
    for (int i=NDIM-1;i>=0;i--)
      {
	ref.coords[i]=(id-cum)/cellStride[type][i];
	cum+=ref.coords[i]*cellStride[type][i];
      }
    id-=cum;
    if ((correct)&&(type==0)) id=swapIndex(id,ref.coords);
    ref.simplex = const_cast<Simplex *>(&allCells[type][id]);
    //ref.simplex = const_cast<Simplex *>(&allCells[type][id-cum]);
  }
  
  int countFacesPerVoxel(int type,int subspace=0, int odd=false)
  {
    int count=0;
    for (unsigned long i=0;i<allCells[type].size();i++)
      {	
	Simplex *s=&allCells[type][i];	
	if (s->getRightBoundaryConfig() != 0) continue;
	int j;
	for (j=0;j<NDIM;j++)
	  {
	    if (!(subspace&(1<<j))) continue;
	    int k;
	    for (k=0;k<s->nVertices;k++)
	      {
		int index=s->vertices[k]->index;
		if ((isAlternate)&&(odd))
		  {
		    if ((~index)&(1<<k))
		      break;
		  }
		else 
		  {
		    if (index&(1<<j))
		      break;
		  }
	      }
	    if (k!=s->nVertices) break;
	  }
	if (j==NDIM) count++;
      }
    return count;
  }

  void computeStrides()
  {
    for (int type=0;type<=NDIM;type++)
      nCellsPerVoxel[type]=countFacesPerVoxel(type);

    voxelStride[0]=1;
    for (int i=0;i<NDIM;i++) 
      voxelStride[i+1]=voxelStride[i]*dims[i];
    nVoxels=voxelStride[NDIM];
    
    //for (int odd=0;odd<=1;odd++)
    for (int type=0;type<=NDIM;type++)
      {
	for (int i=0;i<(1<<NDIM);i++)
	  {
	    facesPerVoxel[type][i] = countFacesPerVoxel(type,i);
	  }
      }
    
    // additional faces on the boundaries
    for (int type=0;type<=NDIM;type++)
      {
	addStride[type][0]=0;
	for (int i=0;i<NDIM;i++) 
	  {
	    if (PERIODIC) 
	      {
		addStride[type][i+1]=0;
		continue;
	      }
	     
	    addStride[type][i+1]=0;
	    for (int j=(1<<i)-1;j>=0;j--)
	      {
		long mask=(1<<i);		 
		//long count=countFacesPerVoxel(type,mask|j);
		long count=facesPerVoxel[type][mask|j];
		int stride=1;
		for (int k=0;k<i;k++) 
		  if (!(j&(1<<k))) stride*=dims[k];

		addStride[type][i+1]+=stride*count;
	      }
	  }
      }
    
    for (int type=0;type<=NDIM;type++)
      {
	cellStride[type][0]=nCellsPerVoxel[type];
	for (int i=0;i<NDIM;i++) 
	  {
	    cellStride[type][i+1] = cellStride[type][i]*dims[i] + addStride[type][i+1];
	  }
	nCells[type]=cellStride[type][NDIM];
      }
  }

  template <typename T>
  long Coords2BoundaryConfig(const T coords[NDIM]) const
  {
    long result=0;
    if (PERIODIC)
      {
	for (int i=0;i<NDIM;i++)
	  {
	    int cur;
	    if (coords[i]==0) cur=0;
	    else if (coords[i]==dims[i]-1) cur=3;
	    else if (coords[i]&1) cur=1;
	    else cur=2;
	    result+=cur<<(2*i);
	  }
      }
    else
      {
	for (int i=0;i<NDIM;i++)
	  {
	    int cur;
	    if (coords[i]==0) cur=0;
	    else if (coords[i]==dims[i]-1) cur=3;
	    else if (coords[i]==dims[i]) return -1;
	    else if (coords[i]&1) cur=1;
	    else cur=2;
	    result+=cur<<(2*i);
	    // we should ahve a case for (coords[i]==dims[i])
	    // but we compute it explicitely for now ...
	  }
      }
    return result;
  }

  long extractBoundaryConfig(long config, int which)
  {
    return (config>>(2*which))&3;
    
  }

  //[boundary configuration][current type][neighbor type]
  std::vector< std::vector<long> > shiftTable[1<<(2*NDIM)][NDIM+1][NDIM+1];

  // Precompute for every possible boundary configurations C (4 by dimension)
  // and for each simplex of type T and index I the shift S(C,N,T,T2) such that
  // its Nth neighbor with type T2 has index I2 = I+S(C,N,T,T2);
  void buildShiftTable()
  {
    long nConf=1<<(2*NDIM);
    int coord[NDIM];
    int confCoord[4];
  
    confCoord[0]=0;
    confCoord[1]=1;
    confCoord[2]=2;
    confCoord[3]=-1;
    
    for (long i=0;i<nConf;i++)
      {
	long j;
	for (j=0;j<NDIM;j++)
	  {
	    coord[j]=confCoord[extractBoundaryConfig(i,j)];
	    if (coord[j]==-1) 
	      coord[j]=(PERIODIC)?(dims[j]-1):dims[j]-1;	      
	    if (coord[j]>dims[j]) break;
	  }
	if (j!=NDIM) continue;
	//if (i==9) {printf("coord = [%d %d]",coord[0],coord[1]);exit(0);}
	for (int type=0;type<=NDIM;type++)
	  {
	    for (int neiType=0;neiType<=NDIM;neiType++)
	      {
		shiftTable[i][type][neiType].resize(allCells[type].size());		
		for (unsigned long k=0;k<allCells[type].size();k++)
		  {
		    std::vector<long> &table = shiftTable[i][type][neiType][k];
		    int k2=k; // FIX ME UGLY
		    if (type==0) k2=swapIndex(k,coord); // ugly !!
		    Simplex *s=&allCells[type][k2]; // k
		    ReferenceCell ref(s,coord);
		    //long refIndex = referenceCellToIndex(ref);
		    /*
		    if ((i==9)&&(type==2)&&(neiType==0))
		      {
			printf("SIMPLEX=%s\n",s->toString().c_str());
			printf("REF=%s with index %ld\n",ref.toString().c_str(),refIndex);
		      }
		    */
		    //long refIndex = simplexToIndex(s,coord);
		    getNeighborsIndex(ref,neiType,std::back_inserter(table));
		    ReferenceCell ref2(&allCells[neiType][0],coord);
		    long refIndex = referenceCellToIndex(ref2);
		    
		    for (unsigned long l=0;l<table.size();l++) 
		      {
			/*
			if ((i==9)&&(type==2)&&(neiType==0))
			printf ("TABLE[%d] = %ld - %ld\n",l,table[l],refIndex);
			*/
			table[l]-=refIndex;
		      }
		    
		  }
		//if ((i==9)&&(type==2)&&(neiType==0)) exit(0);
	      }
	  }
	
	//if (i==9) exit(0);
      }
  }

  void build(TesselationType t = TesselationTypeV::ANY)
  {
    TesselationType actualType=TesselationTypeV::USER_DEFINED;
    if (t!=TesselationTypeV::USER_DEFINED) actualType=setSimplices(t);
    isAlternate = (forceAlternate)?true:checkAlternate();
    if ((isAlternate)&&(!dimsAreEven())&&(PERIODIC))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("All dimensions must be even for alternate tesselations of periodic grids.\n");
	exit(-1);
      }
      
    glb::console->printFlush<LOG_STD>("Building %s %dD %s simplicial grid (type=%s) ... ",
				 (PERIODIC)?"periodic":"non-periodic",
				 NDIM,(isAlternate)?"alternate":"regular",
				 TesselationTypeSelect().getString(actualType).c_str());
    setNeighbors();  
    computeStrides();
    buildShiftTable();
    glb::console->printFlush<LOG_STD>("done.\n");
    /*
    
    for (int type=0;type<=NDIM;type++)
      for (int i=0;i<=NDIM;i++) 
	{
	  glb::console->print<LOG_PEDANTIC>("addStride[t=%d][%d] = %ld\n",type,i,addStride[type][i]);
	}
    
    std::vector<long> result;
    double refCoords[3]={5,5,0};
    int cellIndex=0;
    int cellType=0;
    Simplex *s=&allCells[cellType][cellIndex];
    ReferenceCell ref(s,refCoords);    
    int neiType=NDIM;
    long count=getNeighborsIndex(ref,neiType,std::back_inserter(result));
    glb::console->print<LOG_PEDANTIC>("%s (index=%ld) has %ld neighbors:\n",
				 ref.toString().c_str(),
				 referenceCellToIndex(ref),
				 result.size());
    for (int i=0;i<result.size();i++)
      {
	ReferenceCell c;
	cellToReferenceCell(Cell(neiType,result[i]), c);
	glb::console->print<LOG_PEDANTIC>("   -> %ld: %s\n",result[i],c.toString().c_str());
      }
    
    glb::console->print<LOG_PEDANTIC>("Using table :\n");
    std::vector<Cell> result2;
    Cell cell(ref.simplex->type,referenceCellToIndex(ref));
    getNeighbors(cell,std::back_inserter(result2),neiType);
    glb::console->print<LOG_PEDANTIC>("%s has %ld neighbors:\n",ref.toString().c_str(),result2.size());
    
    
    for (int i=0;i<result2.size();i++)
      {
	ReferenceCell c;
	cellToReferenceCell(result2[i], c);
	glb::console->print<LOG_PEDANTIC>("   -> %ld: %s\n",result2[i].id(),c.toString().c_str());
      }

    exit(0);
    */
    /*
    glb::console->print<LOG_PEDANTIC>("Voxels count : %ld\n",nVoxels);
    for (int type=0;type<NDIM+1;type++)
      {
	glb::console->print<LOG_PEDANTIC>("  created %ld %d-cells (%ld per voxel).\n",
	       nCells[type],type,nCellsPerVoxel[type]);
      }


    for (int type=0;type<=NDIM;type++)
      {
	std::vector<int> count;
	for (int i=1;i<(1<<NDIM);i++)
	  count.push_back(countFacesPerVoxel(type,i));
	glb::console->print<LOG_PEDANTIC>("Counts for type %d: ",type);
	for (int i=0;i<count.size();i++)
	  glb::console->print<LOG_PEDANTIC>("%d ",count[i]);
	glb::console->print<LOG_PEDANTIC>("\n");
      }

    */
    /*
    for (int u=0;u<2;u++)
      for (int v=0;v<2;v++)
	for (int k=0;k<2;k++)
	  {
	    std::vector<long> result;
	    double refCoords[3]={u,v,0};
	    Simplex *s=&allCells[NDIM][k];
	    ReferenceCell ref(s,refCoords);
	    int neiType=NDIM;
	    long count=getNeighborsIndex(ref,neiType,std::back_inserter(result));
	    glb::console->print<LOG_PEDANTIC>("%s (index=%ld) has %ld neighbors:\n",
					 ref.toString().c_str(),
					 referenceCellToIndex(ref),
					 result.size());
	    for (int i=0;i<result.size();i++)
	      {
		ReferenceCell c;
		cellToReferenceCell(Cell(neiType,result[i]), c);
		glb::console->print<LOG_PEDANTIC>("   -> %ld: %s\n",result[i],c.toString().c_str());
	      }
    
	    glb::console->print<LOG_PEDANTIC>("Using table :\n");
	    result.clear();
	    Cell cell(ref.simplex->type,referenceCellToIndex(ref));
	    getNeighbors(cell,std::back_inserter(result),neiType);
	    glb::console->print<LOG_PEDANTIC>("%s has %ld neighbors:\n",ref.toString().c_str(),result.size());
	    for (int i=0;i<result.size();i++)
	      {
		ReferenceCell c;
		cellToReferenceCell(Cell(neiType,result[i]), c);
		glb::console->print<LOG_PEDANTIC>("   -> %ld: %s\n",result[i],c.toString().c_str());
	      }
	  }

    exit(0);
*/
    /*
    for (int type1=NDIM;type1>=0;type1--)
      {
	std::vector<Simplex> &cells =  allCells[type1];
	for (int i=0;i<cells.size();i++)
	  {
	    glb::console->print<LOG_PEDANTIC>("%s has neighbors:\n%s\n----\n",
		   cells[i].toString().c_str(),
		   cells[i].neighborsToString().c_str());
	  }
      }
    */
  }
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

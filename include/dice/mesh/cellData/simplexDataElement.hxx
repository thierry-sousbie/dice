#ifndef __SIMPLEX_DATA_ELEMENT_HXX__
#define __SIMPLEX_DATA_ELEMENT_HXX__

#include "../../mesh/cellData/basicSimplexDataPolicies.hxx"


#include "../../internal/namespace.header"

template <int ID, int N, typename T = double,
	  template <typename,int,class,class,class> class OnInitPolicy = 
	  simplexInitDataPolicy::Copy,
	  template <int,typename,int,class,class,class> class OnRefinePolicy = 
	  simplexRefineDataPolicy::Dummy,
	  template <typename,int> class OnCoarsenPolicy = 
	  simplexCoarsenDataPolicy::Dummy,
	  bool OutputOnDump = true>
class SimplexDataElementT
{
public:  
  typedef T Type;
  typedef SimplexDataElementT<ID,N,T,OnInitPolicy,OnRefinePolicy,OnCoarsenPolicy,OutputOnDump> MyType;
  static const int SIZE = N;
  static const int INDEX = ID;
  static const bool OUTPUT_ON_DUMP=OutputOnDump;
  //static const bool CREATE_FUNCTOR=CreateFunctor;
  //static const size_t OFFSET; // Value is attributed through DECLARE_CELL_DATA_ELEMENT macro
  T value[N];

  template <class M, class S, class DT>
  void init(const M *mesh, const S* simplex, const DT* initVal)
  {
    OnInitPolicy<T,N,M,S,DT>::init(value, mesh, simplex, initVal);
  }

  template <class M, class S>
  void init(const M *mesh, const S* simplex)
  {
    OnInitPolicy<T,N,M,S,T>::init(value, mesh, simplex, static_cast<T*>(NULL));
  }
  
  template <int PASS, class M, class V, class S>
  void refine(const M *mesh, const V* newVertex, const S* simplex, const S* otherSimplex,
	      T* refValue, void *buffer)
  {    
    OnRefinePolicy<PASS,T,N,M,V,S>::refine
      (value,mesh,newVertex,simplex,otherSimplex,refValue,buffer);
  }

  void coarsen(T *rmValue)
  {
    OnCoarsenPolicy<T,N>::coarsen(value,rmValue);
  }

  T* getPointer() {return value;}
  const T* getConstPointer() const {return value;}

  void setValue(T val) {value[0]=val;}
  void setValueAt(T val, int at) {value[at]=val;}
  T  getValue() const {return value[0];}
  T  getValueAt(int at) const {return value[at];} 
  static int getSize() {return SIZE;}
  
  template <class WR>
  static void selfSerialize(const MyType *me, WR *writer)
  {
    writer->write(me->value,N);
  }

  template <class RE>
  static void selfUnSerialize(MyType *me, RE *reader)
  {
    reader->read(me->value,N);
  }
};

/*
template <typename T,
	  template <typename,int,class,class,class> class OnInitPolicy,
	  template <int,typename,int,class,class,class> class OnRefinePolicy,
	  template <typename,int> class OnCoarsenPolicy>
class SimplexDataElementT<0,T,OnInitPolicy,OnRefinePolicy,OnCoarsenPolicy> 
{
public:
  typedef T Type;
  typedef SimplexDataElementT<0,T,OnInitPolicy,OnRefinePolicy,OnCoarsenPolicy> MyType;
  static const int SIZE = 0;
  //static const size_t OFFSET; // Value is attributed through DECLARE_CELL_DATA_ELEMENT macro
  //T value;

  template <class M, class S, class DT>
  void init(const M *mesh, const S* simplex, const DT* initVal)
  {}
  
  template <int PASS, class M, class V, class S>
  void refine(const M *mesh, const V* newVertex, const S* simplex, const S* otherSimplex, T* v)
  {}

  void coarsen()
  {}

  T* getPointer() {return NULL;}
  T  getValue() const {return T();}
  T  getValueAt(int at) const {return T();}
  static int getSize() {return SIZE;}

  template <class WR>
  static void selfSerialize(const MyType *me, WR *writer)
  {}

  template <class RE>
  static void selfUnSerialize(MyType *me, RE *reader)
  {}
};
*/
#include "../../internal/namespace.footer"
#endif

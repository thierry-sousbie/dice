#ifndef __VERTEX_DATA_ELEMENT_HXX__
#define __VERTEX_DATA_ELEMENT_HXX__

#include "../../mesh/cellData/basicVertexDataPolicies.hxx"

#include "../../internal/namespace.header"

template <int ID, int N, typename T = double,
          template <typename, int, class, class, class> class OnInitPolicy =
              vertexInitDataPolicy::Copy,
          template <int, typename, int, class, class, class, class> class OnRefinePolicy =
              vertexRefineDataPolicy::Dummy,
          bool OutputOnDump = true>
class VertexDataElementT
{
public:
  typedef T Type;
  typedef VertexDataElementT<ID, N, T, OnInitPolicy, OnRefinePolicy, OutputOnDump> MyType;
  static const int SIZE = N;
  static const int INDEX = ID;
  static const bool OUTPUT_ON_DUMP = OutputOnDump;
  // static const bool CREATE_FUNCTOR=CreateFunctor;
  // static const size_t OFFSET; // Value is attributed through DECLARE_CELL_DATA_ELEMENT macro
  T value[N];

  template <class M, class DT, class V>
  void init(const M *mesh, const V *vertex, const DT *initVal)
  {
    OnInitPolicy<T, N, M, V, DT>::init(value, mesh, vertex, initVal);
  }

  template <class M, class V>
  void init(const M *mesh, const V *vertex)
  {
    OnInitPolicy<T, N, M, V, T>::init(value, mesh, vertex, static_cast<T *>(NULL));
  }

  // template <class M, class SEG, class V>
  // void refine(const M *mesh, const SEG &seg, const V* vertex, const T* v0, const T* v1)
  template <class M, class SEG, class V, class S>
  void refine(const M *mesh, const SEG &seg, const V *newVertex,
              S *const *simplices, int nSimplices)
  {
    OnRefinePolicy<INDEX, T, N, M, SEG, V, S>::
        refine(value, mesh, seg, newVertex, simplices, nSimplices);
    // OnRefinePolicy<T,N,M,SEG,V>::
    // refine(value,mesh,seg,vertex,v0,v1);
  }

  T *getPointer() { return value; }
  const T *getConstPointer() const { return value; }

  void setValue(T val) { value[0] = val; }
  void setValueAt(T val, int at) { value[at] = val; }
  T getValue() const { return value[0]; }
  T getValueAt(int at) const { return value[at]; }
  static int getSize() { return SIZE; }

  template <class WR>
  static void selfSerialize(const MyType *me, WR *writer)
  {
    writer->write(me->value, N);
  }

  template <class RE>
  static void selfUnSerialize(MyType *me, RE *reader)
  {
    reader->read(me->value, N);
  }
};

/*
// template <typename T,
// 	  template <typename,int,class,class,class> class OnInitPolicy,
// 	  template <typename,int,class,class,class> class OnRefinePolicy>
//class VertexDataElementT<-1,0,T,OnInitPolicy,OnRefinePolicy>
template <>
class VertexDataElementT<-1,0,char,vertexInitPolicy::Dummy,vertexRefinePolicy::Dummy>
{
public:
  typedef char Type;
  typedef Type T;
  typedef VertexDataElementT<-1,0,T,vertexInitPolicy::Dummy,vertexRefinePolicy::Dummy> MyType;
  static const int SIZE = 0;
  static const int INDEX = -1;
  //static const size_t OFFSET; // Value is attributed through DECLARE_CELL_DATA_ELEMENT macro

  template <class M, class V, class DT>
  void init(const M *mesh, const V* vertex, const DT* initVal)
  {}

  template <class M, class V>
  void init(const M *mesh, const V* vertex)
  {}

   template <class M, class SEG, class V>
  void refine(const M *mesh, const SEG &seg, const V* vertex,
        const T* v0, const T* v1)
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

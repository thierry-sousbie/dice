#ifndef __CELL_DATA_MACROS_HXX__
#define __CELL_DATA_MACROS_HXX__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../tools/helpers/helpers_macros.hxx"
#include "../../tools/helpers/helpers.hxx"

#include "../../mesh/cellData/vertexDataElement.hxx"
#include "../../mesh/cellData/simplexDataElement.hxx"
#include "../../mesh/cellData/cellDataFunctors_default.hxx"


#include "../../internal/namespace.header"

#define DEFINE_CELL_DATA_ELEMENT_INFO_STRUCT				\
  struct CellDataElementInfo						\
  {									\
    CellDataElementInfo():						\
      isVoid(true),							\
      id(-1),								\
      offset(0),							\
      name("empty")							\
	{}								\
    									\
    CellDataElementInfo(size_t id_,const std::string &name_, size_t offset_): \
      isVoid(false),							\
      id(id_),								\
      offset(offset_),							\
      name(name_)							\
      {}								\
    									\
    CellDataElementInfo(const CellDataElementInfo &other):		\
      isVoid(other.isVoid),						\
      id(other.id),							\
      offset(other.offset),						\
      name(other.name)							\
	{}								\
    									\
    const bool isVoid;							\
    const size_t id;							\
    const size_t offset;						\
    const std::string name;						\
  };									\


#define DECLARE_CELL_DATA_ELEMENT(DATA_TYPE , TYPE , NAME)		\
  template <>								\
  struct DATA_TYPE ::DeclaredCellData< DATA_TYPE:: TYPE:: INDEX > {	\
  typedef TYPE Result;						        \
static const size_t id = DATA_TYPE:: TYPE:: INDEX;			\
static const size_t offset = OFFSETOF( DATA_TYPE , NAME );		\
static CellDataElementInfo getDataInfo()				\
{return CellDataElementInfo( DATA_TYPE:: TYPE:: INDEX , #NAME , OFFSETOF( DATA_TYPE , NAME ) );} \
};									\
									\

/*
// DO NOT USE THIS EXCEPT WITH C++11
// -> assigning the values in a header file is unsafe in C++03 !
#define FINALIZE_CELL_DATA( DATA_TYPE )					\
  const typename DATA_TYPE ::DataElementInfo DATA_TYPE ::DATA_INFO[] = { \
    DATA_TYPE ::DeclaredCellData<0>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<1>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<2>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<3>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<4>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<5>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<6>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<7>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<8>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<9>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<10>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<11>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<12>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<13>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<14>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<15>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<16>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<17>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<18>::getDataInfo(),			\
    DATA_TYPE ::DeclaredCellData<19>::getDataInfo()};			\
*/								

#define INIT_CELL_DATA( MYTYPE, CELLTYPE )				\
  typedef MYTYPE MyType;						\
  DEFINE_CELL_DATA_ELEMENT_INFO_STRUCT;					\
  template < int ID >							\
  struct DeclaredCellData{						\
    typedef MY_NAMESPACE:: CELLTYPE##DataElementT<-1,0> Result;		\
  };									\
									\
  typedef CellDataElementInfo DataElementInfo;				\
									\
  template <int W> static int getOffset()				\
  {									\
    return MYTYPE ::DeclaredCellData<W>::offset;			\
  }									\
  									\
  template <int W>							\
  typename DeclaredCellData<W>::Result* getDataElementPtr()		\
  {									\
    return static_cast<typename DeclaredCellData<W>::Result*>		\
      ((void*)((size_t)this + MYTYPE ::DeclaredCellData<W>::offset));	\
  }									\
									\
  template <int W>							\
  const typename DeclaredCellData<W>::Result* getDataElementPtr() const	\
  {									\
    return static_cast<const typename DeclaredCellData<W>::Result*>	\
      ((void*)((size_t)this + MYTYPE ::DeclaredCellData<W>::offset));	\
  }									\
									\
  template <class M, int W>						\
  static typename MY_NAMESPACE:: cellDataFunctors::CELLTYPE##DataT<M,W>* \
  newDataElementFunctorPtr(const M* m)					\
  {									\
    return new MY_NAMESPACE:: cellDataFunctors::CELLTYPE##DataT<M,W>(m); \
  }									\
									\
  template <class M,class OutputIterator>				\
  static void createAllDataElementFunctorPtr(const M* m, OutputIterator out) \
  {									\
    MY_NAMESPACE:: cellDataMacroDetails::createAllDataElementFunctorPtr	\
      < DeclaredCellData, MY_NAMESPACE:: cellDataFunctors::CELLTYPE##DataT, M, OutputIterator> \
      (m,out);								\
  }									\
  									\
  static MY_NAMESPACE:: MpiDataType createDataElementsMpiStructType()	\
  {									\
    return MY_NAMESPACE:: cellDataMacroDetails::createDataElementsMpiStructType \
      <MyType,DeclaredCellData>();					\
  }									\
									\
  template <class WR>							\
  static void selfSerialize(const MyType *me, WR *writer)		\
  {									\
    return MY_NAMESPACE:: cellDataMacroDetails::selfSerialize		\
      <MyType,DeclaredCellData,WR>(me,writer);				\
  }									\
									\
  template <class RE>							\
  static void selfUnSerialize(MyType *me, RE *reader)			\
  {									\
    return MY_NAMESPACE:: cellDataMacroDetails::selfUnSerialize		\
      <MyType,DeclaredCellData,RE>(me,reader);				\
  }									\
									\
  template <class M, class SH, class V, class S>			\
  static void onRefineVertex(M *mesh, const SH &seg, V *v, S* const *simplices, int nSimplices,void *buffer) \
  {									\
    MY_NAMESPACE:: cellDataMacroDetails::onRefineVertex			\
      <DeclaredCellData,M,SH,V,S>(mesh,seg,v,simplices,nSimplices,buffer); \
  }									\
  									\
  template <int W, class M, class SH, class V, class S>			\
  static void refineSimplexData(M *mesh, const SH &seg, V *v, S **s0, S **s1, int nSimplices,void *buffer) \
  {									\
    typedef typename MY_NAMESPACE:: hlp::IsTrueT<(W>=0)>::Result Status; \
    MY_NAMESPACE:: cellDataMacroDetails::refineSimplexData		\
      <W,M,SH,V,S>(mesh,seg,v,s0,s1,nSimplices,buffer,Status());	\
  }									\
									\
  template <class M, class SH, class V, class S>			\
  static void onRefineSimplices(M *mesh, const SH &seg, V *v, S **s0, S **s1, int nSimplices,void *buffer) \
  {									\
    MY_NAMESPACE:: cellDataMacroDetails::onRefineSimplices		\
      <DeclaredCellData,M,SH,V,S>(mesh,seg,v,s0,s1,nSimplices,buffer);	\
  }									\
  									\
  template <class M, class S, class V>					\
  static void onCoarsenSimplex(M *mesh, S *keep, int keepVertexIndex, int removeIndex, S *remove, int removeVertexIndex, int keepIndex, const V *removeVertex) \
  {									\
    MY_NAMESPACE:: cellDataMacroDetails::onCoarsenSimplex		\
      <DeclaredCellData,M,S,V>(mesh,keep,keepVertexIndex,removeIndex,remove,removeVertexIndex,keepIndex,removeVertex); \
  }									\
  									\
  
namespace cellDataMacroDetails
{
  // CREATE ALL DATA ELEMENTS FUNCTORS

  template <template <int> class TL, template <class,int> class DT, 
	    class M, int W,class OutputIterator>
  static void createAllDataElementFunctorPtr(const M* m, OutputIterator out, hlp::IsFalse) 
  {}
								
  template <template <int> class TL, template <class,int> class DT, 
	    class M, int W,class OutputIterator>
  static void createAllDataElementFunctorPtr(const M* m, OutputIterator out, hlp::IsTrue) 
  {		    
    int flags=
      (TL<W>::Result::OUTPUT_ON_DUMP)?
      cellDataFunctors::F_NO_FLAG:
      cellDataFunctors::F_SKIP_ON_FILE_DUMP;

    *out=new DT<M,W>(m,flags);
    ++out;   
							
    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status; 
    createAllDataElementFunctorPtr<TL,DT,M,W+1>(m,out,Status());		
  }

  template <template <int> class TL, template <class,int> class DT, 
	    class M,class OutputIterator>
  static void createAllDataElementFunctorPtr(const M* m, OutputIterator out) 
  {									
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0)  >::Result Status; 
    createAllDataElementFunctorPtr<TL,DT,M,0,OutputIterator>(m,out,Status()); 
  }

  // CREATE DATA ELEMENTS MPI_STRUCT_TYPE

  template <class MT, template <int> class TL,int W>
  static void createDataElementsMpiStructType(MpiDataType &mpiDataType,hlp::IsFalse)
  {}

  template <class MT, template <int> class TL,int W>
  static void createDataElementsMpiStructType(MpiDataType &mpiDataType,hlp::IsTrue)
  {
    mpiDataType.push_back<typename TL<W>::Result::Type>
      (TL<W>::Result::SIZE,MT::template DeclaredCellData<W>::offset);
   
    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    createDataElementsMpiStructType<MT,TL,W+1>(mpiDataType,Status());
  }

  template <class MT, template <int> class TL>
  static MpiDataType createDataElementsMpiStructType()
  {
    MpiDataType mpiDataType;
   
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status; 
    createDataElementsMpiStructType<MT,TL,0>(mpiDataType,Status());    
    mpiDataType.commit<MT>();
								
    return mpiDataType;
  }

  // SERIALIZE

  template <class MT, template <int> class TL, class WR, int W>
  static void selfSerialize(const MT *me, WR *writer, hlp::IsFalse)
  {}
  
  template <class MT, template <int> class TL, class WR, int W>
  static void selfSerialize(const MT *me, WR *writer, hlp::IsTrue)
  {    
    TL<W>::Result::selfSerialize(me->template getDataElementPtr<W>(), writer);    

    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    selfSerialize<MT,TL,WR,W+1>(me,writer,Status());
  }
   
  template <class MT, template <int> class TL, class WR>
  static void selfSerialize(const MT *me, WR *writer)
  {
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status; 
    selfSerialize<MT,TL,WR,0>(me,writer,Status());    
  }

  // UNSERIALIZE

  template <class MT, template <int> class TL, class RE, int W>
  static void selfUnSerialize(MT *me, RE *reader, hlp::IsFalse)
  {}
  
  template <class MT, template <int> class TL, class RE, int W>
  static void selfUnSerialize(MT *me, RE *reader, hlp::IsTrue)
  {
    TL<W>::Result::selfUnSerialize(me->template getDataElementPtr<W>(), reader);

    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    selfUnSerialize<MT,TL,RE,W+1>(me,reader,Status());
  }
   
  template <class MT, template <int> class TL, class RE>
  static void selfUnSerialize(MT *me, RE *reader)
  {
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status; 
    selfUnSerialize<MT,TL,RE,0>(me,reader,Status());    
  }
  
  // OnRefineVertex
  template <template <int> class TL, class M, class SH, class V, class S, int W>
  static void onRefineVertex(M *mesh, const SH &seg, V *v, S* const *simplices, 
			     int nSimplices, void *buffer, hlp::IsFalse)
  {}
  
  template <template <int> class TL, class M, class SH, class V, class S, int W>
  static void onRefineVertex(M *mesh, const SH &seg, V *v, S* const *simplices, 
			     int nSimplices, void *buffer, hlp::IsTrue)
  {    
    v->template getDataElementPtr<W>()->refine(mesh,seg,v,simplices,nSimplices);
    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    onRefineVertex<TL,M,SH,V,S,W+1>(mesh,seg,v,simplices,nSimplices,buffer,Status());
  }
   
  template <template <int> class TL, class M, class SH, class V, class S>
  static void onRefineVertex(M *mesh, const SH &seg, V *v, S* const *simplices, 
			     int nSimplices, void *buffer)
  {
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status; 
    v->setRefinedCoords(mesh,seg,simplices,nSimplices,buffer);
    onRefineVertex<TL,M,SH,V,S,0>(mesh,seg,v,simplices,nSimplices,buffer,Status());    
  }

  // OnRefineSimplices
  // Refine data with index W
  template <int W, class M, class SH, class V, class S>
  static void refineSimplexData(M *mesh, const SH &seg, V* v, S** s0, S** s1, 
				int nSimplices, void *buffer, hlp::IsTrue)
  {
    typedef typename S::Data::template DeclaredCellData<W>::Result CellDataW;	
    typedef typename CellDataW::Type Data;
    const long bufferElSize=V::template getSimplexRefineBufferSize<M,SH>();

    char *b = static_cast<char*>(buffer);
    for (int i=0;i<nSimplices;++i)
      {
	// First copy the value in the cell before splitting into refdata
	Data refData[CellDataW::SIZE];
	const Data *dataPtr = s0[i]->template getDataElementPtr<W>()->getPointer();
	std::copy(dataPtr,dataPtr+CellDataW::SIZE,refData);

	// And then compute the new values in the two daughter cells
	s0[i]->template getDataElementPtr<W>()->template refine<0>
	  (mesh,v,s0[i],s1[i],refData,b);
	s1[i]->template getDataElementPtr<W>()->template refine<1>
	  (mesh,v,s1[i],s0[i],refData,b);

	b+=bufferElSize;
      }    
  }

  template <int W, class M, class SH, class V, class S>
  static void refineSimplexData(M *mesh, const SH &seg, V* v, S** s0, S** s1, 
				int nSimplices, void *buffer, hlp::IsFalse)
  {
  }

  template <template <int> class TL, class M, class SH, class V, class S, int W>
  static void onRefineSimplices(M *mesh, const SH &seg, V* v, S** s0, S** s1, 
				int nSimplices, void *buffer, hlp::IsFalse)
  {}
  
  template <template <int> class TL, class M, class SH, class V, class S, int W>
  static void onRefineSimplices(M *mesh, const SH &seg, V* v, S** s0, S** s1, 
				int nSimplices, void *buffer, hlp::IsTrue)
  { 
    refineSimplexData<W,M,SH,V,S>(mesh,seg,v,s0,s1,nSimplices,buffer,hlp::IsTrue());
    /*
    typedef typename S::Data::template DeclaredCellData<W>::Result CellDataW;	
    typedef typename CellDataW::Type Data;
    const long bufferElSize=V::template getSimplexRefineBufferSize<M,SH>();

    char *b = static_cast<char*>(buffer);
    for (int i=0;i<nSimplices;++i)
      {
	// First copy the value in the cell before splitting into refdata
	Data refData[CellDataW::SIZE];
	const Data *dataPtr = s0[i]->template getDataElementPtr<W>()->getPointer();
	std::copy(dataPtr,dataPtr+CellDataW::SIZE,refData);

	// And then compute the new values in the two daughter cells
	s0[i]->template getDataElementPtr<W>()->template refine<0>
	  (mesh,v,s0[i],s1[i],refData,b);
	s1[i]->template getDataElementPtr<W>()->template refine<1>
	  (mesh,v,s1[i],s0[i],refData,b);

	b+=bufferElSize;
      } 
    */
    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    onRefineSimplices<TL,M,SH,V,S,W+1>(mesh,seg,v,s0,s1,nSimplices,buffer,Status());
  }
  
  template <template <int> class TL, class M, class SH, class V, class S>
  static void onRefineSimplices(M *mesh, const SH &seg, V* v, S** s0, S** s1, 
				int nSimplices, void *buffer)
  {
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status;    
    onRefineSimplices<TL,M,SH,V,S,0>(mesh,seg,v,s0,s1,nSimplices,buffer,Status());
				     
  }

  // OnCoarsenSimplex
  template <template <int> class TL, class M, class S, class V, int W>
  static void onCoarsenSimplex(M *mesh, S *keep, int keepVertexIndex, int removeIndex, 
			       S *remove, int removeVertexIndex, int keepIndex, 
			       const V *removeVertex, hlp::IsFalse)
  {}
  
  template <template <int> class TL, class M, class S, class V, int W>
  static void onCoarsenSimplex(M *mesh, S *keep, int keepVertexIndex, int removeIndex, 
			       S *remove, int removeVertexIndex, int keepIndex, 
			       const V *removeVertex, hlp::IsTrue)
  {    
    keep->template getDataElementPtr<W>()->coarsen
      (remove->template getDataElementPtr<W>()->getPointer());
    typedef typename hlp::IsTrueT< (TL<W+1>::Result::SIZE > 0) >::Result Status;
    onCoarsenSimplex<TL,M,S,V,W+1>(mesh,keep,keepVertexIndex,removeIndex,remove,
				   removeVertexIndex,keepIndex,removeVertex,Status());
  }
   
  template <template <int> class TL, class M, class S, class V>
  static void onCoarsenSimplex(M *mesh, S *keep, int keepVertexIndex, int removeIndex, 
			       S *remove, int removeVertexIndex, int keepIndex, 
			       const V *removeVertex)
  {
    typedef typename hlp::IsTrueT< (TL<0>::Result::SIZE > 0) >::Result Status; 
    onCoarsenSimplex<TL,M,S,V,0>(mesh,keep,keepVertexIndex,removeIndex,remove,
				 removeVertexIndex,keepIndex,removeVertex,Status());    
  }


} // namespace

#include "../../internal/namespace.footer"
#endif


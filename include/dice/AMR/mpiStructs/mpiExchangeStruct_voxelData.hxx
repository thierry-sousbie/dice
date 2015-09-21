#ifndef __MPI_EXCHANGE_STRUCT_VOXEL_DATA_HXX__
#define __MPI_EXCHANGE_STRUCT_VOXEL_DATA_HXX__

#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/MPI/mpiDataType.hxx"

#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  
  template <class V>
  struct Mpi_VoxelDataT
  {
    typedef Mpi_VoxelDataT<V> MyType;
    typedef V Voxel;

    typedef typename Voxel::Data Data;
    typedef typename Voxel::ICoord ICoord;
    
    static const int NDIM = V::NDIM;
   
    ICoord index;
    Data   data;
    char   level;
   
    Mpi_VoxelDataT():level(-1)
    {}

    Mpi_VoxelDataT(Voxel *v)
    {
      set(v->getIndex(),v->data,v->getLevel());
    }

    bool isEmpty() const
    {
      return (level<0);
    }
    /*
    unsigned int getBaseCellIndex() const
    {
      return __mpi_base_cell_index;
    }
    */
    void set(ICoord id,Data &d,char l)
    {
      index=id;
      data=d;
      level=l;
    }
   
    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<ICoord>(1,OFFSETOF(MyType,index));
      mpiDataType.push_back<Data>(1,OFFSETOF(MyType,data));
      mpiDataType.push_back<char>(1,OFFSETOF(MyType,level));               

      mpiDataType.commit<MyType>();
      return mpiDataType;          
    }
  };
       
} //namespace

#include "../../internal/namespace.footer"
#endif

#ifndef __MPI_EXCHANGE_STRUCT_CHECK_REFINE_RESULT__
#define __MPI_EXCHANGE_STRUCT_CHECK_REFINE_RESULT__

#include "../../tools/MPI/mpiCommunication.hxx"
#include "../../tools/MPI/mpiDataType.hxx"
#include "../../mesh/mpiStructs/mpiExchangeStruct_base.hxx"


#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  // unfortunately, Mpi_RefineQueryResultT has to be of type POD as it will be 
  // used with MPI, so we cannot inherit and have to specialize the indexed 
  // and the non indexed type for each and every MPI_ type ...  
  struct Mpi_CheckRefineResult
  {    
    static const int IS_INDEXED=0;
    typedef Mpi_CheckRefineResult MyType;
    int index;
    float value;

    Mpi_CheckRefineResult():index(-1)
    {}
   
    bool isEmpty() const
    {
      return (index<0);
    }

    unsigned long getBaseCellIndex()
    {
      PRINT_SRC_INFO(LOG_ERROR);
      glb::console->print<LOG_ERROR>("Wrong type of mpiStruct.\n");
      glb::console->print<LOG_ERROR>("Use a structure with base cell index (template param WITH_INDEX=true)\n");	  
      exit(-1);

      return -1;
    }

    static MpiDataType createMpiStructType()
    {
      MpiDataType mpiDataType;

      mpiDataType.push_back<int>(1,OFFSETOF(MyType,index));      
      mpiDataType.push_back<float>(1,OFFSETOF(MyType,value));
      
      mpiDataType.commit<MyType>();      
      return mpiDataType;
      /*
      MPI_Datatype result;
      std::vector<MPI_Datatype> type;
      std::vector<int> blocklen;
      std::vector<MPI_Aint> disp;
      MyType tmp;      

      type.push_back( MPI_Type<int>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,index));
     
      type.push_back( MPI_Type<float>::get() );
      blocklen.push_back(1);
      disp.push_back(OFFSETOF(MyType,value));

      int n=sizeof(MyType)-(disp.back()+sizeof(value));
      if (n) 
	{	
	  glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",n,sizeof(MyType),(disp.back()+sizeof(value)));
	  for (int i=0;i<n;i++)
	    {
	      type.push_back( MPI_BYTE );
	      blocklen.push_back(1);
	      disp.push_back(disp.back()+1);
	    }
	}

      MPI_Type_create_struct(type.size(), &blocklen[0], &disp[0], &type[0], &result);
      MPI_Type_commit(&result); 
  
      return result;
      */
    }
  };
 
} // namespace

#include "../../internal/namespace.footer"
#endif

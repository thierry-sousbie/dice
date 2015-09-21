#ifndef __MPI_EXCHANGE_STRUCT_BASE_HXX__
#define __MPI_EXCHANGE_STRUCT_BASE_HXX__

#include "../../dice_globals.hxx"
#include "../../tools/MPI/myMpi.hxx"
#include "../../tools/helpers/helpers_macros.hxx"


#include "../../internal/namespace.header"

namespace mpiExchangeStruct
{
  template <bool WITH_INDEX>
  struct Base{
    static const int IS_INDEXED=0;
    bool setBaseCellIndex(unsigned long index)
    {
      return false;
    }  
    template <class MPS>
    bool setMpiBaseCellIndex(MPS &result)
    {
      return false;
    }
  };

  template <>
  struct Base<true>{
    static const int IS_INDEXED=1;
    unsigned long __base_cell_index;
    bool setBaseCellIndex(unsigned long index)
    {
      __base_cell_index=index;
      return true;
    }   
    template <class MPS>
    bool setMpiBaseCellIndex(MPS &result)
    {
      result.__mpi_base_cell_index=__base_cell_index;
      return true;
    }
  };
}

#include "../../internal/namespace.footer"
#endif

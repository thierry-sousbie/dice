#ifndef __GLOBAL_IDENTITY__
#define __GLOBAL_IDENTITY__

#include "../../tools/types/cell.hxx"


#include "../../internal/namespace.header"

// GlobalIdentity is used to give a unique index to a number of elements
// spread over several MPI processes. 
// A global identity is built from a local index which should be unique for 
// each element locally, and a Rank uniquely associated to each 
// MPI process (usually, the rank of the process)
// 
// This class is designed to be memory effecient as the its size is basically
// sizeof(T) (T being the first template argument).

// M is the maximum number of nodes
template <typename T = unsigned long, int M=65534> 
class GlobalIdentityT : public CellT<T, hlp::CountBits<M+1>::value >
{
 public:          
  typedef GlobalIdentityT<T,M> MyT;
  typedef CellT< T , hlp::CountBits<M>::value > Base;   

  static const int MAX_RANKS_COUNT= Base::MAX_TYPE;
  static const T MAX_RANK = Base::MAX_TYPE;  
  static const T EMPTY_RANK = Base::EMPTY_TYPE;  

  typedef typename Base::Value Value;
  typedef typename Base::Type Rank;
  typedef typename Base::Id Id;

  // static const MyT empty;
  // static const MyT max;

  GlobalIdentityT():Base()
    {
      //printf("count : %ld %ld\n",(long)1<<(long)hlp::CountBits<M+1>::value);
    }

  GlobalIdentityT(Base b):Base(b) // implicit conversion 
    {
      //printf("count : %ld %ld\n",(long)1<<(long)hlp::CountBits<M+1>::value);
    }

  GlobalIdentityT(T val):Base(val)
    {
    }
  GlobalIdentityT(long rank, long id):Base(rank,id)
    {
    }
  ~GlobalIdentityT()
    {    
    }
  
  Rank rank() const {return Base::type();}

  void setRank(Rank rank)
  {
    Base::setType(rank);
  }
 

 protected:
  
};

// template <class T, int M> 
// const GlobalIdentityT<T,M> GlobalIdentityT<T,M>::empty = GlobalIdentityT<T,M>::Base::empty;
// template <class T, int M> 
// const GlobalIdentityT<T,M> GlobalIdentityT<T,M>::max = GlobalIdentityT<T,M>::Base::max;

#include "../../internal/namespace.footer"
#endif

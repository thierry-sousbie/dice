#ifndef __CELL_HXX__
#define __CELL_HXX__

#include <iostream>
#include <utility>
#include <cmath>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../../tools/helpers/helpers.hxx"

/**
 * @file 
 * @brief  Defines an struct that can compactly store the identity of a cell as 
 * a pair of integers [t,i] representing its type and identity
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

template <class T = unsigned long, int N = 4>
class CellT {
public:
  
  typedef CellT<T,N> MyT;

  static const T TYPE_NBITS = N;
  static const T TYPE_DEC = (8*sizeof(T)-TYPE_NBITS);
  static const T TYPE_MASK = ( (((T)1)<<TYPE_NBITS) -1 ) << TYPE_DEC;
  static const T ID_MASK = ~TYPE_MASK;
  static const T NOTSET =(((T)1)<<TYPE_DEC)-1;
  static const T EMPTY_ID = NOTSET;
  static const T MAX_ID = EMPTY_ID-1;
  static const T EMPTY_TYPE = (((T)1)<<N)-1;
  static const T MAX_TYPE = EMPTY_TYPE-1;
  static const unsigned long AS_DOUBLE_FACTOR = (10ul)*((MAX_TYPE/10ul)+1);

  typedef T Value;
  typedef typename hlp::MinimalIntegerType< N             , false >::Type Type;
  typedef typename hlp::MinimalIntegerType< 8*sizeof(T)-N , false >::Type Id; 
 
  static const MyT empty;
  static const MyT max;

  CellT(Value val):cell_ID(val)
  {
  }
  
  CellT()
  {
    cell_ID=empty.cell_ID;
  }
  
  CellT(Type type, Id id)
  {
    set(type,id);
  }

  ~CellT() {}

  bool operator<  (const MyT &other) const {return cell_ID <  other.cell_ID;}
  bool operator<= (const MyT &other) const {return cell_ID <= other.cell_ID;}
  bool operator>  (const MyT &other) const {return cell_ID >  other.cell_ID;}
  bool operator>= (const MyT &other) const {return cell_ID >= other.cell_ID;}
  bool operator== (const MyT &other) const {return cell_ID == other.cell_ID;}
  bool operator!= (const MyT &other) const {return cell_ID != other.cell_ID;}

  T operator() (void) const {return cell_ID;}
  MyT &operator=(T const &val) {this->cell_ID=val; return *this;}
  
  Type type() const {return (Type) ((cell_ID&TYPE_MASK) >> TYPE_DEC);}
  Id id() const {return cell_ID&ID_MASK;}
  //bool isNotSet() const {return (id()==NOTSET);}

  bool isEmpty() const {return (*this == empty);}
  bool isMax() const {return (*this == max);}

  bool isMaxId() const {return (id()==MAX_ID);}
  bool isEmptyId() const {return (id()==EMPTY_ID);}

  bool isMaxType() const {return (type()==MAX_TYPE);}
  bool isEmptyType() const {return (type()==EMPTY_TYPE);}
 
  void set(Type type, Id id)
  {  
    cell_ID = make_cell(type,id);
  }

  void setType(Type type)
  {
    cell_ID = make_cell(type,id());
  }

  void setId(Id id)
  {
    cell_ID = make_cell(type(),id);
  }

  void set(Value val)
  {
    cell_ID = val;
  }

  Value get() const 
  {
    return cell_ID;
  }

  std::pair<Type,Id> getAsPair() const 
  {
    return std::make_pair(type(),id());
  }

  double getAsDouble() const 
  {
    
    //if (id()==EMPTY_ID) return -std::numeric_limits<double>::max();
    //else if (id()==MAX_ID) return std::numeric_limits<double>::max();
    //else 
    return static_cast<double>(type())/AS_DOUBLE_FACTOR+static_cast<double>(id());
  }  

  void write(std::ofstream &str) const 
  {
    str.write((const char *) &cell_ID, sizeof(Value));
  }

  void read(std::ifstream &str) const 
  {
    str.read((char *) &cell_ID, sizeof(Value));
  } 
  /*
  static const T emptyValue()
  {
    return empty.cell_ID;
  }
  */
private:
    
  Value cell_ID;
  Value make_cell(Type type, Id id) const 
  {
    return (((Value)type) << TYPE_DEC) + (Value)id;
  }

};

template <class T, int N> 
const CellT<T,N> CellT<T,N>::empty = CellT<T,N>(CellT<T,N>::EMPTY_TYPE,CellT<T,N>::EMPTY_ID);
template <class T, int N> 
const CellT<T,N> CellT<T,N>::max = CellT<T,N>(CellT<T,N>::MAX_TYPE,CellT<T,N>::MAX_ID);

#include "../../internal/namespace.footer"
#endif

#ifndef __MPI_DATA_TYPE_HXX__
#define __MPI_DATA_TYPE_HXX__

#include "../../tools/MPI/myMpi.hxx"


#include "../../internal/namespace.header"

class MpiDataType
{
public:
  typedef MpiDataType MyType;

  MpiDataType():
    commited(false),
    owned(true),
    result(MPI_DATATYPE_NULL)
  {}

  MpiDataType(const MpiDataType &other) // other's ownership is stolen on copy
  {
    commited=false;
    (*this)=other;
  }

  MyType &operator=(const MyType & other)
  {
    if (this != &other) 
      {
	freeResult();
	other.owned=false;
	owned=true;
	commited = other.commited;   
	result = other.result;
	type = other.type;
	blocklen = other.blocklen;
	disp = other.disp;
	typeSize = other.typeSize;
	resultTypeSize = other.resultTypeSize;
	subTypes = other.subTypes;
      }
    return (*this);
  }

  ~MpiDataType()
  {
    freeResult();    
  }

  template <class MPIT>
  void push_back(int blocklen_, MPI_Aint disp_)
  {
    typeSize.push_back(sizeof(MPIT));
    type.push_back(MPI_Type<MPIT>::get());
    blocklen.push_back(blocklen_);
    disp.push_back(disp_);
  }
  
  void push_back(const MpiDataType &type_, int blocklen_, MPI_Aint disp_)
  {          
    subTypes.push_back(type_);
    typeSize.push_back(subTypes.back().resultTypeSize);
    type.push_back(subTypes.back().getType());
    blocklen.push_back(blocklen_);
    disp.push_back(disp_);
  }

  bool isNull() const
  {
    return (result==MPI_DATATYPE_NULL);
  }

  template <class MT>
  void commit()
  {    
    //if (type.size()==0) return;

    if ((commited)&&(owned)&&(!isNull()))
      MPI_Type_free(&result);    

    // This is an empty type, but it should still be at least one byte large ...
    if (type.size() == 0)
      {
	typeSize.push_back(sizeof(char));
	type.push_back(MPI_BYTE);
	blocklen.push_back(1);
	disp.push_back(0);
      }  
    /*
    typeSize.push_back(0);
    type.push_back(MPI_UB);
    blocklen.push_back(1);
    disp.push_back(disp.back()+typeSize.back());
    */
    
    resultTypeSize=sizeof(MT);
    int n=sizeof(MT)-(disp.back()+typeSize.back());
    if (n) 
      {	
	// glb::console->print<LOG_PEDANTIC>("MpiStruct is padded by %d bytes. (%ld / %ld)\n",
	// 				  n,sizeof(MT),(disp.back()+typeSize.back()));
	for (int i=0;i<n;i++)
	  {
	    typeSize.push_back(sizeof(char));
	    type.push_back( MPI_BYTE );
	    blocklen.push_back(1);
	    disp.push_back(disp.back()+1);
	  }
      }
    

    MPI_Type_create_struct(type.size(), &blocklen[0], &disp[0], &type[0], &result);
    MPI_Type_commit(&result); 
    commited=true;
    owned=true;
  }

  MPI_Datatype &getType()
  {
    return result;
  }

private:
  void freeResult()
  {
    if ((commited)&&(owned)&&(!isNull()))
      MPI_Type_free(&result);
    commited=false;
  }
  //MpiDataType(const MyType &other);
  //MyType &operator=(const MyType & other);
  bool commited;
  mutable bool owned;    
  MPI_Datatype result;
  std::vector<MpiDataType> subTypes;

  std::vector<long> typeSize;
  std::vector<MPI_Datatype> type;
  std::vector<int> blocklen;
  std::vector<MPI_Aint> disp;

  long resultTypeSize;
};
 
#include "../../internal/namespace.footer"
#endif

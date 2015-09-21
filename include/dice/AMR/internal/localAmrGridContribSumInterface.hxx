#ifndef __LOCAL_AMR_GRID_CONTRIB_SUM_INTERFACE_HXX__
#define __LOCAL_AMR_GRID_CONTRIB_SUM_INTERFACE_HXX__

#include <limits>
#include <cmath>

#include "../../dice_globals.hxx"
#include "../../tools/helpers/helpers.hxx"
#include "../localAmrGridVisitors.hxx"

#include "../../internal/namespace.header"

namespace internal {
  
  // This is the case where AMR::Data != ST, specialization comes next ...
  // accuracyLevel -> +~10 bits waisted on computational errors
  template <class AMR, class ST, class HT, bool checkAccuracy> 
  class LocalAmrGridContribSumInterfaceT
  {
    typedef typename AMR::Data AmrT;
    typedef typename hlp::IF_<(sizeof(AmrT)<sizeof(unsigned long)),
      unsigned int,unsigned long>::Result Index;
    typedef typename AMR::Voxel Voxel;
    typedef typename hlp::IsPrimitiveType<ST>::Result IsPrimitive;   

    typedef short AccType;
    typedef typename hlp::IsTrueT<checkAccuracy>::Result AccChk;

  public:
    typedef ST SumType;

    LocalAmrGridContribSumInterfaceT(AMR *amr_, double accLevel=checkAccuracy?0.1:0):
      accuracyLevel(0)
    {
      amr=amr_;  
      reprojectionMode=false;
      decimalAccuracyLevel=accLevel;
      unsigned long nLeaves=amr->getNLeaves();
      if (nLeaves >= std::numeric_limits<Index>::max())
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>
	    ("Too many leaves (%ld) in the AMR grid to index them.\n",nLeaves);
	  glb::console->print<LOG_ERROR>("When summation type (ST) differs from AMR::Data, a temporary array index is stored in Voxel::Data.\n");
	  glb::console->print<LOG_ERROR>("Change Amr::Data type from float(%d) to float(%d) or set ST=Amr::Data in the projector.\n",sizeof(AmrT),sizeof(AmrT)*2);
	  exit(-1);
	}
      arr.assign(amr->getNLeaves(),0);
   
      IndexVisitor iv;
      localAmrGridVisitors::LeavesVisitor<AMR,IndexVisitor> liv(iv);
      amr->visitTree(liv,1);   

      if (checkAccuracy)
	{
	  accuracy.assign(amr->getNLeaves(),std::numeric_limits<AccType>::min());
	  frexp(accLevel,&accuracyLevel);
	  // This is the expected computational error (~14 bits waisted)
	  // FIXME: make this a parameter !
	  accuracyLevel = 14-accuracyLevel+1; 
	}
      else
	{
	  accuracyLevel=-1;
	  decimalAccuracyLevel=-1;
	  /*
	  if (accLevel>=0)
	    {
	      PRINT_SRC_INFO(LOG_WARNING);
	      glb::console->print<LOG_WARNING>
		("Template parameter 'checkAccuracy' is false, but an accuracy level was specified as constructor parameter.");
	      glb::console->print<LOG_WARNING>
		("The required accuracy level won't be checked if you do not enable 'checkAccuracy' template parameter.\n");
	    }
	  */
	}	 
    }

    void getAccuracyLevel(int &binaryLevel, double &decimalLevel) const
    {
      binaryLevel=accuracyLevel;
      decimalLevel=decimalAccuracyLevel;
    }
    
    template <class T>
    void add(const Voxel *v,const T &val)
    {    
      Index i=toIndexRef(&v->data);

      if (reprojectionMode)
	{
	  if (checkNeedReprojection(i,AccChk()))
	    {
	      Index j=toIndexRef(&arr[i]);
	      hArr[j] += val;
	    }
	}
      else
	{
	  accCheck(i,val,AccChk());
	  arr[i] += hlp::numericStaticCast<ST>(val);
	}
    }
    
    template <class P>
    void add(const P &p)
    {
      add(p.first,p.second);
    }

    template <class OUT=SumType>
    OUT getValue(const Voxel *v) const
    {
      Index i=toIndexRef(&v->data);
      if ((reprojectionMode)&&(checkNeedReprojection(i,AccChk())))
	{
	  Index j=toIndexRef(&arr[i]);
	  return hlp::numericStaticCast<OUT>(hArr[j]);
	}
      else return hlp::numericStaticCast<OUT>(arr[i]);
    }
    
    template <class P,class OUT=SumType>
    OUT getValue(const P &p)
    {    
      return getValue<OUT>(p.first);
    }    

    void commit(int nThreads)
    { 
      if (reprojectionMode)
	{	  
#pragma omp parallel for num_threads(nThreads)
	  for (unsigned long i=0;i<arr.size();++i)
	    {
	      if (checkNeedReprojection(i,AccChk()))
		arr[i]=hlp::numericStaticCast<ST>(hArr[toIndexRef(&arr[i])]);
	    }
	}

      typedef CommitVisitorT<ST,IsPrimitive::value> CV;
      CV cv(arr);
      localAmrGridVisitors::LeavesVisitor<AMR,CV> lcv(cv);
      amr->visitTree(lcv,nThreads);  

      reprojectionMode=false;
    }

    bool needReprojection(Voxel *v) const
    {  
      return checkNeedReprojection(v,AccChk());
    }
      
    long setReprojectionModeIfNeeded(int nThreads)
    {
      if (checkAccuracy) 
	{
	  typedef CheckAccuracyVisitorT<ST,AccType> CAV;
	  CAV cav(arr,accuracy);
	  localAmrGridVisitors::LeavesVisitor<AMR,CAV> lcav(cav);
	  amr->visitTree(lcav,nThreads);
	  long nFailed=cav.getCount();	  
	  if (nFailed>0) 
	    {
	      reprojectionMode=true;
	      hArr.assign(nFailed,0);
	    }
	  return nFailed;
	}
      return 0;
    }
  
  private:   
    std::vector<ST> arr;
    std::vector<HT> hArr;
    std::vector<AccType> accuracy; // Stores the accuracy of the sum
    AMR *amr;
    int accuracyLevel;
    double decimalAccuracyLevel;
    bool reprojectionMode;
  
    template <class In>
    static Index& toIndexRef(In *ptr)
    {
      return *static_cast<Index*>(static_cast<void*>(ptr));
    }

    template <class In>
    static const Index& toIndexRef(const In *ptr)
    {
      return *static_cast<const Index*>(static_cast<const void*>(ptr));
    }
    
    template <class T>
    void accCheck(Index id, const T &val, hlp::IsTrue) 
    {
      int result;
      // result -> value of the exponent in base b
      frexp(val,&result);
      /*
      std::cout << "received " << val << " N=" << result <<" => "
		<< result+accuracyLevel-std::numeric_limits<T>::digits
		<< std::endl;
      */
      result += accuracyLevel - std::numeric_limits<T>::digits;      
      if (result>accuracy[id]) accuracy[id]=result;
    }
    
    template <class T>
    void accCheck(Index id,const T &val, hlp::IsFalse) 
    {}

    bool checkNeedReprojection(Voxel *v, hlp::IsTrue) const
    {  
      return (accuracy[toIndexRef(&v->data)]==std::numeric_limits<AccType>::max());
    }

    bool checkNeedReprojection(Voxel *v, hlp::IsFalse) const
    {
      return false;
    }

    bool checkNeedReprojection(Index i, hlp::IsTrue) const
    {  
      return (accuracy[i]==std::numeric_limits<AccType>::max());
    }

    bool checkNeedReprojection(Index i, hlp::IsFalse) const
    {
      return false;
    }
    
    // Tree leaves visitors defined from here
    
    class IndexVisitor
    {
      typedef typename AMR::Voxel Voxel;    
    public:
      IndexVisitor(Index initVal=0):counter(initVal)
      {}
     
      void visit(Voxel *voxel)
      { 
	Index *i = static_cast<Index*>(static_cast<void*>(&voxel->data));
	(*i)=counter;
	++counter;	
      }
    private:
      Index counter;
    };

    // primitive type specialization
    template <class ST_, class AT_>
    class CheckAccuracyVisitorT
    {
      typedef typename AMR::Voxel Voxel;    
    public:

      CheckAccuracyVisitorT(std::vector<ST_> &sumArr_,			    
			    std::vector<AT_> &accArr_):
	sumArr(sumArr_),
	accArr(accArr_)	
      {count=0;}

      void visit(Voxel *voxel)
      { 
	Index i = *static_cast<Index*>(static_cast<void*>(&voxel->data));
	int result;      
	frexp(sumArr[i],&result);
	/*	
	if (result<=accArr[i])
	  std::cout << "Voxel " << voxel->getIndex() << "is "<< accArr[i]-result
		    << " digits below threshold with value:"
		    << sumArr[i] << "("<<result<<"/"<<accArr[i]<<")"<<std::endl;
	*/
	if (result<=accArr[i])
	  {
	    Index j;
#pragma omp atomic capture
	    {
	      j=count;
	      count+=1;
	    }
	    toIndexRef(&sumArr[i])=j;
	    accArr[i]=std::numeric_limits<AccType>::max();
	  }
      }
      
      int getCount() const {return count;}
      
    private:
      int count;
      std::vector<ST_> &sumArr;      
      std::vector<AT_> &accArr;
    };

    // primitive type specialization
    template <class ST_, bool BoostMulti>
    class CommitVisitorT
    {
      typedef typename AMR::Voxel Voxel;    
    public:

      CommitVisitorT(const std::vector<ST_> &stArr):
	arr(stArr)
      {}
      
      void visit(Voxel *voxel) const
      { 
	Index i = *static_cast<Index*>(static_cast<void*>(&voxel->data));
	voxel->data = arr[i];	 
      }

    private:
      const std::vector<ST_> &arr;
    };

    // Non primitive type specialization
    template <class ST_>
    class CommitVisitorT<ST_,false>
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::Data Data;
    public:

      CommitVisitorT(const std::vector<ST_> &stArr):
	arr(stArr)
      {}
  
      void visit(Voxel *voxel) const
      { 
	Index i = *static_cast<Index*>(static_cast<void*>(&voxel->data));
	voxel->data = arr[i].template convert_to<Data>();
      }

    private:
      const std::vector<ST_> &arr;
    };

  };
  
  // This is the specialized case where AMR::Data == ST
  // Note that we cannot use it if we are checking accuracy 
  // FIXME: Make yet another specialization ?
  template <typename AMR, class HT>
  class LocalAmrGridContribSumInterfaceT<AMR,typename AMR::Data,HT,0>
  {
    typedef typename AMR::Data AmrT;
    typedef typename AMR::Voxel Voxel;
  public:
    typedef AmrT SumType;

    LocalAmrGridContribSumInterfaceT(AMR *amr_, double accLevel=0){}    

    template <class T>
    static void add(Voxel *v,const T &val)
    {
      v->data+=val;
    }

    template <class P>
    static void add(const P &p)
    {
      p.first->data+=p.second;
    }

    template <class OUT=SumType>
    OUT getValue(const Voxel *v) const
    {
      return hlp::numericStaticCast<OUT>(v->data);
    }

    template <class P,class OUT=SumType>
    OUT getValue(const P &p)
    {    
      return hlp::numericStaticCast<OUT>(p.first->data);
    }      
   
    static void commit(int nThreads)
    {} 

    bool needReprojection(Voxel *v)
    {
      return false;
    }

    int setReprojectionModeIfNeeded(int nThreads)
    {
      return 0;
    }

    void getAccuracyLevel(int &binaryLevel, double &decimalLevel) const
    {
      binaryLevel=0;
      decimalLevel=0;
    }
  };
}

#include "../../internal/namespace.footer"
#endif

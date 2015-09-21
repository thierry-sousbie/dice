#ifndef __SUBGRID_ITERATOR_HXX_
#define __SUBGRID_ITERATOR_HXX_

#include <iterator>
#include "../regularGridNavigation.hxx"
#include "../regularGridFieldLayout.hxx"
#include "../valLocationType.hxx"
#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

namespace internal {

  template <class G, int L = G::LAYOUT> class field_iterator;
  template <class G, int L = G::LAYOUT> class oneField_iterator;
  template <class G, int L = G::LAYOUT> class inOrder_iterator;
  template <class G, int L = G::LAYOUT> class SubgridIteratorT;

  //template <>
  template <class G>
  class SubgridIteratorT<G,regularGridFieldLayout::INTERLEAVED>
  {
  public:
    template <class G2, int I2> friend class SubgridIteratorT;  

    typedef G Grid;
    
    typedef std::forward_iterator_tag iterator_category;
    typedef SubgridIteratorT<Grid> self_type;

    typedef Grid container_type;
    typedef typename container_type::value_type value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long  difference_type;

    //typedef typename container_type::ValLocationType  ValLocationType;
    //typedef typename container_type::ValLocationTypeV ValLocationTypeV;
    typedef typename container_type::GridNav GridNav;
    typedef typename GridNav::Direction Direction;
    
    static const int NDIM = container_type::NDIM;
   
    static const long N_INT_P = (1<<NDIM);
    //static const double N_INT_P_INV = 1.0L/(1<<NDIM);
    static const int IS_INTERLEAVED = G::IS_INTERLEAVED;

  protected:
    long delta[NDIM-1];
    long nFields;
    long i; 
    value_type *data;
    container_type *g;
    int w[NDIM];
    int min[NDIM];
    int max[NDIM];

    void init(bool end)
    {
      nFields = g->getNFields();

      if (end) {i=-1;data=NULL;}
      else 
	{	
	  i=0;
	  for (long ct=0;ct<NDIM;ct++)
	    {	     
	      if (max[ct]<=min[ct])
		{
		  i=-1;data=NULL;
		  return;
		}

	      w[ct] = min[ct];
	      i+=min[ct]*g->getValueStride(ct);
	    }
	  for (long ct=0;ct<NDIM-1;ct++) 
	    delta[ct]=g->getValueStride(ct+1)-
	      (max[ct]-min[ct])*g->getValueStride(ct);
	  
	  data = g->getDataPtr(i);
	  //data+=i*nFields;
	  //i/=nFields;
	}
    }

    void initFromMargin(container_type *g_, Direction region, bool end,
			const int *lowMargin, 
			const int *highMargin)
    {
      if (end) {i=-1;data=NULL;}
      else
	{
	  for (long ct=0;ct<NDIM;ct++)
	    {
	      if (region==GridNav::undefined())
		{
		  min[ct]=0;
		  max[ct]=g->getArrayDim(ct);
		}
	      else if (region&GridNav::dir(ct,-1))
		{
		  min[ct]=0;
		  max[ct]=lowMargin[ct];
		}
	      else if (region&GridNav::dir(ct,+1))
		{
		  min[ct]=max[ct]=g->getArrayDim(ct);
		  min[ct]-=highMargin[ct];
		}
	      else
		{
		  min[ct]=lowMargin[ct];
		  max[ct]=g->getArrayDim(ct);
		  max[ct]-=highMargin[ct];
		}
	    }
	  init(end);
	}
    }

  public:
    
    SubgridIteratorT(container_type *g_,
		     const int *xMin, 
		     const int *xMax,
		     bool end=false):
      g(g_),
      data(g_->getDataPtr())
    {
      for (int i=0;i<NDIM;++i)
	{
	  min[i]=xMin[i];
	  max[i]=xMax[i];
	}
      init(end);
    }
    
    
    SubgridIteratorT(container_type *g_, 
		     Direction region,
		     bool end=false):
      g(g_),
      data(g_->getDataPtr())
    {    
      initFromMargin(region,end,g->getParams().lowMargin,g->getParams().highMargin);     
    }

  
    template <class otherItT>
    SubgridIteratorT(container_type *g_, 
		     const otherItT &it,
		     const int *which=NULL):
      g(g_),
      data(g_->getDataPtr())
    {
      typedef typename std::vector<double>::const_iterator const_iteratorT;
      double C[otherItT::NDIM];
    
      it.valueCoord(C);
      nFields = g->getNFields();

      if (it.data==NULL) {i=-1;data=NULL;}
      else
	{	
	  i=0;
	  for (int ct=0;ct<NDIM;ct++)
	    {
	      int oct=(which!=NULL)?(which[ct]):ct;
	      const_iteratorT pos_it=hlp::findValue(g->getValueCoord(ct),C[oct]);
	    
	      if (pos_it==g->getValueCoord(ct).end())
		{
		  
		  fprintf(stderr,"ERROR in SubgridIteratorT constructor : trying to build an iterator from incompatible grid type.\n");
		  fprintf(stderr," this->coord[%d]=%g != other->Coord[%d]=%g.(delta=%g)\n",ct,*pos_it,oct,C[oct],C[oct]-(*pos_it));
		  exit(-1);
		}
	      
	      w[ct]=min[ct]=std::distance(g->getValueCoord(ct).begin(),pos_it);	    
	      max[ct]=min[ct]+it.max[oct]-it.min[oct];
	      i+=min[ct]*g->getValueStride(ct);
	    }
	  for (int ct=0;ct<NDIM-1;ct++) 
	    delta[ct]=g->getValueStride(ct+1)-
	      (max[ct]-min[ct])*g->getValueStride(ct);
	  data = g->getDataPtr(i);
	  //data+=i*nFields;
	  //i/=nFields;
	}      
    }

    virtual ~SubgridIteratorT(){}
    /*
    void print() const
    {
      int ii;

      char str[2000];
      sprintf(str,"delta = [");
      for (ii=0;ii<NDIM-1;ii++) sprintf(str,"%s%ld ",str,delta[ii]);
      sprintf(str,"%s];nFields=%ld, i=%ld, data=%ld(+%ld); ",str,nFields,i,(long)data,(long)(data-g->getDataPtr()));
      sprintf(str,"%sw = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,w[ii]);sprintf(str,"%s]; ",str);
      sprintf(str,"%smin = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,min[ii]);sprintf(str,"%s]; ",str);
      sprintf(str,"%smax = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,max[ii]);sprintf(str,"%s];",str);
      printf("%s\n",str);
    }
    */
    self_type &operator++()
    {     
      //if (data==NULL) return *this;
      w[0]++;data+=nFields;i++;

      if (w[0]>=max[0])
	{
	  if (NDIM==1) {i=-1;data=NULL;}
	  else
	    {
	      w[0]=min[0];
	      w[1]++;
	      i+=delta[0];data+=nFields*delta[0];
	      if (w[1]>=max[1])
		{
		  w[1]=min[1];
		  int ct=2;
		  while (ct<NDIM)
		    {
		      w[ct]++;i+=delta[ct-1];data+=nFields*delta[ct-1];
		      if (w[ct]>=max[ct])
			{
			  w[ct]=min[ct];
			  ct++;
			}
		      else break;
		    };
		  if (ct==NDIM) {i=-1;data=NULL;}
		}
	    }
	}
     
      return *this;
    }

    self_type &operator+=(long n)
    {
      for (int i=0;(i<n)&&(data!=NULL);++i)
	this->operator++();
      return *this;
    }

    self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    reference operator*() const
    {return *data;}
    
    reference operator[](std::size_t idx) const
    {return data[idx];}   

    long get_i() const {return i;}
    const int *get_w() const {return w;}    
    int get_curField() const {return 0;}

    bool operator==(const self_type& r) const
    {return (data==r.data);}

    bool operator!=(const self_type& r) const
    {return (data!=r.data);}

    bool operator<(const self_type& r) const
    {return (data<r.data);}

    bool operator>(const self_type& r) const
    {return (data>r.data);}	  

    bool operator<=(const self_type& r) const
    {return (data<=r.data);}	  	   

    bool operator>=(const self_type& r) const
    {return (data>=r.data);}
    /*
    

    static long dim(const self_type& it, int d)
    {
      return (d<NDIM)?(it.max[d]-it.min[d]):0;
    }

  
    void cellCoord(double *res) const 
    {
      for (int id=0;id<NDIM;++id) res[id]=g->getCellCoord(id)[w[id]];      
    }

    bool cellSize(double *res) const 
    {
      bool nonNull=true;
      for (int id=0;id<NDIM;++id) 
	if ((res[id]=g->getCellDelta(id)[w[id]])==0) nonNull=false;
      return nonNull;
    }

    double cellVolume() const
    {
      double result=1.;

      for (int id=0;id<NDIM;id++) 
	result*=g->getCellDelta(id)[w[id]];
      
      return result;
    }

    double cellAvg(int s=0) const
    {    
      const std::vector<long> &itPt=g->getIntegrationPointsPtr();      
      double res=0;

      for (int i=0;i<N_INT_P;i+=4)
	res+=data[itPt[i]]+data[itPt[i+1]]+data[itPt[i+2]]+data[itPt[i+3]];
      return res/(1<<NDIM);// *N_INT_P_INV;
    }
    
    template <class T>
    void coords(T res[NDIM]) const 
    {
      for (int id=0;id<NDIM;id++) 
	res[id]=g->getValueCoord(id)[w[id]];
    }

    template <class T>
    void delta(T res[NDIM]) const 
    {
      for (int id=0;id<NDIM;id++) 
	res[id]=g->getValueDelta(id)[w[id]];
    }

    double coord(int dim) const 
    {
      return g->getValueCoord(dim)[w[dim]];
    }

    double delta(int dim) const 
    {
      return g->getValueDelta(dim)[w[dim]];
    }

    double volume(int s=0) const
    {
      double result=1;
      for (int id=0;id<NDIM;id++) 
	result*=g->getVertexDelta(id)[w[id]];     
      return result;    
    }
    */

    void coordAndSize(double coord[NDIM], double delta[NDIM]) const 
    {
      for (int id=0;id<NDIM;id++) 
	{
	  coord[id]=g->getVertexCoord(id)[w[id]];      
	  delta[id]=g->getVertexCoord(id)[w[id]+1]-coord[id];      
	}
    }

    static long dim(const self_type& it, int d)
    {
      return (d<NDIM)?(it.max[d]-it.min[d]):0;
    }

    static long boxSize(const self_type& it)
    {
      long result=1;
      int s;
      if (it.data==NULL) return -1;
      for (s=0;s<NDIM;s++) result *=it.max[s]-it.min[s];
      return result;
    }
    
  };

  //template <>
  template <class G>
  class SubgridIteratorT<G,regularGridFieldLayout::CONSECUTIVE>
  {
  public:
    template <class G2, int I2> friend class SubgridIteratorT;  

    typedef G Grid;
    
    typedef std::forward_iterator_tag iterator_category;
    typedef SubgridIteratorT<Grid> self_type;

    typedef Grid container_type;
    typedef typename container_type::value_type value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long  difference_type;

    //typedef typename container_type::ValLocationType  ValLocationType;
    //typedef typename container_type::ValLocationTypeV ValLocationTypeV;
    typedef typename container_type::GridNav GridNav;
    typedef typename GridNav::Direction Direction;
    
    static const int NDIM = container_type::NDIM;
   
    static const long N_INT_P = (1<<NDIM);
    //static const double N_INT_P_INV = 1.0L/(1<<NDIM);
    static const int IS_INTERLEAVED = G::IS_INTERLEAVED;

  protected:
    long delta[NDIM-1];    
    long i;   
    container_type *g;
    value_type *data; 
    long nFields;
    int w[NDIM];
    int min[NDIM];
    int max[NDIM];
    int curField;

    void init(bool end, int startField=0)
    {
      //nFields = g->getNFields(); // incompatible with onefield_iterator !
      curField = startField;     
      
      if ((end)||(startField>=nFields)) {i=-1;data=NULL;}
      else 
	{	
	  i=0;
	  for (long ct=0;ct<NDIM;ct++)
	    {	     
	      if (max[ct]<=min[ct])
		{
		  i=-1;data=NULL;
		  return;
		}

	      w[ct] = min[ct];
	      i+=min[ct]*g->getValueStride(ct);
	    }
	  for (long ct=0;ct<NDIM-1;ct++) 
	    delta[ct]=g->getValueStride(ct+1)-
	      (max[ct]-min[ct])*g->getValueStride(ct);
	  
	  //data =  g->getDataPtr();
	  data = g->getDataPtr(i,startField); //i + g->getNValues()*startField;
	}
    }

    void initFromMargin(Direction region, bool end,
			const int *lowMargin, 
			const int *highMargin)
    {
      if (end) {i=-1;data=NULL;}
      else
	{
	  for (long ct=0;ct<NDIM;ct++)
	    {
	      if (region==GridNav::undefined())
		{
		  min[ct]=0;
		  max[ct]=g->getArrayDim(ct);
		}
	      else if (region&GridNav::dir(ct,-1))
		{
		  min[ct]=0;
		  max[ct]=lowMargin[ct];
		}
	      else if (region&GridNav::dir(ct,+1))
		{
		  min[ct]=max[ct]=g->getArrayDim(ct);
		  min[ct]-=highMargin[ct];
		}
	      else
		{
		  min[ct]=lowMargin[ct];
		  max[ct]=g->getArrayDim(ct);
		  max[ct]-=highMargin[ct];
		}
	    }
	  init(end);
	}
    }

  public:
    
    SubgridIteratorT(container_type *g_,
		     const int *xMin, 
		     const int *xMax,
		     bool end=false):
      g(g_),
      data(g_->getDataPtr()),
      nFields(g_->getNFields())
    {
      for (int i=0;i<NDIM;++i)
	{
	  min[i]=xMin[i];
	  max[i]=xMax[i];
	}
      init(end);
    }
    
    
    SubgridIteratorT(container_type *g_, 
		     Direction region,
		     bool end=false):
      g(g_),
      data(g_->getDataPtr()),
      nFields(g_->getNFields())
    {    
      initFromMargin(region,end,g->getParams().lowMargin,g->getParams().highMargin);     
    }

  
    template <class otherItT>
    SubgridIteratorT(container_type *g_, 
		     const otherItT &it,
		     const int *which=NULL):
      g(g_),
      data(g_->getDataPtr()),
      nFields(g_->getNFields())
    {
      typedef typename std::vector<double>::const_iterator const_iteratorT;
      double C[otherItT::NDIM];
    
      it.valueCoord(C);
      nFields=g->getNFields();
      curField=0;

      if (it.data==NULL) {i=-1;data=NULL;}
      else
	{	
	  i=0;
	  for (int ct=0;ct<NDIM;ct++)
	    {
	      int oct=(which!=NULL)?(which[ct]):ct;
	      const_iteratorT pos_it=hlp::findValue(g->getValueCoord(ct),C[oct]);
	    
	      if (pos_it==g->getValueCoord(ct).end())
		{
		  
		  fprintf(stderr,"ERROR in SubgridIteratorT constructor : trying to build an iterator from incompatible grid type.\n");
		  fprintf(stderr," this->coord[%d]=%g != other->Coord[%d]=%g.(delta=%g)\n",ct,*pos_it,oct,C[oct],C[oct]-(*pos_it));
		  exit(-1);
		}
	      
	      w[ct]=min[ct]=std::distance(g->getValueCoord(ct).begin(),pos_it);	    
	      max[ct]=min[ct]+it.max[oct]-it.min[oct];
	      i+=min[ct]*g->getValueStride(ct);
	    }

	  for (int ct=0;ct<NDIM-1;ct++) 
	    delta[ct]=g->getValueStride(ct+1)-
	      (max[ct]-min[ct])*g->getValueStride(ct);
	  
	  data = g->getDataPtr(i,0); // should be (i,it.curField)
	  //data += i;
	}      
    }

    virtual ~SubgridIteratorT(){}
    /*
    void print() const
    {
      int ii;

      char str[2000];
      sprintf(str,"delta = [");
      for (ii=0;ii<NDIM-1;ii++) sprintf(str,"%s%ld ",str,delta[ii]);
      sprintf(str,"%s];nFields=%ld, i=%ld, data=%ld(+%ld); ",str,nFields,i,(long)data,(long)(data-g->getDataPtr()));
      sprintf(str,"%sw = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,w[ii]);sprintf(str,"%s]; ",str);
      sprintf(str,"%smin = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,min[ii]);sprintf(str,"%s]; ",str);
      sprintf(str,"%smax = [",str);
      for (ii=0;ii<NDIM;ii++) sprintf(str,"%s%d ",str,max[ii]);sprintf(str,"%s];",str);
      printf("%s\n",str);
    }
    */
    self_type &operator++()
    {     
      //if (data==NULL) return *this;
      w[0]++;data++;i++;

      if (w[0]>=max[0])
	{
	  if (NDIM==1) {i=-1;data=NULL;}
	  else
	    {
	      w[0]=min[0];
	      w[1]++;
	      i+=delta[0];data+=delta[0];
	      if (w[1]>=max[1])
		{
		  w[1]=min[1];
		  int ct=2;
		  while (ct<NDIM)
		    {
		      w[ct]++;i+=delta[ct-1];data+=delta[ct-1];
		      if (w[ct]>=max[ct])
			{
			  w[ct]=min[ct];
			  ct++;
			}
		      else break;
		    };
		  if (ct==NDIM) init(false,curField+1);		    		     
		}
	    }
	}
     
      return *this;
    }

    self_type &operator+=(long n)
    {
      for (int i=0;(i<n)&&(data!=NULL);++i)
	this->operator++();
      return *this;
    }

    self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    reference operator*() const
    {return *data;}
    /*
    reference operator[](std::size_t idx) const
    {return data[nValues * idx];}   
    */
    long get_i() const {return i;}
    const int *get_w() const {return w;}
    int get_curField() const {return curField;}

    bool operator==(const self_type& r) const
    {return (data==r.data);}

    bool operator!=(const self_type& r) const
    {return (data!=r.data);}

    bool operator<(const self_type& r) const
    {return (data<r.data);}

    bool operator>(const self_type& r) const
    {return (data>r.data);}	  

    bool operator<=(const self_type& r) const
    {return (data<=r.data);}	  	   

    bool operator>=(const self_type& r) const
    {return (data>=r.data);}
    /*
    

    static long dim(const self_type& it, int d)
    {
      return (d<NDIM)?(it.max[d]-it.min[d]):0;
    }

  
    void cellCoord(double *res) const 
    {
      for (int id=0;id<NDIM;++id) res[id]=g->getCellCoord(id)[w[id]];      
    }

    bool cellSize(double *res) const 
    {
      bool nonNull=true;
      for (int id=0;id<NDIM;++id) 
	if ((res[id]=g->getCellDelta(id)[w[id]])==0) nonNull=false;
      return nonNull;
    }

    double cellVolume() const
    {
      double result=1.;

      for (int id=0;id<NDIM;id++) 
	result*=g->getCellDelta(id)[w[id]];
      
      return result;
    }

    double cellAvg(int s=0) const
    {    
      const std::vector<long> &itPt=g->getIntegrationPointsPtr();      
      double res=0;

      for (int i=0;i<N_INT_P;i+=4)
	res+=data[itPt[i]]+data[itPt[i+1]]+data[itPt[i+2]]+data[itPt[i+3]];
      return res/(1<<NDIM);// *N_INT_P_INV;
    }
    
    template <class T>
    void coords(T res[NDIM]) const 
    {
      for (int id=0;id<NDIM;id++) 
	res[id]=g->getValueCoord(id)[w[id]];
    }

    double coord(int dim) const 
    {
      return g->getValueCoord(dim)[w[dim]];
    }

    double volume(int s=0) const
    {
      double result=1;
      for (int id=0;id<NDIM;id++) 
	result*=g->getVertexDelta(id)[w[id]];     
      return result;    
    }
    */

    void coordAndSize(double coord[NDIM], double delta[NDIM]) const 
    {
      for (int id=0;id<NDIM;id++) 
	{
	  coord[id]=g->getVertexCoord(id)[w[id]];      
	  delta[id]=g->getVertexCoord(id)[w[id]+1]-coord[id];      
	}
    }

    static long dim(const self_type& it, int d)
    {
      return (d<NDIM)?(it.max[d]-it.min[d]):0;
    }

    static long boxSize(const self_type& it)
    {
      long result=1;
      int s;
      if (it.data==NULL) return -1;
      for (s=0;s<NDIM;s++) result *=it.max[s]-it.min[s];
      return result;
    }

  };

  //template <>
  template <class G>
  class field_iterator<G,regularGridFieldLayout::INTERLEAVED>
    : public SubgridIteratorT<G>
  {
  private:
    typedef SubgridIteratorT<G> Base;
    int fct;
  public:
    typedef field_iterator<G,regularGridFieldLayout::INTERLEAVED> self_type;

    field_iterator(const Base &it):Base(it),fct(0)
    {
    
    }
    
    self_type &operator++()
    {     
      if (Base::data==NULL) return *this;
      ++Base::data;
      if ((++fct)==Base::nFields) {fct=0;Base::w[0]++;Base::i++;}
      else return *this;

      if (Base::w[0]>=Base::max[0])
	{
	  if (Base::NDIM==1) {Base::i=-1;Base::data=NULL;}
	  else
	    {
	      Base::w[0]=Base::min[0];
	      Base::w[1]++;
	      Base::i+=Base::delta[0];Base::data+=Base::nFields*Base::delta[0];
	      if (Base::w[1]>=Base::max[1])
		{
		  Base::w[1]=Base::min[1];
		  int ct=2;
		  while (ct<Base::NDIM)
		    {
		      Base::w[ct]++;Base::i+=Base::delta[ct-1];Base::data+=Base::nFields*Base::delta[ct-1];
		      if (Base::w[ct]>=Base::max[ct])
			{
			  Base::w[ct]=Base::min[ct];
			  ct++;
			}
		      else break;
		    };
		  if (ct==Base::NDIM) {Base::i=-1;Base::data=NULL;}
		}
	    }
	}
     
      return *this;
    }

    self_type &operator+=(long n)
    {
      for (int i=0;(i<n)&&(Base::data!=NULL);++i)
	this->operator++();
      return *this;
    }

    self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    int get_curField() const {return fct;}
  };

  //template <>
  template <class G>
  class field_iterator<G,regularGridFieldLayout::CONSECUTIVE>
    : public SubgridIteratorT<G>
  {
  private:
    typedef SubgridIteratorT<G> Base;
    long fStep;
    long fStep_reset;
    int fct;
  public:
    typedef field_iterator<G,regularGridFieldLayout::CONSECUTIVE> self_type;

    field_iterator(const Base &it):Base(it),fct(0)
    {
      fStep = Base::g->getNValues();
      fStep_reset = -(fStep * Base::nFields);
    }   

    self_type &operator++()
    {     
      if (Base::data==NULL) return *this;
      Base::data+=fStep;
      if ((++fct)==Base::nFields) {fct=0;Base::w[0]++;Base::i++;Base::data+=fStep_reset;}
      else return *this;

      if (Base::w[0]>=Base::max[0])
	{
	  if (Base::NDIM==1) {Base::i=-1;Base::data=NULL;}
	  else
	    {
	      Base::w[0]=Base::min[0];
	      Base::w[1]++;
	      Base::i+=Base::delta[0];Base::data+=Base::delta[0];
	      if (Base::w[1]>=Base::max[1])
		{
		  Base::w[1]=Base::min[1];
		  int ct=2;
		  while (ct<Base::NDIM)
		    {
		      Base::w[ct]++;Base::i+=Base::delta[ct-1];Base::data+=Base::delta[ct-1];
		      if (Base::w[ct]>=Base::max[ct])
			{
			  Base::w[ct]=Base::min[ct];
			  ct++;
			}
		      else break;
		    };
		  if (ct==Base::NDIM) {Base::i=-1;Base::data=NULL;}
		}
	    }
	}
     
      return *this;
    }

    self_type &operator+=(long n)
    {
      for (int i=0;(i<n)&&(Base::data!=NULL);++i)
	this->operator++();
      return *this;
    }

    self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }
  
  }; 

  //template <>
  template <class G>
  class oneField_iterator <G,regularGridFieldLayout::INTERLEAVED>
    : public SubgridIteratorT<G>
  {

  private:
    typedef SubgridIteratorT<G> Base;
  protected :
    int fieldId;
  public:
    typedef oneField_iterator<G,regularGridFieldLayout::INTERLEAVED> self_type;

    oneField_iterator(const Base &it, int id):Base(it)
    {
      if (Base::data!=NULL) 
	Base::data+=id;
      fieldId=id;
    }

    int get_curField() const {return fieldId;}
  }; 

  //template <>
  template <class G>
  class oneField_iterator <G,regularGridFieldLayout::CONSECUTIVE>
    : public SubgridIteratorT<G>
  {

  private:
    typedef SubgridIteratorT<G> Base;
  protected :
    int fieldId;
  public:
    typedef oneField_iterator<G,regularGridFieldLayout::CONSECUTIVE> self_type;

    oneField_iterator(const Base &it, int id):Base(it)
    {
      if (Base::data!=NULL) 
	{
	  Base::data+=Base::g->getNValues()*id;   
	  Base::nFields=1;
	}      
      fieldId=id;
    }

    int get_curField() const {return fieldId;}
  }; 

  
  //template <>
  template <class G>
  class inOrder_iterator<G,regularGridFieldLayout::INTERLEAVED>
    : public field_iterator< G >
  {
    
  private:
    typedef field_iterator<G> Base;
    
  public:
    typedef inOrder_iterator<G,regularGridFieldLayout::INTERLEAVED> self_type;
    
    inOrder_iterator(const Base &it):Base(it)
    {}
  }; 
  


  //template <>
  template <class G>
  class inOrder_iterator <G,regularGridFieldLayout::CONSECUTIVE>
    : public SubgridIteratorT<G>
  {
    
  private:
    typedef SubgridIteratorT<G> Base;
    
  public:
    typedef inOrder_iterator<G,regularGridFieldLayout::CONSECUTIVE> self_type;
    
    inOrder_iterator(const Base &it):Base(it)
    {}
  }; 
  
  

} // internal
#include "../../internal/namespace.footer"
#endif

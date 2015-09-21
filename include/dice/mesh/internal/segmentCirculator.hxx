#ifndef __SEGMENT_CIRCULATOR_HXX__
#define __SEGMENT_CIRCULATOR_HXX__

#include <iterator>


#include "../../internal/namespace.header"


namespace internal {

  template <class H, class S, int NDIM>
  class SegmentCirculatorBaseT;
  
  template <class H, class S> class SegmentCirculatorT : 
    public SegmentCirculatorBaseT<H,S,S::NDIM> 
  {
  public:
    typedef SegmentCirculatorT<H,S> self_type;
    typedef SegmentCirculatorBaseT<H,S,S::NDIM> Base;

    typedef typename Base::Simplex Simplex;
    typedef typename Base::Vertex Vertex;
    typedef typename Base::Handle Handle;
  
    SegmentCirculatorT(const Handle &seg, Simplex *s=NULL):Base(seg,s)
    {

    }

    SegmentCirculatorT(const Handle &seg, Vertex *v, Simplex *s=NULL):Base(seg,v,s)
    {

    }
  };
  
  template <class H, class S>
  class SegmentCirculatorBaseT<H,S,3>
  {
  public:
    typedef H Handle;
    typedef S Simplex;
    typedef typename Simplex::Vertex Vertex;

    typedef SegmentCirculatorBaseT<H,S,3> self_type;

    typedef S* value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long difference_type;

    SegmentCirculatorBaseT(const Handle &seg,const Simplex *s=NULL)
      :status(0)
    {
      if (s==NULL) simplex=seg->getSimplex();    
      else simplex=const_cast<Simplex*>(s);
      vindex=0;
      Vertex *v = s->getVertex(vindex);
      seg->getVertices(vertex);
      if ((v==vertex[0])||(v==vertex[1]))
	{
	  v=s->getVertex(++vindex);
	  if ((v==vertex[0])||(v==vertex[1]))
	    v=s->getVertex(++vindex);    
	}
      simplexRef=const_cast<Simplex*>(simplex);
      vindexRef=vindex;
   
      //if (vertex[0]>vertex[1]) std::swap(vertex[0],vertex[1]);
    }

    SegmentCirculatorBaseT(const Handle &seg,const Vertex *v, const Simplex *s=NULL)
      :status(0)
    {
      if (s==NULL) simplex=seg->getSimplex();     
      else simplex=const_cast<Simplex*>(s);
      vindex=simplex->getVertexIndex(v);    
      simplexRef=const_cast<Simplex*>(simplex);
      vindexRef=vindex;
      seg->getVertices(vertex);
    }

    ~SegmentCirculatorBaseT()
    {
    
    }

    value_type operator->() const
    {
      return simplex;
    }

    value_type operator*() const
    {    
      return simplex;
    }

    self_type &operator++() 
    {  
      Simplex *s = simplex->getNeighbor(vindex);
    
      if (s==NULL) {reverse();return *this;}
      else if (s->isShadow()) {reverse();return *this;}
      else status=0;

      Vertex *v[Simplex::NVERT];
      simplex->getVertices(v);
  
      for (int i=0;i<Simplex::NVERT;i++)
	{
	  if (i==vindex) continue;
	  if ((v[i]==vertex[0])||(v[i]==vertex[1])) continue;	
	  v[0]=v[i];
	  break;
	}    
      simplex = s;
      vindex=simplex->getVertexIndex(v[0]);
    
      return *this;
    }

    const self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    bool operator==(const self_type& r) const
    {return (simplex==r.simplex);}

    bool operator!=(const self_type& r) const
    {return (simplex!=r.simplex);}

    bool crossedBoundary() const
    {
      return status;
    }

    int getNextFacetIndex() const
    {
      return vindex;
    }
 
    /*
      bool operator==(const self_type& r) const
      {return ((simplex==r.simplex)&&(vertex[0]==r.vertex[0])&&(vertex[1]==r.vertex[1]));}

      bool operator!=(const self_type& r) const
      {return ((simplex!=r.simplex)||(vertex[0]!=r.vertex[0])||(vertex[1]!=r.vertex[1]));}
    */

  private:
    Simplex *simplexRef;
    Simplex *simplex;
    Vertex *vertex[2];
    char vindex;
    char vindexRef;
    char status;
    /*
    void reverse()
    {
      //printf("REVERSE");fflush(0);
      
      simplex=simplexRef;
      if (vindexRef<0)
	{
	  vindexRef=-(vindexRef+1);
	  vindex=vindexRef;
	}
      else
	{
	  vindex=vindexRef;
	  vindexRef=-(vindexRef+1);

	  Vertex *v[Simplex::NVERT];
	  simplex->getVertices(v);
	  for (int i=0;i<Simplex::NVERT;i++)
	    {
	      if (i==vindex) continue;
	      if ((v[i]==vertex[0])||(v[i]==vertex[1])) continue;	
	      vindex=i;
	      break;
	    }  	

	  ++(*this);
	}  
	status=1;
    }
    */

    void reverse()
    {            
      simplex=simplexRef;
      vindex=vindexRef;
      
      Vertex *v[Simplex::NVERT];     

      simplex->getVertices(v);
      for (int i=0;i<Simplex::NVERT;i++)
	{
	  if ((i==vindex)||
	      (v[i]==vertex[0])||
	      (v[i]==vertex[1])) continue;	

	  vindex=i;
	  break;
	}  	
      
      while ((simplex->getNeighbor(vindex)!=NULL)&&
	     (!simplex->getNeighbor(vindex)->isShadow()))
	++(*this);
      
      simplex->getVertices(v);
      for (int i=0;i<Simplex::NVERT;i++)
	{
	  if ((i==vindex)||
	      (v[i]==vertex[0])||
	      (v[i]==vertex[1])) continue;	
	 
	  vindex=i;
	  break;
	}  

      status=1;
    }
  };

  template <class H, class S>
  class SegmentCirculatorBaseT<H,S,2>
  {
  public:
    typedef H Handle;
    typedef S Simplex;
    typedef typename Simplex::Vertex Vertex;

    typedef SegmentCirculatorBaseT<H,S,2> self_type;

    typedef S* value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long difference_type;

    SegmentCirculatorBaseT(const Handle &seg, Simplex *s=NULL)
      :status(0)
    {
      if (s==NULL)
	simplex=seg->getSimplex();
      else simplex=s;
      vindex=0;
      Vertex *v = s->getVertex(vindex);
      seg->getVertices(vertex);
      if ((v==vertex[0])||(v==vertex[1]))
	{
	  v=s->getVertex(++vindex);
	  if ((v==vertex[0])||(v==vertex[1]))
	    v=s->getVertex(++vindex);    
	}  
    }

    SegmentCirculatorBaseT(const Handle &seg,const Vertex *v, Simplex *s=NULL)
      :status(0)
    {
      if (s==NULL) simplex=seg->getSimplex();
      else simplex=s;
      vindex=simplex->getVertexIndex(v);
      seg->getVertices(vertex);
    }


    ~SegmentCirculatorBaseT()
    {
    
    }

    value_type operator->() const
    {
      return simplex;
    }

    value_type operator*() const
    {    
      return simplex;
    }

    self_type &operator++() 
    {    
      Vertex *v[Simplex::NVERT];
      /*
	char tmp[255];
	Simplex *nei[Simplex::NNEI];
	simplex->getNeighbors(nei); 
	for (int i=0;i<Simplex::NNEI;i++) if (nei[i]->isShadow()) printf("SHADOW\n");

	sprintf(tmp,"IN : simplex : %ld[%ld %ld %ld], v=[%ld,%ld], i=%d(%ld) ->",(long)simplex,(long)simplex->getVertex(0),(long)simplex->getVertex(1),(long)simplex->getVertex(2),(long)vertex[0],(long)vertex[1],(int)vindex,(long)simplex->getVertex(vindex));
      */

      Simplex *s = simplex->getNeighbor(vindex); 

      if (s==NULL) {status=1;return *this;}
      else if (s->isShadow()) {status=1;return *this;}     
      else status=0;
    
      simplex=s;
      simplex->getVertices(v);
    
      vindex=0;
      if ((v[0]==vertex[0])||(v[0]==vertex[1]))
	{
	  vindex=1;
	  if ((v[1]==vertex[0])||(v[1]==vertex[1]))
	    vindex=2;
	}
      /*
	simplex->getNeighbors(nei); 
	sprintf(tmp,"%s OUT : simplex : %ld[%ld %ld %ld], i=%d(%ld) (%s)\n",tmp,(long)simplex,(long)simplex->getVertex(0),(long)simplex->getVertex(1),(long)simplex->getVertex(2),(int)vindex,(long)v[vindex],simplex->isShadow()?"SHADOW":"NORMAL");
	printf("%s\n",tmp);    
	if (simplex->isShadow()) exit(-1);
      */
      return *this;
    }
    
    const self_type operator++(int)
    {
      self_type it(*this);
      ++(*this);
      return it;
    }

    bool operator==(const self_type& r) const
    {return (simplex==r.simplex);}

    bool operator!=(const self_type& r) const
    {return (simplex!=r.simplex);}

    bool crossedBoundary() const
    {
      return status;
    }

    int getNextFacetIndex() const
    {
      return vindex;
    }
 
    /*
      bool operator==(const self_type& r) const
      {return ((simplex==r.simplex)&&(vertex[0]==r.vertex[0])&&(vertex[1]==r.vertex[1]));}

      bool operator!=(const self_type& r) const
      {return ((simplex!=r.simplex)||(vertex[0]!=r.vertex[0])||(vertex[1]!=r.vertex[1]));}
    */

  private:
    Simplex *simplex;
    Vertex *vertex[2];
    int vindex;
    char status;
  };

  template <class H, class S>
  class SegmentCirculatorBaseT<H,S,1>
  {
  public:
    typedef H Handle;
    typedef S Simplex;
    typedef typename Simplex::Vertex Vertex;

    typedef SegmentCirculatorBaseT<H,S,1> self_type;

    typedef S* value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef long difference_type;

    SegmentCirculatorBaseT(const Handle &seg, Simplex *s=NULL)
    {
      if (s==NULL)
	simplex=seg->getSimplex();
      else 
	simplex=s;
    }

    SegmentCirculatorBaseT(const Handle &seg,const Vertex *v, const Simplex *s=NULL)
    {
      if (s==NULL)
	simplex=seg->getSimplex();
      else 
	simplex=s;
    }
  
    ~SegmentCirculatorBaseT()
    {
    
    }

    value_type operator->() const
    {
      return simplex;
    }

    value_type operator*() const
    {    
      return simplex;
    }

    self_type &operator++() 
    {
      return *this;
    }

    const self_type operator++(int)
    {
      return *this;
    }

    bool operator==(const self_type& r) const
    {return (simplex==r.simplex);}

    bool operator!=(const self_type& r) const
    {return (simplex!=r.simplex);}

    bool crossedBoundary() const
    {
      return 0;
    }

    int getNextFacetIndex() const
    {
      return 0;
    }
 
    /*
      bool operator==(const self_type& r) const
      {return ((simplex==r.simplex)&&(vertex[0]==r.vertex[0])&&(vertex[1]==r.vertex[1]));}

      bool operator!=(const self_type& r) const
      {return ((simplex!=r.simplex)||(vertex[0]!=r.vertex[0])||(vertex[1]!=r.vertex[1]));}
    */
  private:
    Simplex *simplex;
  };

} // namespace internal

#include "../../internal/namespace.footer"
#endif

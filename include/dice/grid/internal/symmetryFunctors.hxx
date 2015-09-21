#ifndef __INTERNAL_REGULAR_GRID_SYMMETRY_FUNCTOR_HXX__
#define __INTERNAL_REGULAR_GRID_SYMMETRY_FUNCTOR_HXX__

#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

namespace internal {
  namespace lrg {
    
    template <class G>
    class SymmetryFunctorInterface
    {
    public:
      typedef typename G::inOrderIterator iterator;
      typedef typename G::Data Data;
      static const int NDIM = G::NDIM;

      SymmetryFunctorInterface (G *g):
	m_grid(g)
      {}

      virtual ~SymmetryFunctorInterface()
      {}
      
      virtual int getSymmetries(iterator &it, Data **ptr) const=0;      
      virtual int getMaxCount() const=0;

    protected:
      G *m_grid;
    };

    template <class G, int I=0>
    class SymmetryFunctorSinglePlane:
      public SymmetryFunctorInterface<G>
    {
      static const int Index = hlp::MinT<G::NDIM-1,I>::value;
    public:
      typedef SymmetryFunctorInterface<G> Base;
      typedef typename Base::Data Data;
      typedef typename Base::iterator iterator;

      SymmetryFunctorSinglePlane(G *g):
	Base(g)
      {
	sz_m1=g->getValueCoord(Index).size();
	halfSz=sz_m1/2;
	sz_m1--;
      }

      ~SymmetryFunctorSinglePlane()
      {}
      
      int getSymmetries(iterator &it, Data **dest) const
      {
	int w[Base::NDIM];
	std::copy_n(it.get_w(),Base::NDIM,w);	

	if (w[Index]<halfSz)
	  {
	    w[Index]=sz_m1-w[Index];
	    dest[0]=Base::m_grid->getDataPtr(w,it.get_curField());
	    return 1;
	  }
	else return 0;
      }
      
      int getMaxCount() const {return 1;}
    private:
      int sz_m1;
      int halfSz;
    };

    template <class G, int I0=0, int I1=1>
    class SymmetryFunctorDualPlane:
      public SymmetryFunctorInterface<G>
    {
      static const int Index0 = hlp::MinT<G::NDIM-1,I0>::value;
      static const int Index1 = hlp::MinT<G::NDIM-1,I1>::value;

    public:
      typedef SymmetryFunctorInterface<G> Base;
      typedef typename Base::Data Data;
      typedef typename Base::iterator iterator;

      SymmetryFunctorDualPlane(G *g):
	Base(g)
      {
	sz0_m1=g->getValueCoord(Index0).size();
	halfSz0=sz0_m1/2;
	sz0_m1--;

	sz1_m1=g->getValueCoord(Index1).size();
	halfSz1=sz1_m1/2;
	sz1_m1--;
      }

      ~SymmetryFunctorDualPlane()
      {}
      
      int getSymmetries(iterator &it, Data **dest) const
      {
	int w[Base::NDIM];
	std::copy_n(it.get_w(),Base::NDIM,w);	

	if ((w[Index0]<halfSz0)&&(w[Index1]<halfSz1))
	  {
	    w[Index0]=sz0_m1-w[Index0];	    
	    dest[0]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index1]=sz1_m1-w[Index1];
	    dest[1]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index0]=sz0_m1-w[Index0];
	    dest[2]=Base::m_grid->getDataPtr(w,it.get_curField());
	    return 3;
	  }
	else return 0;
      }
      
      int getMaxCount() const {return 3;}
    private:
      int sz0_m1;
      int sz1_m1;

      int halfSz0;
      int halfSz1;
    };

    template <class G>
    class SymmetryFunctorTriplePlane:
      public SymmetryFunctorInterface<G>
    {
      static const int Index0 = hlp::MinT<G::NDIM-1,0>::value;
      static const int Index1 = hlp::MinT<G::NDIM-1,1>::value;
      static const int Index2 = hlp::MinT<G::NDIM-1,2>::value;
    public:
      typedef SymmetryFunctorInterface<G> Base;
      typedef typename Base::Data Data;
      typedef typename Base::iterator iterator;

      SymmetryFunctorTriplePlane(G *g):
	Base(g)
      {
	sz0_m1=g->getValueCoord(Index0).size();
	halfSz0=sz0_m1/2;
	sz0_m1--;

	sz1_m1=g->getValueCoord(Index1).size();
	halfSz1=sz1_m1/2;
	sz1_m1--;

	sz2_m1=g->getValueCoord(Index2).size();
	halfSz2=sz2_m1/2;
	sz2_m1--;
      }

      ~SymmetryFunctorTriplePlane()
      {}
      
      int getSymmetries(iterator &it, Data **dest) const
      {
	int w[Base::NDIM];
	std::copy_n(it.get_w(),Base::NDIM,w);	

	if ((w[Index0]<halfSz0)&&(w[Index1]<halfSz1)&&(w[Index2]<halfSz2))
	  {
	    w[Index0]=sz0_m1-w[Index0];	    
	    dest[0]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index1]=sz1_m1-w[Index1];
	    dest[1]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index0]=sz0_m1-w[Index0];
	    dest[2]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index1]=sz1_m1-w[Index1];
	    w[Index2]=sz2_m1-w[Index2];
	    dest[3]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index0]=sz0_m1-w[Index0];	    
	    dest[4]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index1]=sz1_m1-w[Index1];
	    dest[5]=Base::m_grid->getDataPtr(w,it.get_curField());

	    w[Index0]=sz0_m1-w[Index0];
	    dest[6]=Base::m_grid->getDataPtr(w,it.get_curField());
	    return 7;
	  }
	else return 0;
      }
      
      int getMaxCount() const {return 7;}
    protected:
      int sz0_m1;
      int sz1_m1;
      int sz2_m1;
      
      int halfSz0;
      int halfSz1;
      int halfSz2;
    };

    // CONSTANT symmetry
    template <class G, int I>
    class SymmetryFunctorConstant:
      public SymmetryFunctorInterface<G>
    {
      static const int NDIM=G::NDIM;
      static const int Index=hlp::MinT<G::NDIM-1,I>::value;
    public:
      typedef SymmetryFunctorInterface<G> Base;
      typedef typename Base::Data Data;
      typedef typename Base::iterator iterator;

      SymmetryFunctorConstant(G *g):
	Base(g)
      {
	sz=g->getValueCoord(Index).size();	
      }

      ~SymmetryFunctorConstant()
      {}
      
      int getSymmetries(iterator &it, Data **dest) const
      {
	if (it.get_w()[Index] != 0) return 0;

    	int w[Base::NDIM];
	std::copy_n(it.get_w(),Base::NDIM,w);	
	w[Index]++;
	for (int i=0;i<sz-1;i++)
	  {
	    dest[i]=Base::m_grid->getDataPtr(w,it.get_curField());
	    w[Index]++;
	  }
	
	return sz-1;
      }
      
      int getMaxCount() const {return sz-1;} 
      
    protected:     
      int sz;
    };
    
    // CENTRAL symmetry
    template <class G>
    class SymmetryFunctorCentral:
      public SymmetryFunctorInterface<G>
    {
      static const int NDIM=G::NDIM;

    public:
      typedef SymmetryFunctorInterface<G> Base;
      typedef typename Base::Data Data;
      typedef typename Base::iterator iterator;

      SymmetryFunctorCentral(G *g):
	Base(g)
      {
	for (int i=0;i<NDIM;++i)
	  sz_m1[i]=g->getValueCoord(i).size()-1;
      }

      ~SymmetryFunctorCentral()
      {}
      
      int getSymmetries(iterator &it, Data **dest) const
      {
	int w[Base::NDIM];
	std::copy_n(it.get_w(),Base::NDIM,w);	
	
	bool center=(w[0]==0);
	for (int i=1;i<NDIM;++i)
	  center=(center&&(w[i]==0));
	
	if (!center)
	  {
	    for (int i=0;i<NDIM;++i)
	      w[i]=sz_m1[i]-w[i];
	    (*dest)=Base::m_grid->getDataPtr(w,it.get_curField());
	    return 1;
	  }
	else return 0;

      }
      
      int getMaxCount() const {return 1;} 
      
    protected:
      int sz_m1[NDIM];
    };

  }
}

#include "../../internal/namespace.footer"
#endif

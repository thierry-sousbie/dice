#ifndef __MESH_QUADRATURE_FUNCTOR_ADAPTER_HXX__
#define __MESH_QUADRATURE_FUNCTOR_ADAPTER_HXX__

#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

namespace internal {
  namespace mesh {

    template <class FL, class M, class S>
    static void setFunctorListHelper(FL &fl, M *mesh, S *s, hlp::ConstantValue<false>)
    {}
  
    template <class FL, class M, class S>
    static void setFunctorListHelper(FL &fl, M *mesh, S *s, hlp::ConstantValue<true>)
    {
      //typedef typename FL::Next NextList;
      typedef typename hlp::ConstantValue< (FL::INDEX>0) > Continue;

      fl.getObject().set(mesh,s);      
      setFunctorListHelper(fl.getNext(),mesh,s,Continue());
    }
       
    template <class FL, class M, class S>
    void setFunctorList(FL &fl, M *mesh, S *s)
    {
      typedef typename hlp::ConstantValue< (FL::SIZE>0) > Continue;
      setFunctorListHelper(fl,mesh,s,Continue());
    }

  }
}
#include "../../internal/namespace.footer"
#endif

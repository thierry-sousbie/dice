#ifndef __LOOP_UNROLLER_HXX__
#define __LOOP_UNROLLER_HXX__

/**
 * @file 
 * @brief A class to unroll loops for faster iteration
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 *   \{
 */

namespace hlp {
  namespace internal {
    template <long N_UNROLL, long N_CUR=0>
    struct DoUnroll
    {
      typedef DoUnroll<N_UNROLL,N_CUR>    MyType;
      typedef DoUnroll<N_UNROLL,N_CUR+1>  Next;
 
      template <class F>
      static void apply(const F& func)
      {   
	func(N_CUR);
	Next::template apply<F>(func);
      }

      template <class F>
      static void apply(const F& func, long i)
      {    
	func(N_CUR+i);
	Next::template apply<F>(func,i);
      }

      template <class F>
      static void finish(const F& func)
      {
	apply<F>(func);
      }
    };

    template <long N_UNROLL>
    struct DoUnroll<N_UNROLL,N_UNROLL>
    {
      template <class F>
      static void apply(const F& func){}

      template <class F>
      static void apply(const F& func, long i){}

      template <class F>
      static void finish(const F& func){} 
    };
  } // namespace internal;

  template <long N_UNROLL, long N_TOT = N_UNROLL>
  struct Unroller
  {
    static const long N_LEFT = N_TOT%N_UNROLL;

    typedef Unroller<N_UNROLL,N_TOT>       MyType;
    typedef internal::DoUnroll<N_UNROLL,0> Next;
    typedef internal::DoUnroll<N_LEFT,0>   Finish;

    template <class F>
    static void apply(const F &func)
    {   
      Next::template apply<F>(func);
    }

    template <class F>
    static void apply(const F &func, long i)
    {    
      Next::template apply<F>(func,i);
    }

    template <class F>
    static void finish(const F &func)
    {
      Finish::template apply<F>(func,N_TOT-N_LEFT);
    }
  };

} // namespace hlp


/** \}*/
#include "../../internal/namespace.footer"
#endif

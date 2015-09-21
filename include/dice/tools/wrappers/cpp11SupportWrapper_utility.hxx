#ifndef MY_CPP11_SUPPORT_WRAPPER_UTILITY_HXX__
#define MY_CPP11_SUPPORT_WRAPPER_UTILITY_HXX__

/**
 * @file 
 * @brief Imports CPP11 support included in 'utility' to cpp11 namespace using native 
 * CPP11 support or tr1  depending on availability. Also defines declval from boost if
 * needed
 */

#ifdef HAVE_CPP11
#include <utility>
#include "../../internal/namespace.header"

namespace cpp11 {  
  using namespace ::std;
}

#include "../../internal/namespace.footer"

#elif HAVE_TR1

#ifdef HAVE_BOOST
#include <boost/utility/declval.hpp>
#endif

#ifdef HAVE_TR1_HEADER_PREFIX
#include <tr1/utility>
#else //HAVE_TR1_HEADER_PREFIX
#include <utility>
#endif //HAVE_TR1_HEADER_PREFIX

#include "../../internal/namespace.header"

namespace cpp11 {  
#ifdef HAVE_TR1_HEADER_PREFIX
  using namespace ::std::tr1; 
#else
  using namespace ::std;
#endif
#ifdef HAVE_BOOST
  using boost::declval;
#endif
}
#include "../../internal/namespace.footer"

#endif//HAVE_TR1

#endif

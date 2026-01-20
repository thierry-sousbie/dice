#ifndef MY_CPP11_SUPPORT_WRAPPER_NULLPTR_HXX__
#define MY_CPP11_SUPPORT_WRAPPER_NULLPTR_HXX__

/**
 * @file
 * @brief Imports nullptr support or create a mock nullptr
 */

#ifdef HAVE_NO_NULLPTR

const // this is a const object...
    class
{
public:
  template <class T>   // convertible to any type
  operator T *() const // of null non-member
  {
    return 0;
  } // pointer...
  template <class C, class T> // or any type of null
  operator T C::*() const     // member pointer...
  {
    return 0;
  }

private:
  void operator&() const; // whose address can't be taken
} nullptr = {};           // and whose name is nullptr

#endif

#endif

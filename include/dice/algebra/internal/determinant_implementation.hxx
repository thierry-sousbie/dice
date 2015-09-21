#ifndef __MY_DETERMINANT_IMPL_HXX__
#define __MY_DETERMINANT_IMPL_HXX__


#include "../../internal/namespace.header"

namespace internal {

  template <int ND, class T> struct DeterminantT;
  
  template <class T>
  struct DeterminantT<1,T> 
  {
    static const int NDIM = 1;
    typedef T Type;
  
    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      return mat[0][0];  
    }
  };

  template <class T>
  struct DeterminantT<2,T> 
  {
    static const int NDIM = 2;
    typedef T Type;
  
    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      return mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];
    }
  };

  template <class T>
  struct DeterminantT<3,T> 
  {
    static const int NDIM = 3;
    typedef T Type;
  
    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      T det2_12_01 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
      T det2_12_02 = mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
      T det2_12_12 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];

      return 
	mat[0][0] * det2_12_12 - 
	mat[0][1] * det2_12_02 + 
	mat[0][2] * det2_12_01 ;
    }
  };

  template <class T>
  struct DeterminantT<4,T> 
  {
    static const int NDIM = 4;
    typedef T Type;

    // Using Laplace expansion theorem
    // see http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      T s0 = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
      T s1 = mat[0][0] * mat[1][2] - mat[1][0] * mat[0][2];
      T s2 = mat[0][0] * mat[1][3] - mat[1][0] * mat[0][3];
      T s3 = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
      T s4 = mat[0][1] * mat[1][3] - mat[1][1] * mat[0][3];
      T s5 = mat[0][2] * mat[1][3] - mat[1][2] * mat[0][3];

      T c5 = mat[2][2] * mat[3][3] - mat[3][2] * mat[2][3];
      T c4 = mat[2][1] * mat[3][3] - mat[3][1] * mat[2][3];
      T c3 = mat[2][1] * mat[3][2] - mat[3][1] * mat[2][2];
      T c2 = mat[2][0] * mat[3][3] - mat[3][0] * mat[2][3];
      T c1 = mat[2][0] * mat[3][2] - mat[3][0] * mat[2][2];
      T c0 = mat[2][0] * mat[3][1] - mat[3][0] * mat[2][1];
      
      // s1/=10;s2/=10;s3/=10;s4/=10;s5/=100;
      // c1/=10;c2/=10;c3/=10;c4/=10;c5/=100;

      return s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
    }
  };

  template <class T>
  struct DeterminantT<5,T> 
  {
    static const int NDIM = 5;
    typedef T Type;

    // static inline T 
    // compute(const T (&mat)[NDIM][NDIM])
    // {
     
    // }
  };

  // FIXME: No way to do that faster like 4D ? 
  // Would be a pain to write anyway ;)
  template <class T>
  struct DeterminantT<6,T> 
  {
    static const int NDIM = 6;
    typedef T Type;

    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      // det2x2
      const T m01 = mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];
      const T m02 = mat[0][0]*mat[2][1] - mat[2][0]*mat[0][1];
      const T m03 = mat[0][0]*mat[3][1] - mat[3][0]*mat[0][1];
      const T m04 = mat[0][0]*mat[4][1] - mat[4][0]*mat[0][1];
      const T m05 = mat[0][0]*mat[5][1] - mat[5][0]*mat[0][1];
      const T m12 = mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1];
      const T m13 = mat[1][0]*mat[3][1] - mat[3][0]*mat[1][1];
      const T m14 = mat[1][0]*mat[4][1] - mat[4][0]*mat[1][1];
      const T m15 = mat[1][0]*mat[5][1] - mat[5][0]*mat[1][1];
      const T m23 = mat[2][0]*mat[3][1] - mat[3][0]*mat[2][1];
      const T m24 = mat[2][0]*mat[4][1] - mat[4][0]*mat[2][1];
      const T m25 = mat[2][0]*mat[5][1] - mat[5][0]*mat[2][1];
      const T m34 = mat[3][0]*mat[4][1] - mat[4][0]*mat[3][1];
      const T m35 = mat[3][0]*mat[5][1] - mat[5][0]*mat[3][1];
      const T m45 = mat[4][0]*mat[5][1] - mat[5][0]*mat[4][1];
      // minors of rank 3
      const T m012 = m01*mat[2][2] - m02*mat[1][2] + m12*mat[0][2];
      const T m013 = m01*mat[3][2] - m03*mat[1][2] + m13*mat[0][2];
      const T m014 = m01*mat[4][2] - m04*mat[1][2] + m14*mat[0][2];
      const T m015 = m01*mat[5][2] - m05*mat[1][2] + m15*mat[0][2];
      const T m023 = m02*mat[3][2] - m03*mat[2][2] + m23*mat[0][2];
      const T m024 = m02*mat[4][2] - m04*mat[2][2] + m24*mat[0][2];
      const T m025 = m02*mat[5][2] - m05*mat[2][2] + m25*mat[0][2];
      const T m034 = m03*mat[4][2] - m04*mat[3][2] + m34*mat[0][2];
      const T m035 = m03*mat[5][2] - m05*mat[3][2] + m35*mat[0][2];
      const T m045 = m04*mat[5][2] - m05*mat[4][2] + m45*mat[0][2];
      const T m123 = m12*mat[3][2] - m13*mat[2][2] + m23*mat[1][2];
      const T m124 = m12*mat[4][2] - m14*mat[2][2] + m24*mat[1][2];
      const T m125 = m12*mat[5][2] - m15*mat[2][2] + m25*mat[1][2];
      const T m134 = m13*mat[4][2] - m14*mat[3][2] + m34*mat[1][2];
      const T m135 = m13*mat[5][2] - m15*mat[3][2] + m35*mat[1][2];
      const T m145 = m14*mat[5][2] - m15*mat[4][2] + m45*mat[1][2];
      const T m234 = m23*mat[4][2] - m24*mat[3][2] + m34*mat[2][2];
      const T m235 = m23*mat[5][2] - m25*mat[3][2] + m35*mat[2][2];
      const T m245 = m24*mat[5][2] - m25*mat[4][2] + m45*mat[2][2];
      const T m345 = m34*mat[5][2] - m35*mat[4][2] + m45*mat[3][2];
      // minors of rank 4
      const T m0123 = m012*mat[3][3] - m013*mat[2][3] + m023*mat[1][3] - m123*mat[0][3];
      const T m0124 = m012*mat[4][3] - m014*mat[2][3] + m024*mat[1][3] - m124*mat[0][3];
      const T m0125 = m012*mat[5][3] - m015*mat[2][3] + m025*mat[1][3] - m125*mat[0][3];
      const T m0134 = m013*mat[4][3] - m014*mat[3][3] + m034*mat[1][3] - m134*mat[0][3];
      const T m0135 = m013*mat[5][3] - m015*mat[3][3] + m035*mat[1][3] - m135*mat[0][3];
      const T m0145 = m014*mat[5][3] - m015*mat[4][3] + m045*mat[1][3] - m145*mat[0][3];
      const T m0234 = m023*mat[4][3] - m024*mat[3][3] + m034*mat[2][3] - m234*mat[0][3];
      const T m0235 = m023*mat[5][3] - m025*mat[3][3] + m035*mat[2][3] - m235*mat[0][3];
      const T m0245 = m024*mat[5][3] - m025*mat[4][3] + m045*mat[2][3] - m245*mat[0][3];
      const T m0345 = m034*mat[5][3] - m035*mat[4][3] + m045*mat[3][3] - m345*mat[0][3];
      const T m1234 = m123*mat[4][3] - m124*mat[3][3] + m134*mat[2][3] - m234*mat[1][3];
      const T m1235 = m123*mat[5][3] - m125*mat[3][3] + m135*mat[2][3] - m235*mat[1][3];
      const T m1245 = m124*mat[5][3] - m125*mat[4][3] + m145*mat[2][3] - m245*mat[1][3];
      const T m1345 = m134*mat[5][3] - m135*mat[4][3] + m145*mat[3][3] - m345*mat[1][3];
      const T m2345 = m234*mat[5][3] - m235*mat[4][3] + m245*mat[3][3] - m345*mat[2][3];
      // minors of rank 5
      const T m01234 = m0123*mat[4][4] - m0124*mat[3][4] + m0134*mat[2][4] - m0234*mat[1][4] + m1234*mat[0][4];
      const T m01235 = m0123*mat[5][4] - m0125*mat[3][4] + m0135*mat[2][4] - m0235*mat[1][4] + m1235*mat[0][4];
      const T m01245 = m0124*mat[5][4] - m0125*mat[4][4] + m0145*mat[2][4] - m0245*mat[1][4] + m1245*mat[0][4];
      const T m01345 = m0134*mat[5][4] - m0135*mat[4][4] + m0145*mat[3][4] - m0345*mat[1][4] + m1345*mat[0][4];
      const T m02345 = m0234*mat[5][4] - m0235*mat[4][4] + m0245*mat[3][4] - m0345*mat[2][4] + m2345*mat[0][4];
      const T m12345 = m1234*mat[5][4] - m1235*mat[4][4] + m1245*mat[3][4] - m1345*mat[2][4] + m2345*mat[1][4];

      const T m012345 = 
	m01234*mat[5][5] - m01235*mat[4][5] + 
	m01245*mat[3][5] - m01345*mat[2][5] +
	m02345*mat[1][5] - m12345*mat[0][5];

      return m012345;
    }
  };
  /*
  template <class T>
  struct DeterminantT<6,T> 
  {
    static const int NDIM = 6;
    typedef T Type;

    static inline T 
    compute(const T (&mat)[NDIM][NDIM])
    {
      return 0;
    }
  };
  */

}

#include "../../internal/namespace.footer"
#endif

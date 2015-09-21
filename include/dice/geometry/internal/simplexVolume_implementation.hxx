#ifndef __SIMPLEX_VOLUME_IMPLEMENTATION_HXX__
#define __SIMPLEX_VOLUME_IMPLEMENTATION_HXX__

#include "../../algebra/determinant.hxx"
#include "../../tools/helpers/helpers.hxx"

#include "../../internal/namespace.header"

// Compute the volume of a n-dimensional k-simplex
//
// http://www.math.niu.edu/~rusin/known-math/97/volumes.polyh
// Let S have vertices given by the n-dimensional row vectors v_0, v_1, ... v_k.
// Let w_1 = v_1 - v_0, w_2 = v_2 - v_0, ... w_k = v_k - v_0
// and let W be the k by n matrix whose rows are the row vectors w_j,
// 1 \leq j \leq k. Then the GRAM DETERMINANT FORMULA says
//
//                             | det W W^t |^{1/2}
//                 Vol_k (S) = ------------------
//                                     k!
//
// where W^t is the transpose of W. Note that W W^t is a k by k matrix,
// so this formula makes sense.

namespace internal
{

  template <int ND, int NDW, class T>
  struct SimplexVolumeT;

  template <int ND, class T>
  struct SimplexVolumeT<ND,ND,T>
  {
  public:
    static const int NDIM=ND;
    static const int NDIM_W=ND;

    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {         
      static const long factor=hlp::FactorialT<NDIM>::value;
      return DeterminantT<NDIM,T>::compute(vec)/factor;
    }
  };

  /*
  template <class T>
  struct SimplexVolumeT<1,1,T>
  {
  public:
    static const int NDIM=1;
    static const int NDIM_W=1;
  
    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {
      return vec[0][0];
    }  
  };

  template <class T>
  struct SimplexVolumeT<2,2,T>
  {
  public:
    static const int NDIM=2;
    static const int NDIM_W=2;

    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {
      return DeterminantT<NDIM,T>::compute(vec)*0.5F; 
      //return (vec[0][0]*vec[1][1]-vec[1][0]*vec[0][1]) * 0.5;
    }
  };
  
  template <class T>
  struct SimplexVolumeT<3,3,T> 
  {
  public:
    static const int NDIM=3;
    static const int NDIM_W=3;

    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {    
      T det2_12_01 = vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0];
      T det2_12_02 = vec[1][0] * vec[2][2] - vec[1][2] * vec[2][0];
      T det2_12_12 = vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1];

      return (vec[0][0] * det2_12_12 - 
       	      vec[0][1] * det2_12_02 + 
       	      vec[0][2] * det2_12_01) * (1.0F/6.0F); 
      
    }
  };

  template <class T>
  struct SimplexVolumeT<4,4,T> 
  {
  public:
    static const int NDIM=4;
    static const int NDIM_W=4;
  
    // Using Laplace expansion theorem
    // see http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {    
      T s0 = vec[0][0] * vec[1][1] - vec[1][0] * vec[0][1];
      T s1 = vec[0][0] * vec[1][2] - vec[1][0] * vec[0][2];
      T s2 = vec[0][0] * vec[1][3] - vec[1][0] * vec[0][3];
      T s3 = vec[0][1] * vec[1][2] - vec[1][1] * vec[0][2];
      T s4 = vec[0][1] * vec[1][3] - vec[1][1] * vec[0][3];
      T s5 = vec[0][2] * vec[1][3] - vec[1][2] * vec[0][3];

      T c5 = vec[2][2] * vec[3][3] - vec[3][2] * vec[2][3];
      T c4 = vec[2][1] * vec[3][3] - vec[3][1] * vec[2][3];
      T c3 = vec[2][1] * vec[3][2] - vec[3][1] * vec[2][2];
      T c2 = vec[2][0] * vec[3][3] - vec[3][0] * vec[2][3];
      T c1 = vec[2][0] * vec[3][2] - vec[3][0] * vec[2][2];
      T c0 = vec[2][0] * vec[3][1] - vec[3][0] * vec[2][1];
      
      // s1/=10;s2/=10;s3/=10;s4/=10;s5/=100;
      // c1/=10;c2/=10;c3/=10;c4/=10;c5/=100;

      return (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0) * (1.0F/24.0F);
      
    }
  };

  template <class T>
  struct SimplexVolumeT<5,5,T> 
  {
  public:
    static const int NDIM=5;
    static const int NDIM_W=5;
  
    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    {    
      printf("SimplexVolumeT<5,5,T> says : 'IMPLEMENT ME !!!!!' \n");
      printf(" ... and crashes miserably.\n");
      exit(-1);
      return -1;
      
    }
  };

  template <class T>
  struct SimplexVolumeT<6,6,T> 
  {
  public:
    static const int NDIM=6;
    static const int NDIM_W=6;
  
    static inline T compute_S(const T (&vec)[NDIM][NDIM_W])
    { 
      
      printf("SimplexVolumeT<6,6,T> says : 'IMPLEMENT ME !!!!!' \n");
      printf(" ... and crashes miserably.\n");
      exit(-1);
      return -1;
      
    }
  };
  */
  template <class T>
  struct SimplexVolumeT<1,2,T>
  {
  public:
    static const int NDIM=1;
    static const int NDIM_W=2;
  
    //WARNING: NO ORIENTATION ! (use hlp::sgn(vec[0][0]) ??)
    static inline T compute_2S(const T (&vec)[NDIM][NDIM_W])
    {
      return (vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]);
    }
    
  };

  template <class T>
  struct SimplexVolumeT<1,3,T>
  {
  public:
    static const int NDIM=1;
    static const int NDIM_W=3;
  
    //WARNING: NO ORIENTATION ! (use hlp::sgn(vec[0][0]) ??)
    static inline T compute_2S(const T (&vec)[NDIM][NDIM_W])
    {
      return (vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]);
    }
    
  };

  template <class T>
  struct SimplexVolumeT<2,3,T>
  {
  public:
    static const int NDIM=2;
    static const int NDIM_W=3;
  
    static inline T compute_2S(const T (&vec)[NDIM][NDIM_W])
    {
      T u11=vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2];
      T u12=vec[0][0]*vec[1][0]+vec[0][1]*vec[1][1]+vec[0][2]*vec[1][2];
      T u22=vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+vec[1][2]*vec[1][2];
      return (u11*u22-u12*u12)*0.5*0.5;      
    }
    
  };

  template <class T>
  struct SimplexVolumeT<2,4,T>
  {
  public:
    static const int NDIM=2;
    static const int NDIM_W=4;

    static inline T compute_2S(const T (&vec)[NDIM][NDIM_W])
    {
      T u11=vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]+vec[0][3]*vec[0][3];
      T u12=vec[0][0]*vec[1][0]+vec[0][1]*vec[1][1]+vec[0][2]*vec[1][2]+vec[0][3]*vec[1][3];
      T u22=vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+vec[1][2]*vec[1][2]+vec[1][3]*vec[1][3];
      return (u11*u22-u12*u12)*0.5*0.5;    
    }
  };

  template <class T>
  struct SimplexVolumeT<3,6,T> 
  {
  public:
    static const int NDIM=3;
    static const int NDIM_W=6;
  
    static inline T compute_2S(const T (&vec)[NDIM][NDIM_W])
    {    
      T u11=vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]+
	vec[0][3]*vec[0][3]+vec[0][4]*vec[0][4]+vec[0][5]*vec[0][5];

      T u12=vec[0][0]*vec[1][0]+vec[0][1]*vec[1][1]+vec[0][2]*vec[1][2]+
	vec[0][3]*vec[1][3]+vec[0][4]*vec[1][4]+vec[0][5]*vec[1][5];

      T u13=vec[0][0]*vec[2][0]+vec[0][1]*vec[2][1]+vec[0][2]*vec[2][2]+
	vec[0][3]*vec[2][3]+vec[0][4]*vec[2][4]+vec[0][5]*vec[2][5];

      T u22=vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+vec[1][2]*vec[1][2]+
	vec[1][3]*vec[1][3]+vec[1][4]*vec[1][4]+vec[1][5]*vec[1][5];

      T u23=vec[1][0]*vec[2][0]+vec[1][1]*vec[2][1]+vec[1][2]*vec[2][2]+
	vec[1][3]*vec[2][3]+vec[1][4]*vec[2][4]+vec[1][5]*vec[2][5];

      T u33=vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+vec[2][2]*vec[2][2]+
	vec[2][3]*vec[2][3]+vec[2][4]*vec[2][4]+vec[2][5]*vec[2][5];

      T det2_12_01 = u22*u33-u23*u23;
      T det2_12_02 = u12*u33-u23*u13;
      T det2_12_12 = u12*u23-u22*u13;

      return (u11 * det2_12_01 -
	      u12 * det2_12_02 +
	      u13 * det2_12_12) * (1.0F/36.0F);
    }
  };

} // namespace simplexVolumeImp

#include "../../internal/namespace.footer"
#endif

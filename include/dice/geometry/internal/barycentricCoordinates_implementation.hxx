#ifndef __BARYCENTRIC_COORDINATES_IMPL_HXX__
#define __BARYCENTRIC_COORDINATES_IMPL_HXX__

#include "../../algebra/inverseMatrix.hxx"

#include "../../internal/namespace.header"

namespace internal {

  template <int ND, typename CT>
  class BarycentricCoordinatesBaseT
  {
  protected:
    static const int NDIM = ND;

    template <class T>
    bool computeCoeffs(const T (&base)[NDIM][NDIM])
    {      
      return (InverseMatrixT<NDIM,T,CT>::compute(base,matrix)==0);    
    }

    template <class T1, class T2>
    void computeLambda(const T1 *vec, T2 *lambda) const
    {
      lambda[1]=matrix[0][0]*vec[0];
      lambda[0]=1.0F;

      for (int i=1;i<=NDIM;i++)
	{
	  lambda[i]=0;
	  for (int j=0;j<NDIM;j++)
	    lambda[i]+=matrix[i][j]*vec[j];
	  lambda[0]-=lambda[i];
	}     
    }

  protected:
    CT matrix[NDIM][NDIM];
  };
  
  template <typename CT>
  class BarycentricCoordinatesBaseT<1,CT>
  {    
  protected:
    static const int NDIM = 1;
    
    template <class T>
    bool computeCoeffs(const T (&base)[NDIM][NDIM])
    {
      return (InverseMatrixT<NDIM,T,CT>::compute(base,matrix)==0);
    }

    template <class T1, class T2>
    void computeLambda(const T1 *vec, T2 *lambda) const
    {
      lambda[1]=matrix[0][0]*vec[0];
      lambda[0]=1.0F-lambda[1];
    }
    
  protected:
    CT matrix[NDIM][NDIM];    
  };

  template <typename CT>
  class BarycentricCoordinatesBaseT<2,CT>
  {
  protected:
    static const int NDIM = 2;

    template <class T>
    bool computeCoeffs(const T (&base)[NDIM][NDIM])
    {      
      return (InverseMatrixT<NDIM,T,CT>::compute(base,matrix)==0);    
    }
    
    template <class T1, class T2>
    void computeLambda(const T1  vec[], T2 lambda[]) const
    {
      lambda[1] = matrix[0][0]*vec[0] + matrix[0][1]*vec[1];
      lambda[2] = matrix[1][0]*vec[0] + matrix[1][1]*vec[1];
      lambda[0] = 1.0F - lambda[1] - lambda[2];
    }
    
  protected:
    CT matrix[NDIM][NDIM];
  };
  
  template <typename CT>
  class BarycentricCoordinatesBaseT<3,CT>
  {
  protected:
    static const int NDIM = 3;
    
    template <class T>
    bool computeCoeffs(const T (&base)[NDIM][NDIM])
    {
      return (InverseMatrixT<NDIM,T,CT>::compute(base,matrix)==0);          
    }
    
    template <class T1, class T2>
    void computeLambda(const T1 *vec, T2 *lambda) const
    {
      lambda[1] = matrix[0][0]*vec[0] + matrix[0][1]*vec[1] + matrix[0][2]*vec[2];
      lambda[2] = matrix[1][0]*vec[0] + matrix[1][1]*vec[1] + matrix[1][2]*vec[2];
      lambda[3] = matrix[2][0]*vec[0] + matrix[2][1]*vec[1] + matrix[2][2]*vec[2];
      lambda[0] = 1.0F - lambda[1] - lambda[2] - lambda[3];
    }
    
  protected:
    CT matrix[NDIM][NDIM]; 
  };
  
}

#include "../../internal/namespace.footer"
#endif

#ifndef __SURFACE_FITTER_IMPLEMENTATION_HXX__
#define __SURFACE_FITTER_IMPLEMENTATION_HXX__

#include "../../internal/namespace.header"

namespace internal {

  template <int ND, int DEG, int NLEFT, int POWER, int W>
  struct FitterGradient
  {
    enum {isNull = ((W==ND)&&(POWER==0)) || ((W<ND)&&(POWER==NLEFT))};
    enum {isDiff = (W==ND)};
  };

  // Counts the number of coefficients in the taylor expansion
  template <int ND, int DEG, int NLEFT, int POWER>
  struct FitterCoefsCount
  {
    enum {value = 
	  FitterCoefsCount<ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::value +
	  FitterCoefsCount<ND,DEG,NLEFT,POWER-1>::value };
  };

  template <int ND, int DEG, int NLEFT>
  struct FitterCoefsCount<ND,DEG,NLEFT,0>
  {
    enum {value = FitterCoefsCount<ND-1,DEG,NLEFT,NLEFT>::value};
  };

  template <int DEG, int NLEFT, int POWER>
  struct FitterCoefsCount<0,DEG,NLEFT,POWER>
  {
    enum {value = 1};
  };

  template <int DEG>
  struct FitterCoefsCount<0,DEG,0,0>
  {
    enum {value = 1};
  };
  
  template <int ND, int DEG>
  struct FitterCoefsCount<ND,DEG,0,0>
  {
    enum {value = 1};
  };

  template <int ND>
  struct FitterCoefsCount<ND,0,0,0>
  {
    enum {value = 1};
  };

  template <int ND, int DEG>
  struct FitterCoefsCountAll
  {
    enum {value = FitterCoefsCount<ND,DEG,DEG,DEG>::value};
  };
  
  // Compute the coefficients of the taylor expansions (x^iy^j... / i!j!...)
  template <typename T, int ND, int DEG, int NLEFT, int POWER>
  struct FitterCoefs
  {
    static void combine(T *coef, std::vector<T> &result)
    {      
      T tmp = result.back();
      result.back() *= coef[POWER-1];
      /*      
      const char * str[5]={"v","w","z","y","x"};
      char tmpS[255];
      sprintf(tmpS,"* (%s=%g)^%2.2d", str[ND],coef[POWER-1], POWER);
      result.back() += std::string(tmpS);   
      */

      FitterCoefs<T,ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::combine(coef+DEG,result);
      result.push_back(tmp);
      
      FitterCoefs<T,ND,DEG,NLEFT,POWER-1>::combine(coef,result);
    }

    static void combine(T *coef, T *result)
    {      
      T tmp = *result;
      (*result) *= coef[POWER-1];

      FitterCoefs<T,ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::combine(coef+DEG,result);
      result += FitterCoefsCount<ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::value;
      (*result)=tmp;
            
      FitterCoefs<T,ND,DEG,NLEFT,POWER-1>::combine(coef,result);
    }
    
    // Compute the coefficients of the (NDIM-W) component of the gradient of the 
    // taylor expansions, with NDIM the (constant) total number of dimension
    template <int W>
    static void combineGrad(T *coef, T *result)
    { 
      T tmp = *result;
      if (FitterGradient<ND,DEG,NLEFT,POWER,W>::isNull)
	{	 
	  std::fill_n(result,FitterCoefsCount<ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::value,0);
	}
      else 
	{
	  if (FitterGradient<ND,DEG,NLEFT,POWER,W>::isDiff)
	    {
	      if (POWER>1) 
		(*result)*=coef[POWER-2]*(POWER);
	    }	    
	  else (*result) *= coef[POWER-1];
	    

	  FitterCoefs<T,ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::
	    template combineGrad<W>(coef,result);
	}
      
      result += FitterCoefsCount<ND-1,DEG,NLEFT-POWER,NLEFT-POWER>::value;
      (*result)=tmp;     

      FitterCoefs<T,ND,DEG,NLEFT,POWER-1>::
	template combineGrad<W>(coef,result);
    }    
  };

  template <typename T, int ND, int DEG, int NLEFT>
  struct FitterCoefs<T,ND,DEG,NLEFT,0>
  {
    static void combine(T *coef, std::vector<T> &result)
    {      
      /*
	 const char * str[5]={"v","w","z","y","x"};
	 char tmpS[255];
	 sprintf(tmpS,"* %s^%2.2d", str[ND], 0);
	 result.back() += std::string(tmpS);   
      */
      FitterCoefs<T,ND-1,DEG,NLEFT,NLEFT>::combine(coef+DEG,result);
    }

    static void combine(T *coef, T *result)
    {
      
      FitterCoefs<T,ND-1,DEG,NLEFT,NLEFT>::combine(coef+DEG,result);
    }

    template <int W>
    static void combineGrad(T *coef, T *result)
    {
      if (FitterGradient<ND,DEG,NLEFT,0,W>::isNull)
	{
	  std::fill_n(result,FitterCoefsCount<ND-1,DEG,NLEFT,NLEFT>::value,0);
	}
      else
	{
	  FitterCoefs<T,ND-1,DEG,NLEFT,NLEFT>::template combineGrad<W>(coef+DEG,result);
	}
    }
  };

  template <typename T, int DEG, int NLEFT, int POWER>
  struct FitterCoefs<T,0,DEG,NLEFT,POWER>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : A\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  template <typename T, int DEG,  int POWER>
  struct FitterCoefs<T,0,DEG,0,POWER>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : B\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  template <typename T, int DEG, int NLEFT>
  struct FitterCoefs<T,0,DEG,NLEFT,0>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : C\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  template <typename T, int DEG>
  struct FitterCoefs<T,0,DEG,0,0>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : D\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  template <typename T, int ND, int DEG, int POWER>
  struct FitterCoefs<T,ND,DEG,0,POWER>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : E\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  template <typename T, int ND, int DEG>
  struct FitterCoefs<T,ND,DEG,0,0>
  {
    static void combine(T *coef, std::vector<T> &result)
    {/*printf("(%ld) : F\n",result.size()-1);*/}
    static void combine(T *coef, T *result){}
    template <int W> static void combineGrad(T *coef, T *result){}
  };

  // Specializations for DEG<=2
  /*
  template <typename T>
  struct FitterCoefs<T,3,2,2,2>
  {
    static void combine(double *coef, std::vector<T> &result)
    {
      result.resize(10);
      result[0]=1;
      
      result[1]=coef[0];// x
      result[2]=coef[2];// y
      result[3]=coef[4];// z    

      result[4]=coef[1];// x2      
      result[5]=coef[3];// y2
      result[6]=coef[5];// z2

      result[7]=coef[0]*coef[2];// xy
      result[8]=coef[0]*coef[4];// xz
      result[9]=coef[2]*coef[4];// yz
    }    
  };

  template <typename T>
  struct FitterCoefs<T,2,2,2,2>
  {
    static void combine(double *coef, std::vector<T> &result)
    {
      result.resize(6);
      result[0]=1;

      result[1]=coef[0];// x
      result[2]=coef[2];// y

      result[3]=coef[1];// x2      
      result[4]=coef[3];// y2
      result[5]=coef[0]*coef[2];// xy     
    }
  };

  template <typename T>
  struct FitterCoefs<T,1,2,2,2>
  {
    static void combine(double *coef, std::vector<T> &result)
    {     
      result.resize(3);
      result[0]=1;

      result[1]=coef[0];// x
      result[2]=coef[1];// x2      
    }
  };

  template <typename T>
  struct FitterCoefs<T,3,1,1,1>
  {
    static void combine(double *coef, std::vector<T> &result)
    {
      result.resize(4);
      result[0]=1;    

      result[1]=coef[0];// x
      result[2]=coef[1];// y
      result[3]=coef[2];// z    
    }
  };
  
  template <typename T>
  struct FitterCoefs<T,2,1,1,1>
  {
    static void combine(double *coef, std::vector<T> &result)
    {
      result.resize(3);
      result[0]=1;

      result[1]=coef[0];// x
      result[2]=coef[1];// y             
    }
  };

  template <typename T>
  struct FitterCoefs<T,1,1,1,1>
  {
    static void combine(double *coef, std::vector<T> &result)
    {
      result.resize(2);
      result[0]=1;

      result[1]=coef[0];// x            
    }
  };
  */

  template <typename T, int ND, int DEG>
  struct FitterCoefsAll
  {
    typedef FitterCoefs<T,ND,DEG,DEG,DEG> FC;
    
    // Compute the coefficients of all the components of the gradient of the 
    // taylor expansions. Should be called with W==NDIM==ND ...
    template <int W=ND>
    static void combineGrad(T *coef, T *result)
    {
      typedef typename hlp::IsTrueT< (W<=0) >::Result Status;
      combineGrad<W>(coef,result,Status());     
    }

    template <int W>
    static void combineGrad(T *coef, T *result, hlp::IsFalse)
    {      
      (*result)=1;
      FC::template combineGrad<W>(coef,result);

      typedef typename hlp::IsTrueT< ((W-1)<=0) >::Result Status;
      combineGrad<W-1>(coef,result+FitterCoefsCountAll<ND,DEG>::value,Status());
    }

    template <int W>
    static void combineGrad(T *coef, T *result, hlp::IsTrue){}
    

    static void combine(T *coef, std::vector<T> &result)
    {
      result.clear();
      result.push_back(1);
      FC::combine(coef,result);
    }

    static void combine(T *coef, T *result)
    {
      (*result)=1;
      FC::combine(coef,result);
    }
  };

  template <int ND, int DEG>
  struct SurfaceFitterHelperT
  {
    static const int NCOEFS = FitterCoefsCountAll<ND,DEG>::value;

    template <class T,class T2>
    static void computeTaylor(T &coords, std::vector<T2> &result)
    {
      T2 coefs[ND*DEG];            
      int k=0;
      
      // Precompute x^i/i! for x,y,z,...
      for (int i=0;i<ND;++i)
	{
	  coefs[k++]=coords[i];
	  for (int j=1;j<DEG;++j,++k)
	    {
	      coefs[k]= (coefs[k-1]*coords[i]) / (j+1);
	    }	    
	}
      
      FitterCoefsAll<T2,ND,DEG>::combine(coefs,result);
      // for (int i=0;i<result.size();++i)
      // 	printf("(%d): %lg\n",i,result[i]);

      /*
      printf("\n");
      printf("COEFFS COUNT : %ld\n",(long)NCOEFS);
      std::vector<std::string> test;
      test.push_back(std::string("1"));
      FitterCoefs<std::string,3,2,2,2>::combine(coefs,test);
      for (int i=0;i<test.size();++i)
	printf("(%d): %s\n",i,test[i].c_str());
      exit(-1);
      */
    }

    template <class T, class T2>
    static void computeTaylor(T &coords, T2 *result)
    {
      T2 coefs[ND*DEG];      
      //int i,j;
      int k=0;
      
      // Precompute x^i/i! for x,y,z,...
      for (int i=0;i<ND;++i)
	{
	  coefs[k++]=coords[i];
	  for (int j=1;j<DEG;++j,++k)
	    {
	      coefs[k]= (coefs[k-1]*coords[i]) / (j+1);
	    }	    
	}      
      
      FitterCoefsAll<T2,ND,DEG>::combine(coefs,result);
      // for (int i=0;i<NCOEFS;++i)
      // 	printf("(%d): %lg\n",i,result[i]);
      /*
      printf("\n");
      printf("COEFFS COUNT : %ld\n",(long)NCOEFS);
      std::vector<std::string> test;
      test.push_back(std::string("1"));
      FitterCoefs<std::string,3,2,2,2>::combine(coefs,test);
      for (int i=0;i<test.size();++i)
	printf("(%d): %s\n",i,test[i].c_str());
      exit(-1);
      */
    }
    
    // put in result the degree of each term in the taylor expansion
    // -> if ith term is of order n, then result[i]=(1<<n)
    template <class T>
    static void computeTaylorDegree(T *result)
    {
      long coefs[ND*DEG];           
      int k=0;
          
      for (int i=0;i<ND;++i)
	{
	  coefs[k++]=2;
	  for (int j=1;j<DEG;++j,++k)
	    coefs[k]=(coefs[k-1]*2);
	}    
  
      FitterCoefsAll<T,ND,DEG>::combine(coefs,result);
    }

    template <class T, class T2>
    static void computeTaylorGradient(T &coords, T2 result[ND][NCOEFS])
    {
      T2 coefs[ND*DEG];           
      int k=0;
      
      // Precompute x^i/i! for x,y,z,...
      for (int i=0;i<ND;++i)
	{
	  coefs[k++]=coords[i];
	  for (int j=1;j<DEG;++j,++k)
	    {
	      coefs[k]= (coefs[k-1]*coords[i]) / (j+1);
	    }	    
	}      

      FitterCoefsAll<T2,ND,DEG>::combineGrad(coefs,&result[0][0]);
    }    
  };

  template <int NDIM, int DEG, typename T, typename T2>
  inline void removeHighOrderCols(const T& A, T& B, const T2* deg)
  {
    static int NC=SurfaceFitterHelperT<NDIM,DEG>::NCOEFS;
    static int NC_LOW=SurfaceFitterHelperT<NDIM,DEG-1>::NCOEFS;
    
    B.resize(A.n_rows,NC_LOW);

    int j=0;
    for (int i=0;i<NC;++i)
      {
	if ((deg[i]&(1<<DEG))==0)
	  B.col(j++)=A.col(i);
      }
  }

  template <int NDIM, int DEG, typename T, typename T2>
  inline void addHighOrderElements(const T& A, T& B,const T2* deg)
  {
    static int NC=SurfaceFitterHelperT<NDIM,DEG>::NCOEFS;
    B.resize(NC);

    int j=0;
    for (int i=0;i<NC;++i)
      {
	if ((deg[i]&(1<<DEG))==0)
	  B[i]=A[j++];
	else
	  B[i]=0;
      }
  }


} // internal
#include "../../internal/namespace.footer"
#endif

#ifndef __REGULAR_GRID_NAVIGATION_HXX__
#define __REGULAR_GRID_NAVIGATION_HXX__

#include <vector>

#include "../internal/namespace.header"

struct RegularGridNavigation {
  typedef int Direction;
  static const Direction UNDEFINED=(1L<<(sizeof(Direction)-1));
  
  static Direction dir(const char *dim,const char *dir)
  {
    int i=0;
    Direction result=0;

    do {
      int d=2*(dim[i]-'0');
      if (dir[i]=='+') result|= 1<<(d);
      else if (dir[i]=='-') result|= 1<<(d+1);
      else break;
      i++;
    } while(true);

    return result;
  }

  static Direction undefined()
  {
    return UNDEFINED;
  }

  static Direction dir()
  {
    return 0;
  }
  
  static Direction dir(const int dim,const int dir)
  {
    if (!dir) return 0;
    return (dir>0)?(1<<(dim*2)):(1<<(dim*2+1));
  }

  static int getDir(Direction dir, int dim)
  {
    if (dir&(1<<(dim*2))) return 1;
    else if (dir&(1<<(dim*2+1))) return -1;
    return 0;
  }
  
  static void setDim(Direction &which, const int dim, const int dir)
  {
    which &= (~((1<<(2*dim))|(1<<(2*dim+1))));
    if (!dir) return;
    if (dir>0) which |=  (1<<(dim*2));
    else which |=  (1<<(dim*2+1));     
  }

  static Direction dir(const std::vector<int> &dim,const std::vector<int> &dir)
  {  
    Direction tmp=0;
    //long i;
    for (unsigned long i=0;i<dir.size();i++)
      {
	if (!dir[i]) continue;
	if (dir[i]>0) tmp |=  (1<<(dim[i]*2));
	else tmp |=  (1<<(dim[i]*2+1));   
      }
   
    return tmp;
  }

  static Direction reverse(Direction dir)
  {
    if (dir==UNDEFINED) return dir;

    Direction result=0;
    Direction test=3;
    //int i;
    //long bef=(long)dir;
    const long m[4]={0,2,1,0};
    for (unsigned long i=0;i<sizeof(Direction)*8;i+=2)
      result|=(m[(dir&(test<<i))>>i]<<i);
    //printf("%ld -> %ld\n",bef,result);
    return result;
  }

};

#include "../internal/namespace.footer"
#endif

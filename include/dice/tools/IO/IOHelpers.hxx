#ifndef __DICE_MY_IO_HELPERS_HXX__
#define __DICE_MY_IO_HELPERS_HXX__

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../../internal/namespace.header"

namespace myIO {  

  inline bool exists(const char *name) 
  {
    FILE *file= fopen(name, "r");
    if (file) 
      {
	fclose(file);
	return true;
      } 
    else return false;
  }

  inline bool exists(const std::string &name)
  {
    return exists(name.c_str());
  }

  inline int rmFile(const char *name) 
  {
    return remove(name);  
  }

  inline int rmFile(const std::string &name)
  {
    return rmFile(name.c_str());
  }

  int makeDir(const char *name)
  {
    return mkdir(name, S_IRWXU|S_IRWXG);
  }

  int makeDir(const std::string &name)
  {
    return makeDir(name.c_str());
  }
  
  // prototype of generic functions for swapping
  int swapI(int);
  float swapF(float);
  double swapD(double);
  void Dswap2B(void*);
  void Dswap4B(void*);
  void Dswap8B(void*);
  void Dswap2BArr(void*,size_t);
  void Dswap4BArr(void*,size_t);
  void Dswap8BArr(void*,size_t);
    
  // read/write helpers
  size_t fread(void *data,size_t size,size_t nb,FILE *f,int swap);

  template <class S>
  size_t writeEnum(FILE *f, typename S::type val)
  {
    int tmpi=(int)val;
    return fwrite(&tmpi,sizeof(int),1,f);
  }

  template <class S>
  typename S::type readEnum(FILE *f, int swap)
  {
    size_t s;
    int tmpi;
    s=fread(&tmpi,sizeof(int),1,f,swap);
    return (typename S::type) tmpi;
  }
  /*
    std::string genFname(const std::string &fname,const std::string &ext, 
    long rank, long size)
    {
    char res[1024];
    if (size>1)
    sprintf(res,"%s_%6.6ld.%s",fname.c_str(),rank,ext.c_str());
    else
    sprintf(res,"%s.%s",fname.c_str(),ext.c_str());

    return std::string(res);
    }
  */


  /*
    int writeTag(FILE *f, const std::string &str, int size=-1)
    {
    if (size==-1) size=str.length();
    char tag[size];
    memset(tag,0,size*sizeof(char));
    strcpy(tag,str.c_str());
    fwrite(&size,sizeof(int),1,f);
    return fwrite(tag,sizeof(char),size,f);
    }
    
    size_t checkTag(FILE *f, const std::string &ref)
    {    
    size_t size=ref.length();
    char tag[size];
    int rs=-1;
    size_t ret = fread(&rs,sizeof(int),1,f);
    
    if ((static_cast<long>(size)!=rs)&&((static_cast<long>(size)!=swapI(rs))))
    {
    fprintf(stderr,"\nERROR reading file: tag '%s' do not match (size)!!!\n",ref.c_str());
    fclose(f);
    exit(-1);
    }

    ret += fread(tag,sizeof(char),size,f);

    for (unsigned long i=0;i<ref.length();i++) 
    {
    if (tag[i]!=ref[i]) 
    {
    fprintf(stderr,"\nERROR reding file: tags '%s' do not match !!!\n",ref.c_str());
    fclose(f);
    exit(-1);
    }
    }    

    return ret;
    //fprintf(stderr,"tag '%s' is found!!!\n",ref.c_str());
    }
  */
  size_t fread(void *data,size_t size,size_t nb,FILE *f,int swap)
  {
  
    size_t ret = fread(data,size,nb,f);
  
    if (swap)
      {
	switch (size)
	  {
	  case 8: Dswap8BArr(data,nb);break;
	  case 4: Dswap4BArr(data,nb);break;
	  case 2: Dswap2BArr(data,nb);break;
	  }
      }
  
    return ret;
  }
  
  /*
  size_t fwrite_dummy(size_t size, size_t nmemb,FILE *stream)
  {
    char buffer[FWRITE_DUMMY_BFR_SIZE];
    unsigned long bsize=FWRITE_DUMMY_BFR_SIZE/size;
    unsigned long nleft=nmemb;
    size_t ret=0;
    
    while(nleft>0)
      {
	unsigned long nw=(nleft<bsize)?nleft:bsize;
	ret+=fwrite(buffer,size,nw,stream);
	nleft-=nw;
      };

    return ret;
  }
  */
  
  int fread_ului(void *data,size_t sizeOut,size_t nb,size_t sizeIn,FILE *f,int swap)
  {
    int ret;
    if (sizeIn==sizeOut)
      return fread(data,sizeIn,nb,f,swap);
    else if (sizeIn<sizeOut)
      {
	long i;
	ret=fread(data,sizeIn,nb,f,swap);
	unsigned long int *uo = (unsigned long int *)(data);
	unsigned int *ui = (unsigned int *)(data);
      
	for (i=nb-1;i>=0;i--)
	  uo[i]=(unsigned long int) ui[i];
      }
    else
      {      
	if (nb==1)
	  {
	    unsigned long int in;
	    unsigned int *uo = (unsigned int *)(data);
	    ret=fread(&in,sizeIn,nb,f,swap);
	    (*uo)= (unsigned int) in;
	    return ret;
	  }
	unsigned long i;
	data=realloc(data,sizeIn*nb);
	ret=fread(data,sizeIn,nb,f,swap);
	unsigned long int *ui = (unsigned long int *)(data);
	unsigned int *uo = (unsigned int *)(data);
	//assert(0);
	for (i=0;i<nb;i++)
	  uo[i]=(unsigned int) ui[i];
	
	data=realloc(data,sizeOut*nb);
	//assert(0);
      }
    return ret;
  }

  void* fread_fd(void *data,size_t sizeOut,size_t nb,size_t sizeIn,FILE *f,int swap)
  {
    //int ret;  
    if (sizeIn==sizeOut)
      {
	fread(data,sizeIn,nb,f,swap);
	return data;
      }
    else if (sizeIn<sizeOut)
      {
	long i;
	fread(data,sizeIn,nb,f,swap);
	double *uo = (double *)(data);
	float *ui = (float *)(data);
      
	for (i=nb-1;i>=0;i--)
	  uo[i]=(double) ui[i];
      }
    else
      {      
	if (nb==1)
	  {
	    double in;
	    float *uo = (float *)(data);
	    fread(&in,sizeIn,nb,f,swap);
	    (*uo)= (float) in;
	    return data;
	  }
	unsigned long i;
	data=realloc(data,sizeIn*nb);
      
	fread(data,sizeIn,nb,f,swap);
	double *ui = (double *)(data);
	float *uo = (float *)(data);
	//assert(0);
	for (i=0;i<nb;i++)
	  uo[i]=(float) ui[i];
	
	data=realloc(data,sizeOut*nb);
      
	//assert(0);
      } 
    return data;
  }


  size_t freadBE(void *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t res;
    static int isLittle=-1;  
    if (isLittle<0)
      {
	int i=1;
	unsigned char *ic=(unsigned char*)&i;
	if (*ic) isLittle=1;
	else isLittle=0;
      }

    res=fread(ptr,size,nmemb,stream);
    if ((isLittle)&&(size>1))
      {	
	unsigned char a[16];
	unsigned char *cptr=(unsigned char*)ptr;
	for (size_t i=0;i<nmemb*size;i+=size)
	  {
	    for (size_t j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	    for (size_t j=0;j<size;j++) cptr[i+j]=a[j];
	  }
      }
    return res;
  }

  size_t fwriteBE(void *ptr, size_t size, size_t nmemb,FILE *stream)
  {
    size_t res;
    static int isLittle=-1;  
    if (isLittle<0)
      {
	int i=1;
	unsigned char *ic=(unsigned char*)&i;
	if (*ic) isLittle=1;
	else isLittle=0;
      }

    if ((isLittle)&&(size>1))
      {	
	unsigned char a[16];
	unsigned char *cptr=(unsigned char*)ptr;
	for (size_t i=0;i<nmemb*size;i+=size)
	  {
	    for (size_t j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	    for (size_t j=0;j<size;j++) cptr[i+j]=a[j];
	  }
      }

    res=fwrite(ptr,size,nmemb,stream);

    if ((isLittle)&&(size>1))
      {
	unsigned char a[16];
	unsigned char *cptr=(unsigned char*)ptr;
	for (size_t i=0;i<nmemb*size;i+=size)
	  {
	    for (size_t j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	    for (size_t j=0;j<size;j++) cptr[i+j]=a[j];
	  }
      }
    return res;
  }

  int swapI(int val)
  {
    int out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
  }

  float swapF(float val)
  {
    float out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
  }

  double swapD(double val)
  {
    double out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[7]=i[0];
    o[6]=i[1];
    o[5]=i[2];
    o[4]=i[3];
    o[3]=i[4];
    o[2]=i[5];
    o[1]=i[6];
    o[0]=i[7];


    return out; 
  }

  inline void Dswap2B(void *val)
  {
    char *c=(char *)val;
    char a;
    
    a=c[0];c[0]=c[1];c[1]=a; 
  }

  inline void Dswap4B(void *val)
  {
    char *c=(char *)val;
    char a;
    
    a=c[0];c[0]=c[3];c[3]=a;
    a=c[1];c[1]=c[2];c[2]=a; 
 
  }

  inline void Dswap8B(void *val)
  {
    char *c=(char *)val;
    char a;
    
    a=c[0];c[0]=c[7];c[7]=a;
    a=c[1];c[1]=c[6];c[6]=a;
    a=c[2];c[2]=c[5];c[5]=a;
    a=c[3];c[3]=c[4];c[4]=a;
  }

  void Dswap2BArr(void *val,size_t n)
  {
    size_t i;
    char a;

    char *c=(char *)val;

    for (i=0;i<2*n;i+=2)
      {
	a=c[i];
	c[i]=c[i+1];
	c[i+1]=a;
      }

  }

  void Dswap4BArr(void *val,size_t n)
  {
    size_t i;
    char a,b;

    char *c=(char *)val;

    for (i=0;i<4*n;i+=4)
      {
	a=c[i];
	b=c[i+1];
	c[i]=c[i+3];
	c[i+1]=c[i+2];
	c[i+2]=b;
	c[i+3]=a;
      }

  }

  void Dswap8BArr(void *val,size_t n)
  {
    size_t i;
    char a,b,u,v;

    char *c=(char *)val;

    for (i=0;i<8*n;i+=8)
      {
	a=c[i];
	b=c[i+1];
	u=c[i+2];
	v=c[i+3];
	c[i]=c[i+7];
	c[i+1]=c[i+6];
	c[i+2]=c[i+5];
	c[i+3]=c[i+4];
	c[i+4]=v;
	c[i+5]=u;
	c[i+6]=b;
	c[i+7]=a;
      }
  }

  ssize_t
  getDelim (char **lineptr, int *n, int delim, FILE *stream)
  {
    static const int line_size=1024; // Default value for line length.
    int indx = 0;
    int c;
        
    /* Sanity checks.  */
    if (lineptr == NULL || n == NULL || stream == NULL)
      return -1;

    /* Allocate the line the first time.  */
    if (*lineptr == NULL)
      {
	*lineptr = (char *)malloc (line_size);
	if (*lineptr == NULL)
	  return -1;
	*n = line_size;
      }

    /* Clear the line.  */
    memset (*lineptr, '\0', *n);

    while ((c = getc (stream)) != EOF)
      {
	/* Check if more memory is needed.  */
	if (indx >= *n)
	  {
	    *lineptr = (char *)realloc (*lineptr, *n + line_size);
	    if (*lineptr == NULL)
	      {
		return -1;
	      }
	    /* Clear the rest of the line.  */
	    memset(*lineptr + *n, '\0', line_size);
	    *n += line_size;
	  }

	/* Push the result in the line.  */
	(*lineptr)[indx++] = c;

	/* Bail out.  */
	if (c == delim)
	  {
	    break;
	  }
      }
    return (c == EOF) ? -1 : indx;
  }

  ssize_t getLine (char **lineptr, int *n, FILE *stream)
  {
    return getDelim (lineptr, n, '\n', stream);
  }
  
} // myIO

#include "../../internal/namespace.footer"
#endif


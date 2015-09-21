#ifndef __MY_IO_HXX__
#define __MY_IO_HXX__

// Buffered binary read/write with automatic swapping and additional 
// convenience fonctions to serialize objects

#include <stdio.h>
#include <string.h>

#include "../../tools/helpers/helpers.hxx"
#include "IOHelpers.hxx"

#include "../../internal/namespace.header"

namespace myIO {
  /*
  namespace { // anonymous namespace

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
    
  } // end of anonymous namespace 
  */

  static const int endiannessTestValue = 0xabcdef;
  static const unsigned char typeSizeCheckVal[9]={(unsigned char)sizeof(bool),
						  (unsigned char)sizeof(char),
						  (unsigned char)sizeof(short),
						  (unsigned char)sizeof(int),
						  (unsigned char)sizeof(long),
						  (unsigned char)sizeof(long long),
						  (unsigned char)sizeof(float),
						  (unsigned char)sizeof(double),
						  (unsigned char)sizeof(long double)};
  static const char *typeSizeCheckType[] = {"bool",
					    "char",
					    "short",
					    "int",
					    "long",
					    "long long",
					    "float",
					    "double",
					    "long double"};    

  // Default buffer size : 1MB
  template <long bSize=(1L<<20)>
  //template <long bSize>
  class BinaryWriterT
  {
  public:
    static const long bufferSize = bSize;
    static const long paddingSize = 16;

    BinaryWriterT(FILE *f_):
      f(f_),
      curPos(0)
    {	
      // Align the buffer to paddingSize bytes
      memset(buffer_ref,0,bSize+paddingSize);
      buffer = &buffer_ref[0];
      while ((unsigned long)buffer % paddingSize != 0) 
	buffer++;
      //memset(buffer,0,bSize);
      ownFile=false;
    }

    BinaryWriterT(const std::string &fname):      
      curPos(0)
    {	
      // Align the buffer to paddingSize bytes
      memset(buffer_ref,0,bSize+paddingSize);
      buffer = &buffer_ref[0];
      while ((unsigned long)buffer % paddingSize != 0) 
	buffer++;
      //memset(buffer,0,bSize);
      ownFile=true;      
      f=fopen(fname.c_str(),"w");
      if (f==NULL)
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>
	      ("Opening file %s for writing.\n",fname.c_str());
	  /*
	  fprintf(stderr,
		  "FATAL ERROR: could not open file '%s' for writing !\n",fname.c_str());
	  fprintf(stdout,
		  "FATAL ERROR: could not open file '%s' for writing !\n",fname.c_str());
	  */
	  exit(-1);
	}     
    }

    ~BinaryWriterT()
    {
      flush();
      if (ownFile) 
	fclose(f);
    }

    void writeEndiannessCheck(int value = endiannessTestValue)
    {
      write(&value);
    }

    void writeTypeSizeCheck()
    {    
      write(typeSizeCheckVal,9);
    }

    void writeHeader(const std::string &str, float ver=0.0)
    {
      this->writeEndiannessCheck();
      this->write(&str);
      this->write(&ver);
    }

    template <typename TS>
    size_t write(TS* src, size_t N)
    {
      return writeAs<TS,TS>(src,N,typename hlp::SameType<TS,TS>::Result());
    }

    template <typename TS>
    size_t write(TS* src)
    {    
      typedef typename hlp::SameType<TS,std::string>::Result IsString;
      typedef typename hlp::IsVector<TS>::Result IsVector;
      //typedef typename hlp::IsPrimitiveType<TS>::Result IsPrimitiveType;
      typedef hlp::TypePair<IsString,IsVector> Specialization;
      /*
      typedef typename hlp::IsTrueT<(!IsString::value)&&
				    (!IsVector::value)&&
				    (!IsPrimitiveType::value)>::Result TypeError;
      writeTypeError(TypeError());
      */
      return writeOne(src, Specialization());
    }
    
    // Write data after casting to a different type (TS -> TD)
    // TS/TD must be primitive types if TD != TS.
    template <typename TD,typename TS>
    size_t writeAs(const TS* src, size_t N)
    {
      return writeAs<TD,TS>(src,N,typename hlp::SameType<TD,TS>::Result());
    }  

    template <typename TD,typename TS>
    size_t writeAs(const TS* src)
    {      
      return writeAs<TD,TS>(src,typename hlp::SameType<TD,TS>::Result());
    }
 
    // Flush the buffer to file 
    size_t flush()
    {
      size_t res=0;
      if (curPos>0) 
	{
	  //printf("flushing %ld bytes !\n",curPos);
	  res=fwrite(buffer,sizeof(char),curPos,f);    
	  curPos=0;
	}
      return res;
    }

    FILE *getFilePtr()
    {
      flush();
      return f;
    }

  private:     
    // hlp::IsTrue version is not defined so that we get a compile error ...
    void writeTypeError(hlp::IsFalse) {}
    
    template <class TS>
    size_t writeOne(TS *src, hlp::TypePair<hlp::IsTrue,hlp::IsFalse>) // DT == string 
    {      
      int size=src->size();
      write(&size);
      return write(src->c_str(),size+1)+sizeof(int);
    }

    template <class TS>
    size_t writeOne(TS *src, hlp::TypePair<hlp::IsFalse,hlp::IsTrue>) // DT == vector
    {            
      unsigned long size = src->size(); 
      write(&size);
      return write(&src->front(),size)+sizeof(unsigned long);
    }

    template <class TS>
    size_t writeOne(TS *src, hlp::TypePair<hlp::IsFalse,hlp::IsFalse>) // DT == other
    {
      return writeAs<TS,TS>(src,typename hlp::SameType<TS,TS>::Result());      
    }
  
    // TD and TS are different (TD : datatype we want to write)
    // type MUST be primitive
    template <typename TD,typename TS>
    size_t writeAs(const TS* src, size_t N, hlp::IsFalse)
    {
      /*
      for (size_t i=0;i<N;++i)
	{
	  TD tmp=static_cast<TD>(src[i]);
	  writeAs<TD,TD>(&tmp,hlp::IsTrue);
	}
      */
      
      unsigned long delta =sizeof(TD)*N;
      if (curPos+delta>bSize) 
	{
	  flush();
	  if (delta > bSize)
	    {
	      if (sizeof(TD) > bSize)
		{
		  for (size_t i=0;i<N;i++) 
		    {
		      TD res = static_cast<TD>(*src);
		      fwrite(&res,sizeof(TD),1,f);    
		    }
		  return N;
		}
	      else
		{
		  unsigned long flushAt=sizeof(TD) * (unsigned long)(bSize/sizeof(TD));
		  TD *dest = static_cast<TD*>(&buffer[curPos]);
		  for (size_t i=0;i<N;i++) 
		    {
		      (*dest)=static_cast<TD>(src[i]);
		      ++dest;
		      curPos+=sizeof(TD);
		      if (curPos == flushAt) flush(); 
		    }	      
		  return N;
		}
	    }	    
	}

      TD *dest = static_cast<TD*>(&buffer[curPos]);
      for (size_t i=0;i<N;i++) 
	{
	  (*dest)=static_cast<TD>(src[i]);
	  ++dest;
	}
      curPos+=delta;
      return N;
      
    }

    // TD and TS are different (TD : datatype we want to write)
    // type MUST be primitive
    template <typename TD,typename TS>
    size_t writeAs(const TS* src, hlp::IsFalse)
    {   
      TD tmp= static_cast<TD>(*src);
      return writeAs<TD,TD>(&tmp,hlp::IsTrue());
      /*
      if (curPos+sizeof(TD)>bufferSize)
	{
	  flush();
	  if (sizeof(TD) > bufferSize)
	    {
	      TD res = static_cast<TD>(*src);
	      fwrite(&res,sizeof(TD),1,f);    
	      return 1;
	    }
	}
	      
      TD *dest = static_cast<TD *>(static_cast<void*>(&buffer[curPos]));
      (*dest)=static_cast<TD>(*src);
      curPos+=sizeof(TD);
      
      return 1;
      */
    }

    // TD and TS are identical
    template <typename TD,typename TS>
    size_t writeAs(const TS* src, size_t N, hlp::IsTrue)
    {     
      typedef typename hlp::SameType<TS,std::string>::Result IsString;
      typedef typename hlp::IsVector<TS>::Result IsVector;
      
      typedef typename hlp::IsTrueT<(!IsString::value)&&
				    (!IsVector::value)>::Result TypeAsPrimitive;

      return writeSameType(src,N,TypeAsPrimitive());
    }

    // TD and TS are identical
    template <typename TD,typename TS>
    size_t writeAs(const TS* src, hlp::IsTrue)
    {
      typedef typename hlp::SameType<TS,std::string>::Result IsString;
      typedef typename hlp::IsVector<TS>::Result IsVector;
      
      typedef typename hlp::IsTrueT<(!IsString::value)&&
				    (!IsVector::value)>::Result TypeAsPrimitive;
      
      return writeSameType(src,TypeAsPrimitive());
    }
    
    // NON primitive type
    template <typename TS>
    size_t writeSameType(const TS* src, size_t N, hlp::IsFalse)
    {
      //printf("Writing %ld NON primitive types.\n",N);
      size_t count=0;
      for (unsigned long i=0;i<N;++i) count+=write(&src[i]);
      return count;
    }

    // NON primitive type
    template <typename TS>
    size_t writeSameType(const TS* src, hlp::IsFalse)
    {
      return write(src);
    }  
    
    template <typename TS>
    size_t writeSameType(const TS* src, size_t N, hlp::IsTrue)
    {
      unsigned long delta =sizeof(TS)*N;
      // glb::console->print<LOG_STD_ALL>("writing %ld*%ld(=%ldB) el to buffer @cur=%ld / %ld/%ld(%ld).\n",N,sizeof(TS),delta,curPos,bufferSize,bSize,(long)(curPos+delta>bufferSize));
      if (curPos+delta>bufferSize) 
	{
	  flush();
	  if (delta > bufferSize)
	    {
	      fwrite(src,sizeof(TS),N,f); 
	      return N;
	    } 	  
	}
      
      void * __restrict dest_ = buffer+curPos;
      const void * __restrict src_ = src;
      // printf("%ld -> %ld (+delta=%ld)  ref=%ld/%ld)\n",(long)src_,(long)dest_,delta,
      // 	     (long)buffer,buffer_ref);     
      memcpy(dest_,src_,delta);
      curPos+=delta;
      
      return N;
    }

    template <typename TS>
    size_t writeSameType(const TS* src, hlp::IsTrue)
    {     
      // glb::console->print<LOG_STD_ALL>("writing 1: %ld(=%ldB) el to buffer @cur=%ld / %ld(%ld).\n",sizeof(TS),sizeof(TS),curPos,bSize,(long)(curPos+sizeof(TS)>bSize));
      if (curPos+sizeof(TS)>bufferSize)
	{
	  flush();
	  if (sizeof(TS) > bufferSize)
	    {
	      fwrite(src,sizeof(TS),1,f); 
	      return 1;
	    } 
	}
      void * __restrict dest_ = buffer+curPos;
      const void * __restrict src_ = src;
      memcpy(dest_,src_,sizeof(TS));
      curPos+=sizeof(TS);
      return 1;
    }

  private:
    FILE *f; 
    unsigned long curPos; // current position in the buffer  
    char *buffer;
    char buffer_ref[bSize + paddingSize]; // 'paddingSize' spare bytes for data alignement
    bool ownFile;
  };

  template <long bSize=(1L<<20)>
  class BinaryReaderT
  {
  public:
    static const long bufferSize = bSize;   
    static const long paddingSize = 16;
    static const long minBufferSize = 64; // unbuffered read below this threshold   
    typedef BinaryReaderT<bSize> MyType;    
    enum CheckResult {SUCCESS=0, WARNING=1, WARNING_VERSION=2, FAILURE=3, 
		      FAILURE_HEADER=4, FAILURE_ENDIANNESS=5};

    BinaryReaderT(FILE *f_, bool swap_=false, bool autoSwap_=true):
      f(f_),
      curPos(0),      
      needSwap(swap_),
      autoSwap(autoSwap_)
    {
      // Align the buffer to 8 bytes
      memset(buffer_ref,0,bSize+paddingSize);
      buffer=&buffer_ref[0];
      while (((unsigned long)buffer) % paddingSize != 0) 
	buffer++;
      ownFile=false;
      if (bufferSize>=minBufferSize) dummy=fread(buffer,sizeof(char),bufferSize,f);
    }

    BinaryReaderT(const std::string &fname, bool swap_=false, bool autoSwap_=true):      
      curPos(0),    
      needSwap(swap_),
      autoSwap(autoSwap_)
    {
      // Align the buffer to 8 bytes
      memset(buffer_ref,0,bSize+paddingSize);
      buffer=&buffer_ref[0];
      while (((unsigned long)buffer) % paddingSize != 0) 
	buffer++;
      ownFile=true;
      f=fopen(fname.c_str(),"r");
      if (f==NULL)
	{
	  PRINT_SRC_INFO(LOG_ERROR);	  
	  glb::console->print<LOG_ERROR>
	      ("Opening file %s for reading.\n",fname.c_str());
	  /*
	  fprintf(stderr,
		  "FATAL ERROR: could not open file '%s' for reading !\n",fname.c_str());
	  fprintf(stdout,
		  "FATAL ERROR: could not open file '%s' for reading !\n",fname.c_str());
	  */
	  exit(-1);
	}   
      else if (bufferSize>=minBufferSize) dummy=fread(buffer,sizeof(char),bufferSize,f);
    }

    ~BinaryReaderT()
    {
      if (ownFile) fclose(f);
    }
  
    void setNeedSwap(bool b)
    {
      needSwap=b;
    }

    void setAutoSwap(bool b)
    {
      autoSwap=b;
    }

    bool getNeedSwap() const
    {
      return needSwap;
    }

    bool getAutoSwap() const
    {
      return autoSwap;
    }

    static MyType *nullReader()
    {
      return static_cast<MyType*>(NULL);
    }
    
    CheckResult checkEndianness(int value = endiannessTestValue)
    {
      int rv=0;
      this->read(&rv);
      if (rv != value)
	{
	  if (swapI(rv) == value)
	    needSwap = !needSwap;
	  else return FAILURE_ENDIANNESS; // Error
	}
      return SUCCESS;
    }

    template <class LOG_F,class C>
    int checkTypeSize(C *console)
    {
      int result=0;
      unsigned char val[9];
      read(val,9);
      
      for (int i=0;i<9;++i)
	{
	  if (typeSizeCheckVal[i] != val[i])
	    {
	      console->
		template print<LOG_F>("Type '%s' stored in binary file has size %d while its size is %d on current architecture.\n",typeSizeCheckType[i],(int)typeSizeCheckVal[i],(int)val[i]);
	      result |= (1<<i);
	    }
	}
      return result;
    }

    template <class LOG_F,class LOG_W,class C>
    static void reportCheckHeaderFailed(C *c, 
					const std::string &refHeader, float refVersion, 
					const std::string &header, float version,
					CheckResult result,
					bool exitOnFailure=false)
    {
      if (result == FAILURE_HEADER)
	{
	  c->template printNewLine<LOG_W>
	    ("Invalid file format.\n");
	  c->template printNewLine<LOG_F>
	    ("In BinaryReader: headers do not match.\n");
	  c->template print<LOG_F>
	    ("Current header : %s v%.2f.\n",refHeader.c_str(),refVersion);
	  c->template print<LOG_F>
	    ("File header    : %s v%.2f.\n",header.c_str(),version);	  
	  if (exitOnFailure)
	    {
	      PRINT_SRC_INFO(LOG_F);
	      exit(-1);
	    }
	}
      if (result == FAILURE_ENDIANNESS)
	{
	  c->template printNewLine<LOG_W>
	    ("Invalid file format.\n");
	  c->template printNewLine<LOG_W>
	    ("In BinaryReader: endianness test failed before reading header.\n");
	  c->template print<LOG_W>
	    ("Current header : %s v%.2f.\n",refHeader.c_str(),refVersion);
	  c->template print<LOG_W>
	    ("File header    : %s v%.2f.\n",header.c_str(),version);	 
	  if (exitOnFailure)
	    {
	      PRINT_SRC_INFO(LOG_F);
	      exit(-1);
	    }
	}
      if (result==WARNING_VERSION)
	{
	  c->template printNewLine<LOG_W>
	    ("In BinaryReader: class versions do not match.\n");
	  c->template print<LOG_W>
	    ("Current header : %s v%.2f.\n",refHeader.c_str(),refVersion);
	  c->template print<LOG_W>
	    ("File header    : %s v%.2f.\n",header.c_str(),version);
	}
    }
    
    CheckResult checkHeader(const std::string &refHeader, std::string &readHeader)
    {
      float ver;
      return checkHeader(refHeader,readHeader,&ver);
    }
    
    CheckResult checkHeader(const std::string &refHeader, float refVer,  
			    std::string &readHeader, float &readVer)
			    
    {
      if (this->checkEndianness() >= FAILURE) return FAILURE_ENDIANNESS;
      this->read(&readHeader);      
      this->read(&readVer);
      if (readHeader != refHeader) return FAILURE_HEADER;
      if (readVer != refVer) return WARNING_VERSION;
      return SUCCESS;
    }
    /*
    template <class LOG_F, class LOG_W, class T,class C, class R>
    static CheckResult checkHeaderAndReport(C *console, R *reader,
					    float &readVersion)
    //,bool exitOnFailure=true)
    {
      bool exitOnFailure=true;
      return checkHeaderAndReport<LOG_F,LOG_W,T,C,R>
	(console,reader,readVersion,-1);//,exitOnFailure);
    }
*/
    template <class LOG_F, class LOG_W, class T,class C, class R>
    static CheckResult checkHeaderAndReport(C *console, R *reader,
					    float &readVersion,
					    bool exitOnFailure=true)
    {    
      if (reader==NULL)
	{
	  readVersion = T::classVersion();
	  return SUCCESS;
	}
      std::string readHeader;
      CheckResult result = 
	reader->checkHeader(T::classHeader(),T::classVersion(),
			    readHeader,readVersion);

      if (result != CheckResult::SUCCESS)
	{
	  reader->template reportCheckHeaderFailed<LOG_F,LOG_W>
	    (console,T::classHeader(),T::classVersion(),
	     readHeader,readVersion,result,exitOnFailure);

	  if (readVersion<T::compatibleSinceClassVersion())
	    {
	       console->template printNewLine<LOG_F>
		 ("Incompatible versions.\n");
	       if (exitOnFailure) 
		 {
		   PRINT_SRC_INFO(LOG_F);
		   exit(-1);
		 }
	    }
	  else
	    {
	      console->template printNewLine<LOG_W>
		 ("Versions differ but seem compatible, proceeding ...\n");
	    }
	}

      return result;      
    }

    template <class DT>
    size_t read(DT *val)
    {      
      typedef typename hlp::SameType<DT,std::string>::Result IsString;
      typedef typename hlp::IsVector<DT>::Result IsVector;
      typedef hlp::TypePair<IsString,IsVector> Specialization;

      return readOne(val,Specialization());
    }

    template <class DT>
    size_t read(DT *val, long N)
    {
      typedef typename hlp::IsTrueT<(bufferSize >= minBufferSize)>::Result UseBuffer;
      typedef typename hlp::SameType<DT,std::string>::Result IsString;
      typedef typename hlp::IsVector<DT>::Result IsVector;
      typedef typename hlp::IsTrueT<IsString::value || IsVector::value>::Result IsSpecial;
      typedef hlp::TypePair<UseBuffer,IsSpecial> Specialization;

      //typedef typename hlp::IsPrimitiveType<DT>::Result IsPrimitiveType;
      //typedef hlp::TypePair<UseBuffer,IsPrimitiveType> Specialization;
      //typedef hlp::TypeTriplet<UseBuffer,IsString,IsVector> Specialization;
      
      return read(val,N,Specialization());
    }

    template <class DT>
    static bool canAutoSwap()
    {
      return hlp::IsPrimitiveType<DT>::Result::value;
    }

    template <class DT>
    bool willAutoSwap()
    {
      return autoSwap && canAutoSwap<DT>();
    }


    template <class DT>
    static bool swap(DT *val)
    {
      typedef typename hlp::IF_<
	hlp::IsPrimitiveType<DT>::Result::value,
	hlp::ConstantValue<sizeof(DT)>,
	hlp::ConstantValue<0> >::Result SwapSize;
      return swap(val,SwapSize());
    }

    template <class DT>
    static bool swap(DT *val, long N)
    {
      typedef typename hlp::IF_<
	hlp::IsPrimitiveType<DT>::Result::value,
	hlp::ConstantValue<sizeof(DT)>,
	hlp::ConstantValue<0> >::Result SwapSize;

      if (N==1)
	return swap(val,SwapSize());
      else
	return swap(val,N,SwapSize());
    }  

  private:
    template <class DT>
    size_t readOne(DT *val, hlp::TypePair<hlp::IsTrue,hlp::IsFalse>) // DT == string
    {     
      int size;
      size_t N=read(&size)*sizeof(int);
      char buf[size+1];
      N+=read(buf,size+1);
      val->assign(buf);  
      //printf("Reading string '%s'.\n",val->c_str());
      return N;
    }

    template <class DT>
    size_t readOne(DT *val, hlp::TypePair<hlp::IsFalse,hlp::IsTrue>) // DT == vector
    {
      unsigned long size;
      read(&size);
      val->resize(size);
      return read(&val->front(),size);
    }

    template <class DT>
    size_t readOne(DT *val, hlp::TypePair<hlp::IsFalse,hlp::IsFalse>) // DT == other
    {
      return read(val,1);
    }

    // buffered + non special
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsTrue,hlp::IsFalse>) 
    {
      return readFromBuf(val,N);
    }

    // non buffered + non special
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsFalse,hlp::IsFalse>) 
    {
      size_t Nr=fread(val,sizeof(DT),N,f);
      if (needSwap&&autoSwap) swap(val,N);
      return Nr;
    }

    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsFalse,hlp::IsTrue>) 
    {
      return read(val,N,hlp::TypePair<hlp::IsTrue,hlp::IsTrue>());
    }

    // vector or string
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsTrue,hlp::IsTrue>) 
    {
      size_t count=0;
      for (long i=0;i<N;++i) count+=read(&val[i]);
      return count;
    }

    /*
   
    // buffered + NON primitive type
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsFalse,hlp::IsFalse>) 
    {
      // Specialization is commented here ...
      // FIXME: We could use it to swap non primitive type 
      return read(val,N,hlp::TypePair<hlp::IsFalse,hlp::IsTrue>());
      // size_t sum=0;
      // for (long i=0;i<N;++i) sum+=read(&val[i]);
      // return sum;
    }

    // buffered + primitive type
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsFalse,hlp::IsTrue>)
    {
      return readFromBuf(val,N);
    }

    // NON buffered + NON primitive type
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsTrue,hlp::IsFalse>)
    {
      // Specialization is commented here ...
      // FIXME: We could use it to swap non primitive type 
      return read(val,N,hlp::TypePair<hlp::IsTrue,hlp::IsTrue>());
      //return read(val,N,hlp::TypePair<hlp::IsFalse,hlp::IsFalse>());
    }

    // NON buffered + primitive type
    template <class DT>
    size_t read(DT *val, long N, hlp::TypePair<hlp::IsTrue,hlp::IsTrue>)
    {
      size_t Nr=fread(val,sizeof(DT),N,f);
      if (needSwap&&autoSwap) swap(val,N);
      return Nr;
    }
    */
    template <class DT>
    size_t readFromBuf(DT *val, size_t N)
    {
      unsigned long delta = N*sizeof(DT);
      if (curPos+delta > bufferSize)
	{
	  if (delta >= bufferSize) // won't fit in the buffer !
	    {
	      void *val2 = val;	      
	      memcpy(val2,buffer+curPos,bufferSize-curPos);
	      curPos=(bufferSize-curPos);
	      delta-=curPos;
	      dummy=fread(((char*)val2) + curPos,sizeof(char),delta,f);
	      dummy=fread(buffer,sizeof(char),bufferSize,f);
	      curPos=0;
	      if (needSwap&&autoSwap) swap(val,N);
	      return N;
	    }
	  else
	    {
	      void *val2 = (char*)val;
	      memcpy(val2,buffer+curPos,bufferSize-curPos);
	      curPos=(bufferSize-curPos);
	      delta-=curPos;	      
	      dummy=fread(buffer,sizeof(char),bufferSize,f);
	      memcpy(((char*)val2) + curPos,buffer,delta);
	      curPos=delta;
	      if (needSwap&&autoSwap) swap(val,N);
	      return N;

	      // std::copy(buffer+curPos,buffer+bufferSize,buffer);
	      // curPos=bufferSize-curPos;
	      // fread(buffer+curPos,sizeof(char),bufferSize-curPos,f);
	      // curPos=0;
	    }
	}
      
      memcpy(val,&buffer[curPos],delta);
      curPos += delta;
      if (needSwap&&autoSwap) swap(val,N);     
      return N;
    }

  private:

    // the first swap catches any non primitive types
    template <class DT>
    static bool swap(DT *val, long N, hlp::ConstantValue<0>)
    {return false;}
    template <class DT>
    static bool swap(DT *val, long N, hlp::ConstantValue<1>)
    {return true;}
    template <class DT>
    static bool swap(DT *val, long N, hlp::ConstantValue<2>)
    {Dswap2BArr(val,N);return true;}
    template <class DT>
    static bool swap(DT *val, long N, hlp::ConstantValue<4>)
    {Dswap4BArr(val,N);return true;}
    template <class DT>
    static bool swap(DT *val, long N, hlp::ConstantValue<8>)
    {Dswap8BArr(val,N);return true;}
    
    template <class DT>
    static bool swap(DT *val, hlp::ConstantValue<0>)
    {return false;}
    template <class DT>
    static bool swap(DT *val, hlp::ConstantValue<1>)
    {return true;}
    template <class DT>
    static bool swap(DT *val, hlp::ConstantValue<2>)
    {Dswap2B(val);return true;}
    template <class DT>
    static bool swap(DT *val, hlp::ConstantValue<4>)
    {Dswap4B(val);return true;}
    template <class DT>
    static bool swap(DT *val, hlp::ConstantValue<8>)
    {Dswap8B(val);return true;}

  private:    
    FILE *f;
    unsigned long curPos;
    char *buffer; 
    bool needSwap;
    bool autoSwap;
    bool ownFile;
    char buffer_ref[bSize + paddingSize]; // 8 is for data alignement  
    size_t dummy;
  };
  

  /*
    template <long bSize> 
    class BufferedWrite
    {
    public:
    BufferedWrite(FILE *f_=NULL):      
    buffer((unsigned char*)bufferD),
    curID(0)
    {
    setFile(f_);
    }

    void setFile(FILE *f_)
    {      
    f=f_;
    }
    
    // TD : datatype we want to write
    template <typename TD,typename TS>
    int write(TS* src, size_t N)
    {
    unsigned long delta =sizeof(TD)*N;
    // check if (delta>bSize)      
    if (curID+delta>bSize) flush<TD>();

    TD *dest = (TD*)(&buffer[curID]);
    for (size_t i=0;i<N;i++) 
    {
    (*dest)=(TD)src[i];
    ++dest;
    }
    curID+=delta;
    return N;
    }

    // TD : datatype we want to write
    template <typename TD,typename TS>
    int write(TS* src)
    {     
    if (curID+sizeof(TD)>bSize) flush<TD>();
    TD *dest = static_cast<TD *>(static_cast<void*>(&buffer[curID]));
    (*dest)=(TD)(*src);
    curID+=sizeof(TD);
    return 1;
    }

    void reset()
    {
    curID=0;
    }
    
    template <typename TD>
    void flush()
    {
    if (curID>0) 
    {
    //printf("size : %ld %ld\n",sizeof(TD),curID/sizeof(TD));
    fwrite(buffer,sizeof(TD),curID/sizeof(TD),f);    
    }
    curID=0;
    }

    private:
    double bufferD[(bSize/sizeof(double))+1]; // for memory alinement on 8 bytes
    unsigned char *buffer;
    unsigned long curID;
    FILE *f;
    };
  */
};

#include "../../internal/namespace.footer"
#endif

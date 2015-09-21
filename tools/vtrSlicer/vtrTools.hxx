#ifndef VTR_TOOLS_HXX__
#define VTR_TOOLS_HXX__

#include <dice/tools/helpers/helpers.hxx>

namespace vtrSlicer {

  const char *cutName(const char *Name)
  {
    int i,j; 

    for (i=0,j=-1;i<strlen(Name);i++) 
      if (Name[i]=='/') j=i;

    if (j!=-1)
      return &Name[j+1];
    else
      return Name;
  }


  struct Header
  {
    std::vector<std::string> original;
    float version;
    int type;
    int wExtent[NDIM*2];
    int pExtent[NDIM*2];
    long offset[NDIM];

    Header()
    {}

    long getWDim(int index)
    {
      return wExtent[2*index+1]-wExtent[2*index];
    }

    long getPDim(int index)
    {
      return pExtent[2*index+1]-pExtent[2*index];
    }
    
    void setDim(int index, long value)
    {
      if (value<1) value=1;
      wExtent[2*index+1]=wExtent[2*index]+value;
      pExtent[2*index+1]=pExtent[2*index]+value;
    }

    void read(const std::string &fname)
    {   
      std::ifstream stream;
      read(fname,stream);
    }

    std::ifstream &read(const std::string &fname, std::ifstream &stream)
    {
      stream.close();
      stream.open(fname);
      original.clear();
      std::string line;    
      size_t pos=0;

      if (stream.is_open())
	{
	  getline(stream,line);original.push_back(line);
	  pos=check(fname,line,"version");	
	  pos+=line.substr(pos).find("\"")+1;
	
	  version=std::stod(line.substr(pos));	
	  if (version==1.0)
	    {	    
	      pos=line.find("header_type");
	      if (pos!=std::string::npos)
		{
		  if (line.substr(pos).find("UInt64")!=std::string::npos)
		    type=8;
		  else type=4;		  
		}
	      else type=4;
	    }

	  getline (stream,line);original.push_back(line);	
	  pos=check(fname,line,"WholeExtent");	
	  pos+=line.substr(pos).find("\"")+1;
	  for (int i=0;i<NDIM*2;++i)
	    {
	      size_t len;
	      wExtent[i]=std::stoi(line.substr(pos),&len);
	      pos+=len;
	    }

	  getline (stream,line);original.push_back(line);	
	  pos=check(fname,line,"Piece Extent");	
	  pos+=line.substr(pos).find("\"")+1;	
	  for (int i=0;i<NDIM*2;++i)
	    {
	      size_t len;
	      pExtent[i]=std::stoi(line.substr(pos),&len);
	      pos+=len;
	    }
	
	  for (int i=0;i<6;++i)
	    {
	      getline(stream,line);
	      original.push_back(line);
	    }
	  for (int i=0;i<NDIM;++i)
	    {
	      getline(stream,line);original.push_back(line);
	      pos=check(fname,line,"offset");
	      pos+=line.substr(pos).find("\"")+1;
	      offset[i]=std::stol(line.substr(pos),&pos);	    
	    }
	  for (int i=0;i<4;++i)
	    {
	      getline(stream,line);
	      original.push_back(line);
	    }
	
	  while (stream.get()!='_');
	}
      else 
	{
	  std::cout << "Unable to open file "<< fname << std::endl; 
	  exit(-1);
	}
      /*
	for (long i=0;i<original.size();++i)
	std::cout << original[i] << std::endl;
      */
      return stream;
    }

    template <class T>
    void write(const std::string &fname, T *i0, T *di)
    {
      std::ofstream stream(fname);
      write(fname,stream,i0,di);
    }

    void write(const std::string &fname)
    {      
      int i0[NDIM];
      int di[NDIM];
      for (int i=0;i<NDIM;++i)
	{
	  i0[i]=wExtent[2*i];
	  di[i]=wExtent[2*i+1]-wExtent[2*i];
	}
      write(fname,i0,di);
    }

    void write(const std::string &fname, std::ofstream &stream)
    {      
      int i0[NDIM];
      int di[NDIM];
      for (int i=0;i<NDIM;++i)
	{
	  i0[i]=wExtent[2*i];
	  di[i]=wExtent[2*i+1]-wExtent[2*i];
	}
      write(fname,stream,i0,di);
    }
  
    template <class T>
    std::ofstream &write(const std::string &fname, std::ofstream &stream, 
			 T *i0, T *di)
    {
      stream.close();
      stream.open(fname);      

      int index=0;
      T imax[NDIM];
      size_t pos;

      for (int i=0;i<NDIM;++i) imax[i]=i0[i]+di[i];
    
      stream << original[index++] << std::endl;
   
      // Whole extent
      pos=original[index].find("\"")+1;
      stream << original[index].substr(0,pos) << i0[0]<<" "<< imax[0];
      for (int i=1;i<NDIM;++i) stream<<" "<<i0[i]<<" "<< imax[i] ;
      stream << "\">" << std::endl;
      index++;

      // Piece extent
      pos=original[index].find("\"")+1;
      stream << original[index].substr(0,pos) << i0[0]<<" "<< imax[0];
      for (int i=1;i<NDIM;++i) stream <<" "<< i0[i]<<" "<< imax[i] ;
      stream << "\">" << std::endl;     
      index++;

      for (int i=0;i<6;++i)
	stream << original[index++] << std::endl;

      long delta=sizeof(double);
      for (int i=0;i<NDIM;++i) delta*=di[i];
      delta+=type;
        
      for (int i=0;i<NDIM;++i)
	{
	  pos=original[index].find("offset");
	  pos+=original[index].substr(pos).find("\"")+1;
	  stream << original[index].substr(0,pos) << delta <<"\" />"<<std::endl;
	  index++;
	  delta+=type+sizeof(double)*(di[i]+1);
	}
    
      for (int i=0;i<4;++i)
	stream << original[index++] << std::endl;
      stream << "_";  

      return stream;
    }

    size_t check(const std::string &fname, const std::string &line,const std::string &what)
    {
      size_t pos=line.find(what);    
      if (pos==std::string::npos)
	{
	  std::cout << "Invalid VTR file: could not find tag "<< what << " in file " << fname << std::endl;	
	  exit(-1);
	}
    
      return pos;
    }
  };

  struct Indexer
  {
    long x0[NDIM];
    long dx[NDIM];
    long dims[NDIM];
    long stride[NDIM+1];
    long w[NDIM];
    long delta[NDIM];
    long end_;
    long k;
    
    long p0[NDIM];
	    

    Indexer(){}
  
    template <class T=int>
    Indexer(Header &h, T *i0=nullptr, T *di=nullptr)
    {
      init(h,i0,di);
    }

    template <class T=int>
    void init(Header &h, T *i0=nullptr, T *di=nullptr)
    {
      stride[0]=1;
      for (int i=0;i<NDIM;++i)
	{
	  dims[i]=h.pExtent[2*i+1]-h.pExtent[2*i];
	  p0[i]=h.pExtent[2*i];
	  stride[i+1]=stride[i]*dims[i];
	  x0[i]=(i0==nullptr)?0:i0[i];
	  dx[i]=(di==nullptr)?dims[i]:di[i];
	}
      k=0;
      end_=stride[NDIM];
      for (int i=0;i<NDIM;++i)
	{
	  w[i]=0;
	  delta[i]=stride[i+1]-(dx[i]*stride[i]);
	  k+=x0[i]*stride[i];
	}    
    }

    long getNext()
    {
      ++k;
      dice::hlp::getNext<NDIM>(w,k,dx,delta);
      return k;
    }

    long get()
    {
      return k;
    }

    long get(int index)
    {
      return w[index];
    }

    long getWhole(int index)
    {
      return p0[index]+w[index];
    }

    long getDim(int index)
    {
      return dims[index];
    }

    long end()
    {
      return end_;
    }
  };

}

#endif

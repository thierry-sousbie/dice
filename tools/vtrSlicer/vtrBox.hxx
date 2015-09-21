#ifndef VTR_BOX_HXX__
#define VTR_BOX_HXX__

#include <vector>
#include <string>
#include <iostream> 
#include <fstream>
#include <queue>    
#include <stdlib.h>
#include <string.h>
#include <sstream> 

#include "vtrTools.hxx"

namespace vtrSlicer {

  template <class P>
  void box(P &params, 
	   std::vector<std::string> &fName, 
	   std::vector<std::string> &outName)
  {
    int i0[NDIM];
    int di[NDIM];
    std::ostringstream oss;

    for (int i=0;i<NDIM;++i)
      {i0[i]=params.front();params.pop();}
    for (int i=0;i<NDIM;++i)
      {di[i]=params.front();params.pop();}

    oss << ".";
    for (int i=0;i<NDIM;++i)
      oss << i0[i] <<"_";
    for (int i=0;i<NDIM;++i)
      oss << di[i] <<"_";
    oss << "box";

    for (int f=0;f<fName.size();f++)
      {
	Header header;
	std::ifstream ifs;
	std::ofstream ofs;      
	std::string outputName=outName[f]+oss.str()+".vtr";      

	header.read(fName[f],ifs);
	header.write(outputName,ofs,i0,di);

	// switch to binary mode
	size_t ifpos=ifs.tellg();   
	ifs.close();
	ifs.open(fName[f]);
	ifs.seekg(ifpos+header.type,ifs.beg);
      
	ofs.close();
	ofs.open(outputName,std::ofstream::out | std::ofstream::app);

	long nWrite=1;
	for (int j=0;j<NDIM;++j) nWrite*=di[j];
	char tmp[8];
	if (header.type==4) 
	  *reinterpret_cast<unsigned int*>(&tmp)=nWrite*sizeof(double);
	else
	  *reinterpret_cast<unsigned long*>(&tmp)=nWrite*sizeof(double);
	ofs.write(tmp,header.type);

	Indexer indexer(header,i0,di);
	long j0=indexer.get();
	long oldj=j0;
	ifs.seekg(j0*sizeof(double),ifs.cur);
	static const long bufSize=1<<24;

	std::vector<double> vBuffer(bufSize);
	char *buffer=reinterpret_cast<char*>(&vBuffer[0]);

	long nWritten=0;
	for (long j=indexer.getNext();j<indexer.end();j=indexer.getNext())
	  {
	  
	    if (j-oldj == 1) // consecutive in input
	      {	      
		if (j-j0 == bufSize) // buffer is saturated
		  {
		    ifs.read(buffer,sizeof(double)*(j-j0));
		    ofs.write(buffer,sizeof(double)*(j-j0));
		    nWritten+=(j-j0);
		    j0=j;
		  }	      
	      }
	    else // jump in input -> flush
	      {
		ifs.read(buffer,sizeof(double)*(oldj-j0+1));
		ofs.write(buffer,sizeof(double)*(oldj-j0+1));
		nWritten+=(oldj-j0+1);
		//std::cout <<"fmush"<<oldj-j0+1<<" "<<j<<std::endl;
		j0=j;
		ifs.seekg((j-oldj-1)*sizeof(double),ifs.cur);
	      }
	    oldj=j;
	  }
	//std::cout << oldj <<" "<<j0 <<std::endl;
	ifs.read(buffer,sizeof(double)*(oldj-j0+1));
	ofs.write(buffer,sizeof(double)*(oldj-j0+1));
	nWritten+=(oldj-j0+1);
	//std::cout << "wrote " << nWritten << "==" <<nWrite<<std::endl;

	ifs.seekg(ifpos+header.type+indexer.end()*sizeof(double),ifs.beg);
	for (int j=0;j<NDIM;++j)
	  {
	    long delta=(header.wExtent[2*j+1]-header.wExtent[2*j]);
	    ifs.seekg(header.type,ifs.cur);
	    ifs.read(buffer,sizeof(double)*(delta+1));

	    for (int k=0;k<delta;++k) 
	      vBuffer[k+delta+1]=vBuffer[delta]+(vBuffer[k+1]-vBuffer[0]);

	    for (int k=i0[j];k<=i0[j]+di[j];++k)
	      vBuffer[k-i0[j]]=vBuffer[k];

	    char tmp[8];
	    if (header.type==4) 
	      *reinterpret_cast<unsigned int*>(&tmp)=(di[j]+1)*sizeof(double);
	    else
	      *reinterpret_cast<unsigned long*>(&tmp)=(di[j]+1)*sizeof(double);

	    ofs.write(tmp,header.type);
	    ofs.write(buffer,sizeof(double)*(di[j]+1));	  
	  }
      }
  }
}

#endif

#ifndef VTR_PROJECT_HXX__
#define VTR_PROJECT_HXX__

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
  void project(P &params, 
	       std::vector<std::string> &fName, 
	       std::vector<std::string> &outName)
  {
    int projDir=params.front();params.pop();
    Header header0;
    int idx[NDIM]={1,0,0};
    int idy[NDIM]={2,2,1};

    header0.read(fName[0]);

    for (int f=0;f<fName.size();f++)
      {
	Header h;
	h.read(fName[f]);
	for (int i=0;i<NDIM;++i) 
	  {
	    if (h.wExtent[2*i]<header0.wExtent[2*i])
	      header0.wExtent[2*i]=h.wExtent[2*i];
	    if (h.wExtent[2*i+1]>header0.wExtent[2*i+1])
	      header0.wExtent[2*i+1]=h.wExtent[2*i+1];

	    //header0.wExtent[2*i]=0; //works only that way for now ;)
	  }
      }    
    std::copy_n(header0.wExtent,NDIM*2,header0.pExtent);
    
    long dx=header0.getWDim(idx[projDir]);
    std::vector<double> dest(header0.getWDim(idx[projDir])*header0.getWDim(idy[projDir]),0);
    std::cout << "Output dims :" << header0.getWDim(idx[projDir]) <<"x"
	      << header0.getWDim(idy[projDir]) <<std::endl;
    Header header;
    std::ifstream ifs;
    size_t ifpos=0;
    for (int f=0;f<fName.size();f++)
      {	
	ifs.close();
	header.read(fName[f],ifs);
	std::vector<double> buffer(header.getPDim(0)*header.getPDim(1),0);  
  
	ifpos=ifs.tellg();   
	ifs.close();
	ifs.open(fName[f]);
	ifs.seekg(ifpos+header.type,ifs.beg);

	Indexer indexer(header);
	long srcId=-1;
	for (long j=0;j<indexer.end();j=indexer.getNext())
	  {
	    if ((srcId>=buffer.size())||(srcId<0))
	      {
		ifs.read(reinterpret_cast<char*>(&buffer[0]),sizeof(double)*buffer.size());
		srcId=0;
	      }
	    long deltax=(header.wExtent[2*idx[projDir]]-header0.wExtent[2*idx[projDir]]);
	    long deltay=(header.wExtent[2*idy[projDir]]-header0.wExtent[2*idy[projDir]]);
	    long dstId=(indexer.get(idx[projDir])+deltax)+
	      dx*(indexer.get(idy[projDir])+deltay);

	    dest[dstId]+=buffer[srcId++];	    
	  } 
      }    
    
    std::copy_n(header0.wExtent,2*NDIM,header.wExtent);    
    std::vector<double> coords[NDIM];
    
    for (int i=0;i<NDIM;++i)
      {		
	long delta=(header.wExtent[2*i+1]-header.wExtent[2*i]);
	long deltaP=(header.pExtent[2*i+1]-header.pExtent[2*i]);
	long i0P=header.pExtent[2*i]-header.wExtent[2*i];
	//long i0W=header.wExtent[2*i];
	//std::cout<<"delta="<<delta<<std::endl;
	coords[i].resize(delta+1);
	ifs.seekg(header.type,ifs.cur);
	ifs.read(reinterpret_cast<char*>(&coords[i][i0P]),sizeof(double)*(deltaP+1));
	//std::cout <<"read "<<(deltaP+1)<<"@"<<i0P<<"/"<<delta<<std::endl;

	if (i==projDir)
	  {
	    coords[i].resize(2);
	    coords[i][0]=0;
	    coords[i][1]=0;	    
	  }
	else
	  {
	    double di=coords[i][i0P+1]-coords[i][i0P];
	    for (long j=i0P-1;j>=0;--j)
	      coords[i][j]=coords[i][j+1]-di;
	    for (long j=i0P+deltaP+1;j<delta+1;++j)
	      coords[i][j]=coords[i][j-1]+di;
	    // std::cout <<"di="<<di<<std::endl;
	    // for (int j=0;j<coords[i].size();++j) std::cout <<coords[i][j] <<" ";
	    // std::cout << std::endl;
	    
	  }
      }
    
    for (int i=projDir;i<NDIM-1;++i)
      {
	std::swap(coords[i],coords[i+1]);
	std::swap(header0.wExtent[2*i],header0.wExtent[2*(i+1)]);
	std::swap(header0.wExtent[2*i+1],header0.wExtent[2*(i+1)+1]);
	std::swap(header0.pExtent[2*i],header0.pExtent[2*(i+1)]);
	std::swap(header0.pExtent[2*i+1],header0.pExtent[2*(i+1)+1]);
      }
    
    std::ostringstream oss;
    oss << ".Proj" << projDir;    
    std::string outputName=outName[0]+oss.str()+".vtr";
    std::ofstream ofs;    

    header0.setDim(NDIM-1,1);
    header0.write(outputName,ofs);

    ofs.close();
    ofs.open(outputName,std::ofstream::out | std::ofstream::app);

    long nWrite=1;
    for (int j=0;j<NDIM;++j) nWrite*=header0.getWDim(j);
    char tmp[8];
    if (header0.type==4) 
      *reinterpret_cast<unsigned int*>(&tmp)=nWrite*sizeof(double);
    else
      *reinterpret_cast<unsigned long*>(&tmp)=nWrite*sizeof(double);
    ofs.write(tmp,header0.type);
    ofs.write(reinterpret_cast<char*>(&dest[0]),sizeof(double)*dest.size());

    for (int j=0;j<NDIM;++j)
      {
	long delta=coords[j].size()-1;

	char tmp[8];
	if (header0.type==4) 
	  *reinterpret_cast<unsigned int*>(&tmp)=(delta+1)*sizeof(double);
	else
	  *reinterpret_cast<unsigned long*>(&tmp)=(delta+1)*sizeof(double);

	ofs.write(tmp,header0.type);
	ofs.write(reinterpret_cast<char*>(&coords[j][0]),sizeof(double)*coords[j].size());
      }
     
  }
}

#endif

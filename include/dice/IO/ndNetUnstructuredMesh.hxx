#ifndef __DICE_NDNET_UNSTRUCTURED_MESH_WRITER_HXX__
#define __DICE_NDNET_UNSTRUCTURED_MESH_WRITER_HXX__

#include <vector>

#include "../dice_globals.hxx"
#include "../tools/IO/myIO.hxx"

#include "./NDnet/NDNetworkDefines.hxx"
#include "./NDnet/defaultNDnetFilters.hxx"
#include "../tools/helpers/helpers.hxx"

/**
 * @file 
 * @brief  Definition of a class to write unstructured meshes to NDnet format
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup IO
 *   \{
 */

namespace IO {

  enum NDNET_WriterOptions {
    NDNET_Default=0,
    NDNET_WithGhosts=(1<<0), 
    NDNET_WithShadows=(1<<1), 
    NDNET_WithNeighbors=(1<<2),
    NDNET_NoSimplices=(1<<3),
    NDNET_NoSimplexData=(1<<4),
    NDNET_NoVertexData=(1<<5),
    NDNET_TrimExtraVertices=(1<<6)
  };

  template <class M,typename CT=double> 
  class NDnetUnstructuredMeshWriterT 
  {
  public:
    typedef NDnetUnstructuredMeshWriterT<M,CT> MyType;
    static const int NDIM = M::NDIM;
    static const int NDIM_W = M::NDIM_W;
   
    typedef M Mesh;
    typedef CT Float;
    typedef long Int;   

    typedef typename M::Simplex Simplex;
    typedef typename M::GhostSimplex GhostSimplex;
    typedef typename M::ShadowSimplex ShadowSimplex;

    typedef typename M::Facet Facet;
    typedef typename M::Segment Segment;
    typedef typename M::FacetHandle FacetHandle;
    typedef typename M::SegmentHandle SegmentHandle;

    typedef typename M::Vertex Vertex;
    typedef typename M::GhostVertex GhostVertex;
    typedef typename M::ShadowVertex ShadowVertex;

    typedef typename Vertex::LocalIndex LocalIndex;

    NDnetUnstructuredMeshWriterT(M *m, 
				 long opt=NDNET_Default,
				 const char *fileName="mesh",
				 int defaultNThreads=glb::num_omp_threads):
      mesh(m), 
      options(opt),
      useVArr(false),
      useSArr(false),
      useFacetArr(false),
      useSegArr(false),
      defaultNTh(defaultNThreads),
      useVIndexArr(false)
    {
      fName = toFilename(fileName);
      if (opt&NDNET_WithGhosts)
	{
	  gvArr=mesh->getGhostVerticesArray();
	  gsArr=mesh->getGhostSimplicesArray();
	}
      if (opt&NDNET_WithShadows)
	{
	  svArr=mesh->getShadowVerticesArray();
	  ssArr=mesh->getShadowSimplicesArray();
	}

      if (opt&NDNET_NoSimplices)
	{
	  useSArr=true;
	}
      else if (opt&NDNET_WithNeighbors)
	{
	  sArr=mesh->getSimplicesArray();
	  useSArr=true;
	}
     
    }

    void setFileName(const char *fileName)
    {
      fName = toFilename(fileName);
    }

    virtual ~NDnetUnstructuredMeshWriterT()
    {}
   
    static std::string toFilename(const char *name="mesh", bool stripPath=false)
    {
      int start=0;
      if (stripPath)
	{
	  int len=strlen(name);
	  for (int i=0;i<len;++i) 
	    if (name[i]=='/') start=i+1;
	}
      return std::string(&name[start]) + std::string(".NDnet");
    }
    
    template <class F, class CellType=typename F::CellType>
    void filter(const F& functor, int nThreads=-1)
    {    
      if (nThreads<0) nThreads=defaultNTh;
      filter<F>(functor,nThreads,hlp::ConstantType<CellType>());
    }    

    void write(bool quiet=false)
    {
      if (!quiet)
	glb::console->printFlush<LOG_STD>("Dumping %dD unstructured mesh to NDnet file '%s' ... ",NDIM,fName.c_str());
          
      myIO::BinaryWriterT<> writer(fName);
      write(&writer);

      if (!quiet) glb::console->print<LOG_STD>("done.\n");
    }
   
    template <class BW>
    void write(BW *bWriter)
    {      
      if ((options & NDNET_TrimExtraVertices)||(useVArr))
	trimVertices();
      
      std::vector<unsigned long> nCells = getCellsCount();
      double x0[NDIM_W];
      double delta[NDIM_W];
      char comment[80];
    
      sprintf(comment,"Local mesh");
      mesh->getBoundingBox(x0,delta);
      NDnetwork net(NDIM,NDIM_W,x0,delta,&nCells[0],comment,
		    M::BOUNDARY_TYPE==BoundaryType::PERIODIC);

      FILE *f = bWriter->getFilePtr();
      // header
      net.writeHeader(f);

      // vertex coordinates
      unsigned int jj=net.ndims*nCells[0];
      fwrite(&jj,sizeof(unsigned int),1,f);

      if (!useVArr)
	{
	  vArr=mesh->getVerticesArray();
	  useVArr=true;
	}
    
      if (!useVArr)
	{
	  for (auto it=mesh->vertexBegin(); it!=mesh->vertexEnd(); ++it)
	    bWriter->template writeAs<NDNET_FLOAT>((*it)->getCoordsConstPtr(),NDIM_W);
	}
      else
	{
	  for (unsigned long i=0;i<vArr.size();i++)
	    bWriter->template writeAs<NDNET_FLOAT>(vArr[i]->getCoordsConstPtr(),NDIM_W);
	}
    
      for (unsigned long i=0;i<gvArr.size();++i)
	bWriter->template writeAs<NDNET_FLOAT>(gvArr[i]->getCoordsConstPtr(),NDIM_W);
      for (unsigned long i=0;i<svArr.size();++i)
	bWriter->template writeAs<NDNET_FLOAT>(svArr[i]->getCoordsConstPtr(),NDIM_W);
     
      bWriter->flush();
      fwrite(&jj,sizeof(unsigned int),1,f);
     
      // cell count
      jj=(1+net.ndims)*sizeof(NDNET_UINT);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(net.nfaces,sizeof(NDNET_UINT),((size_t)net.ndims+1),f);
      fwrite(&jj,sizeof(unsigned int),1,f);
     
      // defined cells
      jj=(1+net.ndims)*sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(net.haveVertexFromFace,sizeof(int),((size_t)net.ndims+1),f);
      fwrite(&jj,sizeof(unsigned int),1,f);
     
      // segments
      if ((nCells[1]>0)&&(Facet::NVERT!=Segment::NVERT))
	{
	  jj=sizeof(NDNET_UINT)*((size_t)(Segment::NVERT)*nCells[Segment::NVERT-1]);
	  fwrite(&jj,sizeof(unsigned int),1,f);
     
	  for (unsigned long j=0;j<segArr.size();++j)
	    {	    
	      NDNET_UINT index[Segment::NVERT];
	      segArr[j].getVerticesLocalIndex(index);	     
	      if (useVIndexArr)
		for (int i=0;i<Segment::NVERT;++i)
		  index[i]=getIndex(index[i]);
	      bWriter->write(index,Segment::NVERT);	
	    }     
	  bWriter->flush();
	  fwrite(&jj,sizeof(unsigned int),1,f);
	}

      // facets
      if (nCells[NDIM-1]>0)
	{
	  jj=sizeof(NDNET_UINT)*((size_t)(Facet::NVERT)*nCells[Facet::NVERT-1]);
	  fwrite(&jj,sizeof(unsigned int),1,f);
     
	  for (unsigned long j=0;j<facetArr.size();++j)
	    {
	      NDNET_UINT index[Facet::NVERT];
	      facetArr[j].getVerticesLocalIndex(index);	     
	      if (useVIndexArr)
		for (int i=0;i<Facet::NVERT;++i)
		  index[i]=getIndex(index[i]);
	      bWriter->write(index,Facet::NVERT);	     
	    }    

	  if (Facet::NVERT==Segment::NVERT)
	    {
	      for (unsigned long j=0;j<segArr.size();++j)
		{	    
		  NDNET_UINT index[Segment::NVERT];
		  segArr[j].getVerticesLocalIndex(index);	     
		  if (useVIndexArr)
		    for (int i=0;i<Segment::NVERT;++i)
		      index[i]=getIndex(index[i]);
		  bWriter->write(index,Facet::NVERT);	
		}     
	    }
	  bWriter->flush();
	  fwrite(&jj,sizeof(unsigned int),1,f);
	}     
     
      // simplices 
      if (nCells[NDIM]>0)
	{
	  jj=sizeof(NDNET_UINT)*((size_t)(Simplex::NVERT)*nCells[Simplex::NVERT-1]);
	  fwrite(&jj,sizeof(unsigned int),1,f);
       
	  if (!useSArr)
	    {
	      const auto sim_end=mesh->simplexEnd();
	      for (auto it=mesh->simplexBegin();it!=sim_end;++it)
		{	     
		  typename Vertex::LocalIndex tmp[Simplex::NVERT];
		  it->getVerticesLocalIndex(tmp);
		  if (useVIndexArr)
		    for (int i=0;i<Simplex::NVERT;++i)
		      tmp[i]=getIndex(tmp[i]);
		  bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
		}
	    }
	  else
	    {	 	
	      for (unsigned long j=0;j<sArr.size();++j)
		{		     
		  typename Vertex::LocalIndex tmp[Simplex::NVERT];
		  sArr[j]->getVerticesLocalIndex(tmp);
		  if (useVIndexArr)
		    for (int i=0;i<Simplex::NVERT;++i)
		      tmp[i]=getIndex(tmp[i]);
		  bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
		}       
	    }
    
	  long nLocalVertices=(!useVArr)?mesh->getNVertices():vArr.size();
	  for (unsigned long j=0;j<gsArr.size();++j)
	    {	 
	      typename Vertex::LocalIndex tmp[Simplex::NVERT];
	      gsArr[j]->getVerticesLocalIndex(tmp);
	      for (int k=0;k<Simplex::NVERT;k++) 
		{
		  if (gsArr[j]->getVertex(k)->isGhost()) 
		    tmp[k]+=nLocalVertices;
		}
	      bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
	    }   	
	  for (unsigned long j=0;j<ssArr.size();++j)
	    {	
	      typename Vertex::LocalIndex tmp[Simplex::NVERT];
	      ssArr[j]->getVerticesLocalIndex(tmp);
	      for (int k=0;k<Simplex::NVERT;k++) 
		{
		  if (ssArr[j]->getVertex(k)->isGhost()) 
		    tmp[k]+=nLocalVertices;
		  else if (ssArr[j]->getVertex(k)->isShadow()) 
		    tmp[k]+=nLocalVertices+gvArr.size();
		}
	      bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NVERT);	
	    }       
	  bWriter->flush();
        
	  fwrite(&jj,sizeof(unsigned int),1,f);      
	}
     
      // junk
      jj=(1+net.ndims)*sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(net.haveFaceFromVertex,sizeof(int),((size_t)net.ndims+1),f);
      fwrite(&jj,sizeof(unsigned int),1,f);
        
      if (options&NDNET_WithNeighbors) net.haveFaceFromFace[NDIM][NDIM]=1;  

      jj=(1+net.ndims)*(1+net.ndims)*sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      for (long i=0;i<net.ndims+1;i++)
	fwrite(net.haveFaceFromFace[i],sizeof(int),((size_t)net.ndims+1),f);
      fwrite(&jj,sizeof(unsigned int),1,f);

      if (options&NDNET_WithNeighbors)
	{	
	  int i=NDIM;	
	  std::vector<NDNET_IDCUMT> nn(nCells[i]+1);
	  nn[0]=0;
	  for (unsigned long j=1;j<nn.size();++j)
	    nn[j]=nn[j-1]+Simplex::NNEI;
	  jj=sizeof(NDNET_IDCUMT)*nn.size();
	  fwrite(&jj,sizeof(unsigned int),1,f);
	  fwrite(&nn[0],sizeof(NDNET_IDCUMT),((size_t)net.nfaces[i]+1),f);
	  fwrite(&jj,sizeof(unsigned int),1,f);
	  nn.clear();
	
	  jj=sizeof(NDNET_UINT)*(nCells[i]*Simplex::NNEI);
	  fwrite(&jj,sizeof(unsigned int),1,f);
	  const long myRank=mesh->getMpiCom()->rank();
	  for (unsigned long j=0;j<sArr.size();++j)
	    {
	      NDNET_UINT tmp[Simplex::NNEI];
	      Simplex *cur=sArr[j];
	      for (int k=0;k<Simplex::NNEI;++k)
		{
		  //Simplex *cur=sArr[j];
		  Simplex *nei=cur->getNeighbor(k);
		  if (nei==NULL)
		    tmp[k]=j;
		  else if (myRank!=nei->getGlobalIdentity(myRank).rank())
		    tmp[k]=j;
		  else
		    tmp[k]=nei->getLocalIndex();
		}
	      bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
	    }

	  if ((options&NDNET_WithGhosts)||
	      (options&NDNET_WithShadows))
	    {
	      long delta[3];
	      delta[0]=0;
	      delta[1]=sArr.size();
	      delta[2]=delta[1]+gsArr.size();
	      for (unsigned long j=0;j<gsArr.size();++j)
		{
		  NDNET_UINT tmp[Simplex::NNEI];
		  for (int k=0;k<Simplex::NNEI;++k)
		    {
		      Simplex *nei=gsArr[j]->getNeighbor(k);
		      if (nei==NULL)
			tmp[k]=delta[1]+j;
		      else 
			{
			  tmp[k]=nei->getLocalIndex();
			  if (nei->isGhost())
			    tmp[k]+=delta[1];
			  else if (nei->isShadow())
			    tmp[k]+=delta[2]; 	    
			}
		    }
		  bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
		}

	      for (unsigned long j=0;j<ssArr.size();++j)
		{
		  NDNET_UINT tmp[Simplex::NNEI];
		  for (int k=0;k<Simplex::NNEI;++k)
		    {
		      Simplex *nei=ssArr[j]->getNeighbor(k);
		      if (nei==NULL)
			tmp[k]=delta[2]+j;
		      else 
			{
			  tmp[k]=nei->getLocalIndex();
		    
			  if (nei->isGhost())
			    tmp[k]+=delta[1];
			  else if (nei->isShadow())
			    tmp[k]+=delta[2]; 	
			}
		    }
		  bWriter->template writeAs<NDNET_UINT>(tmp,Simplex::NNEI);		    
		}
	    }
	  bWriter->flush();
	  fwrite(&jj,sizeof(unsigned int),1,f);
	}

      bWriter->flush();
      net.haveVFlags=1;
      jj=sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(&net.haveVFlags,sizeof(int),1,f);
      fwrite(&jj,sizeof(unsigned int),1,f);
      
      // vertex flags
      jj=sizeof(unsigned char)*nCells[0];
      fwrite(&jj,sizeof(unsigned int),1,f);
    
      if (!useVArr)
	{
	  for (auto it=mesh->vertexBegin(); it!=mesh->vertexEnd(); ++it)
	    {
	      typename Vertex::Flag flags=(*it)->getFlags();
	      bWriter->template writeAs<unsigned char>(&flags);
	    }
	}
      else
	{
	  for (unsigned long i=0;i<vArr.size();i++)
	    {
	      typename Vertex::Flag flags=vArr[i]->getFlags();
	      bWriter->template writeAs<unsigned char>(&flags);
	    }
	}


      for (unsigned long i=0;i<gvArr.size();i++)
	{
	  typename Vertex::Flag flags=gvArr[i]->getFlags();
	  bWriter->template writeAs<unsigned char>(&flags);
	}
      for (unsigned long i=0;i<svArr.size();i++)
	{
	  typename Vertex::Flag flags=svArr[i]->getFlags();
	  bWriter->template writeAs<unsigned char>(&flags);
	}
    
      bWriter->flush();
      fwrite(&jj,sizeof(unsigned int),1,f);
    
      // cell flags    
      net.haveFFlags[NDIM]=(nCells[NDIM]>0)?1:0;
      jj=sizeof(int)*(net.ndims+1);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(net.haveFFlags,sizeof(int),(net.ndims+1),f);
      fwrite(&jj,sizeof(unsigned int),1,f);

      if (net.haveFFlags[NDIM])
	{
	  //simplex flags
	  jj=sizeof(unsigned char)*nCells[NDIM];
	  fwrite(&jj,sizeof(unsigned int),1,f);
     
	  if (!useSArr)
	    {
	      const auto s_end = mesh->simplexEnd();
	      for (auto it=mesh->simplexBegin(); it!=s_end; ++it)
		{
		  typename Simplex::Flag flags=(*it)->getFlags();
		  bWriter->template writeAs<unsigned char>(&flags);
		}
	    }
	  else
	    {
	      for (unsigned long i=0;i<sArr.size();i++)
		{
		  typename Simplex::Flag flags=sArr[i]->getFlags();
		  bWriter->template writeAs<unsigned char>(&flags);
		}
	    }


	  for (unsigned long i=0;i<gsArr.size();i++)
	    {
	      typename Simplex::Flag flags=gsArr[i]->getFlags();
	      bWriter->template writeAs<unsigned char>(&flags);
	    }
	  for (unsigned long i=0;i<ssArr.size();i++)
	    {
	      typename Simplex::Flag flags=ssArr[i]->getFlags();
	      bWriter->template writeAs<unsigned char>(&flags);
	    }
    
	  bWriter->flush();
	  fwrite(&jj,sizeof(unsigned int),1,f);
	}
      // data
      int nVertexData=mesh->getNVertexFunctor();
      int nSimplexData=(nCells[NDIM]>0)?mesh->getNSimplexFunctor():0;
    
      int nData=0;
      if (!(options&NDNET_NoVertexData))
	for (int i=0;i<nVertexData;i++) 
	  {
	    const auto *vf=mesh->getVertexFunctorPtr(i);
	    if (!(vf->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP))
	      nData+=vf->getSize();      
	  }

      if ((nCells[NDIM]>0)&&(!(options&NDNET_NoSimplexData)))	
	for (int i=0;i<nSimplexData;i++) 
	  {
	    const auto *sf=mesh->getSimplexFunctorPtr(i);
	    if (!(sf->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP))
	      nData+=sf->getSize();
	  }
    
      jj=sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(&nData,sizeof(int),1,f);
      fwrite(&jj,sizeof(unsigned int),1,f);

      if (nData)
	{	
	  int type=0;
	  int count=0;
	  char name[255];
	  int nData = mesh->getNVertexFunctor();
	  if (options&NDNET_NoVertexData) nData=0;

	  for (int i=0;i<nData;i++)
	    {
	      const auto *gvd=mesh->getVertexFunctorPtr(i);
	      count=gvd->getSize();

	      if (gvd->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP)
		continue;

	      for (int ct=0;ct<count;ct++)
		{//printf(name,"%s_%2.2d\n",gsd->getName().c_str(),ct);
		  if (count>1)
		    sprintf(name,"%s_%2.2d",gvd->getName().c_str(),ct);
		  else
		    strcpy(name,gvd->getName().c_str());
		
		  jj=sizeof(int)+255*sizeof(char);
		  fwrite(&jj,sizeof(unsigned int),1,f);
		  fwrite(&type,sizeof(int),1,f);
		  fwrite(name,sizeof(char)*255,1,f);
		  fwrite(&jj,sizeof(unsigned int),1,f);

		  jj=sizeof(double)*nCells[0];
		  fwrite(&jj,sizeof(unsigned int),1,f);
		
		  if (!useVArr)
		    {
		      for (auto it=mesh->vertexBegin(); it!=mesh->vertexEnd(); ++it)
			{
			  double tmp=gvd->get(*it,ct);
			  bWriter->template writeAs<double>(&tmp);	
			}
		    }
		  else
		    {
		      for (unsigned long j=0;j<vArr.size();j++)
			{
			  double tmp=gvd->get(vArr[j],ct);
			  bWriter->template writeAs<double>(&tmp);			
			}
		    }

		  for (unsigned long j=0;j<gvArr.size();j++)
		    {
		      double tmp=gvd->get(gvArr[j],ct);
		      bWriter->template writeAs<double>(&tmp);
		    }
		  for (unsigned long j=0;j<svArr.size();j++)
		    {
		      double tmp=gvd->get(svArr[j],ct);
		      bWriter->template writeAs<double>(&tmp);
		    }		    
	
		  bWriter->flush();
		  fwrite(&jj,sizeof(unsigned int),1,f);
		}		  
	    }
	
	  bWriter->flush();
	  nData = (nCells[NDIM]>0)?mesh->getNSimplexFunctor():0;
	  if (options&NDNET_NoSimplexData) nData=0;
	  type=NDIM;
	
	  for (int i=0;i<nData;i++)
	    {
	      const auto *gsd=mesh->getSimplexFunctorPtr(i);
	      count=gsd->getSize();

	      if (gsd->getFlags()&cellDataFunctors::F_SKIP_ON_FILE_DUMP)
		continue;
	   
	      for (int ct=0;ct<count;ct++)
		{//printf(name,"%s_%2.2d\n",gsd->getName().c_str(),ct);
		  if (count>1)
		    sprintf(name,"%s_%2.2d",gsd->getName().c_str(),ct);
		  else 
		    strcpy(name,gsd->getName().c_str());
		
		  jj=sizeof(int)+255*sizeof(char);
		  fwrite(&jj,sizeof(unsigned int),1,f);
		  fwrite(&type,sizeof(int),1,f);
		  fwrite(name,sizeof(char)*255,1,f);
		  fwrite(&jj,sizeof(unsigned int),1,f);

		  //printf("WRITING field %s (%ld ele).\n",name,);

		  jj=sizeof(double)*nCells[NDIM-1];
		  fwrite(&jj,sizeof(unsigned int),1,f);

		  if (!useSArr)
		    {
		      for (auto it=mesh->simplexBegin(); it!=mesh->simplexEnd(); ++it)
			{
			  double tmp=gsd->get(*it,ct);
			  bWriter->template writeAs<double>(&tmp);	
			}
		    }
		  else
		    {
		      for (unsigned long j=0;j<sArr.size();j++)
			{
			  double tmp=gsd->get(sArr[j],ct);
			  bWriter->template writeAs<double>(&tmp);			
			}
		    }
	
		  for (unsigned long j=0;j<gsArr.size();j++)
		    {
		      double tmp=gsd->get(gsArr[j],ct);
		      bWriter->template writeAs<double>(&tmp);
		    }
		  for (unsigned long j=0;j<ssArr.size();j++)
		    {
		      double tmp=gsd->get(ssArr[j],ct);
		      bWriter->template writeAs<double>(&tmp);
		    }
		 
		  bWriter->flush();
		  fwrite(&jj,sizeof(unsigned int),1,f);
		}		  
	    }
	}
      bWriter->flush();
      //supdata
      jj=sizeof(int);
      fwrite(&jj,sizeof(unsigned int),1,f);
      fwrite(&net.nsupData,sizeof(int),1,f);
      fwrite(&jj,sizeof(unsigned int),1,f);
    }

  private:        
    template <class F>
    void filter(const F &functor, int nThreads, hlp::ConstantType<Simplex>)
    {
      if (!useSArr)
	{
	  if (nThreads==1)
	    {
	      const auto it_end=mesh->simplexEnd();
	      for (auto it=mesh->simplexBegin();it!=it_end;++it)
		{
		  Simplex *s=functor(*it);
		  if (s!=NULL)
		    sArr.push_back(s);
		}
	    }
	  else
	    {
#pragma omp parallel num_threads(nThreads)
	      {
		int th = omp_get_thread_num();
		std::vector<Simplex*> localVec;		
		const auto it_end=mesh->simplexEnd(th,nThreads);
		for (auto it=mesh->simplexBegin(th,nThreads);it!=it_end;++it)
		  {
		    Simplex *s=functor(*it);
		    if (s!=NULL) localVec.push_back(s);		      
		  }
#pragma omp critical
		sArr.insert(sArr.end(),localVec.begin(),localVec.end());
	      }
	    }
	  useSArr=true;	 	  
	}
    }
    
    template <class F>
    void filter(const F &functor, int nThreads, hlp::ConstantType<FacetHandle>)
    {
      if (!useFacetArr)
	{
	  if (nThreads==1)
	    {
	      const auto it_end=mesh->simplexEnd();
	      for (auto it=mesh->simplexBegin();it!=it_end;++it)
		{
		  Simplex *s=(*it);
		  for (int i=0;i<Simplex::NNEI;++i)
		    {
		      if (s->getNeighbor(i)>s)
			{
			  FacetHandle h=functor(s->getFacetHandle(i));
			  if (h->isSet())
			    facetArr.push_back(*h);
			}
		    }
		}
	    }
	  else
	    {
#pragma omp parallel num_threads(nThreads)
	      {
		int th = omp_get_thread_num();
		std::vector<Facet> localVec;		
		const auto it_end=mesh->simplexEnd(th,nThreads);
		for (auto it=mesh->simplexBegin(th,nThreads);it!=it_end;++it)
		  {
		    Simplex *s=(*it);
		    for (int i=0;i<Simplex::NNEI;++i)
		      {
			if (s->getNeighbor(i)>s)
			  {
			    FacetHandle h=functor(s->getFacetHandle(i));
			    if (h->isSet())
			      localVec.push_back(*h);			   
			  }
		      }
		  }
#pragma omp critical
		facetArr.insert(facetArr.end(),localVec.begin(),localVec.end());
	      }
	    }
	  useFacetArr=true;	 	  
	}
    }

    template <class F>
    void filter(const F &functor, int nThreads, hlp::ConstantType<SegmentHandle>)
    {
      if (!useSegArr)
	{
	  if (NDIM==2)
	    {	  
	      if (nThreads==1)
		{
		  const auto it_end=mesh->simplexEnd();
		  for (auto it=mesh->simplexBegin();it!=it_end;++it)
		    {
		      Simplex *s=(*it);
		      for (int i=0;i<Simplex::NSEG;++i)
			{
			  if (s->getNeighbor(i) > s)
			    {
			      SegmentHandle h=functor(s->getSegmentHandle(i));
			      if (h->isSet()) segArr.push_back(*h);
			    }
			}
		    }
		}
	      else // if (nThreads!=1)
		{
#pragma omp parallel num_threads(nThreads)
		  {
		    int th = omp_get_thread_num();
		    std::vector<Segment> localVec;		
		    const auto it_end=mesh->simplexEnd(th,nThreads);
		    for (auto it=mesh->simplexBegin(th,nThreads);it!=it_end;++it)
		      {
			Simplex *s=(*it);
			for (int i=0;i<Simplex::NSEG;++i)
			  {			    
			    if (s->getNeighbor(i)>s)
			      {
				SegmentHandle h=functor(s->getSegmentHandle(i));
				if (h->isSet()) localVec.push_back(*h);			   
			      }
			  }
		      }
#pragma omp critical
		    segArr.insert(segArr.end(),localVec.begin(),localVec.end());
		  }
		}
	      useSegArr=true;	 	  
	    }
	  else // if (NDIM!=2)
	    {
	      if (nThreads==1)
		{
		  const auto it_end=mesh->simplexEnd();
		  for (auto it=mesh->simplexBegin();it!=it_end;++it)
		    {		      
		      SegmentHandle handles[Simplex::NSEG];
		      int nOwned=(*it)->getOwnedLocalSegmentHandles
			(handles,options&NDNET_WithGhosts);
		      
		      for (int i=0;i<nOwned;++i)
			{
			  SegmentHandle h=functor(handles[i]);
			  if (h->isSet()) segArr.push_back(*h);			    
			}		      
		    }
		}
	      else // if (nThreads!=1)
		{
#pragma omp parallel num_threads(nThreads)
		  {
		    int th = omp_get_thread_num();
		    std::vector<Segment> localVec;		
		    const auto it_end=mesh->simplexEnd(th,nThreads);
		    for (auto it=mesh->simplexBegin(th,nThreads);it!=it_end;++it)
		      {			
			SegmentHandle handles[Simplex::NSEG];
			int nOwned=(*it)->getOwnedLocalSegmentHandles
			  (handles,options&NDNET_WithGhosts);
		      
			for (int i=0;i<nOwned;++i)
			  {
			    SegmentHandle h=functor(handles[i]);
			    if (h->isSet()) localVec.push_back(*h);			    
			  }		      			
		      }
#pragma omp critical
		    segArr.insert(segArr.end(),localVec.begin(),localVec.end());
		  }
		}
	      useSegArr=true;	 	  
	    } // NDIM
	} // usesegarr	
    } 
    
    void trimVertices(int nThreads=-1)
    {
      if (nThreads<0) nThreads=defaultNTh;
      if (!useSArr)
	{
	  vArr.clear();
	  useVArr=false;
	}
      else
	{	 	  
	  vIndexArr.resize(mesh->getNVertices(), 
			   std::numeric_limits<LocalIndex>::max());
	  	  
	  if (useVArr) 
	    {
#pragma omp parallel for num_threads(nThreads)
	      for (long i=0;i<vArr.size();++i)
		vIndexArr[vArr[i]->getLocalIndex()]=1;
	    }

#pragma omp parallel for num_threads(nThreads)
	  for (long j=0;j<sArr.size();++j)
	    for (int i=0;i<Simplex::NVERT;++i)
	      vIndexArr[sArr[j]->getVertex(i)->getLocalIndex()]=1;

#pragma omp parallel for num_threads(nThreads)
	  for (long j=0;j<facetArr.size();++j)
	    for (int i=0;i<Facet::NVERT;++i)
	      vIndexArr[facetArr[j].getVertex(i)->getLocalIndex()]=1;

#pragma omp parallel for num_threads(nThreads)
	  for (long j=0;j<segArr.size();++j)	  
	    for (int i=0;i<Segment::NVERT;++i)
	      vIndexArr[segArr[j].getVertex(i)->getLocalIndex()]=1;
	  /*
	  long delta=mesh->getNVertices();
#pragma omp parallel for num_threads(nThreads)
	  for (long j=0;j<gvArr.size();++j)
	    vIndexArr[delta+gvArr[j]->getLocalIndex()]=count++;

	  delta+=mesh->getNGhostVertices();
#pragma omp parallel for num_threads(nThreads)
	  for (long j=0;j<svArr.size();++j)
	    vIndexArr[delta+svArr[j]->getLocalIndex()]=count++;
	  */
	  long count=0;
	  for (long i=0;i<vIndexArr.size();++i)
	    if (vIndexArr[i]!=std::numeric_limits<LocalIndex>::max()) vIndexArr[i]=count++;
	 	  	    
	  vArr.resize(count);
#pragma omp parallel num_threads(nThreads)
	  {
	    int th = omp_get_thread_num();
	    std::vector<Segment> localVec;		
	    const auto it_end=mesh->vertexEnd(th,nThreads);
	    for (auto it=mesh->vertexBegin(th,nThreads);it!=it_end;++it)	  
	      {
		LocalIndex id=(*it)->getLocalIndex();
		if (vIndexArr[id]!=std::numeric_limits<LocalIndex>::max())
		  vArr[vIndexArr[id]]=(*it);
	      }
	  }
	  
	  useVIndexArr=true;
	  useVArr=true;	
	}
    }
    
    template <typename T>
    NDNET_UINT getIndex(T val)
    {
      if (useVIndexArr)
	return static_cast<NDNET_UINT>(vIndexArr[val]);
      else
	return static_cast<NDNET_UINT>(val);
    }

    template <typename V>
    NDNET_UINT getIndex(V *vertex)
    {
      if (useVIndexArr)
	return static_cast<NDNET_UINT>(vIndexArr[vertex->getLocalIndex()]);
      else
	return static_cast<NDNET_UINT>(vertex->getLocalIndex());
    }
    
  protected:
  
    std::vector<unsigned long> getCellsCount()
    {
      std::vector<unsigned long> nCells(NDIM+1,0);
    
      if (!useVArr)
	nCells[0]=mesh->getNVertices();
      else 
	nCells[0]=vArr.size();

      nCells[0]+=gvArr.size()+svArr.size();

      if (!useSArr)
	nCells[NDIM]=mesh->getNSimplices();
      else 
	nCells[NDIM]=sArr.size();

      nCells[NDIM]+=gsArr.size()+ssArr.size();
    
      nCells[NDIM-1]=facetArr.size();     
      nCells[1]+=segArr.size(); // merged with facets if (NDIM==2)            

      return nCells;
    }   

    Mesh *mesh;
    std::string fName;
    long options;

    bool useVArr;
    bool useSArr;
    bool useFacetArr;
    bool useSegArr;

    int defaultNTh;

    std::vector<Vertex*> vArr;
    std::vector<GhostVertex*> gvArr;
    std::vector<ShadowVertex*> svArr;

    std::vector<Simplex*> sArr;  
    std::vector<GhostSimplex*> gsArr;  
    std::vector<ShadowSimplex*> ssArr;  

    std::vector<Segment> segArr;
    std::vector<Facet> facetArr;
    
    std::vector<LocalIndex> vIndexArr;
    bool useVIndexArr;
  };

} // namespace IO
    
/** \}*/
#include "../internal/namespace.footer"
#endif

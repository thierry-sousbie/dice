#ifndef __DICE_VTK_AMR_WRITER_HXX__
#define __DICE_VTK_AMR_WRITER_HXX__

#include "../dice_globals.hxx"
#include "../tools/IO/myIO.hxx"
#include "../tools/helpers/helpers.hxx"

/**
 * @file 
 * @brief  Definition of a class to write AMR grids to VTK format
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup IO
 *   \{
 */

namespace IO {
  /**
   * \class VtkAmrGridWriterT
   * \brief a class to write AMR grids to VTK format
   * \tparam AMR The class of the AMR grid object to write in VTK format
   */
  template <class AMR, class CT=double, bool ForceHeaderType32=false> 
  class VtkAmrWriterT
  {
  public:
    typedef VtkAmrWriterT<AMR,CT,ForceHeaderType32> MyType;
    
    static const int NDIM=AMR::NDIM;
    typedef typename AMR::Voxel Voxel;
    typedef CT Float;
    typedef typename Voxel::Data Data;
    typedef long Int;

    static const int HEADER_TYPE= ((sizeof(unsigned long)==4)||(ForceHeaderType32))?32:64;
    typedef typename hlp::MinimalIntegerType<HEADER_TYPE,false>::Type HInt;

    VtkAmrWriterT(AMR *amr_, const char *fileName="vtkAmr"):
      amr(amr_)
    {      
      fName = toFilename(fileName);
      vtk_cellT=(NDIM==2)?8:11; // pixel or voxel     
    }
    
    virtual ~VtkAmrWriterT()
    {}

    void write(bool fast=true, bool quiet=false)
    {
      //if (mpiCom->size()==1)
      if (!quiet)
	glb::console->printFlush<LOG_STD>("Dumping %dD AMR grid to VTK file '%s' ... ",
					  NDIM,fName.c_str());
      // else
      // 	glb::console->printFlush<LOG_STD>("Dumping %dD AMR grid to VTK file '%s' (%d files) ... ",
      // 					  NDIM,fName.c_str(),mpiCom->size());

      //printf("Writing to file %s ... ",fName.c_str());fflush(0);
      FILE *f=fopen(fName.c_str(),"w");      
      myIO::BinaryWriterT<> writer(f);
      write<>(&writer,fast);
      fclose(f);
      
      if (!quiet) glb::console->print<LOG_STD>("done.\n");
    }

    template <class BW>
    void write(BW *bWriter, bool fast=true)
    {      
      //amr->assignVerticesToLeaves();      
      writeHeader(bWriter->getFilePtr(),true,fast);
      
      //unsigned long nVertex=amr->getUniqueVerticesCount();
      unsigned long nCells=amr->getNLeaves();

      fprintf(bWriter->getFilePtr(),"<AppendedData encoding=\"raw\">\n_");
      
      // coordinates
      bWriter->write(&size[0]);
      amr->visitTree(RetrieveStructureFastVisitorT<AMR,BW,0>(amr,bWriter),1);

      // connectivity
      bWriter->write(&size[1]); 
      for (Int i=0;i<(nCells<<NDIM);++i)
	bWriter->write(&i);

      // offsets
      bWriter->write(&size[2]);
      for (Int i=(1<<NDIM);i<(nCells+1)*(1<<NDIM);i+=(1<<NDIM))       
	bWriter->write(&i);

      // type
      bWriter->write(&size[3]);
      for (Int i=0;i<nCells;i++)
	bWriter->write(&vtk_cellT);
      
      HInt s=sizeof(Data)*nCells;
      bWriter->write(&s);
      amr->visitTree(RetrieveStructureFastVisitorT<AMR,BW,1>(amr,bWriter),1);
      
      fprintf(bWriter->getFilePtr(),"</AppendedData>\n");
      fprintf(bWriter->getFilePtr(),"</VTKFile>\n");

      /*
      std::vector<Float> coords(nVertex*3,0);
      std::vector<Voxel*> voxelId(nVertex*(1<<NDIM));
      std::vector<unsigned char> vertexIndex(nVertex*(1<<NDIM)); 
      std::vector<double> data(nCells);

      // coords.assign(nVertex*3,0);
      // voxelId.resize(nVertex*(1<<NDIM));
      // vertexIndex.resize(nVertex*(1<<NDIM));
	
      visitor.init(amr,&coords[0],&voxelId[0],&vertexIndex[0]);
      amr->visitTree(visitor,1);
     
      std::vector<Int> connectivity(nCells<<NDIM);
      
      for (long i=0;i<connectivity.size();++i)
	connectivity[i]=i%nVertex;

      fprintf(bWriter->getFilePtr(),"<AppendedData encoding=\"raw\">\n_");
      bWriter->write(&size[0]);
      bWriter->write(&coords[0],coords.size());
      bWriter->write(&size[1]); 
      bWriter->write(&connectivity[0],connectivity.size());
      bWriter->write(&size[2]);
      for (Int i=(1<<NDIM);i<(nCells+1)*(1<<NDIM);i+=(1<<NDIM))       
	bWriter->write(&i);
      bWriter->write(&size[3]);
      for (Int i=0;i<nCells;i++)
	bWriter->write(&vtk_cellT);
      
      unsigned int s=sizeof(double)*nCells;
      bWriter->write(&s);
      bWriter->write(&data[0],data.size());
      
      fprintf(bWriter->getFilePtr(),"</AppendedData>\n");
      fprintf(bWriter->getFilePtr(),"</VTKFile>\n");
      */
    }

  protected:

    std::string toFilename(const char *name="vtkAmr")
    {
      return std::string(name) + std::string(".amr.vtu");
    }
    
    void writeHeader(FILE *f, bool binary, bool fast)
    {
      char format[255];      
      long offset=0;
      long sizeindex=0;    
 
      unsigned long nVertex=amr->getUniqueVerticesCount();
      unsigned long nCells=amr->getNLeaves();

      if (fast) nVertex=(nCells<<NDIM);

      sprintf(format,"%s",(binary)?"appended":"ascii");
      fprintf(f,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt%d\">\n",HEADER_TYPE);
      fprintf(f,"  <UnstructuredGrid>\n");
      fprintf(f,"    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n",
	      (long)nVertex,(long)nCells);

      fprintf(f,"<Points>\n");
      if (binary)
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\" />\n",
		(long)3,"coords",sizeof(Float)*8,format,offset);  
      else
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" >\n",
		(long)3,"coords",sizeof(Float)*8,format);  
      
      if (binary)
	{	
	  size[sizeindex]=sizeof(Float)*3*nVertex;
	  offset+=size[sizeindex++]+sizeof(HInt);
	  fprintf(f,"</Points>\n");
	}
      else
	{
	  //Print array here ...
	  //fprintf(f,"%g %g %g\n",v[0],v[1],v[2]);
	  for (int i=0;i<nVertex;++i) fprintf(f,"%g %g %g\n",
					      1.0+((5*i)%6),
					      1.0+((2*i+1)%5),
					      1.0+((3*i+2)%7));
	  fprintf(f,"</DataArray>\n</Points>\n");
	}
      
      fprintf(f,"<Cells>\n");
      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",sizeof(Int)*8,"connectivity",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",sizeof(Int)*8,"connectivity",format);

      if (binary) 
	{
	  size[sizeindex]=sizeof(Int)*(1<<NDIM)*nCells;
	  offset+=size[sizeindex++]+sizeof(HInt);	 
	}
      else
	{
	  for (long i=0;i<nCells;i++) 
	    {
	      fprintf(f,"0 1 2 3\n");
	      //fprintf(f," %ld",(long)net->f_vertexIndex[cell_type][i]);
	      //if ((i+1)%(1<<NDIM) == 0) fprintf(f,"\n");
	    }
	  fprintf(f,"</DataArray>\n");
	}
      
      if (binary)
	 fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",sizeof(Int)*8,"offsets",format,offset);
       else
	 fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",sizeof(Int)*8,"offsets",format);

      if (binary) 
	{
	  size[sizeindex]=sizeof(Int)*(nCells);
	  offset+=size[sizeindex++]+sizeof(HInt);
	}
      else
	{
	  for (long i=(1<<NDIM);i<(nCells+1)*(1<<NDIM);i+=(1<<NDIM))       
	    fprintf(f,"%ld\n",i);
	  // for (ct=(cell_type+1);ct<(cell_type+1)*(ncells+1);ct+=(cell_type+1)) 
	  //   fprintf(f,"%ld\n",(long)ct);
	  fprintf(f,"</DataArray>\n");
	}
      
       if (binary)
	fprintf(f,"<DataArray type=\"UInt8\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n","types",format,offset);
      else
	fprintf(f,"<DataArray type=\"UInt8\" Name=\"%s\" format=\"%s\">\n","types",format);
      
      if (binary) 
	{
	  size[sizeindex]=sizeof(unsigned char)*nCells;
	  offset+=size[sizeindex++]+sizeof(HInt);
	}
      else
	{
	  for (long i=0;i<nCells;i++) fprintf(f,"%u\n",(unsigned int)vtk_cellT);
	  //for (i=0;i<ncells;i++) fprintf(f,"%u\n",(unsigned int)vtk_cellT);
	  fprintf(f,"</DataArray>\n");
	}
      fprintf(f,"</Cells>\n");

      // CELL DATA            
      char name[]="value";
      fprintf(f,"<CellData>\n");
      if (binary)
	{		  
	  fprintf(f,"<DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\"/>\n",name,sizeof(Data)*8,format,offset);
	  size[sizeindex]=sizeof(Data)*nCells;
	  offset+=size[sizeindex++]+sizeof(HInt);
	  //dataSize[nData]=net->nfaces[cell_type];
	  //dataPtr[nData++]=(void*)d;		  
	}
      else
	{
	  fprintf(f,"<DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" >\n",name,sizeof(Data)*8,format);
		  
	  //for (long j=0;j<nCells;j++) fprintf(f,"%g\n",);	
	  fprintf(f,"</DataArray>\n");
	}
     
      fprintf(f,"</CellData>\n"); 
      
      fprintf(f,"   </Piece>\n");
      fprintf(f,"  </UnstructuredGrid>\n");
      if (!binary) fprintf(f,"</VTKFile>\n");
    }

  private:

    template <class AMRV, class W, int WD>
    class RetrieveStructureFastVisitorT
    {     
      typedef typename AMRV::Voxel Voxel;
      static const long NDIM=AMRV::NDIM;
    public:
      RetrieveStructureFastVisitorT(AMRV *amr_, W *bWriter_):amr(amr_),bWriter(bWriter_)
      {}

      static void initialize(Voxel *rootVoxel) {}
  
      bool visit(Voxel *voxel) const
      { 
	if (voxel->isLeaf())
	  {
	    if (WD)
	      bWriter->write(&(voxel->data));
	    else
	      {
		static const Float zero=0;
		Float corners[2][NDIM];
		amr->index2CornerCoordsAndOpp(voxel->getIndex(),voxel->getLevel(),
					      &corners[0][0],&corners[1][0]);	   
		for (int i=0;i<(1L<<NDIM);++i)
		  {
		    for (int j=0;j<NDIM;++j)
		      bWriter->write(&(corners[(i>>j)&1][j]));
		    for (int j=NDIM;j<3;++j)
		      bWriter->write(&zero);
		  }
	      }	  
	    return false;
	  }
	return true;      
      }      

      static bool visit(Voxel *voxel,int i) 
      {
	return true;
      }

      static void visited(Voxel *voxel) 
      {
      }

      static void visited(Voxel *voxel, int i) 
      {
      }

      void init(AMRV *amr_, W *bWriter_) 
      {
	amr=amr_;
	bWriter=bWriter_;
      }
      
    private:
      mutable AMRV *amr;
      mutable W *bWriter;
    };
    /*
    template <class AMRV>
    class RetrieveStructureVisitorT
    {
      typedef typename AMRV::Voxel Voxel;
      static const long NDIM=AMRV::NDIM;
    public:
      RetrieveStructureVisitorT():amr(NULL) {}

      static void initialize(Voxel *rootVoxel) {}

      bool visit(Voxel *voxel)
      { 
	if (voxel->isLeaf())
	  {	    	    
	    int dir[1<<NDIM][NDIM];	   
	    for (int i=0;i<(1L<<NDIM);++i)
	      {
		if (voxel->getVertexFlag(i))
		  {		    
		    amr->getVertexNeighbors(voxel,i,voxelId);	
		    amr->index2CornerCoords(voxel->getIndex(),voxel->getLevel(),coords,i);
		    amr->getVertexNeighborsDirection(i,dir);
		     
		    for (int j=0;j<(1<<NDIM);++j)
		      {
			(*vertexIndex)=0;
			for (int k=0;k<NDIM;++k)
			  (*vertexIndex) |= ((1-dir[j][k])>>1)<<k;
			++vertexIndex;
		      }

		    coords += 3;
		    voxelId += (1<<NDIM);		    
		    //data++;
		  }
	      }
	    return false;
	  }
	return true;      
      }

      static void visited(Voxel *voxel) {}

      void init(AMRV *amr_, Float *coords_, Voxel **voxelId_, unsigned char *vertexIndex_) 
      {
	amr=amr_;
	coords=coords_;
	voxelId=voxelId_;
	vertexIndex = vertexIndex_;
      }
      
    private:
      mutable AMRV *amr;
      
      Float *coords;
      Voxel **voxelId;
      unsigned char *vertexIndex; 
    };
    */
    AMR *amr;
    std::string fName;
    //RetrieveStructureVisitorT<AMR> visitor;    

    // file data
    unsigned char vtk_cellT;
    //long offset;
    HInt size[50];
    //int sizeindex;
  };

}


/** \}*/
#include "../internal/namespace.footer"
#endif

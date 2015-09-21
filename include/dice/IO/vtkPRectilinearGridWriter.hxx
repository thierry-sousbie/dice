#ifndef __DICE_VTK_P_RECTILINEAR_GRID_WRITER_HXX__
#define __DICE_VTK_P_RECTILINEAR_GRID_WRITER_HXX__

#include "../dice_globals.hxx"
#include "../tools/MPI/mpiCommunication.hxx"

#include "../tools/IO/myIO.hxx"
#include "./vtkRectilinearGridWriter.hxx"

#include "../grid/valLocationType.hxx"

#include "../tools/helpers/helpers.hxx"

/**
 * @file 
 * @brief  Definition of a class to write rectilinear grids grids to VTK format
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup IO
 *   \{
 */

namespace IO {
  /**
   * \class VtkPRectilinearGridWriterT
   * \brief a class to write parallel rectilinear grids to VTK format
   * \tparam AMR The class of the AMR grid object to write in VTK format
   */
  template <class G,typename CT=double, bool ForceHeaderType32=false> 
  class VtkPRectilinearGridWriterT
  {
  public:
    typedef VtkPRectilinearGridWriterT<G,CT,ForceHeaderType32> MyType;
    typedef typename G::LocalGrid LocalGrid;
    typedef VtkRectilinearGridWriterT<LocalGrid,CT,ForceHeaderType32>  LocalWriter;
    static const int NDIM=G::NDIM;  

    typedef G Grid;
    typedef CT Float;
    typedef typename G::Data Data;
    typedef long Int;

    static const int HEADER_TYPE= ((sizeof(unsigned long)==4)||(ForceHeaderType32))?32:64;
    typedef typename hlp::MinimalIntegerType<HEADER_TYPE,false>::Type HInt;

    //typedef typename G::ValLocationType  ValLocationType;
    //typedef typename G::ValLocationTypeV ValLocationTypeV;

    VtkPRectilinearGridWriterT(Grid *g, MpiCommunication *com, 
			       const char *globalFileName,
			       const char *fileNameFormat):
      grid(g),
      mpiCom(com),      
      globalFName(globalFileName),
      fNameFormat(fileNameFormat)
    {     
      
      if (mpiCom->size()>1)
	{
	  char tmp[255];	 
	  sprintf(tmp,fileNameFormat,mpiCom->rank());
	  localFName = tmp;
	}
      else localFName = globalFName;
      
      fName = toFilename(globalFileName);      
    }

    virtual ~VtkPRectilinearGridWriterT()
    {}

    
    void write(bool quiet=false)
    {      
      if (mpiCom->size()<=1)
	{
	  LocalWriter localWriter(grid->getLocalGrid(),localFName.c_str());
	  localWriter.write(quiet);
	}
      else
	{
	  if (!quiet)
	    glb::console->printFlush<LOG_STD>("Dumping %dD parallel rectilinear grid to VTK file '%s' (%d chunks) ... ",NDIM,fName.c_str(),mpiCom->size());
	  
	  LocalWriter localWriter(grid->getLocalGrid(),localFName.c_str());
	  localWriter.write(true);
            
	  if (mpiCom->rank()==0)
	    {
	      FILE *f=fopen(fName.c_str(),"w");      
	      myIO::BinaryWriterT<> writer(f);
	      write<>(&writer);
	      fclose(f);
	    }
      
	  if (!quiet) glb::console->print<LOG_STD>("done.\n");
	}
    }
   
  protected:

    template <class BW>
    void write(BW *bWriter)
    {
      setGridInfo();
      writeHeader(bWriter->getFilePtr(),true);
    }

    std::string toFilename(const char *name="vtkPRGrid")
    {
      return std::string(name) + std::string(".pvtr");
    }

    void writeHeader(FILE *f, bool binary)
    {
      char format[255];      
      //long offset=0;
      //long sizeindex=0;      

      //LocalGrid *lGrid = grid->getLocalGrid();
                     
      sprintf(format,"%s",(binary)?"appended":"ascii");
      fprintf(f,"<VTKFile type=\"PRectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt%d\">\n",HEADER_TYPE);
      fprintf(f," <PRectilinearGrid WholeExtent = \"%d %d %d %d %d %d\" GhostLevel=\"0\">\n",
	      wholeExtent[0][0],wholeExtent[0][1],
	      wholeExtent[1][0],wholeExtent[1][1],
	      wholeExtent[2][0],wholeExtent[2][1]);

      
      char dataType[255];
      if (onPoints)
	{	  
	  sprintf(dataType,"PPointData");
	}
      else
	{
	  fprintf(f,"   <PPointData>\n");
	  fprintf(f,"   </PPointData>\n");
	  sprintf(dataType,"PCellData");
	}
      char name[255];
      fprintf(f,"   <%s>\n",dataType);
      for (int i=0;i<nFields;++i)
	{	  
	  if (nFields==1) sprintf(name,"%s",grid->getName().c_str());
	  else sprintf(name,"%s_%2.2d",grid->getName().c_str(),i);
	  fprintf(f,"    <PDataArray Name=\"%s\" type=\"Float%2.2ld\"/>\n",name,sizeof(Data)*8);
	}

      fprintf(f,"   </%s>\n",dataType);

      if (onPoints)
	{
	  fprintf(f,"   <PCellData>\n");
	  fprintf(f,"   </PCellData>\n");
	}
      
      fprintf(f,"   <PCoordinates>\n");
      fprintf(f,"    <PDataArray Name=\"%s\" type=\"Float%2.2ld\" />\n",
	      "coordsX",sizeof(Float)*8);  
      fprintf(f,"    <PDataArray Name=\"%s\" type=\"Float%2.2ld\" />\n",
	      "coordsY",sizeof(Float)*8);  
      fprintf(f,"    <PDataArray Name=\"%s\" type=\"Float%2.2ld\" />\n",
	      "coordsZ",sizeof(Float)*8); 
      fprintf(f,"   </PCoordinates>\n");

      for (int i=0;i<mpiCom->size();++i)
	{
	  char tmp[1024];
	  setExtent(i);
	  sprintf(tmp,fNameFormat.c_str(),i);
	  std::string name = LocalWriter::toFilename(tmp,true);
	  fprintf(f,"  <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>\n",
		  extent[0][0],extent[0][1],
		  extent[1][0],extent[1][1],
		  extent[2][0],extent[2][1],
		  name.c_str());
	}     
      fprintf(f," </PRectilinearGrid>\n");
      fprintf(f,"</VTKFile>\n");      
    }

    void setGridInfo()
    {
      onPoints= (grid->getValLocation(0)==ValLocationTypeV::VERTEX); 
      nFields = grid->getNFields();
      grid->getBoundingBox(bBox);  
      for (int i=0;i<NDIM;++i) 
	{
	  wholeExtent[i][0]=0;
	  wholeExtent[i][1]=wholeExtent[i][0]+grid->getResolution(i);
	}
      for (int i=NDIM;i<3;++i) 
	{
	  wholeExtent[i][0]=0;
	  wholeExtent[i][1]=1;
	}
    }

    void setExtent(int i)
    {
      typedef typename Grid::Params Params;
      const Params &params = grid->getRemoteGridParams(i);
      for (int i=0;i<NDIM;++i) 
	{
	  extent[i][0]=params.position[i];
	  extent[i][1]=extent[i][0]+params.resolution[i];
	}
      for (int i=NDIM;i<3;++i)
	{	 
	  extent[i][0]=0;extent[i][1]=1;
	}
    }

    Grid *grid;
    MpiCommunication *mpiCom;

    std::string fName;
    std::string localFName;
    std::string globalFName;
    std::string fNameFormat;

    double bBox[3][2];    

    int wholeExtent[3][2];
    int extent[3][2];

    //long nCells;
    //long nVertices;
    //long nValues;

    //double cellDims[3];
    //double vertexDims[3];
    //double valueDims[3];

    int nFields;
    bool onPoints;
  };

} // namespace IO

/** \}*/
#include "../internal/namespace.footer"
#endif

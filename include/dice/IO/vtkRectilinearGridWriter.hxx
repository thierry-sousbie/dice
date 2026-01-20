#ifndef __DICE_VTK_RECTILINEAR_GRID_WRITER_HXX__
#define __DICE_VTK_RECTILINEAR_GRID_WRITER_HXX__

#include "../dice_globals.hxx"
#include "../tools/IO/myIO.hxx"

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

namespace IO
{

  /**
   * \class VtkRectilinearGridWriterT
   * \brief a class to write Rectilinear grids to VTK format
   * \tparam AMR The class of the AMR grid object to write in VTK format
   */
  template <class G, typename CT = double, bool ForceHeaderType32 = false>
  class VtkRectilinearGridWriterT
  {
  public:
    typedef VtkRectilinearGridWriterT<G, CT, ForceHeaderType32> MyType;

    static const int NDIM = G::NDIM;

    typedef G Grid;
    typedef CT Float;
    typedef typename G::Data Data;
    typedef long Int;

    static const int HEADER_TYPE = ((sizeof(unsigned long) == 4) || (ForceHeaderType32)) ? 32 : 64;
    typedef typename hlp::MinimalIntegerType<HEADER_TYPE, false>::Type HInt;

    // typedef typename G::ValLocationType  ValLocationType;
    // typedef typename G::ValLocationTypeV ValLocationTypeV;

    VtkRectilinearGridWriterT(Grid *g, const char *fileName = "vtkRGrid") : grid(g)
    {
      fName = toFilename(fileName);
    }

    virtual ~VtkRectilinearGridWriterT()
    {
    }

    void write(bool quiet = false)
    {
      if (!quiet)
        glb::console->printFlush<LOG_STD>("Dumping %dD rectilinear grid to VTK file '%s' ... ", NDIM, fName.c_str());

      FILE *f = fopen(fName.c_str(), "w");
      myIO::BinaryWriterT<> writer(f);
      write<>(&writer);
      fclose(f);

      if (!quiet)
        glb::console->print<LOG_STD>("done.\n");
    }

    template <class BW>
    void write(BW *bWriter)
    {
      setGridInfo();
      writeHeader(bWriter->getFilePtr(), true);

      fprintf(bWriter->getFilePtr(), "<AppendedData encoding=\"raw\">\n_");

      Data *data = grid->getDataPtr();
      int sizeIndex = 0;

      // data
      if (nFields == 1)
      {
        bWriter->write(&size[sizeIndex++]);
        bWriter->write(data, nValues);
      }
      else if (Grid::IS_INTERLEAVED)
      {
        for (int j = 0; j < nFields; ++j)
        {
          bWriter->write(&size[sizeIndex++]);

          long index = j;
          for (long i = 0; i < nValues; i++)
          {
            bWriter->write(&data[index]);
            index += nFields;
          }
        }
      }
      else
      {
        for (long j = 0; j < nFields; ++j)
        {
          bWriter->write(&size[sizeIndex++]);

          Data *curData = data + nValues * j;
          bWriter->write(curData, nValues);
        }
      }

      // coordinates
      bWriter->write(&size[sizeIndex++]);
      const std::vector<Float> &coord0 = grid->getVertexCoord(0);
      bWriter->write(&coord0[0], coord0.size());
      // char tmp[10000];
      // for (int i=0;i<coord0.size();++i) sprintf("")
      bWriter->write(&size[sizeIndex++]);
      const std::vector<Float> &coord1 = grid->getVertexCoord(1);
      bWriter->write(&coord1[0], coord1.size());

      bWriter->write(&size[sizeIndex++]);
      if (NDIM > 2)
      {
        const std::vector<Float> &coord2 = grid->getVertexCoord(2);
        bWriter->write(&coord2[0], coord2.size());
      }
      else
      {
        std::vector<Float> coord2(2, 0);
        bWriter->write(&coord2[0], coord2.size());
      }

      fprintf(bWriter->getFilePtr(), "</AppendedData>\n");
      fprintf(bWriter->getFilePtr(), "</VTKFile>\n");
    }

    static std::string toFilename(const char *name = "vtkRGrid", bool stripPath = false)
    {
      int start = 0;
      if (stripPath)
      {
        int len = strlen(name);
        for (int i = 0; i < len; ++i)
          if (name[i] == '/')
            start = i + 1;
      }
      return std::string(&name[start]) + std::string(".vtr");
    }

  protected:
    void writeHeader(FILE *f, bool binary)
    {
      char format[255];
      long offset = 0;
      long sizeindex = 0;

      sprintf(format, "%s", (binary) ? "appended" : "ascii");
      fprintf(f, "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt%d\">\n", HEADER_TYPE);
      fprintf(f, " <RectilinearGrid WholeExtent = \"%d %d %d %d %d %d\">\n",
              extent[0][0], extent[0][1],
              extent[1][0], extent[1][1],
              extent[2][0], extent[2][1]);
      fprintf(f, "  <Piece Extent=\"%d %d %d %d %d %d\">\n",
              extent[0][0], extent[0][1],
              extent[1][0], extent[1][1],
              extent[2][0], extent[2][1]);

      char dataType[255];
      if (onPoints)
      {
        sprintf(dataType, "PointData");
      }
      else
      {
        fprintf(f, "   <PointData>\n");
        fprintf(f, "   </PointData>\n");
        sprintf(dataType, "CellData");
      }

      char name[255];
      fprintf(f, "   <%s>\n", dataType);
      for (int i = 0; i < nFields; ++i)
      {
        if (nFields == 1)
          sprintf(name, "%s", grid->getName().c_str());
        else
          sprintf(name, "%s_%2.2d", grid->getName().c_str(), i);

        if (binary)
        {
          fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\"/>\n", name, sizeof(Data) * 8, format, offset);
          size[sizeindex] = sizeof(Data) * nValues;
          offset += size[sizeindex++] + sizeof(HInt);
        }
        else
        {
          fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" >\n", name, sizeof(Data) * 8, format);
          // for (long j=0;j<nCells;j++) fprintf(f,"%g\n",);
          fprintf(f, "    </DataArray>\n");
        }
      }
      fprintf(f, "   </%s>\n", dataType);

      if (onPoints)
      {
        fprintf(f, "   <CellData>\n");
        fprintf(f, "   </CellData>\n");
      }

      fprintf(f, "   <Coordinates>\n");
      if (binary)
      {
        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\" />\n", "coordsX", sizeof(Float) * 8, format, offset);
        size[sizeindex] = sizeof(Float) * vertexDims[0];
        offset += size[sizeindex++] + sizeof(HInt);

        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\" />\n", "coordsY", sizeof(Float) * 8, format, offset);
        size[sizeindex] = sizeof(Float) * vertexDims[1];
        offset += size[sizeindex++] + sizeof(HInt);

        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" offset=\"%ld\" />\n", "coordsZ", sizeof(Float) * 8, format, offset);
        size[sizeindex] = sizeof(Float) * vertexDims[2];
        offset += size[sizeindex++] + sizeof(HInt);
      }
      else
      {
        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" />\n", "coordsX", sizeof(Float) * 8, format);
        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" />\n", "coordsY", sizeof(Float) * 8, format);
        fprintf(f, "    <DataArray Name=\"%s\" type=\"Float%2.2ld\" format=\"%s\" />\n", "coordsZ", sizeof(Float) * 8, format);
      }
      fprintf(f, "   </Coordinates>\n");

      fprintf(f, "  </Piece>\n");
      fprintf(f, " </RectilinearGrid>\n");
      if (!binary)
        fprintf(f, "</VTKFile>\n");
    }

  protected:
    void setGridInfo()
    {
      grid->getBoundingBox(bBox);
      grid->getGlobalBoundingBox(gBBox);

      onPoints = (grid->getValLocation(0) == ValLocationTypeV::VERTEX);
      nFields = grid->getNFields();

      nCells = 1;
      nVertices = 1;
      nValues = 1;
      for (int i = 0; i < NDIM; ++i)
      {
        cellDims[i] = grid->getCellCoord(i).size();
        vertexDims[i] = cellDims[i] + 1;
        valueDims[i] = grid->getValueCoord(i).size();
        nCells *= cellDims[i];
        nVertices *= vertexDims[i];
        nValues *= (onPoints) ? vertexDims[i] : cellDims[i];
      }
      if (NDIM < 3)
      {
        cellDims[2] = 1;
        vertexDims[2] = 2;
        valueDims[2] = (onPoints) ? 2 : 1;
        nVertices *= 2;
        if (onPoints)
          nValues *= 2;
      }

      for (int i = 0; i < NDIM; ++i)
      {
        wholeExtent[i][0] = 0;
        wholeExtent[i][1] = wholeExtent[i][0] + grid->getGlobalResolution(i);

        extent[i][0] = grid->getPosition(i);
        extent[i][1] = extent[i][0] + grid->getResolution(i);
      }
      for (int i = NDIM; i < 3; ++i)
      {
        wholeExtent[i][0] = 0;
        wholeExtent[i][1] = 1;
        extent[i][0] = 0;
        extent[i][1] = 1;
      }
    }

    Grid *grid;
    std::string fName;
    HInt size[50];

    double bBox[3][2];
    double gBBox[3][2];

    int wholeExtent[3][2];
    int extent[3][2];

    long nCells;
    long nVertices;
    long nValues;

    double cellDims[3];
    double vertexDims[3];
    double valueDims[3];

    int nFields;
    bool onPoints;
  };

} // namespace IO

/** \}*/
#include "../internal/namespace.footer"
#endif

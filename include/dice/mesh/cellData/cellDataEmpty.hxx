#ifndef __CELL_DATA_EMPTY_HXX__
#define __CELL_DATA_EMPTY_HXX__

#include "../../mesh/cellData/cellDataMacros.hxx"

#include "../../internal/namespace.header"

class VertexDataEmpty
{
public:
  INIT_CELL_DATA(VertexDataEmpty, Vertex)
};
// FINALIZE_CELL_DATA( VertexDataEmpty )

class SimplexDataEmpty
{
public:
  INIT_CELL_DATA(SimplexDataEmpty, Simplex)
};
// FINALIZE_CELL_DATA( SimplexDataEmpty )

#include "../../internal/namespace.footer"
#endif

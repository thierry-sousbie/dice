#ifndef __REGULAR_GRID_SLICER_BASE_HXX__
#define __REGULAR_GRID_SLICER_BASE_HXX__

#include "../valLocationType.hxx"
#include "../../internal/namespace.header"

namespace internal {

  template <class G>
  class RegularGridSlicerBaseT {
  public:
    typedef RegularGridSlicerBaseT<G> MyType;

    static const int NDIM = G::NDIM;

    typedef G Grid;
    typedef typename Grid::Params Params;
    typedef typename Grid::Scale Scale;
    //typedef typename Grid::ValLocationType ValLocationType;
    //typedef typename Grid::ValLocationTypeV ValLocationTypeV;
    typedef typename Scale::ScaleTypeV ScaleTypeV;
    typedef typename Scale::ScaleType  ScaleType;
 
  protected:
  
    static void getChunkIndices(const Params &gp, int index, int size, 
				int dir, int &n0, int &n1)
    {
      double deltaD;
      //long deltaI;
      int gp_res=gp.resolution[dir];
      ValLocationType vl = gp.valLocation[dir];

      if (vl == ValLocationTypeV::VERTEX)
	gp_res++;
    
      deltaD=(double(gp_res)/size);      
      //deltaI=long(deltaD*(index+1))-long(deltaD*index);
      n0 = long(deltaD*index);
      n1 = long(deltaD*(index+1));

      if (index==0) n0=0;
      if (index==size-1) n1=gp_res;
    }

    static Params divide(const Params &gp, long n0, long n1,
			 int lowMargin, int highMargin, bool periodic, 
			 int dir)
    {
      Params out = gp;    
      ValLocationType vl = gp.valLocation[dir];
    
      int &out_lowMargin=out.lowMargin[dir];
      int &out_highMargin=out.highMargin[dir];    
      int &out_res=out.resolution[dir];

      int gp_res=gp.resolution[dir];
      const ScaleType &gp_scale=gp.scale[dir];

      if (vl == ValLocationTypeV::VERTEX)
	gp_res++;

      if (periodic)
	{
	  out_lowMargin=lowMargin;
	  out_highMargin=highMargin;
	}
      else
	{
	  if (n0>=lowMargin) out_lowMargin=lowMargin;
	  else out_lowMargin=n0;

	  if (gp_res-n1>=highMargin) out_highMargin=highMargin;
	  else out_highMargin=gp_res-n1;
	}    

      n0 -= out_lowMargin;
      n1 += out_highMargin;
    
      out_res = n1-n0;    

      if (vl == ValLocationTypeV::VERTEX)
	out_res--;

      out.x0[dir]=Scale::valueAt
	(gp.x0[dir],gp.x0[dir]+gp.delta[dir],gp.resolution[dir],
	 n0,gp_scale,ValLocationTypeV::VERTEX);

      out.delta[dir]=Scale::valueAt
	(gp.x0[dir],gp.x0[dir]+gp.delta[dir],gp.resolution[dir],
	 n0+out_res,gp_scale,ValLocationTypeV::VERTEX);

      out.delta[dir]-=out.x0[dir];
    
      if (n0==0) out.x0[dir]=gp.x0[dir];
      if (n0+out_res==gp.resolution[dir]) 
	{
	  out.delta[dir] = gp.x0[dir]+gp.delta[dir] - out.x0[dir];
	}
        
      out.position[dir]=n0;
      
      if (glb::debug)
	glb::console->print<LOG_DEBUG>
	  (" Regular slicer base: index[%ld,%ld] :(%e,%e,%d)==> [%g,%g] res=%d\n",
	   n0,n1,gp.x0[dir],gp.x0[dir]+gp.delta[dir],gp.resolution[dir],
	   out.x0[dir],out.x0[dir]+out.delta[dir],out.resolution[dir]);
      
      return out;
    }    
  };

} // internal

#include "../../internal/namespace.footer"
#endif

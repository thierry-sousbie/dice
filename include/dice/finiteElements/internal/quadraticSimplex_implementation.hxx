#ifndef __QUADRATIC_SIMPLEX_IMPLEMENTATION_HXX__
#define __QUADRATIC_SIMPLEX_IMPLEMENTATION_HXX__

#include <Eigen/Dense>

#include "../../internal/namespace.header"

namespace internal {
  namespace quadSimplex {    
         
    template <int ND, int ORDER> struct GetIndexHelper;
    template <int ND, bool HIERACHICAL, int ORDER> struct ShapeFunctionsHelper;
    template <int ND, int NDW, bool HIERACHICAL, int ORDER> class TransformHelper;

    template <int NDIM, class M> 
    void standardToHierarchical(M &cM)
    {
      static const int NSEG = (NDIM*(NDIM+1))/2;
      static const int NVERT = NDIM+1;
      int v0=0;
      int v1=1;
      for (int i=0;i<NSEG;++i)
	{ 
	  cM.col(NVERT+i)-=
	    (cM.col(v0)+cM.col(v1))*0.5;
	  v1++;
	  if (v1==NVERT) 
	    {
	      v0++;
	      v1=v0+1;
	    }
	}
    }

    template <int NDIM, class M> 
    void standardToHierarchical(M &cM, int seg)
    {
      static const int NVERT = NDIM+1;
      int v0=0;
      int v1=seg;
      int delta=NDIM;

      for (;v1>delta;v1-=delta,v0++)
	{delta--;}
      
      cM.col(NVERT+seg)-=
	(cM.col(v0)+cM.col(v1))*0.5;     
    }

    template <int NDIM,class M> 
    void hierarchicalToStandard(M &cM)
    {
      static const int NSEG = (NDIM*(NDIM+1))/2;
      static const int NVERT = NDIM+1;
      int v0=0;
      int v1=1;
      for (int i=0;i<NSEG;++i)
	{ 
	  cM.col(NVERT+i)+=
	    (cM.col(v0)+cM.col(v1))*0.5;
	  v1++;
	  if (v1==NVERT) 
	    {
	      v0++;
	      v1=v0+1;
	    }
	}
    }

    template <int NDIM,class M> 
    void hierarchicalToStandard(M &cM, int seg)
    {
      static const int NVERT = NDIM+1;
      int v0=0;
      int v1=1;
      int delta=NDIM;
      for (v1=seg;v1>delta;v1-=delta,v0++)
	{delta--;}
      
      cM.col(NVERT+seg)+=
	(cM.col(v0)+cM.col(v1))*0.5;     
    }

    template <int NDIM>
    struct GetIndexHelper<NDIM,2>
    {
      static int get(int vi, int si)
      {
	if (si<0) return vi;
	else if (vi>si) std::swap(vi,si);

	return NDIM*(vi+1) - (vi*(vi-1))/2 + si-vi;    
      }      
    };

    template <int NDIM>
    struct ShapeFunctionsHelper<NDIM,false,2>
    {
      template <class I, class O>
      static void evalShapeFunctions(const I in, O &out)
      {
	for (int i=0;i<NDIM+1;++i)
	  out[i]=in[i]*(in[i]*2.0-1.0);

	int n=NDIM+1;
	for (int i=0;i<NDIM;++i)
	  for (int j=i+1;j<NDIM+1;++j,++n)
	    out[n]=in[i]*in[j]*4.0;
      }
    };

    template <int NDIM>
    struct ShapeFunctionsHelper<NDIM,true,2>
    {
      template <class I, class O>
      static void evalShapeFunctions(const I in, O &out)
      {
	for (int i=0;i<NDIM+1;++i)
	  out[i]=in[i];

	int n=NDIM+1;
	for (int i=0;i<NDIM;++i)
	  for (int j=i+1;j<NDIM+1;++j,++n)
	    out[n]=in[i]*in[j]*4.0;
      }
    };

    // Shape functions
    template <>
    struct ShapeFunctionsHelper<3,false,2>
    {
      template <class I, class O>
      static void evalShapeFunctions(const I in, O &out)
      {
	out[0]=in[0]*(in[0]*2.0-1.0);
	out[1]=in[1]*(in[1]*2.0-1.0);
	out[2]=in[2]*(in[2]*2.0-1.0);
	out[3]=in[3]*(in[3]*2.0-1.0);
	out[4]=in[0]*in[1]*4.0;
	out[5]=in[0]*in[2]*4.0;
	out[6]=in[0]*in[3]*4.0;
	out[7]=in[1]*in[2]*4.0;
	out[8]=in[1]*in[3]*4.0;
	out[9]=in[2]*in[3]*4.0;
      }

      template <class I, class OM>
      static void evalDiffShapeFunctions(const I in, OM &out)
      {
	out << 
	  in[0]*4-1, 0, 0, 0,
	  0, in[1]*4-1, 0, 0,
	  0, 0, in[2]*4-1, 0,
	  0, 0, 0, in[3]*4-1,
	  in[1]*4.0, in[0]*4.0, 0, 0,
	  in[2]*4.0, 0, in[0]*4.0, 0,
	  in[3]*4.0, 0, 0, in[0]*4.0,
	  0, in[2]*4.0, in[1]*4.0, 0,
	  0, in[3]*4.0, 0, in[1]*4.0,
	  0, 0, in[3]*4.0, in[2]*4.0;
      }
    };
    
    template <>
    struct ShapeFunctionsHelper<2,false,2>
    {
      template <class I, class O>
      static void evalShapeFunctions(const I in, O &out)
      {
	out[0]=in[0]*(in[0]*2.0-1.0);
	out[1]=in[1]*(in[1]*2.0-1.0);
	out[2]=in[2]*(in[2]*2.0-1.0);
	out[3]=in[0]*in[1]*4.0;
	out[4]=in[0]*in[2]*4.0;
	out[5]=in[1]*in[2]*4.0;
      }

      template <class I, class OM>
      static void evalDiffShapeFunctions(const I in, OM &out)
      {
	out << 
	  in[0]*4-1, 0, 0,
	  0, in[1]*4-1, 0,
	  0, 0, in[2]*4-1,
	  in[1]*4.0, in[0]*4.0, 0,
	  in[2]*4.0, 0, in[0]*4.0,
	  0, in[2]*4.0, in[1]*4.0;
      }
    };

    // Hierachical
    template <int NDIM_W>
    class TransformHelper<3,NDIM_W,true,2>
    {    
    public:
      static const int NDIM=3;

      template <class I, class CM, class OM>
      static void evalQuarterJacobian(const I in, const CM& cM, OM &out)
      {	
	out.row(0) << 0.25, 0.25, 0.25, 0.25;
	for (int dim=0;dim<NDIM_W;++dim)
	  out.row(dim+1) << 
	    cM(dim,0)*0.25+cM(dim,4)*in[1]+cM(dim,5)*in[2]+cM(dim,6)*in[3],
	    cM(dim,1)*0.25+cM(dim,4)*in[0]+cM(dim,7)*in[2]+cM(dim,8)*in[3],
	    cM(dim,2)*0.25+cM(dim,5)*in[0]+cM(dim,7)*in[1]+cM(dim,9)*in[3],
	    cM(dim,3)*0.25+cM(dim,6)*in[0]+cM(dim,8)*in[1]+cM(dim,9)*in[2];
      }     
    };

    // Standard
    template <int NDIM_W>
    class TransformHelper<3,NDIM_W,false,2>
    {    
    public:
      static const int NDIM=3;

      template <class I, class CM, class OM>
      static void evalQuarterJacobian(const I in, const CM& cM, OM &out)
      {	
	out.row(0) << 0.25, 0.25, 0.25, 0.25;
	for (int dim=0;dim<NDIM_W;++dim)
	  out.row(dim+1) << 
	    cM(dim,0)*(in[0]-0.25)+cM(dim,4)*in[1]+cM(dim,5)*in[2]+cM(dim,6)*in[3],
	    cM(dim,1)*(in[1]-0.25)+cM(dim,4)*in[0]+cM(dim,7)*in[2]+cM(dim,8)*in[3],
	    cM(dim,2)*(in[2]-0.25)+cM(dim,5)*in[0]+cM(dim,7)*in[1]+cM(dim,9)*in[3],
	    cM(dim,3)*(in[3]-0.25)+cM(dim,6)*in[0]+cM(dim,8)*in[1]+cM(dim,9)*in[2];
      }     
    
      // TODO: NEED TO IMPLEMENT A FAST VERSION OF THIS (AS IN 2D!)
      /*
      template <class CM, class O>
      static void evalJacobianDetAtVertices(const CM& cM, O out[4]);
      
      // slow version  ;)      
      template <class CM, class O>
      static void evalJacobianDetAtVertices(const CM& cM, O out[4])
      {
	for (int i=0;i< NDIM+1;++i)
	  {
	    O c[NDIM+1]={0};
	    c[i]=1;
	    out[i]=evalJacobianDet(c);
	  }
      }
      */
      
    };

    // Hierarchical
    template <int NDIM_W>
    class TransformHelper<2,NDIM_W,true,2>
    {
    public:
      static const int NDIM=2;
      
      template <class I, class CM, class OM>
      static void evalQuarterJacobian(const I in, const CM& coordsMat, OM &out)
      {	
	out.row(0) << 0.25, 0.25, 0.25;
	for (int dim=0;dim<NDIM_W;++dim)
	  out.row(dim+1) << 
	    coordsMat(dim,0)*0.25+coordsMat(dim,3)*in[1]+coordsMat(dim,4)*in[2],
	    coordsMat(dim,1)*0.25+coordsMat(dim,3)*in[0]+coordsMat(dim,5)*in[2],
	    coordsMat(dim,2)*0.25+coordsMat(dim,4)*in[0]+coordsMat(dim,5)*in[1];
      }
    };

    // Standard
    template <int NDIM_W>
    class TransformHelper<2,NDIM_W,false,2>
    {
    public:
      static const int NDIM=2;
      
      template <class I, class CM, class OM>
      static void evalQuarterJacobian(const I in, const CM& coordsMat, OM &out)
      {	
	out.row(0) << 0.25, 0.25, 0.25;
	for (int dim=0;dim<NDIM_W;++dim)
	  out.row(dim+1) <<
	    coordsMat(dim,0)*(in[0]-0.25)+coordsMat(dim,3)*in[1]+coordsMat(dim,4)*in[2],
	    coordsMat(dim,1)*(in[1]-0.25)+coordsMat(dim,3)*in[0]+coordsMat(dim,5)*in[2],
	    coordsMat(dim,2)*(in[2]-0.25)+coordsMat(dim,4)*in[0]+coordsMat(dim,5)*in[1];
      }
     
      // TODO: NEED TO IMPLEMENT A FAST VERSION OF THIS (AS IN 2D!)
      /*
      template <class CM, class O>
      static void evalJacobianDetAtVertices(const CM& cM, O out[4]);
      
      // slow version  ;)      
      template <class CM, class O>
      static void evalJacobianDetAtVertices(const CM& cM, O out[4])
      {
	for (int i=0;i< NDIM+1;++i)
	  {
	    O c[NDIM+1]={0};
	    c[i]=1;
	    out[i]=evalJacobianDet(c);
	  }
      }
      */
      /*
      template <class CM, class O>
      static void evalJacobianDetAtVertices(const CM& cM, O out[3])
      {
	O a=cM(0,3)-cM(0,1)/4;
	O b=cM(1,4)-cM(1,2)/4;
	O c=cM(1,3)-cM(1,1)/4;
	O d=cM(0,4)-cM(0,2)/4;
	out[0]=(a*b-c*d)+(cM(0,0)*(c-b) + cM(1,0)*(d-a))*0.75;

	// 0->1, 1->2, 2->0, 3->5, 4->3
	a=cM(0,5)-cM(0,2)/4;
	b=cM(1,3)-cM(1,0)/4;
	c=cM(1,5)-cM(1,2)/4;
	d=cM(0,3)-cM(0,0)/4;
	out[1]=(a*b-c*d)+(cM(0,1)*(c-b) + cM(1,1)*(d-a))*0.75;
	
	// 0->2, 1->0, 2->1, 3->4, 4->5
	a=cM(0,4)-cM(0,0)/4;
	b=cM(1,5)-cM(1,1)/4;
	c=cM(1,4)-cM(1,0)/4;
	d=cM(0,5)-cM(0,1)/4;
	out[2]=(a*b-c*d)+(cM(0,2)*(c-b) + cM(1,2)*(d-a))*0.75;
      }
      */

      /*
      template <typename I, class CM, class OT=double>
      static OT evalJacobianDet(const I in, const CM& coordsMat)
      {

	static const OT factor = 
	      static_cast<OT>(hlp::IntPower<4,NDIM+1>::value)/
	      static_cast<OT>(hlp::FactorialT<NDIM>::value);
	Eigen::Matrix<OT,NDIM_W+1,NDIM+1> J;
	evalQuarterJacobian(in,coordsMat,J);
	if (NDIM==NDIM_W)
	  {
	    // This is faster ...
	    return J.determinant()*factor;
	  }
	else
	  {
	    OT det2=(J.transpose()*J).determinant();	   

	    if (det2<0)
	      return -sqrt(-det2)*factor;
	    else
	      return sqrt(det2)*factor;	   
	  }
      }
      */
    };
 
  } 
}

#include "../../internal/namespace.footer"
#endif

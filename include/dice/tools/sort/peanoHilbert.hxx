#ifndef __PEANO_HILBERT_HXX__
#define __PEANO_HILBERT_HXX__

/**
 * @file 
 * @brief  A Peano-Hilbert N-D coordinates to 1-D distance mapping class
 * @author Thierry Sousbie (from J.K. Lawder, 2000)
 */

#include "../../internal/namespace.header"
/** \addtogroup TOOLS 
 *   \{
 */

/**
 * \class PeanoHilbertT
 * \brief  A static class used to convert coordinates to distance over a peano hilbert curve
 * or reciprocally. It work in any number of dimensions and to any order ...
 * The core of the code is from J.K. Lawder, 2000.
 * \tparam ND The number of dimensions
 * \tparam maxOrder the maximum order of the PH curve: it may have a resolution of (1<<maxOrder) at most along each dimension.
 */

template <int ND, int maxOrder=32>
class PeanoHilbertT
{
public:
  
  typedef typename hlp::MinimalIntegerType<maxOrder,false>::Type Index; 
  
  /** 
   *  \struct HCode 
   *  A structure used to store distance along the PH curve as a set of 
   *  NDIM unsigned integers. The lower bits (i.e least important) are stored in hcode[0],
   *  while the higher ones are in hcode[NDIM-1].
   */
  struct HCode {
    static const int NDIM = ND;
    Index hcode[NDIM];

    /** \brief Converts the HCode to a floating point distance (precision may be lost)
     *  \tparam FT the floating point type to convert to
     */
    template <typename FT>
    FT to() const
    {
      FT result=static_cast<FT>(hcode[0]);
      FT fac=std::numeric_limits<Index>::max();
      for (int i=1;i<NDIM;++i)
	{
	  result += static_cast<FT>(hcode[i]) * fac;
	  fac *= std::numeric_limits<Index>::max();
	}
      return result/fac;
    }
    
    bool operator== (const HCode &rhs) const
     {      
      for (int i=0;i<NDIM;++i) 
	if (hcode[i]!=rhs.hcode[i]) return false;
      
      return true;
    }

    bool operator< (const HCode &rhs) const
    {
      for (int i=NDIM-1;i>0;--i) 
	{
	  if (hcode[i]<rhs.hcode[i]) return true;      
	  if (hcode[i]>rhs.hcode[i]) return false;
	}
      if (hcode[0]<rhs.hcode[0]) return true;
      return false;      
    }
  };

  static const int NDIM = ND;
  
  /** \brief Convert coordinates to distance along the PH curve
   *  \param[in]  coords coordinates
   *  \param[out] length distance along the curve
   *  \param      x0     Coordinates of the lower left corner of the bounding box
   *  \param      delta_inv inverse of the size of the bounding box along each dimension
   *  \param      order  order of the PH curve (i.e. its precision is 1<<order along each dimension). order should be lower than or equal to maxOrder.
   */
  template <typename CT, typename DT>
  static void coordsToLength(CT coords[NDIM], HCode &length, 
			     DT x0[NDIM], DT delta_inv[NDIM],
			     int order=maxOrder)
  {
    Point p;
    for (int i=0;i<NDIM;++i)
      {
	p.hcode[i] = static_cast<Index>
	  (((coords[i]-x0[i])*delta_inv[i]) * std::numeric_limits<Index>::max());
      }
    length = H_encode(p,order);
  }

  /** \brief Convert coordinates to distance along the PH curve
   *  \param[in]  coords coordinates
   *  \param[out] length distance along the curve
   *  \param      order  order of the PH curve (i.e. its precision is 1<<order along each dimension). order should be lower than or equal to maxOrder.
   */
  template <typename CT>
  static void coordsToLength(CT coords[NDIM], HCode &length, int order=maxOrder)
  {
    Point p;
    for (int i=0;i<NDIM;++i)
      {
	p.hcode[i] = static_cast<HCT>
	  (coords[i] * std::numeric_limits<Index>::max());
      }
    length = H_encode(p,order);
  }

  /** \brief Convert a length along the PH curve to coordinates.   
   *  \param[in] length distance along the curve
   *  \param[out]  coords coordinates  
   *  \param      order  order of the PH curve (i.e. its precision is 1<<order along each dimension). order should be lower than or equal to maxOrder.
   */
  template <typename CT>
  static void lengthToCoords(const HCode &length, CT coords[NDIM], int order=maxOrder)
  {
    HCode point = H_decode(length,order);
    for (int i=0;i<NDIM;++i)
      {
	coords[i]=static_cast<CT>(point.hcode[i]) / std::numeric_limits<Index>::max();
      }
  }

  /** \brief Convert a length along the PH curve to coordinates.   
   *  \param[in] length distance along the curve
   *  \param[out]  coords coordinates
   *  \param      x0     Coordinates of the lower left corner of the bounding box
   *  \param      delta  The size of the bounding box along each dimension
   *  \param      order  order of the PH curve (i.e. its precision is 1<<order along each dimension). order should be lower than or equal to maxOrder.
   */
  template <typename CT, typename DT>
  static void lengthToCoords(const HCode &length, CT coords[NDIM], 
			    DT x0[NDIM], DT delta[NDIM],
			    int order=maxOrder)
  {
    HCode point = H_decode(length,order);
    for (int i=0;i<NDIM;++i)
      {
	coords[i]=static_cast<CT>(point.hcode[i])/std::numeric_limits<Index>::max();
	coords[i]= x0[i] + (coords[i] * delta[i]);
      }
  }
  
private: 
typedef Index HCT;
typedef HCode Point;
  static HCT g_mask(int i) {return HCT(1)<<(NDIM-1-i);}

  // From here, this is an adapted version of J.K. Lawder, 2000
  /*
    This code assumes the following:
    The macro ORDER corresponds to the order of curve and is 32,
    thus coordinates of points are 32 bit values.
    A HCT should be a 32 bit unsigned integer.
    The macro NDIM corresponds to the number of dimensions in a
    space.
    The derived-key of a Point is stored in an HCode which is an
    array of HCT. The bottom bit of a derived-key is held in the
    bottom bit of the hcode[0] element of an HCode and the top bit
    of a derived-key is held in the top bit of the hcode[NDIM-1]
    element of and HCode.
    g_mask is a global array of masks which helps simplify some
    calculations - it has NDIM elements. In each element, only one
    bit is zeo valued - the top bit in element no. 0 and the bottom
    bit in element no. (NDIM - 1). eg.
    #if NDIM == 5 const HCT g_mask[] = {16, 8, 4, 2, 1}; #endif
    #if NDIM == 6 const HCT g_mask[] = {32, 16, 8, 4, 2, 1}; #endif
    etc...
  */

  
  /*===========================================================*/
  /* calc_P */
  /*===========================================================*/
  static HCT calc_P (int i, HCode &H, int order=maxOrder)
  {
    int element;
    HCT P, temp1, temp2;
    element = i / order;
    P = H.hcode[element];
    if (i % order > order - NDIM)
      {
	temp1 = H.hcode[element + 1];
	P >>= i % order;
	temp1 <<= order - i % order;
	P |= temp1;
      }
    else
      P >>= i % order; /* P is a NDIM bit hcode */
    /* the & masks out spurious highbit values */

    if (NDIM<order)
      P &= (1 << NDIM) -1;

    return P;
  }

  /*===========================================================*/
  /* calc_P2 */
  /*===========================================================*/
  static HCT calc_P2(HCT S)
  {
    int i;
    HCT P;
    P = S & g_mask(0);
    for (i = 1; i < NDIM; i++)
      if( (S & g_mask(i)) ^ ((P >> 1) & g_mask(i)))
	P |= g_mask(i);
    return P;
  }

  /*===========================================================*/
  /* calc_J */
  /*===========================================================*/
  static HCT calc_J (HCT P)
  {
    int i;
    HCT J;
    J = NDIM;
    for (i = 1; i < NDIM; i++)
      if ((P >> i & 1) == (P & 1))
	continue;
      else
	break;
    if (i != NDIM)
      J -= i;
    return J;
  }

  /*===========================================================*/
  /* calc_T */
  /*===========================================================*/
  static HCT calc_T (HCT P)
  {
    if (P < 3)
      return 0;
    if (P % 2)
      return (P - 1) ^ (P - 1) / 2;
    return (P - 2) ^ (P - 2) / 2;
  }

  /*===========================================================*/
  /* calc_tS_tT */
  /*===========================================================*/
  static HCT calc_tS_tT(HCT xJ, HCT val)
  {
    HCT retval, temp1, temp2;
    retval = val;
    if (xJ % NDIM != 0)
      {
	temp1 = val >> xJ % NDIM;
	temp2 = val << (NDIM - xJ % NDIM);
	retval = temp1 | temp2;
	retval &= (HCT(1) << NDIM) - 1;
      }
    return retval;
  }

  /*===========================================================*/
  /* H_decode */
  /*===========================================================*/
  /* For mapping from one dimension to NDIM dimensions */
  static Point H_decode (HCode &H, int order=maxOrder)
  {
    HCT mask = HCT(1) << order - 1,
      A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    Point pt = {0};
    int i = order * NDIM - NDIM, j;
    P = calc_P(i, H);
    J = calc_J(P);
    xJ = J - 1;
    A = S = tS = P ^ P / 2;
    T = calc_T(P);
    tT = T;
    /*--- distrib bits to coords ---*/
    for (j = NDIM - 1; P > 0; P >>=1, j--)
      if (P & 1)
	pt.hcode[j] |= mask;
    for (i -= NDIM, mask >>= 1; i >=0; i -= NDIM, mask >>= 1)
      {
	P = calc_P(i, H);
	S = P ^ P / 2;
	tS = calc_tS_tT(xJ, S);
	W ^= tT;
	A = W ^ tS;
	/*--- distrib bits to coords ---*/
	for (j = NDIM - 1; A > 0; A >>=1, j--)
	  if (A & 1)
	    pt.hcode[j] |= mask;
	if (i > 0)
	  {
	    T = calc_T(P);
	    tT = calc_tS_tT(xJ, T);
	    J = calc_J(P);
	    xJ += J - 1;
	  }
      }
    return pt;
  }

  /*===========================================================*/
  /* H_encode */
  /*===========================================================*/
  /* For mapping from NDIM dimensions to one dimension */
  static HCode H_encode(Point &pt, int order=maxOrder)
  {
    HCT mask = (HCT)1 << (order - 1), element,
      A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    HCode h = {{0}};
    int i = order * NDIM - NDIM, j;
    for (j = A = 0; j < NDIM; j++)
      if (pt.hcode[j] & mask)
	A |= g_mask(j);
    S = tS = A;
    P = calc_P2(S);
    /* add in NDIM bits to hcode */
    element = i / order;
    if (i % order > order - NDIM)
      {
	h.hcode[element] |= P << (i % order);
	h.hcode[element + 1] |= P >> (order - (i % order));
      }
    else
      h.hcode[element] |= P << (i - element * order);
    J = calc_J(P);
    xJ = J - 1;
    T = calc_T(P);
    tT = T;
    for (i -= NDIM, mask >>= 1; i >=0; i -= NDIM, mask >>= 1)
      {
	for (j = A = 0; j < NDIM; j++)
	  if (pt.hcode[j] & mask)
	    A |= g_mask(j);
	W ^= tT;
	tS = A ^ W;
	S = calc_tS_tT(xJ, tS);
	P = calc_P2(S);
	/* add in NDIM bits to hcode */
	element = i / order;
	if (i % order > order - NDIM)
	  {
	    h.hcode[element] |= P << i % order;
	    h.hcode[element + 1] |= P >> (order - i % order);
	  }
	else
	  h.hcode[element] |= P << (i - element * order);
	if (i > 0)
	  {
	    T = calc_T(P);
	    tT = calc_tS_tT(xJ, T);
	    J = calc_J(P);
	    xJ += J - 1;
	  }
      }
    return h;
  }

  
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

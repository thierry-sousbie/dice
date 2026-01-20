#ifndef __LOCAL_AMR_GRID_DEFAULT_OVERLAP_VISITORS_HXX__
#define __LOCAL_AMR_GRID_DEFAULT_OVERLAP_VISITORS_HXX__

#include <stdio.h>

#include "../../internal/namespace.header"

namespace internal
{

  namespace localAmrGridVisitor
  {
    template <class AMR, class V>
    class BoxOverlapHandlerT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM = AMR::NDIM;

    public:
      BoxOverlapHandlerT(V &visitor_) : visitor(visitor_)
      {
      }

      void initialize_overlap_handler(const ICoord iBox_[2][NDIM],
                                      ICoord voxIDim_[NDIM],
                                      ICoord iLen_)
      {
        iBox = &iBox_[0][0];
        voxIDim = voxIDim_;
        iLen = iLen_;
        /**/
        iLenRef = 0; // AMR::BBOX_ILEN;
        /**/
      }

      static void initialize(Voxel *root) {}

      bool visit(Voxel *voxel)
      {
        return visitor.visit(voxel);
      }

      bool visit(Voxel *voxel, int i)
      {
        /**/
        if (iLenRef > iLen)
        {
          iLen >>= 1;
          return true;
        }
        else
          iLenRef = 0;
        /**/

        bool inBound = true;
        // update voxIDim and check bounding box
        for (ICoord j = 0; j < NDIM; ++j)
        {
          voxIDim[j] += ((i & (1L << j)) >> j) * iLen;

          if ((voxIDim[j] > iBox[j + NDIM]) ||
              (voxIDim[j] + iLen <= iBox[j]))
            inBound = false; // out of bound!
          /**/
          else if ((voxIDim[j] <= iBox[j]) &&
                   (voxIDim[j] + iLen > iBox[j + NDIM]))
            iLenRef = iLen;
          /**/
        }
        iLen >>= 1;
        return inBound;
      }

      void visited(Voxel *voxel)
      {
        return visitor.visited(voxel);
      }

      void visited(Voxel *voxel, int i)
      {
        iLen <<= 1;
        /**/
        if (iLen >= iLenRef)
          /**/
          for (ICoord j = 0; j < NDIM; ++j)
            voxIDim[j] -= ((i & (1L << j)) >> j) * iLen;
      }

    protected:
      const ICoord *iBox;
      ICoord *voxIDim;
      ICoord iLen;
      ICoord iLenRef;
      V &visitor;
    };
  }

  // These visitors are supposed to be used as parameters to
  // localAmrGridVisitor::BoxOverlapHandlerT
  namespace localAmrGridOverlapVisitor
  {

    // To be used through BoxOverlapHandlerT
    template <class AMR, class OutputIterator>
    class GetLeavesBoxOverlapT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM = AMR::NDIM;

    public:
      GetLeavesBoxOverlapT(OutputIterator out_, bool nonNullOnly_ = false) : out(out_),
                                                                             nonNullOnly(nonNullOnly_)
      {
      }

      bool visit(Voxel *voxel)
      {
        if (voxel->isLeaf())
        {
          if ((!nonNullOnly) || (voxel->data != 0))
          {
            (*out) = voxel;
            ++out;
          }
          return false;
        }
        return true;
      }

      static void visited(Voxel *voxel)
      {
      }

    private:
      OutputIterator out;
      const bool nonNullOnly;
    };

    // To be used through BoxOverlapHandlerT
    template <class AMR>
    class RefineT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM = AMR::NDIM;

    public:
      RefineT(AMR *amr_, int level_) : amr(amr_),
                                       level(level_)
      {
      }

      bool visit(Voxel *voxel)
      {
        if (voxel->getLevel() < level)
        {
          if (voxel->isLeaf())
            voxel->refine(amr);
          return true;
        }
        else
          return false;
        /*
        if (voxel->isLeaf())
          {
            if (voxel->getLevel()<level)
              voxel->refine(amr);
            else return false;
          }
        return true;
        */
      }

      static void visited(Voxel *voxel)
      {
      }

    private:
      AMR *amr;
      int level;
    };

    // To be used through BoxOverlapHandlerT
    template <class AMR>
    class ThreadedRefineT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM = AMR::NDIM;

    public:
      ThreadedRefineT(AMR *amr_, int level_) : amr(amr_),
                                               level(level_)
      {
      }

      bool visit(Voxel *voxel)
      {
        if (voxel->getLevel() < level)
        {
          if (voxel->isLeaf())
            voxel->refineCritical(amr);
          return true;
        }
        return false;
      }

      static void visited(Voxel *voxel)
      {
      }

    private:
      AMR *amr;
      int level;
    };

    // To be used through BoxOverlapHandlerT
    template <class AMR>
    class SetRootLevelT
    {
      typedef typename AMR::Voxel Voxel;
      typedef typename AMR::ICoord ICoord;

      static const long NDIM = AMR::NDIM;

    public:
      SetRootLevelT(AMR *amr_, int level_, const Voxel *rootVoxels,
                    std::vector<unsigned char> &rootVoxelsLevel) : amr(amr_),
                                                                   level(level_),
                                                                   ref(rootVoxels),
                                                                   rootLevel(rootVoxelsLevel)
      {
      }

      bool visit(const Voxel *voxel)
      {
        long i = std::distance(ref, voxel);
        if (rootLevel[i] < level)
          rootLevel[i] = level;

        return false;
      }

      static void visited(Voxel *voxel)
      {
      }

      void setLevel(int lvl)
      {
        level = lvl;
      }

    private:
      AMR *amr;
      int level;
      const Voxel *ref;
      std::vector<unsigned char> &rootLevel;
    };

  }
}

#include "../../internal/namespace.footer"

#endif

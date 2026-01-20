#ifndef __REGULAR_GRID_SLICER_FFTW_HXX__
#define __REGULAR_GRID_SLICER_FFTW_HXX__

#include <fftw3.h>

#include "../fftw3-mpi_dummy.h"

#include "../../dice_globals.hxx"
#include "../../grid/internal/regularGridSlicerBase.hxx"

#include "../../tools/MPI/mpiCommunication.hxx"

#include "../../grid/valLocationType.hxx"

#include "../../internal/namespace.header"

namespace internal
{

  template <class G>
  class RegularGridSlicerFFTWT : public RegularGridSlicerBaseT<G>
  {
  public:
    typedef RegularGridSlicerBaseT<G> Base;
    typedef RegularGridSlicerSimpleT<G> MyType;

    static const int NDIM = G::NDIM;

    typedef G Grid;
    typedef typename Grid::Params Params;
    typedef typename Grid::Scale Scale;
    // typedef typename Grid::ValLocationType ValLocationType;
    // typedef typename Grid::ValLocationTypeV ValLocationTypeV;
    typedef typename Scale::ScaleTypeV ScaleTypeV;
    typedef typename Scale::ScaleType ScaleType;
    typedef typename Grid::GridNav GridNav;
    typedef typename Grid::Direction Direction;

    RegularGridSlicerFFTWT(const Params &gp_,
                           MpiCommunication *com_,
                           int nThreads_,
                           bool periodic_,
                           long allocFactor_ = 1) : allocFactor(allocFactor_)
    {
      initialize(gp_, com_, nThreads_, periodic_);
    }

    RegularGridSlicerFFTWT()
    {
    }

    void slice(int index)
    {
      neighbors.clear();
      long nSlices = nChunks;
      globalGridParams = gp;
      for (int i = 0; i < NDIM; ++i)
        sliceCount[i] = 1;
      sliceCount[NDIM - 1] = nSlices;
      myIndex = index;
      for (int i = 0; i < NDIM; ++i)
        slicePos[i] = 0;
      slicePos[NDIM - 1] = index;

      gridParams = this->divide(gp,
                                local_0_start[index],
                                local_0_start[index] + local_n0[index],
                                0, 0, periodic, NDIM - 1);

      if (slicePos[NDIM - 1] + 1 >= sliceCount[NDIM - 1])
      {
        if (periodic)
          neighbors.push_back(Neighbor(NDIM - 1, 1, 0));
      }
      else
        neighbors.push_back(Neighbor(NDIM - 1, 1, myIndex + 1));

      if (slicePos[NDIM - 1] <= 0)
      {
        if (periodic)
          neighbors.push_back(Neighbor(NDIM - 1, -1, sliceCount[0] - 1));
      }
      else
        neighbors.push_back(Neighbor(NDIM - 1, -1, myIndex - 1));

      gridParams.haveParentGrid = true;
      gridParams.parentNDim = NDIM;
      gridParams.minElementsCount = 0;

      for (int i = 0; i < NDIM; ++i)
        gridParams.position[i] = 0;
      gridParams.position[NDIM - 1] = local_0_start[index];

      for (int i = 0; i < NDIM; ++i)
      {
        gridParams.parentResolution[i] = gp.resolution[i];
        gridParams.parentX0[i] = gp.x0[i];
        gridParams.parentDelta[i] = gp.delta[i];
      }
      /*
      printf("nalloc[%d] = %ld * %ld (%d %d %d)\n",
       index,nAlloc[index],allocFactor,
       gridParams.resolution[0],
       gridParams.resolution[1],
       gridParams.resolution[2]);
      */
      gridParams.minElementsCount = nAlloc[index] * allocFactor;
    }

    const Params &getLocalGridParams() const
    {
      return gridParams;
    }

    const Params &getGlobalGridParams() const
    {
      return globalGridParams;
    }

    long neighborsCount() const
    {
      return neighbors.size();
    }

    template <class T, class T2>
    void getNeighborInfo(int which, T &neiDim, T &neiDir, T2 &neiIndex) const
    {
      neiDim = neighbors[which].dim;
      neiDir = neighbors[which].dir;
      neiIndex = neighbors[which].index;
    }

    long getSliceCount(int dim) const
    {
      return sliceCount[dim];
    }

    long getSlicePos(int dim) const
    {
      return slicePos[dim];
    }

    MpiCommunication *getMpiCom() const
    {
      return mpiCom;
    }

    long getNChunks()
    {
      return nChunks;
    }

    long getNThreads()
    {
      return nChunks;
    }

    bool needFacade() const
    {
      if (periodic)
        return false;
      return true;
    }

    bool toFacade(const Params &gp_)
    {
      if (!needFacade())
        return false;
      initialize(gp_, mpiCom, numThreads, false, true);

      return true;
    }

    template <class MPC>
    static void rearrangeFromFacadeToGrid(Grid *grid, Grid *gridFacade, MPC *com)
    {
      typedef typename Grid::LocalGrid::inOrderIterator IOT;
      typedef typename Grid::Data Data;

      if (grid == gridFacade)
        return;

      auto *localFacade = gridFacade->getLocalGrid();
      auto *localGrid = grid->getLocalGrid();

      if (com->size() == 1)
      {
        // No need for buffer here, we can do that inplace
        int iMin[NDIM] = {0};
        int iMax[NDIM];

        for (int i = 0; i < NDIM; ++i)
          iMax[i] = localFacade->getResolution(i);

        auto *src = localFacade->getDataPtr();
        auto *dest = localGrid->getDataPtr() + localGrid->getNValues() / 2;
        std::copy(src, src + localFacade->getNValues(), dest);

        // so that lower part becomes 0 where appropriate
        std::fill_n(localGrid->getDataPtr(), localGrid->getNValues() / 2, 0);

        auto git = IOT(localGrid->subbox_begin(iMin, iMax));
        std::copy(dest, dest + localFacade->getNValues(), git);

        // so that upper part becomes 0 where appropriate
        std::fill_n(dest, localGrid->getNValues() / 2, 0);

        return;
      }

      int smallFullRes[NDIM];
      int fullRes[NDIM];
      int smallRes[NDIM];
      int res[NDIM];
      int lastSmallRes[NDIM];
      int lastRes[NDIM];
      int curRes[NDIM];
      int curSmallRes[NDIM];
      // long smallCount=1;
      // long count=1;

      for (int i = 0; i < NDIM; ++i)
      {
        smallFullRes[i] = gridFacade->getResolution(i);
        fullRes[i] = grid->getResolution(i);
        smallRes[i] = localFacade->getResolution(i);
        res[i] = localGrid->getResolution(i);
        lastSmallRes[i] = localFacade->getResolution(i);
        lastRes[i] = localGrid->getResolution(i);
        // smallCount*=smallRes[i];
        // count*=res[i];
      }

      if (com->rank() == com->size() - 1)
      {
        res[NDIM - 1] = (fullRes[NDIM - 1] - lastRes[NDIM - 1]) / (com->size() - 1);
        smallRes[NDIM - 1] = (smallFullRes[NDIM - 1] - lastSmallRes[NDIM - 1]) / (com->size() - 1);
        std::copy(lastRes, lastRes + NDIM, curRes);
        std::copy(lastSmallRes, lastSmallRes + NDIM, curSmallRes);
      }
      else
      {
        lastRes[NDIM - 1] = fullRes[NDIM - 1] - (com->size() - 1) * res[NDIM - 1];
        lastSmallRes[NDIM - 1] = smallFullRes[NDIM - 1] - (com->size() - 1) * smallRes[NDIM - 1];
        std::copy(res, res + NDIM, curRes);
        std::copy(smallRes, smallRes + NDIM, curSmallRes);
      }

      int sendCount[2 + 1] = {0};
      int receiveCount[3 + 1] = {0};
      int dest[2];
      int src[3];
      int sendOffset[2];
      int receiveOffset[3];

      dest[0] = (com->rank() * smallRes[NDIM - 1]) / res[NDIM - 1];
      sendOffset[0] = 0;
      sendCount[0] = res[NDIM - 1] * (dest[0] + 1) - smallRes[NDIM - 1] * com->rank();
      if (sendCount[0] < curSmallRes[NDIM - 1])
      {
        dest[1] = dest[0] + 1;
        sendOffset[1] = sendCount[0];
        sendCount[1] = curSmallRes[NDIM - 1] - sendCount[0];
      }
      else
      {
        sendCount[0] = curSmallRes[NDIM - 1];
      }

      int idx = -1;
      int oldSrc = (com->rank() * res[NDIM - 1]) / smallRes[NDIM - 1];
      int oldOffset = 0;
      int oldCount = 0;

      if (oldSrc < com->size())
        do
        {
          idx++;
          src[idx] = oldSrc++;

          receiveOffset[idx] = oldOffset + oldCount;
          oldOffset = receiveOffset[idx];

          int offset = (com->rank() * res[NDIM - 1]) + oldOffset - src[idx] * smallRes[NDIM - 1];
          int nMaxRcv = (src[idx] == com->size() - 1) ? lastSmallRes[NDIM - 1] : smallRes[NDIM - 1];
          nMaxRcv -= offset;
          if (nMaxRcv < 0)
            break;

          receiveCount[idx] = std::min(curRes[NDIM - 1] - receiveOffset[idx], nMaxRcv);
          oldCount = receiveCount[idx];
        } while ((receiveOffset[idx] + receiveCount[idx] < curRes[NDIM - 1]) &&
                 (src[idx] != com->size() - 1));

      /*
      // Print debug info
      glb::console->print<LOG_STD_ALL>
  ("process %d sends to [%d: %d @%d] [%d: %d @%d]\n",
   com->rank(),
   dest[0],sendCount[0],sendOffset[0],
   dest[1],sendCount[1],sendOffset[1]);

      glb::console->print<LOG_STD_ALL>
  ("process %d [%d] receive from [%d: %d @%d] [%d: %d @%d] [%d: %d @%d]\n",
   com->rank(),curRes[NDIM-1],
   src[0],receiveCount[0],receiveOffset[0],
   src[1],receiveCount[1],receiveOffset[1],
   src[2],receiveCount[2],receiveOffset[2]);

      com->barrier();
      */

      MPI_Request req[5];
      int nReq = 0;
      // long sliceCount=smallCount/smallRes[NDIM-1];
      long sliceCount = 1;
      for (int i = 0; i < NDIM - 1; ++i)
        sliceCount *= smallRes[i];

      for (int i = 0; sendCount[i]; ++i)
      {
        com->Isend(localFacade->getDataPtr() + sendOffset[i] * sliceCount,
                   sendCount[i] * sliceCount, dest[i], &req[nReq++]);
      }

      // Where we will receive the data
      int count = 1;
      for (int i = 0; i < NDIM; ++i)
        count *= curRes[i];
      std::vector<Data> buffer(count >> (NDIM - 1));

      for (int i = 0; receiveCount[i]; ++i)
      {
        com->Irecv(&buffer[receiveOffset[i] * sliceCount],
                   receiveCount[i] * sliceCount, src[i], &req[nReq++]);
      }

      com->Waitall(nReq, req);
      localGrid->erase();

      int iMin[NDIM] = {0};
      int iMax[NDIM];
      for (int i = 0; i < NDIM - 1; ++i)
        iMax[i] = iMin[i] + curSmallRes[i];
      iMax[NDIM - 1] = iMin[NDIM - 1] + curRes[NDIM - 1];

      auto git = IOT(localGrid->subbox_begin(iMin, iMax));
      std::copy(buffer.begin(), buffer.end(), git);

      // REMOVE
      /*
      com->barrier();
      grid->toVtk("Grid","Grid_%d");
      gridFacade->toVtk("GridFacade","GridFacade_%d");
      com->barrier();
      exit(-1);
      */
      /*

      int posLow[NDIM];
      int posHigh[NDIM];
      int resolution[NDIM];
      int subResolution[NDIM];
      int nSendRecv=1;
      for (int i=0;i<NDIM;++i)
  {
    resolution[i]=localGrid->getResolution(i);
    posHigh[i]=posLow[i]=iMin[i];
    subResolution[i]=iMax[i]-iMin[i];
    nSendRecv *= subResolution[i];
  }
      posHigh[NDIM-1]+=subResolution[NDIM-1];

      // Create data types
      MPI_Datatype lowerPartType;
      MPI_Datatype upperPartType;
      glb::console->print<LOG_STD_ALL>
  ("sub[%d] : PL=[%d %d] PH=[%d %d] SR=[%d %d] R=[%d %d]\n",com->rank(),
   posLow[0],posLow[1],posHigh[0],posHigh[1],
   subResolution[0],subResolution[1],
   resolution[0],resolution[1]);

      MPI_Type_create_subarray(NDIM,resolution,subResolution,posLow,MPI_ORDER_FORTRAN,
             MPI_Type<Data>::get(),&lowerPartType);
      MPI_Type_commit(&lowerPartType);

      MPI_Type_create_subarray(NDIM,resolution,subResolution,posHigh,MPI_ORDER_FORTRAN,
             MPI_Type<Data>::get(),&upperPartType);
      MPI_Type_commit(&upperPartType);

      //remove
      com->barrier();

      MPI_Request sendReq;
      bool sending=false;
      // odd ranks transfer their lower part to destination upper part
      if (com->rank()&1)
  {
    glb::console->print<LOG_STD_ALL>
      ("%d Isending to %d\n",com->rank(),com->rank()/2);
    com->IsendMpiType
      (localGrid->getDataPtr(),1,lowerPartType,com->rank()/2,&sendReq,0);
    sending=true;
  }

      if (com->rank()< com->size()/2)
  {
    // No possible conflict here
    glb::console->print<LOG_STD_ALL>
      ("%d receiving from %d\n",com->rank(),com->rank()*2+1);
    com->RecvMpiType(localGrid->getDataPtr(),1,upperPartType,com->rank()*2+1,0);
  }
      if (sending) com->Wait(&sendReq);

      sending=false;
      // even ranks transfer their lower part to destination lower part
      if ((!(com->rank()&1))&&(com->rank()>0))
  {
    glb::console->print<LOG_STD_ALL>("%d Isending to %d\n",com->rank(),com->rank()/2);
    sending=true;
    com->IsendMpiType
      (localGrid->getDataPtr(),1,lowerPartType,com->rank()/2,&sendReq,1);
  }

      if ((com->rank()<(com->size()+1)/2)&&(com->rank()>0))
  {

    if (sending)
      {
         glb::console->print<LOG_STD_ALL>
     ("%d Receiving(+send) from %d\n",com->rank(),com->rank()*2);

        // We will erase data we are sending if we are not carefull
        std::vector<Data> buffer(nSendRecv);
        com->Recv(&buffer[0],buffer.size(),com->rank()*2,1);
        com->Wait(&sendReq);
        auto git = IOT(localGrid->subbox_begin(iMin,iMax));
        std::copy(buffer.begin(),buffer.end(),git);
      }
    else
      {
        glb::console->print<LOG_STD_ALL>
     ("%d Receiving from %d\n",com->rank(),com->rank()*2);
        com->RecvMpiType(localGrid->getDataPtr(),nSendRecv,
             lowerPartType,com->rank()*2,1);
      }


    // std::vector<Data> buffer(nSendRecv);
    // com->Recv(&buffer[0],buffer.size(),com->rank()*2,1);
    // if (sending) com->Wait(&sendReq);
    // auto git = IOT(localGrid->subbox_begin(iMin,iMax));
          // std::copy(buffer.begin(),buffer.end(),git);

  }
      else if (sending) com->Wait(&sendReq);
      */
      /*
      glb::console->print<LOG_STD_ALL>("%d passed",com->rank());
      // erase where needed
      if (com->rank()>=(com->size()+1)/2)
  grid->erase();

      MPI_Type_free(&lowerPartType);
      MPI_Type_free(&upperPartType);

      grid->toVtk("Grid","Grid_%d");
      gridFacade->toVtk("GridFacade","GridFacade_%d");
      com->barrier();
      exit(-1);
      */
    }

    template <class MPC>
    static void rearrangeFromGridToFacade(Grid *grid, Grid *gridFacade, MPC *com)
    {
      typedef typename Grid::LocalGrid::inOrderIterator IOT;
      typedef typename Grid::Data Data;

      if (grid == gridFacade)
        return;

      auto *localFacade = gridFacade->getLocalGrid();
      auto *localGrid = grid->getLocalGrid();

      if (com->size() == 1)
      {
        // No need for buffer here, we can do that inplace
        int iMin[NDIM] = {0};
        int iMax[NDIM];

        for (int i = 0; i < NDIM; ++i)
          iMax[i] = gridFacade->getResolution(i);

        auto git = IOT(localGrid->subbox_begin(iMin, iMax));
        auto git_end = IOT(localGrid->subbox_end(iMin, iMax));
        auto fit = IOT(localFacade->begin());

        std::copy(git, git_end, fit);
        return;
      }

      int smallFullRes[NDIM];
      int fullRes[NDIM];
      int smallRes[NDIM];
      int res[NDIM];
      int lastSmallRes[NDIM];
      int lastRes[NDIM];
      int curRes[NDIM];
      int curSmallRes[NDIM];
      // long smallCount=1;
      // long count=1;

      for (int i = 0; i < NDIM; ++i)
      {
        smallFullRes[i] = gridFacade->getResolution(i);
        fullRes[i] = grid->getResolution(i);
        smallRes[i] = localFacade->getResolution(i);
        res[i] = localGrid->getResolution(i);
        lastSmallRes[i] = localFacade->getResolution(i);
        lastRes[i] = localGrid->getResolution(i);
      }

      if (com->rank() == com->size() - 1)
      {
        res[NDIM - 1] = (fullRes[NDIM - 1] - lastRes[NDIM - 1]) / (com->size() - 1);
        smallRes[NDIM - 1] = (smallFullRes[NDIM - 1] - lastSmallRes[NDIM - 1]) / (com->size() - 1);
        std::copy(lastRes, lastRes + NDIM, curRes);
        std::copy(lastSmallRes, lastSmallRes + NDIM, curSmallRes);
      }
      else
      {
        lastRes[NDIM - 1] = fullRes[NDIM - 1] - (com->size() - 1) * res[NDIM - 1];
        lastSmallRes[NDIM - 1] = smallFullRes[NDIM - 1] - (com->size() - 1) * smallRes[NDIM - 1];
        std::copy(res, res + NDIM, curRes);
        std::copy(smallRes, smallRes + NDIM, curSmallRes);
      }
      /*
      long smallCount=1;
      long count=1;
      for (int i=0;i<NDIM;++i)
  {
    smallCount*=smallRes[i];
    count*=res[i];
  }
      */
      int sendCount[3 + 1] = {0};
      int receiveCount[2 + 1] = {0};
      int dest[3];
      int src[2];
      int sendOffset[3];
      int receiveOffset[2];

      src[0] = (com->rank() * smallRes[NDIM - 1]) / res[NDIM - 1];
      receiveOffset[0] = 0;
      receiveCount[0] = res[NDIM - 1] * (src[0] + 1) - smallRes[NDIM - 1] * com->rank();
      if (receiveCount[0] < curSmallRes[NDIM - 1])
      {
        src[1] = src[0] + 1;
        receiveOffset[1] = receiveCount[0];
        receiveCount[1] = curSmallRes[NDIM - 1] - receiveCount[0];
      }
      else
      {
        receiveCount[0] = curSmallRes[NDIM - 1];
      }

      int idx = -1;
      int oldDest = (com->rank() * res[NDIM - 1]) / smallRes[NDIM - 1];
      int oldOffset = 0;
      int oldCount = 0;

      if (oldDest < com->size())
        do
        {
          idx++;
          dest[idx] = oldDest++;

          sendOffset[idx] = oldOffset + oldCount;
          oldOffset = sendOffset[idx];

          int offset = (com->rank() * res[NDIM - 1]) + oldOffset - dest[idx] * smallRes[NDIM - 1];
          int nMaxRcv = (dest[idx] == com->size() - 1) ? lastSmallRes[NDIM - 1] : smallRes[NDIM - 1];
          nMaxRcv -= offset;
          if (nMaxRcv < 0)
            break;

          sendCount[idx] = std::min(curRes[NDIM - 1] - sendOffset[idx], nMaxRcv);
          oldCount = sendCount[idx];
        } while ((sendOffset[idx] + sendCount[idx] < curRes[NDIM - 1]) &&
                 (dest[idx] != com->size() - 1));
      /*
      glb::console->print<LOG_STD_ALL>
  ("process %d sends to [%d: %d @%d] [%d: %d @%d] [%d: %d @%d]\n",
   com->rank(),
   dest[0],sendCount[0],sendOffset[0],
   dest[1],sendCount[1],sendOffset[1],
   dest[2],sendCount[2],sendOffset[2]);

      glb::console->print<LOG_STD_ALL>
  ("process %d [%d] receive from [%d: %d @%d] [%d: %d @%d]\n",
   com->rank(),curSmallRes[NDIM-1],
   src[0],receiveCount[0],receiveOffset[0],
   src[1],receiveCount[1],receiveOffset[1]);

      com->barrier();
      //exit(-1);
      */

      MPI_Request req[5];
      int nReq = 0;
      int subPos[NDIM] = {0};
      int subRes[NDIM];

      for (int i = 0; i < NDIM - 1; ++i)
        subRes[i] = smallRes[i];

      for (int i = 0; sendCount[i]; ++i)
      {
        subPos[NDIM - 1] = sendOffset[i];
        subRes[NDIM - 1] = sendCount[i];
        // FIXME: I do not know exaclty how fast this is ...
        MPI_Datatype subType;
        MPI_Type_create_subarray(NDIM, res, subRes, subPos, MPI_ORDER_FORTRAN,
                                 MPI_Type<Data>::get(), &subType);
        MPI_Type_commit(&subType);

        com->IsendMpiType(localGrid->getDataPtr(), 1, subType, dest[i], &req[nReq++]);

        MPI_Type_free(&subType);
      }

      int sliceCount = 1;
      for (int i = 0; i < NDIM - 1; ++i)
        sliceCount *= smallRes[i];
      std::vector<Data> buffer(sliceCount * res[NDIM - 1]);

      for (int i = 0; receiveCount[i]; ++i)
      {
        com->Irecv(&buffer[receiveOffset[i] * sliceCount],
                   receiveCount[i] * sliceCount, src[i], &req[nReq++]);
      }

      com->Waitall(nReq, req);
      std::copy(buffer.begin(), buffer.end(), localFacade->getDataPtr());
      // com->barrier();
      // gridFacade->toVtk("GridFacade","GridFacade_%d");
    }

  private:
    void initialize(const Params &gp_,
                    MpiCommunication *com_,
                    int nThreads_,
                    bool periodic_,
                    bool facade = false)
    {
      gp = gp_;
      mpiCom = com_;
      numThreads = nThreads_;
      periodic = periodic_;
      nChunks = mpiCom->size();

      // We'll use padding in the non-periodic case
      if ((!periodic) && (!facade))
      {
        for (int i = 0; i < NDIM; ++i)
        {
          gp.delta[i] *= 2;
          gp.resolution[i] *= 2;
        }
      }

      local_n0.assign(nChunks, 0);
      local_0_start.assign(nChunks, 0);
      nAlloc.assign(nChunks, 0);

      int rank = mpiCom->rank();
      ptrdiff_t N[NDIM];
      ptrdiff_t Nr[NDIM];

      // We will only perform r2c and c2r transforms, so we do not need to allocate
      // full arrays
      std::copy(gp.resolution, gp.resolution + NDIM, N);

      N[0] = N[0] / 2 + 1;
      for (int i = 0; i < NDIM; ++i)
        Nr[i] = N[NDIM - i - 1];

      nAlloc[rank] = fftw_mpi_local_size(NDIM, Nr, mpiCom->getCom(),
                                         &local_n0[rank],
                                         &local_0_start[rank]);
      nAlloc[rank] *= 2; // because we allocate double, not fftw_complex elements...

      mpiCom->Allgather_inplace(local_n0);
      mpiCom->Allgather_inplace(local_0_start);
      mpiCom->Allgather_inplace(nAlloc);

      slice(mpiCom->rank());
    }

    struct Neighbor
    {
      Neighbor(int dm, int dr, int id) : dim(dm), dir(dr), index(id)
      {
      }

      void set(int dm, int dr, int id)
      {
        dim = dm;
        dir = dr;
        index = id;
      }

      int dim;
      int dir;
      int index;
    };

    Params globalGridParams;
    Params gridParams;

    std::vector<Neighbor> neighbors;

    int sliceCount[NDIM]; // int slicerSize[DIMS];  // number of slices in each dim
    int myIndex;          // int myID;
    int slicePos[NDIM];   // int myCoords[DIMS]; // coords within the slices
    // int localGridPos[NDIM];//int myFullGridPos[DIMS]; // coords within the fullgrid
    bool periodic;

    int numThreads;
    MpiCommunication *mpiCom;
    Params gp;
    long nChunks;
    long allocFactor;

    std::vector<ptrdiff_t> local_n0;
    std::vector<ptrdiff_t> local_0_start;
    std::vector<ptrdiff_t> nAlloc;
  };
} // internal

#include "../../internal/namespace.footer"
#endif

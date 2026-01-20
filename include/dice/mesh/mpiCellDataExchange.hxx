#ifndef __MPI_CELL_DATA_EXCHANGE_HXX__
#define __MPI_CELL_DATA_EXCHANGE_HXX__

#include <vector>

#include "../dice_globals.hxx"

#include "../tools/MPI/mpiCommunication.hxx"
#include "../tools/MPI/mpiDataType.hxx"

/**
 * @file
 * @brief A helper class used to communicate data between cells of a mesh
 * and their shadow/ghosts images on remote nodes.
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup MESH
 *   \{
 */

/**
 * \class MpiCellDataExchangeT
 * \brief This helper structure is used to transfer data from cells to their ghosts/shadows
 * ghosts/shadows copies/existing on other processes
 *
 * Exchange functions have a template argument class DT that defines the type of
 * the class used to stores the actual data to transfer. This class should define a
 * DT::MpiStruct type to which it can be converted and that can be used for MPI transfer
 * (i.e. a POD structure). Similarly, DT::IS_INDEXED indicates whether the simplices
 * to transfer are indexed (i.e. a subset of the cells designated by their index can be
 * selected for the transfer) or if they should all be transfered. Using indexed transfer
 * is preferable only if a small number of cell need to tranfer data as the additional
 * index must also bee added to the transfered data ...
 *
 * Simplices send[i] are sent to process with rank sendRank[i].
 * Exchange simplices receive[i] are received from processes with rank receiveRank[i]
 *
 * data exchange can be achieved using exchange(...) functions, or by
 * explicitly sending/receiving data with the isend/ireceive and waiting for
 * the exchange to finish with wait(...).
 * The later version may be faster in some cases when one needs to cache
 * computations while the transfer occurs ...
 * \tparam C  The type of local cell (e.g. Simplex or Vertex)
 * \tparam EC The type of the exchange cell (e.g. the image of a remote cell, such as a
 *            GhostSimplex or ShadowSimplex)
 */

template <class C, class EC>
struct MpiCellDataExchangeT
{
  typedef C Cell;
  typedef EC ExchangeCell;

  typedef MpiCellDataExchangeT<C, EC> MyType;

  static std::string classHeader() { return "mpi_cell_exchange"; }
  static float classVersion() { return 0.10; }
  static float compatibleSinceClassVersion() { return 0.10; }

  const MpiCommunication *mpiCom;
  const int defaultTag;

  std::vector<int> sendRank;    /**< the rank of the processes to which data needs to be sent*/
  std::vector<int> receiveRank; /**< the rank of processes from which data is received*/
  std::vector<long> nSendCum;
  std::vector<long> nReceiveCum;

  /** An array of the simplices on the local node that have images on another node.
   *  send[i] contains local cells that have a ghost / shadow image on the node
   *  with rank sendRank[i]
   */
  std::vector<std::vector<Cell *>> send;

  /** An array of the shared/shadow simplices on the local. receive[i] contains cells
   *  whose local image is in send[i] on node with rank receiveRank[i].
   */
  std::vector<std::vector<ExchangeCell *>> receive;

  MpiCellDataExchangeT(MpiCommunication *mpiCom_ = glb::mpiComWorld) : mpiCom(mpiCom_),
                                                                       defaultTag(mpiCom_->reserveTags(1)),
                                                                       send(mpiCom_->size()),
                                                                       receive(mpiCom_->size())
  {
    mpiCom->barrier();
  }

  template <class W>
  void serialize(W *writer)
  {
    writer->writeHeader(classHeader(), classVersion());

    writer->write(&send);
    writer->write(&receive);
    writer->write(&sendRank);
    writer->write(&receiveRank);
    writer->write(&nSendCum);
    writer->write(&nReceiveCum);
  }

  template <class CU, class ECU>
  void updatePointers(const CU &cellUpdater, const ECU &exchangeCellUpdater)
  {
    for (unsigned long i = 0; i < send.size(); ++i)
      for (unsigned long j = 0; j < send[i].size(); ++j)
      {
        send[i][j] = cellUpdater(send[i][j]);
        if (send[i][j] == NULL)
        {
          PRINT_SRC_INFO(LOG_ERROR);
          glb::console->print<LOG_ERROR>("when unserializing mpiCellExchange: could no update 'send' element.\n");
          exit(-1);
        }
      }

    for (unsigned long i = 0; i < receive.size(); ++i)
      for (unsigned long j = 0; j < receive[i].size(); ++j)
      {
        receive[i][j] = exchangeCellUpdater(receive[i][j]);
        if (receive[i][j] == NULL)
        {
          PRINT_SRC_INFO(LOG_ERROR);
          glb::console->print<LOG_ERROR>("when unserializing mpiCellExchange: could no update 'receive' element.\n");
          exit(-1);
        }
      }
  }

  template <class R, class CU, class ECU>
  void unSerialize(R *reader, const CU &cellUpdater, const ECU &exchangeCellUpdater)
  {
    float version;
    R::template checkHeaderAndReport<LOG_ERROR, LOG_WARNING, MyType>(glb::console, reader, version, true);

    reader->read(&send);
    reader->read(&receive);
    reader->read(&sendRank);
    reader->read(&receiveRank);
    reader->read(&nSendCum);
    reader->read(&nReceiveCum);

    updatePointers(cellUpdater, exchangeCellUpdater);
    /*
    for (unsigned long i=0;i<send.size();++i)
      for (unsigned long j=0;j<send[i].size();++j)
  {
    send[i][j]=cellUpdater(send[i][j]);
    if (send[i][j]==NULL)
      {
        PRINT_SRC_INFO(LOG_ERROR);
        glb::console->print<LOG_ERROR>("when unserializing mpiCellExchange: could no update 'send' element.\n");
        exit(-1);
      }
  }

    for (unsigned long i=0;i<receive.size();++i)
      for (unsigned long j=0;j<receive[i].size();++j)
  {
    receive[i][j]=exchangeCellUpdater(receive[i][j]);
    if (receive[i][j]==NULL)
      {
        PRINT_SRC_INFO(LOG_ERROR);
        glb::console->print<LOG_ERROR>("when unserializing mpiCellExchange: could no update 'receive' element.\n");
        exit(-1);
      }
  }
    */
  }

  template <class CT>
  void reclaimUnusedMemory(CT &data)
  {
    CT tmp;
    data.swap(tmp);
  }

  /** \brief Deallocates the memory of all cleared vector.
   * That's equivalent to calling clear() followed by a shrink_to_fit() from C++11 specs
   */
  void reclaimUnusedMemory()
  {
    for (unsigned long i = 0; i < send.size(); ++i)
    {
      if ((send[i].size() == 0) && (send[i].capacity() != 0))
        reclaimUnusedMemory(send[i]);
    }
    for (unsigned long i = 0; i < receive.size(); ++i)
    {
      if ((receive[i].size() == 0) && (receive[i].capacity() != 0))
        reclaimUnusedMemory(receive[i]);
    }
  }

  std::vector<long> getNSendCum()
  {
    std::vector<long> result(sendRank.size() + 1);
    result[0] = 0;
    for (unsigned long i = 0; i < sendRank.size(); i++)
      result[i + 1] = result[i] + send[sendRank[i]].size();
    return result;
  }
  void updateNSendCum()
  {
    nSendCum = getNSendCum();
  }

  std::vector<long> getNReceiveCum()
  {
    std::vector<long> result(receiveRank.size() + 1);
    result[0] = 0;
    for (unsigned long i = 0; i < receiveRank.size(); i++)
      result[i + 1] = result[i] + receive[receiveRank[i]].size();
    return result;
  }
  void updateNReceiveCum()
  {
    nReceiveCum = getNReceiveCum();
  }
  void updateNCum()
  {
    updateNSendCum();
    updateNReceiveCum();
  }

  template <class L>
  void print(const char *what)
  {
    if (!glb::console->willPrint<L>())
      return;
    int myRank = mpiCom->rank();
    int nParts = mpiCom->size();
    mpiCom->barrier();
    glb::console->setRankPrint(false);
    for (int i = 0; i < nParts; i++)
    {
      if (i == myRank)
      {
        if (receiveRank.size() == 0)
          glb::console->printToBuffer<L>("Process %d receives %s from 0 neighbors.\n", i, what);
        else
        {
          glb::console->printToBuffer<L>("Process %d receives %s from %d neighbors: %ld(%ld)",
                                         i, what, receiveRank.size(),
                                         receiveRank[0],
                                         receive[receiveRank[0]].size());
          for (unsigned long i = 1; i < receiveRank.size(); i++)
            glb::console->printToBuffer<L>(",%ld(%ld)", receiveRank[i],
                                           receive[receiveRank[i]].size());
          glb::console->printToBuffer<L>(".\n");

          // for (int i=0;i<receiveRank.size();i++)
          //   {
          //     glb::console->printToBuffer<L>(" RECEIVE RANK %ld -> %ld: ",(long)receiveRank[i],myRank);
          //   for (int j=0;j<receive[receiveRank[i]].size();++j)
          //     glb::console->printToBuffer<L>("(%ld,%ld),",
          // 				   receive[receiveRank[i]][j]->getGlobalIdentity().rank(),
          // 				   receive[receiveRank[i]][j]->getGlobalIdentity().id());
          //   glb::console->printToBuffer<L>(".\n");
          //   }

          glb::console->flushBuffer<L>();
        }

        if (sendRank.size() == 0)
          glb::console->printToBuffer<L>("Process %d sends %s to 0 neighbors.\n", i, what);
        else
        {
          glb::console->printToBuffer<L>("Process %d sends %s to %d neighbors: %ld(%ld)",
                                         i, what, sendRank.size(),
                                         sendRank[0],
                                         send[sendRank[0]].size());
          for (unsigned long i = 1; i < sendRank.size(); i++)
            glb::console->printToBuffer<L>(",%ld(%ld)", sendRank[i],
                                           send[sendRank[i]].size());

          glb::console->printToBuffer<L>(".\n");
          // for (int i=0;i<sendRank.size();i++)
          //   {
          //     glb::console->printToBuffer<L>(" SEND RANK %ld -> %ld: ",(long)sendRank[i],myRank);
          //   for (int j=0;j<send[sendRank[i]].size();++j)
          //     glb::console->printToBuffer<L>("(%ld,%ld),",
          // 				   send[sendRank[i]][j]->getGlobalIdentity(myRank).rank(),
          // 				   send[sendRank[i]][j]->getGlobalIdentity(myRank).id());
          //   glb::console->printToBuffer<L>(".\n");
          //   }

          glb::console->flushBuffer<L>();
        }
      }
    }
    glb::console->setRankPrint(true);
  }
  /*
  template <class DT>
  void exchange(std::vector< std::vector<DT> > &toSend,
    std::vector<std::vector<typename DT::MpiStruct> > &receiveData,
    MPI_Datatype mpiStructType,
    MpiCommunication *mpiCom)
  {
    const long nSend=sendRank.size();
    const long nReceive=receiveRank.size();
    const long nTot = nSend+nReceive;

    std::vector< std::vector<typename DT::MpiStruct> > sendData(nSend);
    MPI_Request req[nTot];

    //for (long i=0;i<nSend;i++)
    //sendData[i].resize(send[sendRank[i]].size());
    for (long i=0;i<nSend;i++)
      sendData[i].resize(toSend[i].size());

    receiveData.resize(nReceive);
    for (long i=0;i<nReceive;i++)
      receiveData[i].resize(receive[receiveRank[i]].size());

    for (long i=0;i<nReceive;i++)
      {
  mpiCom->Irecv(&receiveData[i][0],receiveData[i].size(),
          mpiStructType,receiveRank[i],&req[i]);
      }

    for (int i=0;i<nSend;i++)
      {
  //DT *toSendCur = &toSend[nSendCum[i]];
#pragma omp parallel for
  for (long j=0;j<send[sendRank[i]].size();j++)
    sendData[i][j]=toSend[i][j].getMpiStruct();
  mpiCom->Isend(&sendData[i][0],sendData[i].size(),
          mpiStructType,sendRank[i],&req[i+nReceive]);

      }

    for (int i=0;i<nTot;i++) mpiCom->wait(&req[i]);
  }
  */
  /*
  template <class DT>
  void exchange(std::vector< std::vector<DT> >&toSend,
    std::vector<DT> &receiveData,
    MPI_Datatype mpiStructType,
    MpiCommunication *mpiCom)
  {
    std::vector< std::vector<typename DT::MpiStruct> > mpiReceive;
    exchange(toSend,mpiReceive,mpiStructType,mpiCom);

    receiveData.resize(nReceiveCum.back());
    for (int i=0;i<mpiReceive.size();i++)
      {
  std::vector<ExchangeCell*> &receiveCell =
    receive[receiveRank[i]];
#pragma omp parallel for
  for (long j=0;j<mpiReceive[i].size();j++)
    receiveData[receiveCell[j]->getLocalIndex()].set(receiveCell[j],mpiReceive[i][j]);
      }
  }
  */

  // send a mpi request for all the simplices that need to be received
  template <class DT_MPI>
  void iReceive(std::vector<std::vector<DT_MPI>> &receiveData,
                std::vector<MPI_Request> &reqs,
                MPI_Datatype mpiStructType, int tag = 0)
  {
    const long nSend = sendRank.size();
    const long nReceive = receiveRank.size();

    receiveData.resize(nReceive);
    for (long i = 0; i < nReceive; i++)
      receiveData[i].resize(receive[receiveRank[i]].size());

    reqs.resize(nSend + nReceive);

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->IrecvMpiType(&receiveData[i][0], receiveData[i].size(),
                           mpiStructType, receiveRank[i], &reqs[i], tag);
    }
  }

  // send a mpi request for all the simplices that need to be sent
  // to proces with rank sendRank[which]
  template <class DT>
  void iSend(int which, std::vector<DT> &toSend,
             std::vector<typename DT::MpiStruct> &sendData,
             std::vector<MPI_Request> &reqs,
             MPI_Datatype mpiStructType, int tag = 0)
  {
    const long nReceive = receiveRank.size();
    sendData.resize(toSend.size());

#pragma omp parallel for
    for (unsigned long j = 0; j < toSend.size(); j++)
      sendData[j] = toSend[j].getMpiStruct();

    mpiCom->IsendMpiType(&sendData[0], sendData.size(),
                         mpiStructType, sendRank[which], &reqs[nReceive + which], tag);
  }

  // send a mpi request for all the simplices that need to be sent
  // to proces with rank sendRank[which]
  template <class DT>
  void iSend(int which, std::vector<DT> &toSend,
             std::vector<MPI_Request> &reqs,
             MPI_Datatype mpiStructType, int tag = 0)
  {
    const long nReceive = receiveRank.size();

    mpiCom->IsendMpiType(&toSend[0], toSend.size(),
                         mpiStructType, sendRank[which], &reqs[nReceive + which], tag);
  }

  // send a mpi request for all the simplices that need to be sent
  // to proces with rank sendRank[which]
  template <class DT>
  void iSend(std::vector<std::vector<DT>> &toSend,
             std::vector<MPI_Request> &reqs,
             MPI_Datatype mpiStructType, int tag = 0)
  {
    const long nReceive = receiveRank.size();
    const long nSend = sendRank.size();

    for (long i = 0; i < nSend; i++)
      mpiCom->IsendMpiType(&toSend[i][0], toSend[i].size(),
                           mpiStructType, sendRank[i], &reqs[nReceive + i], tag);
  }

  // wait for all the interprocess comm staterted with iSend/iReceive to finish
  // and retieve data in the MPI-friendly format DT::MpiStruct used for transfer
  template <class DT>
  void wait(std::vector<std::vector<typename DT::MpiStruct>> &mpiReceive,
            std::vector<MPI_Request> &reqs,
            MPI_Datatype mpiStructType,
            MpiCommunication *mpiCom)
  {
    const long nReceive = receiveRank.size();

    for (unsigned long i = 0; i < reqs.size(); i++)
    {
      MPI_Status status;
      // FIXME: Use wait_all here ?
      mpiCom->Wait(&reqs[i], &status);
      if ((DT::IS_INDEXED) && (i < nReceive))
        mpiReceive[i].resize(mpiCom->getCount(&status, mpiStructType));
    }
  }

  // wait for all the interprocess comm staterted with iSend/iReceive to finish
  // and retieve data in the regular data type DT
  template <class DT>
  void wait(std::vector<std::vector<typename DT::MpiStruct>> &mpiReceive,
            std::vector<DT> &receiveData,
            std::vector<MPI_Request> &reqs,
            MPI_Datatype mpiStructType)
  {
    const long nReceive = receiveRank.size();

    for (unsigned long i = 0; i < reqs.size(); i++)
    {
      MPI_Status status;
      // FIXME: Use wait_all here ?
      mpiCom->Wait(&reqs[i], &status);
      if ((DT::IS_INDEXED) && (i < nReceive))
        mpiReceive[i].resize(mpiCom->getCount(&status, mpiStructType));
    }
    convert(mpiReceive, receiveData);
  }

  // convert MPI-friendly data DT::MpiStruct used for transfer to regular data type DT
  template <class DT>
  void convert(std::vector<std::vector<typename DT::MpiStruct>> &mpiReceive,
               std::vector<DT> &receiveData)
  {
    if (DT::IS_INDEXED)
    {
      receiveData.clear();
      for (unsigned long i = 0; i < mpiReceive.size(); i++)
      {
        std::vector<ExchangeCell *> &receiveCell = receive[receiveRank[i]];
        for (unsigned long j = 0; j < mpiReceive[i].size(); j++)
        {
          DT data;
          data.set(receiveCell[mpiReceive[i][j].getBaseCellIndex()], mpiReceive[i][j]);
          receiveData.push_back(data);
        }
      }
    }
    else
    {
      receiveData.resize(nReceiveCum.back());
      for (unsigned long i = 0; i < mpiReceive.size(); i++)
      {
        std::vector<ExchangeCell *> &receiveCell =
            receive[receiveRank[i]];
#pragma omp parallel for
        for (unsigned long j = 0; j < mpiReceive[i].size(); j++)
          receiveData[receiveCell[j]->getLocalIndex()].set(receiveCell[j], mpiReceive[i][j]);
      }
    }
  }

  // convert MPI-friendly data DT::MpiStruct used for transfer to regular data type DT
  template <class DT>
  void reversedConvert(std::vector<std::vector<typename DT::MpiStruct>> &mpiReceive,
                       std::vector<DT> &receiveData)
  {
    if (DT::IS_INDEXED)
    {
      receiveData.clear();
      for (unsigned long i = 0; i < mpiReceive.size(); i++)
      {
        std::vector<Cell *> &receiveCell = send[sendRank[i]];
        for (unsigned long j = 0; j < mpiReceive[i].size(); j++)
        {
          DT data;
          data.set(receiveCell[mpiReceive[i][j].getBaseCellIndex()], mpiReceive[i][j]);
          receiveData.push_back(data);
        }
      }
    }
    else
    {
      receiveData.resize(nReceiveCum.back());
      for (unsigned long i = 0; i < mpiReceive.size(); i++)
      {
        std::vector<Cell *> &receiveCell = send[sendRank[i]];
#pragma omp parallel for
        for (unsigned long j = 0; j < mpiReceive[i].size(); j++)
          receiveData[receiveCell[j]->getLocalIndex()].set(receiveCell[j], mpiReceive[i][j]);
      }
    }
  }

  // exchange all the data, output is in the MPI friendly format DT::MpiStruct used for transfer
  template <class DT>
  void exchangeStruct(std::vector<std::vector<DT>> &toSend,
                      std::vector<std::vector<typename DT::MpiStruct>> &receiveData,
                      MPI_Datatype mpiStructType, bool autoResize = true, int tag = 0)
  {
    const long nSend = sendRank.size();
    const long nReceive = receiveRank.size();
    const long nTot = nSend + nReceive;

    std::vector<std::vector<typename DT::MpiStruct>> sendData(nSend);
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    for (long i = 0; i < nSend; i++)
      sendData[i].resize(toSend[i].size());

    receiveData.resize(nReceive);

    if (autoResize)
    {
      for (long i = 0; i < nReceive; i++)
        receiveData[i].resize(receive[receiveRank[i]].size());
    }

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->IrecvMpiType(&receiveData[i][0], receiveData[i].size(),
                           mpiStructType, receiveRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
#pragma omp parallel for
      for (unsigned long j = 0; j < toSend[i].size(); j++)
        sendData[i][j] = toSend[i][j].getMpiStruct();

      mpiCom->IsendMpiType(&sendData[i][0], sendData[i].size(),
                           mpiStructType, sendRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
    if (DT::IS_INDEXED)
    {
      for (long i = 0; i < nReceive; ++i)
        receiveData[i].resize(mpiCom->getCount(&status[i], mpiStructType));
    }
    /*
    for (int i=0;i<nTot;i++)
      {
  //MPI_Status status;
  mpiCom->Wait(&req[i],&status[i]);
  if ((DT::IS_INDEXED)&&(i<nReceive))
    receiveData[i].resize(mpiCom->getCount(&status[i],mpiStructType));
      }
    */
  }

  // exchange all the data, output is in the MPI friendly format DT::MpiStruct used
  // for transfer
  template <class DT>
  void exchangeStruct(std::vector<std::vector<DT>> &toSend,
                      std::vector<std::vector<DT>> &receiveData,
                      MPI_Datatype mpiStructType, bool autoResize = true, int tag = 0)
  {
    const long nSend = sendRank.size();
    const long nReceive = receiveRank.size();
    const long nTot = nSend + nReceive;
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    receiveData.resize(nReceive);

    if (autoResize)
    {
      for (long i = 0; i < nReceive; i++)
        receiveData[i].resize(receive[receiveRank[i]].size());
    }

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->IrecvMpiType(&receiveData[i][0], receiveData[i].size(),
                           mpiStructType, receiveRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
      mpiCom->IsendMpiType(&toSend[i][0], toSend[i].size(),
                           mpiStructType, sendRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
    if (DT::IS_INDEXED)
    {
      for (long i = 0; i < nReceive; ++i)
        receiveData[i].resize(mpiCom->getCount(&status[i], mpiStructType));
    }
    /*
    for (long i=0;i<nTot;++i)
      {
  MPI_Status status;
  int index;
  mpiCom->WaitAny(&req[i],index,&status);
  if ((DT::IS_INDEXED)&&(index<nReceive))
    receiveData[index].resize(mpiCom->getCount(&status,mpiStructType));
      }
    */
    /*
    for (long i=0;i<nTot;i++)
      {
  MPI_Status status;
  //FIXME: Use wait_all here , use probe ?
  mpiCom->Wait(&req[i],&status);
  if ((DT::IS_INDEXED)&&(i<nReceive))
    receiveData[i].resize(mpiCom->getCount(&status,mpiStructType));
      }
    */
  }

  // exchange all the data, DT is a basic type (no mpiStruct needed)
  template <class DT>
  void exchange(std::vector<std::vector<DT>> &toSend,
                std::vector<std::vector<DT>> &receiveData,
                int tag = 0)
  {
    const long nSend = sendRank.size();
    const long nReceive = receiveRank.size();
    const long nTot = nSend + nReceive;
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    receiveData.resize(nReceive);

    for (long i = 0; i < nReceive; i++)
      receiveData[i].resize(receive[receiveRank[i]].size());

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->Irecv(&receiveData[i][0], receiveData[i].size(),
                    receiveRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
      mpiCom->Isend(&toSend[i][0], toSend[i].size(),
                    sendRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
  }

  // exchange the data but send and receive are reversed.
  // DT is a basic type (no mpiStruct needed)
  template <class DT>
  void reversedExchange(std::vector<std::vector<DT>> &toSend,
                        std::vector<std::vector<DT>> &receiveData,
                        int tag = 0)
  {
    const long nSend = receiveRank.size();
    const long nReceive = sendRank.size();
    const long nTot = nSend + nReceive;
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    receiveData.resize(nReceive);

    for (long i = 0; i < nReceive; i++)
      receiveData[i].resize(send[sendRank[i]].size());

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->Irecv(&receiveData[i][0], receiveData[i].size(),
                    sendRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
      mpiCom->Isend(&toSend[i][0], toSend[i].size(),
                    receiveRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
  }

  // exchange the data but send and receive are reversed ...
  template <class DT>
  void reversedExchangeStruct(std::vector<std::vector<DT>> &toSend,
                              std::vector<std::vector<DT>> &receiveData,
                              MPI_Datatype mpiStructType, bool autoResize = true, int tag = 0)
  {
    const long nSend = receiveRank.size();
    const long nReceive = sendRank.size();
    const long nTot = nSend + nReceive;
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    receiveData.resize(nReceive);

    if (autoResize)
    {
      for (long i = 0; i < nReceive; i++)
        receiveData[i].resize(send[sendRank[i]].size());
    }

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->IrecvMpiType(&receiveData[i][0], receiveData[i].size(),
                           mpiStructType, sendRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
      mpiCom->IsendMpiType(&toSend[i][0], toSend[i].size(),
                           mpiStructType, receiveRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
    if (DT::IS_INDEXED)
    {
      for (long i = 0; i < nReceive; ++i)
        receiveData[i].resize(mpiCom->getCount(&status[i], mpiStructType));
    }
    /*
    for (int i=0;i<nTot;i++)
      {
  MPI_Status status;
  //FIXME: Use wait_all here , use probe ?
  mpiCom->Wait(&req[i],&status[i]);
  if ((DT::IS_INDEXED)&&(i<nReceive))
    receiveData[i].resize(mpiCom->getCount(&status[i],mpiStructType));
      }
    */
  }

  // exchange all the data, output is in the MPI friendly format DT::MpiStruct used for transfer
  template <class DT>
  void reversedExchangeStruct(std::vector<std::vector<DT>> &toSend,
                              std::vector<std::vector<typename DT::MpiStruct>> &receiveData,
                              MPI_Datatype mpiStructType, bool autoResize = true, int tag = 0)
  {
    const long nSend = receiveRank.size();
    const long nReceive = sendRank.size();
    const long nTot = nSend + nReceive;

    std::vector<std::vector<typename DT::MpiStruct>> sendData(nSend);
    MPI_Request req[nTot];
    MPI_Status status[nTot];

    for (long i = 0; i < nSend; i++)
      sendData[i].resize(toSend[i].size());

    receiveData.resize(nReceive);

    if (autoResize)
    {
      for (long i = 0; i < nReceive; i++)
        receiveData[i].resize(send[sendRank[i]].size());
    }

    for (long i = 0; i < nReceive; i++)
    {
      mpiCom->IrecvMpiType(&receiveData[i][0], receiveData[i].size(),
                           mpiStructType, sendRank[i], &req[i], tag);
    }

    for (long i = 0; i < nSend; i++)
    {
#pragma omp parallel for
      for (unsigned long j = 0; j < toSend[i].size(); j++)
        sendData[i][j] = toSend[i][j].getMpiStruct();

      mpiCom->IsendMpiType(&sendData[i][0], sendData[i].size(),
                           mpiStructType, receiveRank[i], &req[i + nReceive], tag);
    }

    mpiCom->Waitall(nTot, req, status);
    if (DT::IS_INDEXED)
    {
      for (long i = 0; i < nReceive; ++i)
        receiveData[i].resize(mpiCom->getCount(&status[i], mpiStructType));
    }
    /*
    for (int i=0;i<nTot;i++)
      {
  //MPI_Status status;
  mpiCom->Wait(&req[i],&status[i]);
  if ((DT::IS_INDEXED)&&(i<nReceive))
    receiveData[i].resize(mpiCom->getCount(&status[i],mpiStructType));
      }
    */
  }

  // exchange and convert MPI friendly output DT::MpiStruct used for transfer
  // to the input data type DT
  template <class DT>
  void exchangeStruct(std::vector<std::vector<DT>> &toSend,
                      std::vector<DT> &receiveData,
                      MPI_Datatype mpiStructType, int tag = 0)
  {
    std::vector<std::vector<typename DT::MpiStruct>> mpiReceive;
    exchangeStruct(toSend, mpiReceive, mpiStructType, true, tag);
    convert(mpiReceive, receiveData);
  }

  // exchange and convert MPI friendly output DT::MpiStruct used for transfer
  // to the input data type DT
  template <class DT>
  void reversedExchangeStruct(std::vector<std::vector<DT>> &toSend,
                              std::vector<DT> &receiveData,
                              MPI_Datatype mpiStructType, int tag = 0)
  {
    std::vector<std::vector<typename DT::MpiStruct>> mpiReceive;
    reversedExchangeStruct(toSend, mpiReceive, mpiStructType, true, tag);
    reversedConvert(mpiReceive, receiveData);
  }
};

/** \}*/
#include "../internal/namespace.footer"
#endif

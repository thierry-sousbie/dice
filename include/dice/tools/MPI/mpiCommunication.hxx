#ifndef __MPI_COM_HXX__
#define __MPI_COM_HXX__

#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <list>
#include <algorithm>

#ifndef USE_MPI
#include <ctime>
#endif

#include <string.h>

#include "../../dice_globals.hxx"

#include "../../tools/OMP/openMP_interface.hxx"
#include "./myMpi.hxx"
#include "./mpiDataType.hxx"

#include "../../internal/namespace.header"

#pragma GCC diagnostic ignored "-Wunused-private-field"

// internal namespace for local  stuff
namespace internal
{

  template <bool disabled = false>
  class MpiOmpLockerT
  {
  public:
    MpiOmpLockerT()
    {
      omp_init_lock(&lock);
    }

    ~MpiOmpLockerT()
    {
      omp_destroy_lock(&lock);
    }

    void set()
    {
      omp_set_lock(&lock);
    }

    void unset()
    {
      omp_unset_lock(&lock);
    }

  private:
    omp_lock_t lock;
  };

  template <>
  class MpiOmpLockerT<true>
  {
  public:
    static void set() {}
    static void unset() {}
  };

  template <bool disabled = false>
  class MpiOmpLockChekerT
  {
  public:
    MpiOmpLockChekerT(MpiOmpLockerT<disabled> &locker, int level) : locker_(locker),
                                                                    level_(level)
    {
      if (level_ == 2)
        locker_.set();                                       // SERIALIZED
      else if ((level_ == 1) && (omp_get_thread_num() != 0)) // FUNNELED
      {
        fprintf(stderr, "ERROR: trying to make MPI call from thread %d.\n",
                omp_get_thread_num());
        fprintf(stderr, "       MPI is now configured in MPI_THREAD_FUNNELED mode, meaning only the main thread should make MPI calls.\n");
        fprintf(stderr, "       Increase threading support to MPI_THREAD_SERIALIZED or higher.\n");
        printf("ERROR: trying to make MPI call from thread %d.\n",
               omp_get_thread_num());
        printf("       MPI is now configured in MPI_THREAD_FUNNELED mode, meaning only the main thread should make MPI calls.\n");
        printf("       Increase threading support to MPI_THREAD_SERIALIZED or higher.\n");
        exit(-1);
      }
    }

    ~MpiOmpLockChekerT()
    {
      if (level_ == 2)
        locker_.unset(); // SERIALIZED
    }

  private:
    MpiOmpLockerT<disabled> &locker_;
    const int level_;
  };

  template <>
  class MpiOmpLockChekerT<true>
  {
  public:
    MpiOmpLockChekerT(MpiOmpLockerT<true> &locker, int level)
    {
    }
  };

  template <bool disabled = false>
  class MpiCallTimerT
  {
  public:
    MpiCallTimerT(TimerPool::Timer *t, int th = 0) : timer(t),
                                                     thread(th)
    {
      if (omp_get_thread_num() == thread)
        timer->start();
    }

    ~MpiCallTimerT()
    {
      if (omp_get_thread_num() == thread)
        timer->stop();
    }

  private:
    TimerPool::Timer *timer;
    int thread;
  };

  template <>
  class MpiCallTimerT<true>
  {
  public:
    MpiCallTimerT(TimerPool::Timer *t, int th = 0) {}
  };
}

class MpiCommunication
{

private:
  static const int maxReservedTags = (MPID_TAG_UB >> 1);

  typedef internal::MpiCallTimerT<false> MpiCallTimer;

  int myRank;
  int nProcs;
  MPI_Comm com;
  bool finalize;
  bool initialized;
  bool deleteCom;
  int reservedTagsStart;
  int originThreadId;
  int threadSupport;

  TimerPool::Timer *localComTimer;
  TimerPool::Timer *globalComTimer;
  TimerPool::Timer *barrierTimer;

  std::list<MPI_Request> pendingRequests;
#ifdef USE_MPI

  typedef internal::MpiOmpLockerT<false> MpiOmpLocker;
  typedef internal::MpiOmpLockChekerT<false> MpiOmpLockChecker;
  mutable MpiOmpLocker locker;

private:
  void initCom(MPI_Comm myCom)
  {
    com = myCom;
    MPI_Comm_rank(com, &myRank);
    MPI_Comm_size(com, &nProcs);
  }

  void init(int *argc, char ***argv)
  {
    const char *threadLevel[4] = {"MPI_THREAD_SINGLE",
                                  "MPI_THREAD_FUNNELED",
                                  "MPI_THREAD_SERIALIZED",
                                  "MPI_THREAD_MULTIPLE"};

    if (!initialized)
    {
      finalize = true;
      deleteCom = false;

#ifdef USE_OPENMP
      int required = MPI_THREAD_MULTIPLE;
      int provided;
      MPI_Init_thread(argc, argv, required, &provided);
      printf("MPI initialized with the following thread support: %s\n",
             threadLevel[provided]);
      if (provided == 0)
      {
        fprintf(stderr, "ERROR: MPI initialisation error !\n");
        fprintf(stderr, "       This MPI implementation does not support any type of multithreading.\n");
        fprintf(stderr, "       Disable multiThreading or change your MPI library.\n");

        printf("ERROR: MPI initialisation error !\n");
        printf("       This MPI implementation does not support multithreading.\n");
        printf("       Disable multiThreading or change your MPI library.\n");
        exit(-1);
      }
      threadSupport = provided;
#else
      threadSupport = 0;
      printf("MPI initialized with the following thread support: %s\n",
             threadLevel[threadSupport]);
      MPI_Init(argc, argv);
#endif
    }
    else
    {
      fprintf(stderr, "ERROR: trying to initialize MPI more than once!\n");
      fprintf(stderr, "       init() may only be called once !\n");
      printf("ERROR: trying to initialize MPI more than once!\n");
      printf("       init() may only be called once !\n");
      exit(-1);
    }

    initialized = true;
    initCom(MPI_COMM_WORLD);
  }

  /*
  MpiCommunication(MPI_Comm myCom):
    finalize(false),
    initialized(false),
    deleteCom(true),
    reservedTagsStart(100)
  {
    initCom(myCom);
    originThreadId = omp_get_thread_num();
  }
  */
public:
  MpiCommunication(MPI_Comm com_ = MPI_COMM_WORLD) : finalize(false),
                                                     initialized(false),
                                                     deleteCom(false),
                                                     reservedTagsStart(100)
  {
    originThreadId = omp_get_thread_num();
    finalize = false;
    deleteCom = false;
    setCom(com_);
    globalComTimer = glb::timerPool->pop("MPI_globalCom");
    localComTimer = glb::timerPool->pop("MPI_localCom");
    barrierTimer = glb::timerPool->pop("MPI_barrier");
  }

  MpiCommunication(int *argc, char ***argv) : finalize(true),
                                              initialized(false),
                                              deleteCom(false),
                                              reservedTagsStart(100)
  {
    originThreadId = omp_get_thread_num();
    init(argc, argv);
    globalComTimer = glb::timerPool->pop("MPI_globalCom", false);
    localComTimer = glb::timerPool->pop("MPI_localCom", false);
    barrierTimer = glb::timerPool->pop("MPI_barrier", false);
  }

  MPI_Comm getCom() const
  {
    return com;
  }

  ~MpiCommunication()
  {
    if (deleteCom)
      MPI_Comm_free(&com);
    if (finalize)
      MPI_Finalize();
  }

  void setCom(MPI_Comm com_)
  {
    // finalize=false;
    // deleteCom=false;
#ifdef USE_OPENMP
    MPI_Query_thread(&threadSupport);
#else
    threadSupport = 0;
#endif
    initialized = true;
    initCom(com_);
  }

  /*
  MpiCommunication split(int color, int key)
  {
    MPI_Comm myCom;
    MPI_Comm_split(com,color,key,&myCom);
    return MpiCommunication(myCom);
  }
  */
  void debug(int n = -1) const
  {
    int wait = 1;
    char hostname[256];
    MPI_Barrier(MPI_COMM_WORLD);
    gethostname(hostname, sizeof(hostname));
    printf("Process %d on %s has PID: %d . ('set wait=0' to continue)\n", rank(), hostname, getpid());
    fflush(stdout);
    int nPass = 0;
    while (wait && (nPass != n))
    {
      sleep(1); // will most probably break here ...
      // communicate so that it's enough to 'set waiting = 0' on one process
      MPI_Allreduce(MPI_IN_PLACE, &wait, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      nPass++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  template <typename T, bool UseAsBarrier = false>
  std::pair<T, T> minMax(const T &val) const
  {
    T tmp[2];
    tmp[0] = val;
    tmp[1] = -val;

    MpiCallTimer timer(UseAsBarrier ? barrierTimer : globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, tmp, 2, MPI_Type<T>::get(), MPI_MIN, com);

    return std::make_pair(tmp[0], -tmp[1]);
  }

  template <typename T, bool UseAsBarrier = false>
  std::pair<std::pair<T, T>, T> minMaxSum(const T &val) const
  {
    MpiCallTimer timer(UseAsBarrier ? barrierTimer : globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Op op;
    MPI_Op_create((MPI_User_function *)MPI_Min_Max_Sum_impl<T>, 1, &op);

    T tmp[3];
    tmp[0] = val;
    tmp[1] = val;
    tmp[2] = val;

    MPI_Allreduce(MPI_IN_PLACE, tmp, 3, MPI_Type<T>::get(), op, com);
    MPI_Op_free(&op);

    return std::make_pair(std::make_pair(tmp[0], tmp[1]), tmp[2]);
  }

  template <typename T, bool UseAsBarrier = false>
  T max(const T &val) const
  {
    T tmp = val;
    MpiCallTimer timer(UseAsBarrier ? barrierTimer : globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, &tmp, 1, MPI_Type<T>::get(), MPI_MAX, com);

    return tmp;
  }

  template <typename T>
  void max(T *val, long N) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, val, N, MPI_Type<T>::get(), MPI_MAX, com);
  }

  template <typename T>
  T min(const T &val) const
  {
    T tmp = val;
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, &tmp, 1, MPI_Type<T>::get(), MPI_MIN, com);

    return tmp;
  }

  template <typename T>
  void min(T *val, long N) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, val, N, MPI_Type<T>::get(), MPI_MIN, com);
  }

  template <typename T>
  T sum(const T &val) const
  {
    T tmp = val;
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, &tmp, 1, MPI_Type<T>::get(), MPI_SUM, com);

    return tmp;
  }

  template <typename T>
  void sum(T *val, long N) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Allreduce(MPI_IN_PLACE, val, N, MPI_Type<T>::get(), MPI_SUM, com);
  }

  void barrier() const
  {
    MpiCallTimer timer(barrierTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Barrier(com);
  }

  void exit(int val = -1) const
  {
    MPI_Barrier(com);
    exit(val);
  }

  template <class Container, bool UseAsBarrier = false>
  int Bcast(Container &buffer, int root = 0, long count = -1) const
  {
    MpiCallTimer timer(UseAsBarrier ? barrierTimer : globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    // MpiCallTimer timer(globalComTimer,originThreadId);
    if (count < 0)
      count = buffer.size();
    int res = MPI_Bcast(&buffer[0], count, MPI_Type<typename Container::value_type>::get(),
                        root, com);

    return res;
  }

  template <typename T, bool UseAsBarrier = false>
  int Bcast(T *buffer, int root = 0, long count = 1) const
  {
    MpiCallTimer timer(UseAsBarrier ? barrierTimer : globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    // MpiCallTimer timer(globalComTimer,originThreadId);
    int res = MPI_Bcast(buffer, count, MPI_Type<T>::get(), root, com);

    return res;
  }

  template <typename T>
  int Gather(T *sendBuffer, long sendCount, T *rcvBuffer, long rcvCount, long root) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    return MPI_Gather(sendBuffer, sendCount, MPI_Type<T>::get(),
                      rcvBuffer, rcvCount, MPI_Type<T>::get(),
                      root, com);
  }

  /*
  template <typename T>
  int Gather_inplace(T* buffer, long count, long root)
  {
    return MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
       buffer,count,MPI_Type<T>::get(),root,com);
  }

  template <class Container>
  int Gather(Container &sendBuffer, Container &rcvBuffer, long root)
  {
    rcvBuffer.resize(sendBuffer.size()*size());
    return MPI_Gather(&sendBuffer[0],sendBuffer.size(),MPI_Type<typename Container::value_type>::get(),
          &rcvBuffer[0],sendBuffer.size(),MPI_Type<typename Container::value_type>::get(),
          root,com);
  }
  */

  // blocking send/receive
  template <typename T>
  int Recv(T *rcv, long count, long node, int tag = MPI_ANY_TAG) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Recv(rcv, count, MPI_Type<T>::get(), node, tag, com, &status);

    return res;
  }

  template <typename T>
  int Send(T *snd, long count, long node, int tag = 0) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Send(snd, count, MPI_Type<T>::get(), node, tag, com);

    return res;
  }
  /*
  template <class Container>
  int Recv(Container &buffer, long node, int tag=MPI_ANY_TAG)
  {
    MPI_Status status;
    return MPI_Recv(buffer,buffer.size(),
        MPI_Type<typename Container::value_type>::get(),node,tag,com,&status);
  }

  template <class Container>
  int Send(Container &buffer, long node, int tag=0)
  {
    return  MPI_Send(buffer,buffer.size(),
         MPI_Type<typename Container::value_type>::get(),node,tag,com);
  }
  */
  // non blocking send/receive
  template <typename T>
  int Irecv(T *rcv, long count, long node, MPI_Request *req, int tag = MPI_ANY_TAG) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    return MPI_Irecv(rcv, count, MPI_Type<T>::get(), node, tag, com, req);
  }

  template <typename T>
  int Isend(T *snd, long count, long node, MPI_Request *req, int tag = 0) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    return MPI_Isend(snd, count, MPI_Type<T>::get(), node, tag, com, req);
  }

  // Structured data send / receive
  template <typename T>
  int RecvMpiType(T *rcv, long count, MPI_Datatype dataType, long node,
                  int tag = MPI_ANY_TAG) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Recv(rcv, count, dataType, node, tag, com, &status);

    return res;
  }

  template <typename T>
  int SendMpiType(T *snd, long count, MPI_Datatype dataType, long node, int tag = 0) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Send(snd, count, dataType, node, tag, com);

    return res;
  }

  template <typename T>
  int IrecvMpiType(T *rcv, long count, MPI_Datatype dataType, long node,
                   MPI_Request *req, int tag = MPI_ANY_TAG) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Irecv(rcv, count, dataType, node, tag, com, req);
    return res;
  }

  template <typename T>
  int IsendMpiType(T *snd, long count, MPI_Datatype dataType, long node,
                   MPI_Request *req, int tag = 0) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Isend(snd, count, dataType, node, tag, com, req);
    return res;
  }
  /*
  template <class Container>
  int Irecv(Container &buffer, MPI_Datatype dataType, long node, MPI_Request *req, int tag=MPI_ANY_TAG)
  {
    return MPI_Recv(buffer,buffer.size(),dataType,node,tag,com,req);
  }

  template <class Container>
  int Isend(Container &buffer, MPI_Datatype dataType, long node, MPI_Request *req, int tag=0)
  {
    return MPI_Send(buffer,buffer.size(),dataType,node,tag,com,req);
  }
  */

  int Wait(MPI_Request *req) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Wait(req, &status);

    return res;
  }

  int Wait(MPI_Request *req, MPI_Status *status) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Wait(req, status);

    return res;
  }

  int Waitall(int count, MPI_Request *req) const
  {
    MPI_Status status[count];
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitall(count, req, status);

    return res;
  }

  int Waitall(int count, MPI_Request *req, MPI_Status *status) const
  {
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitall(count, req, status);

    return res;
  }

  int Waitall(std::vector<MPI_Request> &req) const
  {
    MPI_Status status[req.size()];
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitall(req.size(), &req[0], status);

    return res;
  }

  int Waitall(std::vector<MPI_Request> &req, std::vector<MPI_Status> &status) const
  {
    status.resize(req.size());
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitall(req.size(), &req[0], &status[0]);

    return res;
  }

  int Waitany(int count, MPI_Request *req, int *index, MPI_Status *status) const
  {
    // int index;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitany(count, req, index, status);

    return res;
  }

  int Waitany(std::vector<MPI_Request> &req, int *index, MPI_Status *status) const
  {
    // int index;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Waitany(req.size(), &req[0], index, status);

    return res;
  }

  int getCount(MPI_Status *status, MPI_Datatype &datatype) const
  {
    int count;
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    MPI_Get_count(status, datatype, &count);
    return count;
  }

  // NOTE: this does not work ...
  // bug occurs only in very obscure circumstances ...
  bool AllGatherIsBugged()
  {
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    std::vector<long> count(nProcs * nProcs, 0);
    for (long i = 0; i < nProcs * nProcs; i++)
    {
      if (myRank == (i / nProcs))
        count[i] = myRank;
      else
        count[i] = 0;
    }

    Allgather_inplace(count);
    bool bug = false;
    for (long i = 0; i < nProcs; i++)
      for (long j = 0; j < nProcs; j++)
      {
        if (count[i * nProcs + j] != i)
          bug = true;
      }
    return bug;
    //(bool)max(bug);
  }

  // NOTE :: there is a stranhge rare bug in MPI_Allgather of openMPI v1.4
  // use allreduce instead if you can and performence is not critical
  template <typename T>
  int Allgather_inplace(T *buffer, long count) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                            buffer, count, MPI_Type<T>::get(), com);

    return res;
  }

  template <class Container>
  int Allgather_inplace(Container &buffer) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    int res = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                            &buffer[0], buffer.size() / size(),
                            MPI_Type<typename Container::value_type>::get(), com);

    return res;
  }

  template <class Container>
  int Alltoall(Container &bufferIn, Container &bufferOut) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->start();
    int sz = bufferIn.size() / size();
    int res = MPI_Alltoall(&bufferIn[0], sz, MPI_Type<typename Container::value_type>::get(),
                           &bufferOut[0], sz, MPI_Type<typename Container::value_type>::get(),
                           com);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->stop();
    return res;
  }

  template <class T>
  int Alltoall(T *bufferIn, T *bufferOut, int sz) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    // if (omp_get_thread_num()==originThreadId) globalComTimer->start();
    int res = MPI_Alltoall(&bufferIn[0], sz, MPI_Type<T>::get(),
                           &bufferOut[0], sz, MPI_Type<T>::get(),
                           com);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->stop();
    return res;
  }

  template <class T>
  int Alltoallv(T *sendBuf, int *sendCount, int *sendDisp,
                MPI_Datatype sendType,
                T *receiveBuf, int *receiveCount, int *receiveDisp,
                MPI_Datatype receiveType) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    // if (omp_get_thread_num()==originThreadId) globalComTimer->start();
    int res = MPI_Alltoallv(sendBuf, sendCount, sendDisp, sendType,
                            receiveBuf, receiveCount, receiveDisp, receiveType,
                            com);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->stop();
    return res;
  }

  template <class T>
  int Alltoallv(T *sendBuf, int *sendCount, int *sendDisp,
                T *receiveBuf, int *receiveCount, int *receiveDisp) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    // if (omp_get_thread_num()==originThreadId) globalComTimer->start();
    int res = MPI_Alltoallv(sendBuf, sendCount, sendDisp, MPI_Type<T>::get(),
                            receiveBuf, receiveCount, receiveDisp, MPI_Type<T>::get(),
                            com);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->stop();
    return res;
  }

  template <class T>
  int Alltoallw(T *sendBuf, int *sendCounts, int *sendDisps,
                MPI_Datatype *sendTypes,
                T *receiveBuf, int *receiveCounts, int *receiveDisps,
                MPI_Datatype *receiveTypes) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    // if (omp_get_thread_num()==originThreadId) globalComTimer->start();
    int res = MPI_Alltoallw(sendBuf, sendCounts, sendDisps, sendTypes,
                            receiveBuf, receiveCounts, receiveDisps, receiveTypes, com);
    // if (omp_get_thread_num()==originThreadId) globalComTimer->stop();
    return res;
  }

  template <typename T>
  int Allreduce_inplace(T *buffer, long count, MPI_Op op) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    int res = MPI_Allreduce(MPI_IN_PLACE, buffer, count, MPI_Type<T>::get(), op, com);

    return res;
  }
  template <class Container>
  int Allreduce_inplace(Container &buffer, MPI_Op op) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    int res = MPI_Allreduce(MPI_IN_PLACE, &buffer[0], buffer.size(),
                            MPI_Type<typename Container::value_type>::get(), op, com);

    return res;
  }

  template <typename T>
  int Reduce_inplace(T *buffer, int root, long count, MPI_Op op) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    int res;
    if (myRank == root)
    {
      res = MPI_Reduce(MPI_IN_PLACE, buffer, count, MPI_Type<T>::get(), op, root, com);
    }
    else
    {
      res = MPI_Reduce(buffer, NULL, count, MPI_Type<T>::get(), op, root, com);
    }

    return res;
  }
  template <class Container>
  int Reduce_inplace(Container &buffer, int root, MPI_Op op) const
  {
    MpiCallTimer timer(globalComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    int res;
    if (myRank == root)
    {
      res = MPI_Reduce(MPI_IN_PLACE, &buffer[0], buffer.size(),
                       MPI_Type<typename Container::value_type>::get(), op, root, com);
    }
    else
    {
      res = MPI_Reduce(&buffer[0], NULL, buffer.size(),
                       MPI_Type<typename Container::value_type>::get(), op, root, com);
    }

    return res;
  }
  /*
  bool Iprobe(int count)
  {
    int count;
    int flag;
    MPI_Status status;
    MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, com, &flag, &status );
    if (!flag) return false;
    MPI_Get_count( &status, MPI_INT, &count );
    return true;
  }
  */

  template <class TT>
  int ProbeCount(int &count, int tag = MPI_ANY_TAG) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Probe(MPI_ANY_SOURCE, tag, com, &status);
    MPI_Get_count(&status, MPI_Type<TT>::get(), &count);

    return status.MPI_SOURCE;
  }

  int ProbeCount(int &count, MPI_Datatype &dataType, int tag = MPI_ANY_TAG) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Probe(MPI_ANY_SOURCE, tag, com, &status);
    MPI_Get_count(&status, dataType, &count);

    return status.MPI_SOURCE;
  }

  int Probe(int tag = MPI_ANY_TAG) const
  {
    MPI_Status status;
    MpiCallTimer timer(localComTimer, originThreadId);
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    MPI_Probe(MPI_ANY_SOURCE, tag, com, &status);

    return status.MPI_SOURCE;
  }
  /*
  void waitMyTurn(int nSim, int tag=maxReservedTags+1) const
  {
    int tmp=0;
    if (myRank>=nSim) Recv(&tmp,0,myRank-nSim,tag);
  }

  void finishedMyTurn(int nSim, int tag=maxReservedTags+1) const
  {
    int tmp=0;
    if (myRank+nSim < nProcs)
      Send(&tmp,0,myRank+nSim,tag);
  }
  */
  int reserveTags(int N)
  {

    if (reservedTagsStart + N > maxReservedTags)
    {
      fprintf(stderr, "ERROR in mpiCommunication::reserveTags: no more tags available (%d max).\n", maxReservedTags);
      printf("ERROR in mpiCommunication::reserveTags: no more tags available (%d max).\n", maxReservedTags);
      exit(-1);
      return -1;
    }
    int res = reservedTagsStart;
    reservedTagsStart += N;
    return res;
  }

  int abort(int errorCode = 0)
  {
    MpiOmpLockChecker lockChecker(locker, threadSupport);

    return MPI_Abort(com, errorCode);
  }

#else // we do not have MPI
private:
  // time_t refTime;
public:
  MpiCommunication(int *argc, char ***argv, MPI_Comm com_ = NULL)
  {
    myRank = 0;
    nProcs = 1;
    reservedTagsStart = 0;
    MPI_Wtime(); // This needs to be initialized ...
    // time(&refTime);
    globalComTimer = glb::timerPool->pop("MPI_globalCom");
    localComTimer = glb::timerPool->pop("MPI_localCom");
    barrierTimer = glb::timerPool->pop("MPI_barrier");
  }

  MpiCommunication(MPI_Comm com_ = NULL)
  {
    myRank = 0;
    nProcs = 1;
    reservedTagsStart = 0;
    MPI_Wtime();
    // time(&refTime);
    globalComTimer = glb::timerPool->pop("MPI_globalCom");
    localComTimer = glb::timerPool->pop("MPI_localCom");
    barrierTimer = glb::timerPool->pop("MPI_barrier");
  }

  ~MpiCommunication()
  {
  }

  void setCom(MPI_Comm com_)
  {
  }

  int Wait(MPI_Request *req) const
  {
    return 0;
  }

  int Wait(MPI_Request *req, MPI_Status *status) const
  {
    return 0;
  }

  int Waitall(int count, MPI_Request *req) const
  {
    return 0;
  }

  int Waitall(int count, MPI_Request *req, MPI_Status *status) const
  {
    return 0;
  }

  int Waitall(std::vector<MPI_Request> &req) const
  {
    return 0;
  }

  int Waitall(std::vector<MPI_Request> &req, std::vector<MPI_Status> &status) const
  {
    return 0;
  }

  int Waitany(int count, MPI_Request *req, int *index, MPI_Status *status) const
  {
    return 0;
  }

  int Waitany(std::vector<MPI_Request> &req, int *index, MPI_Status *status) const
  {
    return 0;
  }

  int getCount(MPI_Status *status, MPI_Datatype &datatype) const
  {
    return 0;
  }

  bool AllGatherIsBugged() const { return false; }

  MpiCommunication split(int color, int key) const
  {
    return MpiCommunication(0, NULL);
  }

  template <class Container>
  int Alltoall(Container &bufferIn, Container &bufferOut) const
  {
    std::copy(bufferIn.begin(), bufferIn.end(), bufferOut.begin());
    // bufferOut = bufferIn;
    return 0;
  }

  template <class T>
  int Alltoall(T *bufferIn, T *bufferOut, int sz) const
  {
    std::copy_n(bufferIn, sz, bufferOut);
    return 0;
  }

  template <class T>
  int Alltoallv(T *sendBuf, int *sendCount, int *sendDisp,
                MPI_Datatype sendType,
                T *receiveBuf, int *receiveCount, int *receiveDisp,
                MPI_Datatype receiveType) const
  {
    // memcpy(receiveBuf+receiveDisp[0],sendBuf+sendDisp[0],sizeof(T)*(*sendCount));
    return 0;
  }

  template <class T>
  int Alltoallv(T *sendBuf, int *sendCount, int *sendDisp,
                T *receiveBuf, int *receiveCount, int *receiveDisp) const
  {
    // memcpy(receiveBuf+receiveDisp[0],sendBuf+sendDisp[0],sizeof(T)*(*sendCount));
    return 0;
  }

  template <class T>
  int Alltoallw(T *sendBuf, int *sendCounts, int *sendDisps,
                MPI_Datatype *sendTypes,
                T *receiveBuf, int *receiveCounts, int *receiveDisps,
                MPI_Datatype *receiveTypes) const
  {
    // memcpy(receiveBuf+receiveDisps[0],sendBuf+sendDisps[0],sizeof(T)*(*sendCounts));
    return 0;
  }

  void debug(int n = -1) const {}

  template <typename T, bool UseAsBarrier = false>
  std::pair<T, T> minMax(const T &val) const
  {
    return std::make_pair(val, val);
  }

  template <typename T, bool UseAsBarrier = false>
  std::pair<std::pair<T, T>, T> minMaxSum(const T &val) const
  {
    return std::make_pair(std::make_pair(val, val), val);
  }

  template <typename T, bool UseAsBarrier = false>
  T max(const T &val) const { return val; }
  template <typename T>
  void max(T *val, long N) const {}
  template <typename T>
  T min(const T &val) const { return val; }
  template <typename T>
  void min(T *val, long N) const {}
  template <typename T>
  T sum(const T &val) const { return val; }
  template <typename T>
  void sum(T *val, long N) const {}

  void barrier() const {}
  void exit(int val = -1) const { ::exit(val); }

  MPI_Comm getCom() const
  {
    return NULL;
  }

  template <typename T>
  int Recv(T *rcv, long count, long node, int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  template <typename T>
  int Send(T *snd, long count, long node, int tag = 0) const
  {
    return 0;
  }

  template <typename T>
  int Allgather_inplace(T *buffer, long count) const
  {
    return 0;
  }
  template <class Container>
  int Allgather_inplace(Container &buffer) const
  {
    return 0;
  }

  template <typename T>
  int Allreduce_inplace(T *buffer, long count, MPI_Op op) const
  {
    return 0;
  }
  template <class Container>
  int Allreduce_inplace(Container &buffer, MPI_Op op) const
  {
    return 0;
  }

  template <class Container, bool UseAsBarrier = false>
  int Bcast(Container &buffer, int root = 0, long count = -1) const
  {
    return 0;
  }

  template <typename T, bool UseAsBarrier = false>
  int Bcast(T *buffer, int root = 0, long count = 1) const
  {
    return 0;
  }

  template <typename T>
  int Gather(T *sendBuffer, long sendCount, T *rcvBuffer, long rcvCount, long root)
  {
    std::copy(sendBuffer, sendBuffer + sendCount, rcvBuffer);
    return 0;
  }

  template <typename T>
  int Irecv(T *rcv, long count, long node, MPI_Request *req, int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  template <typename T>
  int Isend(T *snd, long count, long node, MPI_Request *req, int tag = 0) const
  {
    return 0;
  }

  template <typename T>
  int Irecv(T *rcv, long count, MPI_Datatype dataType, long node,
            MPI_Request *req, int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  template <typename T>
  int Isend(T *snd, long count, MPI_Datatype dataType, long node,
            MPI_Request *req, int tag = 0) const
  {
    return 0;
  }

  template <typename T>
  int RecvMpiType(T *rcv, long count, MPI_Datatype dataType,
                  long node, int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  template <typename T>
  int SendMpiType(T *snd, long count, MPI_Datatype dataType, long node, int tag = 0) const
  {
    return 0;
  }

  template <typename T>
  int IrecvMpiType(T *rcv, long count, MPI_Datatype dataType, long node,
                   MPI_Request *req, int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  template <typename T>
  int IsendMpiType(T *snd, long count, MPI_Datatype dataType, long node,
                   MPI_Request *req, int tag = 0) const
  {
    return 0;
  }

  template <typename T>
  int Reduce_inplace(T *buffer, int root, long count, MPI_Op op) const
  {
    return 0;
  }

  template <class Container>
  int Reduce_inplace(Container &buffer, int root, MPI_Op op) const
  {
    return 0;
  }

  template <class TT>
  int ProbeCount(int &count, int tag = MPI_ANY_TAG) const
  {
    count = 0;
    return 0;
  }

  int ProbeCount(int &count, MPI_Datatype &dataType, int tag = MPI_ANY_TAG) const
  {
    count = 0;
    return 0;
  }

  int Probe(int tag = MPI_ANY_TAG) const
  {
    return 0;
  }

  void waitMyTurn(int nSim, int tag = maxReservedTags + 1) const
  {
  }

  void finishedMyTurn(int nSim, int tag = maxReservedTags + 1)
  {
  }

  int reserveTags(int N)
  {
    if (reservedTagsStart + N > maxReservedTags)
      return -1;
    int res = reservedTagsStart;
    reservedTagsStart += N;
    return res;
  }

  int abort(int errorCode = 0)
  {
    exit(errorCode);
    return errorCode;
  }
  //};

#endif // HAVE_MPI
       /*
         static double Wtime()
         {
           return MPI_Wtime();
         }
       */
  double reportTimeSpent(double &globalCom,
                         double &localCom,
                         double &barrierCom) const
  {
    globalCom = globalComTimer->totalSpent();
    localCom = localComTimer->totalSpent();
    barrierCom = barrierTimer->totalSpent();

    return globalCom + localCom + barrierCom;
  }

  void reportTimeSpent(double globalCom[3],
                       double localCom[3],
                       double barrierCom[3],
                       double total[3], int rk[2]) const

  {
    double arr[3][4];
    arr[0][0] = globalComTimer->totalSpent();
    arr[0][1] = localComTimer->totalSpent();
    arr[0][2] = barrierTimer->totalSpent();
    arr[0][3] = arr[0][0] + arr[0][1] + arr[0][2];
    double localTotalSpent = arr[0][3];

    for (int i = 0; i < 4; i++)
      arr[1][i] = arr[0][i];
    for (int i = 0; i < 4; i++)
      arr[2][i] = arr[0][i] / nProcs;

    min(arr[0], 4);
    max(arr[1], 4);
    sum(arr[2], 4);

    for (int i = 0; i < 3; i++)
    {
      globalCom[i] = arr[i][0];
      localCom[i] = arr[i][1];
      barrierCom[i] = arr[i][2];
      total[i] = arr[i][3];
    }

    if (total[0] == localTotalSpent)
      rk[0] = rank();
    else
      rk[0] = -1;
    if (total[1] == localTotalSpent)
      rk[1] = rank();
    else
      rk[1] = -1;

    max(rk, 2);
    /*
    Allreduce_inplace(arr,8,MPI_MIN);
    Allreduce_inplace(&arr[8],4,MPI_SUM);

    globalCom[0]=arr[0];globalCom[1]=-arr[4];globalCom[2]=arr[8]/nProcs;
    localCom[0]=arr[1];localCom[1]=-arr[5];localCom[2]=arr[9]/nProcs;
    barrierCom[0]=arr[2];barrierCom[1]=-arr[6];barrierCom[2]=arr[10]/nProcs;

    for (int i=0;i<3;i++)
      total[i]=arr[3+4*i];
    total[1]=-total[1];

    if (total[0] == localTotalSpent) rk[0]=rank();
    else rk[0]=-1;
    if (total[1] == localTotalSpent) rk[1]=rank();
    else rk[1]=-1;

    max(rk,2);
    */
  }

  int rank() const { return myRank; }
  int size() const { return nProcs; }

private:
  template <typename T>
  static void MPI_Min_Max_Sum_impl(T *invec, T *inoutvec, int *len, MPI_Datatype *datatype)
  {
    for (int i = 0; i < (*len); i += 3)
    {
      if (invec[i] < inoutvec[i])
        inoutvec[i] = invec[i];
      if (invec[i + 1] > inoutvec[i + 1])
        inoutvec[i + 1] = invec[i + 1];
      inoutvec[i + 2] += invec[i + 2];
    }
  }
};

#include "../../internal/namespace.footer"
#endif

#ifndef __TIMER_POOL_HXX__
#define __TIMER_POOL_HXX__

#include <vector>
#include <list>
#include <set>
#include <map>
#include <sstream>

#include "../../dice_globals.hxx"

#include "../../tools/MPI/myMpi.hxx"

/**
 * @file
 * @brief  Implementation of a timer pool, used to pop and manage timers
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 *   \{
 */
/**
 * \brief  A pool of timers that can be used to measure performences of the code
 */
class TimerPool
{
public:
  TimerPool()
  {
  }

  ~TimerPool()
  {
  }

  /** \brief A timer that can measure lapses of time and store them
   */
  class Timer
  {
  public:
    /** \brief constructor
     * \param name_ The name of the timer
     * \param keepLog If true (default), keep a log of every time the timer was stopped and
     * restarted. If false, only the cumulated time is saved.
     */
    Timer(const std::string name_ = "dummy", bool keepLog = true) : name(name_),
                                                                    state(false)
    {
      refTime = MPI_Wtime();
      keep = keepLog;
      cumulated = 0;
      if (!keep)
      {
        record.resize(2);
        record[0].first = refTime;
        record[0].second = 0;
        record[1].first = refTime;
        record[1].second = refTime;
      }
    }

    ~Timer()
    {
    }

    /** \brief start measuring time
     * \return The current time
     */
    double start()
    {
      if (state != false)
      {
        printf("WARNING: Timer '%s' desynchronized : calling start() before stopping !\n",
               name.c_str());
        fprintf(stderr,
                "WARNING: Timer '%s' desynchronized : calling start() before stopping !\n",
                name.c_str());
      }
      if (keep)
        record.push_back(std::pair<double, double>(0, 0));

      record.back().first = MPI_Wtime();
      state = true;
      return record.back().first - refTime;
    }

    /** \brief stop measuring time
     * \return The time spent since start() was last called
     */
    double stop()
    {
      // #ifdef USE_OPENMP
      //       if (omp_get_thread_num()!=0) return -1;
      // #endif
      if (state != true)
      {
        printf("WARNING: Timer '%s' desynchronized : calling stop() before starting !\n",
               name.c_str());
        fprintf(stderr,
                "WARNING: Timer '%s' desynchronized : calling stop() before starting !\n",
                name.c_str());
      }
      record.back().second = MPI_Wtime();
      state = false;
      // if (!keep) record[0].second += record.back().second-record.back().first;
      cumulated += record.back().second - record.back().first;
      return record.back().second - record.back().first;
    }

    /** \brief Check time spent since last call to start(), but do not stop the timer and
     * do not record anything.
     * \return time spent since last call to start()
     */
    double check()
    {
      if (state == true)
        return MPI_Wtime() - record.back().first;

      return 0;
    }

    /** \brief returns the time spent between the last calls to start() and stop()
     */
    double lastSpent()
    {
      if (record.size() == 0)
        return 0;
      return record.back().second - record.back().first;
    }

    /** \brief returns the total cumlated time spent
     */
    double totalSpent()
    {
      // #ifdef USE_OPENMP
      //       if (omp_get_thread_num()!=0) return -1;
      // #endif
      double result = 0;
      return cumulated + check();
      /*
      if (keep)
  {
    for (unsigned long i=0;i<record.size();++i)
      result += record[i].second-record[i].first;
  }
      else result = record[0].second;
      */
      return result;
    }

    const std::string &getName() const { return name; }

  private:
    double cumulated;
    bool keep;
    double refTime;
    std::vector<std::pair<double, double>> record;
    const std::string name;
    bool state;
  };

  /** \brief Creates a new timer and returns a pointer to it.
   * \param name The name of the timer
   * \param keepLog If true (default), keep a log of every time the timer was stopped and
   * \return a pointer to the timer. If a timer with the same name exists, it is returned,
   * else anew one is allocated.
   */
  Timer *pop(const std::string name, bool keepLog = true)
  {
    // #ifdef USE_OPENMP
    //     if (omp_get_thread_num()!=0) return NULL;
    // #endif
    timersItMap_iterator it = timersItMap.find(name);

    if (it == timersItMap.end())
    {
      timersList.push_back(Timer(name, keepLog));
      it = timersItMap.insert(std::make_pair(name, &timersList.back())).first;
    }

    return it->second;
    /*
    return &((*it).second);
      return &(*timersSet.insert(Timer(name,keepLog)).first);
    */

    /*
  timersList.push_back(Timer(name,keepLog));
  return &timersList.back();
    */
  }

  /** \brief return a string representing the cumulated total time spent within each timer
   */
  std::string getTimings()
  {
    static std::ostringstream oss;

    oss.str(std::string(""));
    oss.clear();

    if (timersList.size() != 0)
    {
      std::list<Timer>::iterator it = timersList.begin();
      oss << (*it).totalSpent();
      for (++it; it != timersList.end(); ++it)
      {
        oss << " " << (*it).totalSpent();
      }
    }

    return oss.str();
  }

  /** \brief return a string representing the poped timers names
   */
  std::string getTimersName()
  {
    static std::ostringstream oss;

    oss.str(std::string("#"));
    oss.clear();

    if (timersList.size() != 0)
    {
      std::list<Timer>::iterator it = timersList.begin();
      oss << (*it).getName();
      for (++it; it != timersList.end(); ++it)
      {
        oss << " " << (*it).getName();
      }
    }

    return oss.str();
  }

private:
  typedef std::map<std::string, Timer *>::iterator timersItMap_iterator;
  std::map<std::string, Timer *> timersItMap;
  std::list<Timer> timersList;
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

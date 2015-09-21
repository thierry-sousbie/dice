#ifndef __EVENT_MANAGER_HXX__
#define __EVENT_MANAGER_HXX__

#include <set>

#include "../strings/stringTokenizer.hxx"

#include "../../dice_globals.hxx"

#include "../../internal/namespace.header"

class EventManager
{
public:
  EventManager()
  {}

  explicit EventManager(std::string &cmd)
  {
    addEvents(cmd);
  }

  EventManager(double start, double every)
  {
    addEvents(start,every);
  }

  EventManager(double start, double stop, double every)
  {
    addEvents(start,stop,every);
  }

  // If cmd starts with a ":", a digit ot '-', it may be of the form "start:stop:every", "start::every" or "value",
  // else it is interpreted as a filename containing a list of events timings
  // Different commands can be combined, separating them with any of " ={}[],;\t"
  void addEvents(std::string &cmd)
  {
    std::vector<std::string> tokens;    
    StringTokenizer::split(tokens,cmd.c_str(),std::string(" ={}[],;\t"));
    
    for (unsigned long i=0;i<tokens.size();++i)
      {
	//printf("Token[%ld] = '%s'\n",i,tokens[i].c_str());

	if ((tokens[i].at(0)==':')||(std::isdigit(tokens[i].at(0)))||(tokens[i].at(0)=='-'))
	  {
	    std::vector<std::string> t;
	    std::vector<double> val;

	    StringTokenizer::split(t,tokens[i].c_str(),std::string(":"));
	    val.resize(t.size());
	    for (unsigned long j=0;j<t.size();++j)
	      StringTokenizer::from_string(val[j],t[j]);
	     
	    if (t.size()==1) addSingleEvent(val[0]);
	    else if (t.size()==2) addEvents(val[0],val[1]);
	    else if (t.size()==3) addEvents(val[0],val[1],val[2]);
	    else 
	      {
		PRINT_SRC_INFO(LOG_ERROR);
		glb::console->print<LOG_ERROR>("Processing string '%s'\n",cmd.c_str());
		glb::console->print<LOG_ERROR>("Unknown token '%s'\n",tokens[i].c_str());
		exit(-1);
	      }
	  }
	else addEventsFromFile(tokens[i]);	
      }
  }
    
  void addEventsFromFile(std::string &fname)
  {
    PRINT_SRC_INFO(LOG_ERROR);
    glb::console->print<LOG_ERROR>("Cannot add events from file '%s'.\n",fname.c_str());
    glb::console->print<LOG_ERROR>("Events from files NOT IMPLEMENTED YET\n");
    exit(-1);
  }

  void addSingleEvent(double at)
  {
    events.insert(at);
  }

  void addEvents(double start, double every)
  {
    addEvents(start,time_max(),every);
  }

 
  void addEvents(double start, double stop, double every)
  {
    //printf("Adding event : [%g : %g : %g]\n",start,stop,every);
    eventTriplets.insert(Triplet(start,stop,every));
  }

  void reset()
  {
    events.clear();
  }

  bool checkEventTimeFrame(double t0, double t1, double epsilon=0) const
  {
    std::set<double>::iterator sit = 
      hlp::findFirstHigherOrEqual(events.begin(),events.end(),t0);
    std::set<double>::iterator sit2 = 
      hlp::findFirstHigherOrEqual(events.begin(),events.end(),t1);

    if (sit!=sit2) return true;

    std::set<Triplet>::iterator it;
    for (it=eventTriplets.begin(); it!=eventTriplets.end(); ++it)
      {
	if (t1<=(*it).start) continue;

	long count0 = static_cast<long>((t0-(*it).start)/(*it).every);
	long count1 = static_cast<long>((t1-(*it).start)/(*it).every);	

	double time0 = (*it).start + ((*it).every)*(count0);
	double time1 = (*it).start + ((*it).every)*(count1);
	
	if (time0 != t0) time0 = (*it).start + ((*it).every)*(count0+1);
	if (time1 != t1) time1 = (*it).start + ((*it).every)*(count1+1);
	
	if (time1!=time0) 
	  {
	    if (time0<=(*it).stop) 
	      return true;
	  }
      }

    return false;
  }

  double getNextEventTime(double curTime, double epsilon=0) const
  {
    double result;
    
    std::set<double>::iterator sit = 
      hlp::findFirstHigher(events.begin(),events.end(),curTime+epsilon);

    if (sit!=events.end())
      result=(*sit);
    else result = time_max();

    std::set<Triplet>::iterator it;
    for (it=eventTriplets.begin(); it!=eventTriplets.end(); ++it)
      {
	long count = static_cast<long>((curTime-(*it).start)/(*it).every);
	double time = (*it).start + ((*it).every)*(count+1);
	
	if ((time-curTime)<=epsilon)
	  time=(*it).start + ((*it).every)*(count+2);
		
	if (time<=(*it).stop)
	  result = std::min(result,time);
	else if ((*it).stop>curTime)	
	  result = std::min(result,(*it).stop);
      }

    return result;
  }

  bool checkEvent(double curTime, double epsilon=0) const
  {    
    std::set<double>::iterator sit = 
      hlp::findFirstHigher(events.begin(),events.end(),curTime);
    std::set<double>::iterator b_sit=events.end();
    //double result;

    if ((sit!=events.begin())&&(events.size()>0))
      {
	b_sit=sit;
	b_sit--;
      }

    if (sit!=events.end())
      if (fabs(*sit-curTime) <= epsilon) return true;

    if (b_sit!=events.end())
      if (fabs(*b_sit-curTime) <= epsilon) return true;

    std::set<Triplet>::iterator it;
    for (it=eventTriplets.begin(); it!=eventTriplets.end(); ++it)
      {
	long count = static_cast<long>((curTime-(*it).start)/(*it).every);

	double time = (*it).start + ((*it).every)*count;
	time = std::min(time,(*it).stop);
	if (fabs(time-curTime)<=epsilon) return true;
	
	time=(*it).start + ((*it).every)*(count+1);
	time = std::min(time,(*it).stop);
	if (fabs(time-curTime)<=epsilon) return true;	
      }

    return false;
  }

  
  /*
  double setEpsilon(double e)
  {
    epsilon=e;
  }
  */
private:
  /*
  static double getEpsilon(double t)
  {
    return fabs(t*20.0L*std::numeric_limits<double>::epsilon());
  }
  */

  static double time_min()
  {
    return -std::numeric_limits<double>::max();
  }

  static double time_max()
  {
    return std::numeric_limits<double>::max();
  }

  struct Triplet
  {
    Triplet(double f, double s, double t)
    {
      start=f;
      stop=s;
      every=t;
    }
    
    double start;
    double stop;
    double every; 
    
    bool operator<(const Triplet &rhs) const
    {
      if (every<rhs.every) return true;
      else if (every==rhs.every)
	{
	  if (stop<rhs.stop) return true;
	  else if (stop==rhs.stop)
	    {
	      if (start<rhs.start) return true;
	    }
	}	
      return false;
    }
  };

  std::set<double> events;  
  std::set<Triplet> eventTriplets;
  //double epsilon;
};

#include "../../internal/namespace.footer"
#endif

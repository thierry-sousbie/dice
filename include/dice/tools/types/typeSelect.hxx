#ifndef __TYPE_SELECT_HXX__
#define __TYPE_SELECT_HXX__

#include <map>
#include <string>

#include "../../dice_globals.hxx"


#include "../../internal/namespace.header"

template <class S>
struct TypeSelectT {
private:
  typedef typename S::Type T;
  typedef typename std::map<std::string,T> Map;
  typedef typename Map::iterator MapIt;
  typedef typename std::map<T,std::string> MapInv;
  typedef typename MapInv::iterator MapInvIt;
  
  Map m_;
  MapInv minv_;

public:
  TypeSelectT() {}
  virtual ~TypeSelectT() {}

  void insert(const std::string &str,const T &val)
  {
    m_.insert(std::make_pair(str,val));
    minv_.insert(std::make_pair(val,str));
  }

  T getVal(const std::string &str, bool errorIfUndefined=true) 
  { 
    MapIt it = m_.find(str);
    if (it==m_.end()) 
      {
	if (errorIfUndefined)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("unknown %s value: '%s'\n",name().c_str(),str.c_str());
	    exit(-1);	    
	  }
	return S::UNDEFINED;
      }
    else return it->second;
  }

  std::string getString(const T &val, bool errorIfUndefined=true) 
  { 
    MapInvIt it = minv_.find(val);
    if (it==minv_.end())
      {
	if (errorIfUndefined)
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("No string associated to type '%s'\n",name().c_str());
	    exit(-1);
	  }
	return "undefined";
      }
    else return it->second;
  }

  std::string getAllString(const char *format,const char* sep=",")
  {
    std::string result;
    std::string sep_str(sep);
    MapIt it=m_.begin();
    result = (*it).first;
    for (++it;it!=m_.end();++it)
      result += sep_str+(*it).first;
    
    char tmp[512];
    sprintf(tmp,format,result.c_str());

    return std::string(tmp);
  }

  virtual std::string name()=0;

};

#include "../../internal/namespace.footer"
#endif

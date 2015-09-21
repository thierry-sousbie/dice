#ifndef __PARAMS_MANAGER_HXX__
#define __PARAMS_MANAGER_HXX__

#include <string>
#include <sstream> 
#include <cstddef>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>

#include <stdlib.h>
#include <stdio.h>

#include "../helpers/helpers.hxx"

/**
 * @file 
 * @brief  Implementation of a parameter manager (used for parameters IO)
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 *   \{
 */
/**
 * \class ParamsManagerT
 * \brief  A parameter manager designed to be coupled with a parameter parser and used to
 * serialize / unserialize / manage command line (or paramter file) parameters. 
 * This class is especially designed to make it easy to assign value of parameters that 
 * that can be both set from command line and whose value was stored from a previous run.
  *
 * \tparam P a parameter parser (typically a ParamsParser)
 */
template <class P>//, class C, class LOG>
class ParamsManagerT {  
  //typedef C Console;

public:
  typedef P Parser;

  /** \brief Define respective priorities of assigned parameter values between parser and
   *  unserialized file */
  enum Priority {FILE_FIRST=0, PARSER_FIRST=1, IGNORE_FILE=2};

  ParamsManagerT(Parser *p)://, Console *c):
    parser(p)
  {
    //console = c;
  }
  
  ~ParamsManagerT()
  {
    
  }
 
  /** \brief Get a multivalued parameter value from the parameter parser or a 
   *  serialized object (through reader). 
   *  Every time a parameter is read, we store a structure with its 
   *  information so that it can be serialized automatically. Priority defines whether 
   *  a value read from the serialized object should override the parsed value or not.
   *
   *  \param what The parameter name
   *  \param cat The category name (parameter format is cat.what)
   *  \param def The default value
   *  \param index The index of the this parameter value, >0 when the parameter takes
   *  more than one value
   *  \param reader The binary reader from where the parameters are unserialized (may
   *  be NULL)
   *  \param p The priority of the assigned value compared to the unserialized one 
   *  (see Priority enum)  
   *  \param comment a comment describing the parameter
   *  \param useReader if false, skip reading the parameter from "reader". This is usefull
   *  when one wants to manage a parameter that is not present in older versions of the 
   *  serialized file.
   *  \tparam T The type of the parameter to read
   *  \tparam R The class of the binary reader (typically a BinaryReaderT)
   */  
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def, int index,
	R *reader, Priority p, const std::string &comment, bool useReader=true)//, bool store=true)
  {       
    bool store=true;
    T result=def;
    Parameter prm;
     
    switch (p) {

    case FILE_FIRST:
      if ((reader!=NULL)&&(store))
	{
	  //reader->read(&result);
	  if (useReader) 
	    result = prm.template read<T>(reader,what,cat);
	  
	  T tmp=parser->template get<T>(what,cat,result,comment,index);
	  
	  if (!useReader) result=tmp;
	  
	  if (tmp != result)
	    {
	      std::ostringstream oss;
	      oss << result;
	      glb::console->
		template print<LOG_STD>("Parsed value overridden from file: %s.%s = %s\n",
				    cat.c_str(),what.c_str(),oss.str().c_str());
	    }
	  break;
	}

    case PARSER_FIRST:
      if ((reader!=NULL)&&(store))
	{
	  if (useReader) result = prm.template read<T>(reader,what,cat); //reader->read(&result);
	}
      result = parser->template get<T>(what,cat,result,comment,index);
      break;

    case IGNORE_FILE:    
      if ((reader!=NULL)&&(store)) 
	if (useReader) prm.template read<T>(reader,what,cat);
      
      result = parser->template get<T>(what,cat,result,comment,index);
      break;
    }

    if (store) allParameters.push_back(Parameter(what,cat,result));    
    
    return result;
  }

  template <typename T, class R>
  void discard(const std::string &what, const std::string &cat, R *reader, 
	       bool condition=true)
  {
    if (condition)
      {
	Parameter prm;
	prm.template read<T>(reader,what,cat);
      }
  }

  /** \brief Get a parameter value from the parameter parser or a serialized object 
   *  (through reader). Every time a parameter is read, we store a structure with its 
   *  information so that it can be serialized automatically. Priority defines whether 
   *  a value read from the serialized object should override the parsed value or not.
   *
   *  \param what The parameter name
   *  \param cat The category name (parameter format is cat.what)
   *  \param def The default value   
   *  \param reader The binary reader from where the parameters are unserialized (may
   *  be NULL)
   *  \param p The priority of the assigned value compared to the unserialized one 
   *  (see Priority enum)
   *  \param useReader if false, skip reading the parameter from "reader". This is usefull
   *  when one wants to manage a parameter that is not present in older versions of the 
   *  serialized file.
   *  \tparam T The type of the parameter to read
   *  \tparam R The class of the binary reader (typically a BinaryReaderT)
   */  
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def,
	R *reader, Priority p, bool useReader=true)//, bool store=true)
  {       
    return get<T,R>(what,cat,def,-1,reader,p,"No description available",useReader);//,store);
  }

  /** \brief Get a parameter value from the parameter parser or a serialized object 
   *  (through reader). Every time a parameter is read, we store a structure with its 
   *  information so that it can be serialized automatically. Priority defines whether 
   *  a value read from the serialized object should override the parsed value or not.
   *
   *  \param what The parameter name
   *  \param cat The category name (parameter format is cat.what)
   *  \param def The default value   
   *  \param index The index of the this parameter value, >0 when the parameter takes
   *  more than one value
   *  \param reader The binary reader from where the parameters are unserialized (may
   *  be NULL)
   *  \param p The priority of the assigned value compared to the unserialized one 
   *  (see Priority enum)
   *  \param useReader if false, skip reading the parameter from "reader". This is usefull
   *  when one wants to manage a parameter that is not present in older versions of the 
   *  serialized file.
   *  \tparam T The type of the parameter to read
   *  \tparam R The class of the binary reader (typically a BinaryReaderT)
   */  
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def, int index,
	R *reader, Priority p, bool useReader=true)//, bool store=true)
  {       
    return get<T,R>(what,cat,def,index,reader,p,"No description available",useReader);//,store);
  }

  /** \brief Get a multivalued parameter value from the parameter parser or a 
   *  serialized object (through reader). 
   *  Every time a parameter is read, we store a structure with its 
   *  information so that it can be serialized automatically. Priority defines whether 
   *  a value read from the serialized object should override the parsed value or not.
   *
   *  \param what The parameter name
   *  \param cat The category name (parameter format is cat.what)
   *  \param def The default value
   *  \param reader The binary reader from where the parameters are unserialized (may
   *  be NULL)
   *  \param p The priority of the assigned value compared to the unserialized one 
   *  (see Priority enum)  
   *  \param comment a comment describing the parameter
   *  \param useReader if false, skip reading the parameter from "reader". This is usefull
   *  when one wants to manage a parameter that is not present in older versions of the 
   *  serialized file.
   *  \tparam T The type of the parameter to read
   *  \tparam R The class of the binary reader (typically a BinaryReaderT)
   */  
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def, 
	R *reader, Priority p, const std::string &comment, bool useReader=true)//, bool store=true)
  {       
    return get<T,R>(what,cat,def,-1,reader,p,comment,useReader);//,store);
  }
  
  /** \brief Get a multivalued parameter value from the parameter parser or a 
   *  serialized object (through reader).  See other versions for more info
   */
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def, 
	R *reader, Priority p, const char *comment, bool useReader=true)//, bool store=true)
  {       
    return get<T,R>(what,cat,def,-1,reader,p,std::string(comment),useReader);//,store);
  }

  /** \brief Get a multivalued parameter value from the parameter parser or a 
   *  serialized object (through reader).  See other versions for more info
   */
  template <typename T, class R>
  T get(const std::string &what,const std::string &cat, const T &def, int index,
	R *reader, Priority p, const char *comment, bool useReader=true)
  {
    return get<T,R>(what,cat,def,index,reader,p,std::string(comment),useReader);//,store);
  }

  /** \brief Serialize managed parameters that were parsed through get member
   *  \tparam W The class of the binary writer (typically a BinaryWriterT)  
   */
  template <class W>
  void write(W *writer) const
  {
    for (unsigned long i=0;i<allParameters.size();++i)
      {
	allParameters[i].write(writer);
	/*
	const Parameter &p = allParameters[i];
	if (p.size<0)
	  writer->write(&p.str);
	else
	  writer->write(p.value,p.size);
	*/
      }
  }
  
  /** \brief Print a report to console stating the values of parsed parameters, including
   *  the non stored ones.
   *  \param defaultOnly If true (default), print only the parameters that were assigned
   *  the default value
   *  \tparam LT the Log level of the console (typically LOG_STD)
   */
  template <class LT>
  void reportParsed(bool defaultOnly=true)
  {     
    parser->template report<LT>(defaultOnly);    
  }

  /** \brief Print a report to console stating the values of managed parameters
   *  \tparam LT the Log level of the console (typically LOG_STD)
   */
  template <class LT>
  void report() const
  {        
    if (allParameters.size()==0) 
      {
	glb::console->print<LT>("Parameter manager: nothing to report !\n");
	return;
      }

    std::vector<Parameter> all(allParameters.begin(),allParameters.end());
    std::stable_sort(all.begin(),all.end());

    bool first=true;
    std::string curCat = std::string();
    std::string curWhat = std::string();

    for (unsigned long i=0;i<all.size();++i)
      {
	const Parameter &p = all[i];
	if (curCat != p.cat)
	  {
	    if (!first) {glb::console->print<LT>("\n");}
	    else first=false;
	    glb::console->print<LT>(" Category: '%s'\n",p.cat.c_str());
	    glb::console->print<LT>("   %s =",p.what.c_str());
	    curCat=p.cat;
	  }
	else if (curWhat != p.what)
	  {
	    glb::console->print<LT>("\n");
	    glb::console->print<LT>("   %s =",p.what.c_str());
	    curWhat = p.what;
	  }
	  
	glb::console->print<LT>(" %s",p.str.c_str());
      }
    glb::console->print<LT>("\n");
  }

private:

  // Store information about the parameters
  struct Parameter
  {
    Parameter() 
    {}

    template <class T>
    Parameter(const std::string &w, const std::string &c,const T &v)
    {
      init(w,c,v,typename hlp::SameType<T,std::string>::Result());      
    }

    // Type is a string 
    template <class T>
    void init(const std::string &w, const std::string &c,const T &v, hlp::IsTrue)
    {
      what=w;
      cat=c;
      str=v;
      size=-1;
    }

    // Type is a base type (i.e. int, float, ...)
    template <class T>
    void init(const std::string &w, const std::string &c,const T &v, hlp::IsFalse)
    {
      what=w;
      cat=c;
      std::ostringstream oss;
      oss << v;
      str = oss.str();
      memcpy(value,&v,sizeof(T));
      size=sizeof(T);
    }

    std::string cat;
    std::string what;

    int size;
    std::string str;    
    char value[8];
    
    bool operator<(const Parameter &other) const 
    {
      if (cat<other.cat) return true;
      if (cat == other.cat)
	{
	  if (what<other.what) return true;
	  return false;
	}
      return false;	
    }

    template <class W>
    void write(W *writer) const
    {
      writer->write(&cat);
      writer->write(&what);
      writer->write(&size);
      
      if (size<0) writer->write(&str);
      else writer->write(value,size);
    }

    template <class T, class R>
    T read(R *reader, const std::string &tWhat, const std::string &tCat)
    {
      T result;
      int typeSize = int((hlp::SameType<T,std::string>::Result::value)?-1:sizeof(T));

      reader->read(&cat);
      reader->read(&what);
      reader->read(&size);

      if ((cat != tCat)||(what!=tWhat))
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Reading from file: parameter name mismatch.\n");
	  glb::console->print<LOG_ERROR>("From file : %s.%s (size %d).\n",cat.c_str(),what.c_str(),size);
	  glb::console->print<LOG_ERROR>("Required  : %s.%s (size %d).\n",tCat.c_str(),tWhat.c_str(),typeSize);
	  exit(-1);//return T();
	}
      else if (typeSize != size)
	{
	  PRINT_SRC_INFO(LOG_ERROR);
	  glb::console->print<LOG_ERROR>("Parameters sizes do not match.\n");	  
	  glb::console->print<LOG_ERROR>("From file : %s.%s (size %d).\n",cat.c_str(),what.c_str(),size);
	  glb::console->print<LOG_ERROR>("Required  : %s.%s (size %d).\n",tCat.c_str(),tWhat.c_str(),typeSize);
	  exit(-1);//return T();
	}

      reader->read(&result);

      if (hlp::SameType<T,std::string>::Result::value)
	str = result;
      else
	memcpy(value,&result,size);      
           
      return result;
    }   
    

  };
  
  std::vector< Parameter > allParameters;

public:
  /** \brief returns a pointer to the parameter parser. This is usefull to parse parameters
   * without actually recording them within the manager.
   */
  Parser *getParser()
  {
    return parser;
  }

private:
  Parser *parser;
  //Console *console;
};

/** \}*/
#include "../../internal/namespace.footer"
#endif

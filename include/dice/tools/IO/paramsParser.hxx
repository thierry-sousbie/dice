#ifndef __PARAMS_PARSER_HXX__
#define __PARAMS_PARSER_HXX__

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


#include "../../dice_globals.hxx"

#include "../strings/stringTokenizer.hxx"

/**
 * @file 
 * @brief  Implementation of a parameter parser (command line and parameter file reader)
 * @author Thierry Sousbie
 */

#include "../../internal/namespace.header"

/** \addtogroup TOOLS
 *   \{
 */

class ParamsParser {
public:
  static const std::string defaultCategory() {return "default";}

private:
 
  struct paramValueT {
    paramValueT() {}
    std::string comment;
    std::vector< std::string > val;
    std::vector< std::string > fname;
    std::vector< int > line;
    std::vector< int > status; // 0 => not used, 1=> ignored, 2=> used
    //paramValue():wasRead(false){}
  };

  typedef std::map< std::string, paramValueT > paramsT;
  typedef paramsT::iterator paramsT_it;
  typedef paramsT::const_iterator paramsT_cit;

  typedef std::list< paramsT > paramsListT;
  typedef paramsListT::iterator paramsListT_it;
  typedef paramsListT::const_iterator paramsListT_cit;

  typedef std::map<std::string,paramsListT_it> categoryT;
  typedef std::map<std::string,paramsListT_it>::iterator categoryT_it;
  typedef std::map<std::string,paramsListT_it>::const_iterator categoryT_cit;

  std::string execFileName;
  std::string paramFileName;

  paramsListT params;
  categoryT category;
  

  std::pair<std::string,std::string> getKey(const std::string &rawKey)
  {
    std::pair<std::string,std::string> result;
    size_t p;
    p=rawKey.find_first_of(std::string("."));
    //printf("key=%s\n",rawKey.c_str());
    if (p==std::string::npos)
      return std::make_pair(defaultCategory(),rawKey);

    if ((p!=rawKey.find_last_of(std::string(".")))||(p==0)||(p==rawKey.length()-1))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("Ill formated keyword '%s'.\n",rawKey.c_str());
	exit(-1);
      }

    return std::make_pair(rawKey.substr(0,p),rawKey.substr(p+1));
  }

  void insertItem(const std::string &rawKey,const std::string &value, 
		  const std::string &comment, int status=0)
  {
    std::pair<std::string,std::string> key=getKey(rawKey);
    insertItem(key.first,key.second,value,comment,status);
  }

  void insertItem(const std::string &cat, const std::string &what, 
		  const std::string &value, const std::string &comment, 
		  int status=0)
  {
    categoryT_it cit=category.find(cat);
    paramsListT_it lit;
    
    if (cit==category.end())
      {
	//printf("ParamsParser: new category '%s'\n",key.first.c_str());
	lit=params.insert(params.begin(),paramsT());
	category.insert(std::make_pair(cat,lit));
      }
    else lit=cit->second;

    paramsT_it it = lit->find(what);

    if (it==lit->end())
      {
	it = lit->insert(std::make_pair(what,paramValueT())).first;	
      }
    
    it->second.val.push_back(value);
    it->second.status.push_back(status);
    it->second.comment=comment;
    //it->second.val.push_back(fname);
    //it->second.val.push_back(line);

    // empty value means 'true'
    //if (value.size()==0) it->second.val.back().push_back(std::string("1"));


    /*
    std::pair<paramsT_it,bool> lit->insert(std::make_pair(key.second,value));

    if (! lit->insert(std::make_pair(key.second,value)).second)
      {
	//printf("WARNING: multiple definitions of '%s' ignored.\n",rawKey.c_str());
      }
    */
  }

  bool parse(std::string fname)
  {
    std::ifstream t(fname.c_str());

    if (!t) return false;

    std::string line;
    std::vector<std::string> tokens;

    //std::string key;
    //std::vector<std::string> value;

    while (!t.eof())
      {
	tokens.clear();
	std::getline(t,line);

	//printf("line:  %s",line.c_str());
	if (line.find_first_of("#")!=std::string::npos)
	  line=line.substr(0,line.find_first_of("#"));	  
	//printf("-> %s\n",line.c_str());

	StringTokenizer::split(tokens,line,std::string(" =:{}[],\t"));

	if (tokens.size())
	  {
	    std::string key=tokens[0];
	    if (tokens.size()==1)
	      insertItem(key,std::string("1"),"unknown parameter (from file)");
	    for (unsigned int i=1;i<tokens.size();++i)
	      insertItem(key,tokens[i],"unknown parameter (from file)");
	    //std::vector<std::string> value(tokens.begin()+1,tokens.end());
	    //insertItem(key,value);
	  }
      }
    return true;
  }

  bool parse(int argc, char **argv)
  {
    int n=0;
    std::string key;
    //std::vector<std::string> value;

    while (n<argc)
      {
	if (argv[n][0]!='-')
	  {
	    PRINT_SRC_INFO(LOG_ERROR);
	    glb::console->print<LOG_ERROR>("parsing command line at %s\n",argv[n]);
	    exit(-1);
	  }
	
	//value.clear();
	key=std::string(&argv[n++][1]);
	int nInserted=0;
	while ((n<argc)&&((argv[n][0]!='-')||(std::isdigit(argv[n][1])))) 
	  {
	    //value.push_back(argv[n++]);
	    insertItem(key,std::string(argv[n++]),"unknown parameter (from cmd line)");
	    nInserted++;
	    //value.clear();
	  }
	if (nInserted==0) 
	  insertItem(key,std::string("1"),"unknown parameter (from cmd line)");
	/*
	while ((n<argc)&&((argv[n][0]!='-')||(isdigit(argv[n][1])))) value.push_back(argv[n++]);
	
	insertItem(key,value);
	*/
      }

    return true;
  }

public:
  ParamsParser()
  {
    
  }
  
  ParamsParser(int argc, char **argv, std::string fname=std::string("params.ini"))
  {
    init(argc,argv,fname);
  }
  
  ~ParamsParser()
  {

  }

  void init(int argc=0, char **argv=NULL, std::string fname=std::string("params.ini"))
  {
    int i0=0;

    if (argc>i0)
      execFileName=std::string(argv[i0++]);

    if ((argc>i0)&&(fname!=std::string("")))
      {
	if (argv[i0][0]!='-') 
	  fname=std::string(argv[i0++]);
      }
    //printf("Parsing file : %s\n",fname.c_str());
    parse(argc-i0, &argv[i0]);
    if (parse(fname)) 
      paramFileName=fname;
    else if (fname!=std::string("params.ini"))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("could not read parameter file '%s'.\n",fname.c_str());
	exit(-1);
      }
    else paramFileName=std::string();
  }

  template <class LT>
  void report(bool defaultOnly=true)
  {     
    bool comment=glb::console->getVerboseLevel()>=3;
    
    if (defaultOnly)
      glb::console->print<LT>("The following parameters have been assigned default value:\n");
    else
      glb::console->print<LT>("The following parameters have been parsed:\n");

    glb::console->print<LT>(" Input file: '%s'\n",paramFileName.c_str());
    for (categoryT_it cit=category.begin();cit!=category.end();cit++)
      {
	glb::console->print<LT>(" Category: '%s'\n",cit->first.c_str());
	std::string prefix("");
	if (cit->first!=std::string("default"))
	  prefix=cit->first+std::string(".");
	for (paramsT_it pit=cit->second->begin();pit!=cit->second->end();pit++)
	  {
	    if ((pit->second.status[0]>1)&&(defaultOnly)) continue;

	    if (comment) glb::console->print<LT>(" -> ");	    
	    //glb::console->print<LT>("   %s = ",pit->first.c_str());
	    for (unsigned long j=0;j<pit->second.val.size();j++)
	      {
		if (pit->second.val.size()>1)
		  {
		    if (comment)
		      {
			if (j>0)
			  glb::console->print<LT>("     %s%s[%ld] = ",prefix.c_str(),pit->first.c_str(),j);
			else
			  glb::console->print<LT>(" %s%s[%ld] = ",prefix.c_str(),pit->first.c_str(),j);
		      }
		    else
		      {
			if (j>0)
			  glb::console->print<LT>("     %s%s[%ld] = ",prefix.c_str(),pit->first.c_str(),j);
			else
			  glb::console->print<LT>("  -> %s%s[%ld] = ",prefix.c_str(),pit->first.c_str(),j);
		      }
		  }
		else
		  {
		    if (comment)
		      glb::console->print<LT>(" %s%s = ",prefix.c_str(),pit->first.c_str());
		    else
		      glb::console->print<LT>("  -> %s%s = ",prefix.c_str(),pit->first.c_str());
		  }

		//if (j!=0) glb::console->print<LT>("  @[%d] = ",j);
		//unsigned long i;
		
		if (pit->second.status[j]==0)
		  glb::console->print<LT>("%s (unknown parameter)\n",pit->second.val[j].c_str());
		else if (pit->second.status[j]==1)
		  glb::console->print<LT>("%s\n",pit->second.val[j].c_str());
		else
		  glb::console->print<LT>("%s\n",pit->second.val[j].c_str());		
	
		/*
		for (i=0;i<pit->second.val[j].size();i++)
		  {
		    if (pit->second.status[j]==0)
		      glb::console->print<LT>("%s (unknown parameter)",pit->second.val[j][i].c_str());
		    else if (pit->second.status[j]==1)
		      glb::console->print<LT>("%s ",pit->second.val[j][i].c_str());
		    else
		      glb::console->print<LT>("%s ",pit->second.val[j][i].c_str());
		  }
		if (i==0) 
		  {
		    if (pit->second.status[j]==0)
		      glb::console->print<LT>("true (unknown parameter)\n");
		    else 
		      glb::console->print<LT>("true \n");
		  }
		else glb::console->print<LT>("\n");
		*/
	      }
	    if (comment) glb::console->print<LT>("   %s\n\n",pit->second.comment.c_str());
	  }
      }

    //glb::console->flushBuffer<LT>();
  }
  
  template <class LT>
  bool reportUnused()
  {
    bool found=false;    
    
    for (categoryT_it cit=category.begin();cit!=category.end();cit++)
      {
	//printf("cit->first.c_str() = %s\n",cit->first.c_str());
	for (paramsT_it pit=cit->second->begin();pit!=cit->second->end();pit++)
	  {
	    //printf("pit->second.wasRead.size() = %ld\n",pit->second.wasRead.size());
	    
	    for (unsigned long j=0;j<pit->second.status.size();j++)
	      {
		//printf("   '%s'.'%s' (%d)\n ",cit->first.c_str(),pit->first.c_str(),(int)pit->second.wasRead[j]);
		if (!pit->second.status[j])
		  {		    
		    if (!found)
		      {
			glb::console->printToBuffer<LT>("The following parameter(s) definitions were ignored:\n");
			found=true;
		      }
		    if (pit->second.status.size()>1)
		      {
			glb::console->printToBuffer<LT>("   %s.%s[%ld] = ",cit->first.c_str(),pit->first.c_str(),j);
			glb::console->printToBuffer<LT>("%s\n",pit->second.val[j].c_str());
		      }
		    else
		      {
			glb::console->printToBuffer<LT>("   %s.%s = ",cit->first.c_str(),pit->first.c_str());
			glb::console->printToBuffer<LT>("%s\n",pit->second.val[j].c_str());
		      }
		    /*
		    unsigned long i;
		    for (i=0;i<pit->second.val[j].size();i++)
		      {
			glb::console->printToBuffer<LT>("%s ",pit->second.val[j][i].c_str());
		      }
		    if (i==0) glb::console->printToBuffer<LT>("true\n");
		    else glb::console->printToBuffer<LT>("\n");
		    */
		    
		  }
	      }	    
	  }
      }

    glb::console->flushBufferNewLine<LT>();
    return found;
  }

  template <typename T>
  T getOrDie(const std::string &what,const std::string &cat, int index = -1) const
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("category '%s' not found for keyword '%s'.\n",cat.c_str(),what.c_str());
	exit(-1);
      }

    paramsT &p=*(cit->second);
    paramsT_it it=p.find(what);

    if (it==p.end())
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("keyword '%s' not found in category '%s'.\n",what.c_str(),cat.c_str());
	exit(-1);
      }

    if ((index>=0)&&(it->second.val.size()<=index))
      {
	PRINT_SRC_INFO(LOG_ERROR);
	glb::console->print<LOG_ERROR>("'%dth' definition of keyword '%s.%s' not available (%ld total).\n",index,cat.c_str(),what.c_str(),it->second.val.size());
	exit(-1);
      }

    int id=(index<0)?0:index;   
    T result;
    StringTokenizer::from_string<T>(result,it->second.val[id]);
    it->second.status[id]=2;    
      
    return result;   
  }

  template <typename T>
  T get(const std::string &what, const std::string &cat, const T &def,const std::string &comment, int index=-1)
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end())
      {
	std::string tmp;
	StringTokenizer::to_string(def,tmp);
	insertItem(cat,what,tmp,comment,1);

	return def;
      }

    paramsT &p=*(cit->second);
    paramsT_it it=p.find(what);

    if (it==p.end())
      {	
	std::string tmp;
	StringTokenizer::to_string(def,tmp);
	insertItem(cat,what,tmp,comment,1);

	return def;
      }

    if ((index>=0)&&(it->second.val.size()<=static_cast<unsigned int>(index)))
      {
	std::string tmp;
	StringTokenizer::to_string(def,tmp);
	insertItem(cat,what,tmp,comment,1);

	return def;
      }

    int id = (index<0)?0:index;   
    T result;
    StringTokenizer::from_string<T>(result,it->second.val[id]);    
    it->second.status[id]=2;
    it->second.comment = comment;

    return result;
  }

  template <typename T>
  T get(const std::string &what, const std::string &cat, const T &def, int index=-1)
  {
    return get(what,cat,def,"No description available",index);
  }

  long isDefined(const std::string &what, const std::string &cat) const
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end()) return 0;
   
    paramsT &p=*(cit->second);
    paramsT_cit it=p.find(what);

    if (it==p.end()) return 0;

    return it->second.val.size();
  }


};

/** \}*/
#include "../../internal/namespace.footer"
#endif

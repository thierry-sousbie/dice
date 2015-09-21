#ifndef __CONSOLE_HXX__
#define __CONSOLE_HXX__

#include <cstdarg>
#include <string.h>
#include <sstream> 

#include "../../dice_globals.hxx"

// Included here so that files using the console automatically have the LOG_TRAITS
#include "../../tools/IO/log_traits.hxx" 


#include "../../internal/namespace.header"

class Console
{
public: 

  class Out : public std::ostringstream {
  private:
    friend class Console;
    typedef std::ostringstream Base;
    template <class T>
    void print(Console *c,MpiCommunication *com=glb::mpiComWorld)
    {
      c->print<T>(com,"%s",Base::str().c_str());
      Base::str(std::string());
      Base::clear();      
    }
  };

  Out out;
  char buffer[1<<14]; // 16ko
 
  Console(int verbose_=2, int level_=0):
    verbose(verbose_),
    level(level_),
    startingLine(true),
    outputFile(NULL)
  {
    setLevel(level);

    memset(before,0,1024);
    strcpy(before,"");

    memset(after,0,1024);
    strcpy(after,"");

    memset(buffer,0,(1<<14));
    strcpy(buffer,"");

    allowRankPrint=true;
  }

  ~Console()
  {
    if (outputFile != NULL)
      fclose(outputFile);
  }

  void closeOutputToFile()
  {
    if (outputFile!=NULL)
      fclose(outputFile);
    outputFile=NULL;
  }

  void outputToFile(const char *fname)
  {
    closeOutputToFile();
    outputFile=fopen(fname,"w");
  }

  void setStartingLine(bool val=true)
  {
    startingLine=val;
  }
  
  void indent()
  {
    setLevel(level+2);
  }

  void unIndent()
  {
    setLevel(level-2);
  }

  bool getRankPrint()
  {
    return allowRankPrint;
  }

  void setRankPrint(bool val=true)
  {
    allowRankPrint=val;
  }

  void setLevel(int level_)
  {   
    char tmp[1024];
    level=(level_<0)?0:level_;

    if (level==0)
      {
	strcpy(before,"");
	return;
      }

    strcpy(tmp," ");
    for (int i=1;i<level;i++)
      {
	sprintf(before," %s",tmp);
	strcpy(tmp,before);
      }
  }
    
  void setVerboseLevel(int verbose_) {verbose=verbose_;}
  int getVerboseLevel() {return verbose;}

  template <class T>
  bool willPrint() const
  {
    return (verbose>=T::verboseLevel);
  }

  template <class T>
  void print(MpiCommunication *com, const char *s, ...)
  {  
    if (verbose<T::verboseLevel) return;
    va_list fmtargs;
    va_start(fmtargs,s);
    return print_<T>(com,s,fmtargs);
  }

  template <class T>
  void print(const char *s, ...)
  {  
    if (verbose<T::verboseLevel) return;
    va_list fmtargs;
    va_start(fmtargs,s);
    return print_<T>(glb::mpiComWorld,s,fmtargs);
  }

  template <class T>
  void printNewLine(MpiCommunication *com, const char *s, ...)
  {  
    if (verbose<T::verboseLevel) return;
    va_list fmtargs;
    va_start(fmtargs,s);
    startingLine=true;
    //char tmp[255];
    //sprintf(tmp,"\n%s",s);
    return print_<T>(com,s,fmtargs);
  }

  template <class T>
  void printNewLine(const char *s, ...)
  {  
    if (verbose<T::verboseLevel) return;
    va_list fmtargs;
    va_start(fmtargs,s);
    startingLine=true;
    //char tmp[255];
    //sprintf(tmp,"\n%s",s);
    return print_<T>(glb::mpiComWorld,s,fmtargs);
  }  

  template <class T>
  void printFlush(MpiCommunication *com, const char *s, ...)
  {      
    va_list fmtargs;
    va_start(fmtargs,s);
    print_<T>(com,s,fmtargs);
    fflush(0);
  }

  template <class T>
  void printFlush(const char *s, ...)
  {   
    va_list fmtargs;
    va_start(fmtargs,s);
    print_<T>(glb::mpiComWorld,s,fmtargs);
    fflush(0);
  }

  template <class T>
  void printToBuffer(const char *s, ...)
  {  
    if (verbose<T::verboseLevel) return;
    va_list fmtargs;
    va_start(fmtargs,s);
    char msg[2048];
    vsprintf(msg,s,fmtargs);
    strcat(buffer,msg);
  } 
  
  template <class T>
  void flushBuffer(MpiCommunication *com=glb::mpiComWorld)
  {
    if (!strcmp(buffer,"")) return;
    long len=strlen(buffer);
    if (buffer[len-1]=='\n')
      {
	buffer[len-1]='\0';
	print<T>(com,"%s\n",buffer);
      }
    else print<T>(com,"%s",buffer);
    strcpy(buffer,"");
  }

  template <class T>
  void flushBufferNewLine(MpiCommunication *com=glb::mpiComWorld)
  {
    if (!strcmp(buffer,"")) return;
    startingLine=true;
    long len=strlen(buffer);
    if (buffer[len-1]=='\n')
      {
	buffer[len-1]='\0';
	print<T>(com,"%s\n",buffer);
      }
    else print<T>(com,"%s",buffer);
    strcpy(buffer,"");
  }

  template <class T>
  void printSrcInfo(const char *pfunc,const char *func, const char *file, const int line,MpiCommunication *com=glb::mpiComWorld)
  {
    if (verbose<T::verboseLevel) return;
    print<T>(com,"\nIn '%s' line %d.\n",file,line);
    print<T>(com,"Function '%s':\n",func);
    print<T>(com,"  %s.\n",pfunc);
  } 

  
  template <class T>
  void flushStream(MpiCommunication *com=glb::mpiComWorld)
  {
    out.print<T>(this,com);
  }

private:
  int verbose;
  int level;
  bool startingLine;
  bool allowRankPrint;

  char before[1024];
  char after[1024];

  FILE *outputFile;
  
  template <class T>
  void print_(MpiCommunication *com, const char *s, va_list &fmtargs)
  {  
    if (verbose<T::verboseLevel) return;
    bool notMainRank=((T::mpiMainRankOnly)&&(com->rank() != 0));
    //if ((T::mpiMainRankOnly)&&(com->rank() != 0)) return;
   
    char msg[2048];
    char newS[2048];    
    long len=strlen(s);

    if (startingLine) 
      sprintf(msg,"%s%s%s%s",T::header,before,s,after);
    else
      sprintf(msg,"%s%s",s,after);
    /*
    if (startingLine) 
      {
	if ((!T::mpiMainRankOnly)&&(allowRankPrint))
	  sprintf(newS,"[%d/%d]%s",com->rank(),com->size(),T::header);
	else 
	  strcpy(newS,T::header);
	    
	sprintf(msg,"%s%s%s%s",newS,before,s,after);
      }
    else sprintf(msg,"%s%s",s,after);
    */

    vsprintf(newS,msg,fmtargs);
       
    if (outputFile!=NULL)
      fprintf(outputFile,"%s",newS);

    if ((startingLine)&&(!T::mpiMainRankOnly)&&(allowRankPrint))
      {
	sprintf(msg,"[%d/%d] %s",com->rank(),com->size(),newS);
	strcpy(newS,msg);
      }

    if (s[len-1]=='\n')
      startingLine=true;
    else
      startingLine=false;

    if (notMainRank) return;

    for (int i=0;i<T::nOutput;i++)
      {
	  fprintf(T::output[i],"%s",newS);
      }

  } 
  
};

#include "../../internal/namespace.footer"
#endif

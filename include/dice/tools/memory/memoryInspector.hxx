#ifndef __MEMORY_INSPECTOR_HXX__
#define __MEMORY_INSPECTOR_HXX__

#include <sys/types.h>
#include <sys/time.h>
//#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../dice_globals.hxx"


#include "../../internal/namespace.header"

class MemoryInspector
{
public:
 
  MemoryInspector()
  {
    pid=getpid();
    sprintf(path, "/proc/%d/status", pid);
    FILE *f=fopen(path,"r");
    if (f==NULL)
      {
	PRINT_SRC_INFO(LOG_WARNING);
	glb::console->print<LOG_WARNING>("Could not open '%s' for memory inspection.\n  MemoryInspector will be disabled.\n",path);
	pid=0;
      }
    else fclose(f);
  }

  ~MemoryInspector()
  {
    
  }

  template <class L>
  bool reportAllProcesses(MpiCommunication *com=glb::mpiComWorld)
  {   
    double result[5];
    if (!getStat(result)) return false;
    for (int i=1;i<5;i++) result[i] /= (1<<20);
    if (com->size()>1) 
      {
	double resultMax[5];
	for (int i=1;i<5;i++) resultMax[i] = result[i];
        com->max(resultMax,5);
	com->min(result,5);
	glb::console->print<L>("Memory status (min/max) [size, peak]: [%lg/%lg, %lg/%lg] Go.\n",result[1],resultMax[1],result[2],resultMax[2]);
      }
    else
      glb::console->print<L>("Memory status [size, peak]: [%lg, %lg] Go.\n",result[1],result[2]);
    return true;
  }

  template <class L>
  bool report()
  {
    //double resultMax[5];
    double result[5];
    if (!getStat(result)) return false;
    for (int i=1;i<5;i++) result[i] /= (1<<20);
    
    glb::console->print<L>("Memory status [size, peak]: [%lg, %lg] Go.\n",result[1],result[2]);
    return true;
  }

  // returns 5 floating point values :  [time (s), size (kB), peak (kB), rss (kB), hwm (kB)]
  template <class OutputIterator>
  bool getStat(OutputIterator out) 
  {
    if (pid==0)
      {
	for (int i=0;i<5;i++) (*out)=0;
	return false;
      }
    struct timeval time_ref,time;
    static long ncall=0;
    static double t_ref=0;
    double t=0;
    char *vmsize=NULL;
    char *vmpeak=NULL;
    char *vmrss=NULL;
    char *vmhwm=NULL;
    size_t len=128;
    char *line=new char[len];//(char*)malloc(len);

    if (!ncall)
      {
	gettimeofday(&time_ref, NULL);
	t_ref=(double)time_ref.tv_sec + ((double)time_ref.tv_usec)*1.E-6;
	t=0;
	ncall++;
      }
    else
      {
	gettimeofday(&time, NULL);
	t  = (double)time.tv_sec + ((double)time.tv_usec)*1.E-6;
	t -= t_ref;
	ncall++;
      }

    FILE *f = fopen(path, "r");
    while (!vmsize || !vmpeak || !vmrss || !vmhwm)
      {
	if (getline(&line, &len, f) == -1)
	  {
	    // Some of the information is missing
	    PRINT_SRC_INFO(LOG_WARNING);
	    glb::console->print<LOG_WARNING>("Memory information not available.\n  MemoryInspector will be disabled.\n",path);
	    pid=0;
	    for (int i=0;i<5;i++) (*out)=0;
	    free(vmpeak);free(vmsize);free(vmrss);free(vmhwm);
	    free(line);
	    fclose(f);
	    return false;
	  }
	if (!strncmp(line, "VmPeak:", 7))
	  {
	    vmpeak = strdup(&line[7]);
	  }
	else if (!strncmp(line, "VmSize:", 7))
	  {
	    vmsize = strdup(&line[7]);
	  }
	else if (!strncmp(line, "VmRSS:", 6))
	  {
	    vmrss = strdup(&line[7]);
	  }
	else if (!strncmp(line, "VmHWM:", 6))
	  {
	    vmhwm = strdup(&line[7]);
	  }
      }
    fclose(f);
    //free(line);
    delete[] line; 

    /* Get rid of " kB\n"*/
    /*
    len = strlen(vmsize);
    vmsize[len - 4] = 0;
    len = strlen(vmpeak);
    vmpeak[len - 4] = 0;
    len = strlen(vmrss);
    vmrss[len - 4] = 0;
    len = strlen(vmhwm);
    vmhwm[len - 4] = 0;
    */
    *out=t;++out;
    *out=atof(vmsize);++out;
    *out=atof(vmpeak);++out;   
    *out=atof(vmrss);++out;
    *out=atof(vmhwm);++out;

    free(vmpeak);
    free(vmsize);
    free(vmrss);
    free(vmhwm);

    return true;
  }

private:
  char path[256];
  pid_t pid;
};

#include "../../internal/namespace.footer"
#endif

#include <vector>
#include <string>
#include <iostream> 
#include <fstream>
#include <queue>    
#include <stdlib.h>
#include <string.h>
#include <sstream> 

static const int NDIM=3;

#include <dice/tools/helpers/helpers.hxx>
#include "vtrTools.hxx"
#include "vtrBox.hxx"
#include "vtrProject.hxx"


void usage(char *exe)
{
  printf("\n");
  printf("This program is used to extract slices from 2D and 3D vtr files.\n");
  printf("Usage : %s <filenames> [-box <i0> <j0> <k0> <di> <dj> <dk>] [-project <dir>]\n",exe);
  printf("Notes:\n");
  printf(" - Use k0=0 and dk=1 for 2D images.\n");
  printf(" - Several 'filenames' and '-box' options can be used at the same times\n");
}

int main(int argc, char **argv)
{
  std::vector<std::string> fName;
  std::vector<std::string> outName;
  enum Action { act_center, act_slice, act_box, act_project };
  std::vector<Action> actions;
  //std::vector<Header> headers;
  std::queue<double> params;

  for (int i=1;i<argc;)
    {
      std::string arg=argv[i++];
      if (arg[0]!='-')
	{
	  fName.push_back(arg);
	  outName.push_back(vtrSlicer::cutName(arg.c_str()));
	  //headers.resize(headers.size()+1);
	  //headers.back().read(arg);
	}
      /*
      else if (arg==std::string("-center"))
	{
	  actions.push_back(act_center);
	  for (int j=0;j<NDIM;++j)
	    {
	      double val=strtod(argv[i++],NULL);
	      params.push(val);
	    }
	}
      else if (arg==std::string("-slice"))
	{
	  actions.push_back(act_slice);
	  for (int j=0;j<2;++j)
	    {
	      double val=strtod(argv[i++],NULL);
	      params.push(val);
	    }
	}
      */
      else if (arg==std::string("-box"))
	{
	  actions.push_back(act_slice);
	  for (int j=0;j<2*NDIM;++j)
	    {
	      double val=strtod(argv[i++],NULL);
	      params.push(val);
	    }
	}      
      else if (arg==std::string("-project"))
	{
	  actions.push_back(act_project);
	  double val=strtod(argv[i++],NULL);
	  params.push(val);
	}    
      else
	{
	  std::cout << "Unknown argument : " << arg <<std::endl;
	  usage(argv[0]);
	  exit(-1);
	  exit(-1);
	}
    }

  if (fName.size()==0)
    {
      std::cout << "At least one file name must be provided !" <<std::endl;
      usage(argv[0]);
      exit(-1);
    }

  for (int a=0;a<actions.size();a++)
    {
      if (actions[a]==act_box)
	vtrSlicer::box(params,fName,outName);
      else if (actions[a]==act_project)
	vtrSlicer::project(params,fName,outName);
    }

  return 0;
}

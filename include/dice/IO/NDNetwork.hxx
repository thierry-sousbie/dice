#ifndef __ND_NETWORK_HXX__
#define __ND_NETWORK_HXX__

// This is ugly because it is legacy code, sorry ;(

#include <stdlib.h>
#include <stdio.h>

#include "../tools/types/types.hxx"
#include "../tools/IO/IOHelpers.hxx"
#include "../tools/IO/console.hxx"

#include "../internal/namespace.header"

// #define NDNET_UINT                UINT
// #define NDNET_INT                 INT
// #define NDNET_FLOAT               FLOAT
// #define NDNET_IDCUMT              long

#define NDNETWORK_DATA_STR_SIZE  16
#define NDNETWORK_TAG            "NDNETWORK"
#define NDNETWORK_ASCII_TAG      "ANDNET"
#define NDNETWORK_ASCII_ADDTAG   "[ADDITIONAL_DATA]"
#define NDNETWORK_V_FLAG_TAG     "internal_vertex_flags"
#define NDNETWORK_F_FLAG_TAG     "internal_face_flags"

namespace IO {

  struct NDnetwork_Data
  {
    int type; // the cell-type
    char name[255];  // name of the field
    double *data;  // value for each of the nfaces[type] type-faces
  };

  struct NDnetwork_SupData
  {
    int type; //0 for vertex, t for faces of type (i.e. dimension) t.
    char name[255];
    int datasize;
    char datatype[255];// a string to identity how data should be casted
    void *data;
  };

  template <class F = FLOAT, class UI = UINT, class I = INT, class IDC = long>
  class NDnetwork
  {
  public:

    typedef UI NDNET_UINT;
    typedef I NDNET_INT;
    typedef F NDNET_FLOAT;
    typedef IDC NDNET_IDCUMT;

    char comment[80];
    int periodicity;
    int ndims; // the number of spatial dimensions
    int ndims_net; // number of dimension of the network itself (e.g. 2 for a sphere embedded in 3D)
    int isSimpComplex;  // 1 if network is a simplicial complex (always true in disperse)
    double *x0;  // origin of the bounding box
    double *delta;  // size of the bounding box
    int indexSize; // size of NDNET_UINT type in Bytes
    int cumIndexSize; // size of NDNET_IDCUMT type in Bytes
    int floatSize; // size of NDNET_IDCUMT type in Bytes
    char dummy[160-4*2-4]; // dummy data reserved for future extensions
  
    NDNET_UINT nvertex;  // total number of vertices
    NDNET_FLOAT *v_coord; //vertices coodinates (X_0,Y_0,Z_0,X_1,Y_1,...,Z_nvertex-1)
    
    NDNET_UINT *nfaces; // number of faces of a given type t is given by nfaces[t]
  
    int *haveVertexFromFace; // haveVertexFromFace[n] is 1 if we have an explicit definition of the n-faces (at least one type of face must be defined).
    NDNET_IDCUMT **f_numVertexIndexCum;// cumulative number of vertice in the faces f of type t, NULL for simplicial faces
    NDNET_UINT **f_vertexIndex; // list of vertices defining the n-faces is stored in f_vertexIndex[n], all vertices being enumerated for each face (the indices of the vertices in the kth n-face start at f_vertexIndex[n][(n+1)*k] )
    // see also macro  NUM_VERTEX_IN_FACE(net,type,face) and VERTEX_IN_FACE(net,type,face)
  
   
    int *haveFaceFromVertex; // haveFaceFromVertex[n] is 1 if we have an explicit list of all the n-faces that contain each vertex (used to navigate within the network)
    NDNET_IDCUMT **v_numFaceIndexCum; // cumulative number of faces of type t a vertex v belongs to
    NDNET_UINT **v_faceIndex; // indices of the faces of type t a vertex v belongs to ( the list of n-faces of vertex k starts at v_faceIndex[n][net->v_numFaceIndexCum[n][k]] and ends at v_faceIndex[n][net->v_numFaceIndexCum[n][k+1]] )
    // see also macro  NUM_FACE_IN_VERTEX(net,type,vertex) and  FACE_IN_VERTEX(net,type,vertex)
        
    int **haveFaceFromFace; // haveFaceFromFace[k][n] is 1 if we have an explicit list of all the n-faces that have a boundary/co-boundary relation with each k-face (used to navigate within the network)
    NDNET_IDCUMT ***f_numFaceIndexCum; //  cumulative number of n-faces having a boundary / co-boundary relation with each k-face f_numFaceIndexCum[k][n]
    NDNET_UINT ***f_faceIndex; // indices of the faces (similar to v_faceIndex)
    // see also macro NUM_FACE_IN_FACE(net,ref_type,ref_face,type) and FACE_IN_FACE(net,ref_type,ref_face,type)
    
    int haveVFlags;  // do we have flags associated to each vertex ?
    int *haveFFlags;  // do we have flags associated to each n-face ?
    unsigned char *v_flag; // nvertex flag values (1 for each vertex) or NULL 
    unsigned char **f_flag; // nfaces[n] flag values (1 of each n-face) or NULL
 
    int ndata; // number of additional data fields.
    NDnetwork_Data *data; // array of all additionnal data (data in total)
    
    int nsupData;
    NDnetwork_SupData *supData;

    bool fileLoaded;
    
    template <class T1,class T2,class T3,class T4,class T5>
    NDnetwork(T1 nDims, T2 nDimsW, T3 x0_[], T4 delta_[], T5 nCells[], 
	      const char *cmt, bool periodic)
    {
      //memset(&(*this),0,sizeof(NDnetwork));
      setNULL();
      create(nDimsW,nDims,nCells);
      strcpy(comment,cmt);
      if (periodic) periodicity=((1<<nDims)-1);
      for (int i=0;i<nDims;i++) 
	{
	  x0[i]=x0_[i];
	  delta[i]=delta_[i];
	}   
      for (int i=nDims;i<nDimsW;i++) 
	{
	  x0[i]=0;
	  delta[i]=0;
	}
      fileLoaded=false;
    }

    NDnetwork()
    {
      setNULL();
      fileLoaded=false;
      /*
      if (!IsNDnetwork(filename))
	{
	  dice::glb::console->printFlush<dice::LOG_ERROR>
	  ("Invalid initial conditions type: '%s'.\n",icType.c_str());
	  fprintf(stderr,"ERROR loading file")
	  isLoaded=false;
	}
      else
	{
	  read(fname);
	  isLoaded=true;
	}
      */
    }

    ~NDnetwork()
    {
      freeData();
    }

    bool isLoaded() const
    {
      return fileLoaded;
    }

    static int IsNDnetwork(const char *filename)
    {
      int i;
      char tag[NDNETWORK_DATA_STR_SIZE+2];
      FILE *f;
      int swap=0;
  
      memset(tag,0,16*sizeof(char));
   
      if(!(f = fopen(filename,"r")))
	{
	  glb::console->printFlush<LOG_WARNING>
	  ("File %s does not exist.\n",filename);
	  //fprintf(stderr,"File %s does not exist.\n",filename);
	  return 0;
	}

      myIO::fread(&i,sizeof(int),1,f,swap);
      if (i!=NDNETWORK_DATA_STR_SIZE) swap=1-swap;
      myIO::fread(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f,swap);
      myIO::fread(&i,sizeof(int),1,f,swap);
  
      fclose(f); 
      tag[NDNETWORK_DATA_STR_SIZE+1]='\0';

      if (strcmp(tag,NDNETWORK_TAG)) 
	{
	  f = fopen(filename,"r");
	  char *line=NULL;int n;
	  myIO::getLine(&line,&n,f);      
	  fclose(f);
      
	  line[strlen(line)-1]='\0';
	  if (!strcmp(line,NDNETWORK_ASCII_TAG)) 
	    {
	      free(line);
	      return 2;
	    }

	  free(line);
	  return 0;
	}

      return 1;
    }
    
    void writeHeader(FILE *f)
    {
      long i,j;
      char tag[NDNETWORK_DATA_STR_SIZE];
       
      memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
      strcpy(tag,NDNETWORK_TAG);
      i=NDNETWORK_DATA_STR_SIZE;
     
      fwrite(&i,sizeof(int),1,f);
      fwrite(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f);
      fwrite(&i,sizeof(int),1,f);

      j=sizeof(int)*2;
      fwrite(&j,sizeof(int),1,f);
      fwrite(&ndims,sizeof(int),1,f);
      fwrite(&ndims_net,sizeof(int),1,f);
      fwrite(&j,sizeof(int),1,f);  
      j=sizeof(NDNET_UINT)+sizeof(int)*4+80*sizeof(char)+ndims*sizeof(double)*2+160*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(comment,sizeof(char),80,f);
      fwrite(&periodicity,sizeof(int),1,f);
      fwrite(&isSimpComplex,sizeof(int),1,f);
      fwrite(x0,sizeof(double),ndims,f);
      fwrite(delta,sizeof(double),ndims,f);
      fwrite(&indexSize,sizeof(int),1,f);
      fwrite(&cumIndexSize,sizeof(int),1,f);
      fwrite(&floatSize,sizeof(int),1,f);
      fwrite(dummy,sizeof(char),160-sizeof(int)*3,f);
      fwrite(&nvertex,sizeof(NDNET_UINT),1,f);
      fwrite(&j,sizeof(int),1,f);

    }

    int read(const char *filename, int verbose=1)
    {
      long i,j,k;
      int dummy;
      char tag[NDNETWORK_DATA_STR_SIZE];
      FILE *f;
      int swap=0;
          
      //if (IsNDnetwork(filename) == 2) return Load_NDnetwork_ASCII(filename);
      if (IsNDnetwork(filename) != 1)
	{
	  glb::console->printFlush<LOG_WARNING>
	    ("File %s is not a valid NDnetwork.\n",filename);
	  return -1;
	}

      freeData();
      //net = calloc(1,sizeof(NDnetwork));
      memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
  
      f=fopen(filename,"r");
  
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      if (dummy!=NDNETWORK_DATA_STR_SIZE) swap = 1-swap;
      myIO::fread(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
  
      if (strcmp(tag,NDNETWORK_TAG))
	{
	  fclose(f);
	  fprintf (stderr,"File %s has an unknown format.\n",filename);
	  return -2;
	}
      
  
      j=sizeof(int);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(&ndims,sizeof(int),1,f,swap);
      myIO::fread(&ndims_net,sizeof(int),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      if (verbose>0)
	 glb::console->printFlush<LOG_STD>
	   ("Loading %dD network from file \"%s\" ...",ndims,filename);

      x0=(double*)malloc(ndims*sizeof(double));
      delta=(double*)malloc(ndims*sizeof(double));

      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(comment,sizeof(char),80,f,swap);
      myIO::fread(&periodicity,sizeof(int),1,f,swap);
      myIO::fread(&isSimpComplex,sizeof(int),1,f,swap);
      myIO::fread(x0,sizeof(double),ndims,f,swap);
      myIO::fread(delta,sizeof(double),ndims,f,swap);
      myIO::fread(&indexSize,sizeof(int),1,f,swap);
      
      if (indexSize != 8) indexSize=4;
      myIO::fread(&cumIndexSize,sizeof(int),1,f,swap);
      if (cumIndexSize != 8) cumIndexSize=4;
      myIO::fread(&floatSize,sizeof(int),1,f,swap);
      if (floatSize != 8) floatSize=4;
      char trash[256];
      myIO::fread(trash,sizeof(char),160-3*sizeof(int),f,swap);
      
      //myIO::fread(&nvertex,sizeof(NDNET_UINT),1,f,swap);
      myIO::fread_ului(&nvertex,sizeof(NDNET_UINT),1,indexSize,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      // FIXME : NO SUPPORT FOR DOUBLE COORDS IMPLEMENTED, JUST CONVERTING NOW
      // should be:
      //myIO::fread_fd(v_coord,sizeof(NDNET_FLOAT),(size_t)nvertex*(size_t)ndims,floatSize,f,swap);
      //j=ndims*nvertex;
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      v_coord=(NDNET_FLOAT*) malloc(sizeof(NDNET_FLOAT)*(size_t)ndims*(size_t)nvertex);
      //printf("->%ld\n",(unsigned long)v_coord);
      //v_coord=myIO::fread_fd(v_coord,sizeof(float),(size_t)nvertex*(size_t)ndims,floatSize,f,swap);
      myIO::fread_fd(v_coord,sizeof(NDNET_FLOAT),(size_t)nvertex*(size_t)ndims,floatSize,f,swap);
      // whatever the input format, it is always converted single precision !
      floatSize = sizeof(float);
      //printf("=->%ld\n",(unsigned long)v_coord);
      //myIO::fread(v_coord,sizeof(float)*(size_t)ndims*(size_t)nvertex,1,f,swap);   
      myIO::fread(&dummy,sizeof(int),1,f,swap);
     
      //j=(1+ndims)*sizeof(uint);
      nfaces=(NDNET_UINT*)malloc(sizeof(NDNET_UINT)*((size_t)ndims+1));
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread_ului(nfaces,sizeof(NDNET_UINT),((size_t)ndims+1),indexSize,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      //j=(1+ndims)*sizeof(int);
      haveVertexFromFace=(int*)malloc(sizeof(int)*((size_t)ndims+1));
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(haveVertexFromFace,sizeof(int)*((size_t)ndims+1),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      f_vertexIndex = (NDNET_UINT**)calloc(sizeof(NDNET_UINT*),(1+ndims));
      f_numVertexIndexCum = (NDNET_IDCUMT**)calloc(sizeof(NDNET_IDCUMT*),(1+ndims));

      for (i=0;i<1+ndims;i++)
	{    
      
	  if (haveVertexFromFace[i])
	    {
	      if (!isSimpComplex)
		{
	    
		  //j=sizeof(uint)*((size_t)nfaces[i]+1);
		  f_numVertexIndexCum[i]=
		    (NDNET_IDCUMT*)malloc(sizeof(NDNET_IDCUMT)*((size_t)nfaces[i]+1));
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		  myIO::fread_ului(f_numVertexIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)nfaces[i]+1),cumIndexSize,f,swap);
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
	    
		  //j=sizeof(uint)*((size_t)f_numVertexIndexCum[i][nfaces[i]]);
		  f_vertexIndex[i]=
		    (NDNET_UINT*)malloc(sizeof(NDNET_UINT)*((size_t)f_numVertexIndexCum[i][nfaces[i]]));
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		  myIO::fread_ului(f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)f_numVertexIndexCum[i][nfaces[i]]),indexSize,f,swap);
		  myIO::fread(&dummy,sizeof(int),1,f,swap); 
		}
	      else
		{
		  //j=sizeof(uint)*((size_t)(i+1)*nfaces[i]);
		  f_vertexIndex[i]=(NDNET_UINT*)malloc(sizeof(NDNET_UINT)*((size_t)(i+1)*nfaces[i]));
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		  myIO::fread_ului(f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)(i+1)*nfaces[i]),indexSize,f,swap);
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		}
	  
	    }
	}
  
      j=(1+ndims)*sizeof(int);
      haveFaceFromVertex=(int*)malloc(sizeof(int)*((size_t)ndims+1));
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(haveFaceFromVertex,sizeof(int)*((size_t)ndims+1),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      v_faceIndex = (NDNET_UINT**)calloc(sizeof(NDNET_UINT*),(1+ndims));
      v_numFaceIndexCum = (NDNET_IDCUMT**)calloc(sizeof(NDNET_IDCUMT*),(1+ndims));
      //printf("%ld == %d\n",sizeof(NDNET_IDCUMT),cumIndexSize);
      for (i=0;i<1+ndims;i++)
	{     
	  //printf("HELLO %ld\n",i);
	  if (haveFaceFromVertex[i])
	    {
	      //j=sizeof(uint)*((size_t)nvertex+1);
	
	      v_numFaceIndexCum[i]=(NDNET_IDCUMT*)malloc(sizeof(NDNET_IDCUMT)*((size_t)nvertex+1));
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread_ului(v_numFaceIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)nvertex+1),cumIndexSize,f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	  
	      //j=sizeof(uint)*((size_t)v_numFaceIndexCum[i][nvertex]);
	      v_faceIndex[i]=(NDNET_UINT*)malloc(sizeof(NDNET_UINT)*((size_t)v_numFaceIndexCum[i][nvertex]));
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread_ului(v_faceIndex[i],sizeof(NDNET_UINT),((size_t)v_numFaceIndexCum[i][nvertex]),indexSize,f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);	  
	    }
	}
  
      haveFaceFromFace =(int **)calloc(sizeof(int *),(1+ndims));
      f_faceIndex = (NDNET_UINT***)calloc(sizeof(NDNET_UINT**),(1+ndims));
      f_numFaceIndexCum = (NDNET_IDCUMT***)calloc(sizeof(NDNET_IDCUMT**),(1+ndims));

      myIO::fread(&dummy,sizeof(int),1,f,swap);
      for (k=0;k<ndims+1;k++)
	{
	  f_faceIndex[k] = (NDNET_UINT**)calloc(sizeof(NDNET_UINT*),(1+ndims));
	  f_numFaceIndexCum[k] = (NDNET_IDCUMT**)calloc(sizeof(NDNET_IDCUMT*),(1+ndims));

	  haveFaceFromFace[k]=(int*)malloc(sizeof(int)*((size_t)ndims+1));
	  myIO::fread(haveFaceFromFace[k],sizeof(int)*((size_t)ndims+1),1,f,swap);
	}   
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      for (i=0;i<1+ndims;i++)
	{
	  for (k=0;k<1+ndims;k++)
	    {
	  
	      if (haveFaceFromFace[i][k])
		{
		  //j=sizeof(uint)*(1+(size_t)nfaces[i]);
		  f_numFaceIndexCum[i][k]=(NDNET_IDCUMT*)malloc(sizeof(NDNET_IDCUMT)*((size_t)nfaces[i]+1));
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		  myIO::fread_ului(f_numFaceIndexCum[i][k],sizeof(NDNET_IDCUMT),((size_t)nfaces[i]+1),cumIndexSize,f,swap);
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
	      
		  //j=sizeof(uint)*((size_t)f_numFaceIndexCum[i][k][nfaces[i]]);
		  f_faceIndex[i][k]=(NDNET_UINT*)malloc(sizeof(NDNET_UINT)*((size_t)f_numFaceIndexCum[i][k][nfaces[i]]));
		  myIO::fread(&dummy,sizeof(int),1,f,swap);
		  myIO::fread_ului(f_faceIndex[i][k],sizeof(NDNET_UINT),((size_t)f_numFaceIndexCum[i][k][nfaces[i]]),indexSize,f,swap);
		  myIO::fread(&dummy,sizeof(int),1,f,swap);	  
		}
	    }
	}

      j=sizeof(int);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(&haveVFlags,sizeof(int),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);

      if (haveVFlags)
	{
	  j=sizeof(unsigned char)*nvertex;
	  v_flag=(unsigned char*)malloc(sizeof(unsigned char)*(size_t)nvertex);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);
	  myIO::fread(v_flag,sizeof(unsigned char),(size_t)nvertex,f,swap);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);
	}

      haveFFlags=(int*)calloc(ndims+1,sizeof(int));
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(haveFFlags,sizeof(int),(ndims+1),f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
  
      f_flag=(unsigned char **)calloc(ndims+1,sizeof(unsigned char *));
  
      for (i=0;i<1+ndims;i++)
	if (haveFFlags[i])
	  {	  
	    j=sizeof(unsigned char)*nfaces[i];
	 
	    if (j) f_flag[i]=(unsigned char*)malloc(sizeof(unsigned char)*(size_t)nfaces[i]);
	    myIO::fread(&dummy,sizeof(int),1,f,swap);
	    if (j) myIO::fread(f_flag[i],sizeof(unsigned char),(size_t)nfaces[i],f,swap);
	    myIO::fread(&dummy,sizeof(int),1,f,swap);
	  }

  
      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(&ndata,sizeof(int),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
  

      data=(NDnetwork_Data*)malloc(sizeof(NDnetwork_Data)*ndata);
      for (i=0;i<ndata;i++)
	{     
	  //j=sizeof(int)+255*sizeof(char);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);
	  myIO::fread(&(data[i].type),sizeof(int),1,f,swap);
	  myIO::fread(data[i].name,sizeof(char)*255,1,f,swap);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);

	  if (data[i].type==0)
	    {
	      //j=sizeof(double)*nvertex;
	      data[i].data=(double*)malloc(sizeof(double)*(size_t)nvertex);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread(data[i].data,sizeof(double),(size_t)nvertex,f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	    }
	  else
	    {
	      //j=sizeof(double)*nfaces[data[i].type];
	      data[i].data=(double*)malloc(sizeof(double)*(size_t)nfaces[data[i].type]);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread(data[i].data,sizeof(double),(size_t)nfaces[data[i].type],f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	    }  
	}

      myIO::fread(&dummy,sizeof(int),1,f,swap);
      myIO::fread(&nsupData,sizeof(int),1,f,swap);
      myIO::fread(&dummy,sizeof(int),1,f,swap);
  
      supData=(NDnetwork_SupData*)malloc(sizeof(NDnetwork_SupData)*nsupData);
      for (i=0;i<nsupData;i++)
	{     
	  j=sizeof(int)+255*sizeof(char);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);
	  myIO::fread(&(supData[i].type),sizeof(int),1,f,swap);
	  myIO::fread(supData[i].name,sizeof(char)*255,1,f,swap);
	  myIO::fread(&(supData[i].datasize),sizeof(int),1,f,swap);
	  myIO::fread(supData[i].datatype,sizeof(char)*255,1,f,swap);
	  myIO::fread(&dummy,sizeof(int),1,f,swap);

	  if (supData[i].type==0)
	    {
	      //j=(size_t)supData[i].datasize*nvertex;
	      supData[i].data=(NDnetwork_SupData*)malloc((size_t)supData[i].datasize*(size_t)nvertex);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread(supData[i].data,(size_t)supData[i].datasize,(size_t)nvertex,f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	    }
	  else
	    {
	      //j=(size_t)supData[i].datasize*nfaces[supData[i].type];
	      supData[i].data=(NDnetwork_SupData*)malloc((size_t)supData[i].datasize*(size_t)nfaces[supData[i].type]);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	      myIO::fread(supData[i].data,(size_t)supData[i].datasize,(size_t)nfaces[supData[i].type],f,swap);
	      myIO::fread(&dummy,sizeof(int),1,f,swap);
	    }  
	}

      fclose(f);
      fileLoaded=true;

      floatSize=sizeof(NDNET_FLOAT);
      indexSize=sizeof(NDNET_UINT);
      cumIndexSize=sizeof(NDNET_IDCUMT);
 
      if (verbose>0) glb::console->print<LOG_STD>(" done.\n");
      return 0;
    }

    int write(const char *filename, int verbose=0)
    {
      int j;
      long i,k;
      char tag[NDNETWORK_DATA_STR_SIZE];
      FILE *f;
  
      if (verbose>1) printf ("Saving %dD network to file %s ...",ndims,filename);fflush(0);

      floatSize=sizeof(NDNET_FLOAT);
      indexSize=sizeof(NDNET_UINT);
      cumIndexSize=sizeof(NDNET_IDCUMT);
      
      memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
      strcpy(tag,NDNETWORK_TAG);
      i=NDNETWORK_DATA_STR_SIZE;
  
      f=fopen(filename,"w");
      fwrite(&i,sizeof(int),1,f);
      fwrite(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f);
      fwrite(&i,sizeof(int),1,f);

      j=sizeof(int)*2;
      fwrite(&j,sizeof(int),1,f);
      fwrite(&ndims,sizeof(int),1,f);
      fwrite(&ndims_net,sizeof(int),1,f);
      fwrite(&j,sizeof(int),1,f);  
      j=sizeof(NDNET_UINT)+sizeof(int)*4+80*sizeof(char)+ndims*sizeof(double)*2+152*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(comment,sizeof(char),80,f);
      fwrite(&periodicity,sizeof(int),1,f);
      fwrite(&isSimpComplex,sizeof(int),1,f);
      fwrite(x0,sizeof(double),ndims,f);
      fwrite(delta,sizeof(double),ndims,f);
      fwrite(&indexSize,sizeof(int),1,f);
      fwrite(&cumIndexSize,sizeof(int),1,f);
      fwrite(&floatSize,sizeof(int),1,f);      
      fwrite(dummy,sizeof(char),160-3*sizeof(int),f);
      fwrite(&nvertex,sizeof(NDNET_UINT),1,f);
      fwrite(&j,sizeof(int),1,f);
  
      j=ndims*nvertex;
      fwrite(&j,sizeof(int),1,f);
      fwrite(v_coord,sizeof(NDNET_FLOAT),(size_t)ndims*(size_t)nvertex,f);
      fwrite(&j,sizeof(int),1,f);
  
      j=(1+ndims)*sizeof(NDNET_UINT);
      fwrite(&j,sizeof(int),1,f);
      fwrite(nfaces,sizeof(NDNET_UINT),((size_t)ndims+1),f);
      fwrite(&j,sizeof(int),1,f);

      j=(1+ndims)*sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      fwrite(haveVertexFromFace,sizeof(int),((size_t)ndims+1),f);
      fwrite(&j,sizeof(int),1,f);

  
      for (i=0;i<1+ndims;i++)
	{
      
	  if (haveVertexFromFace[i])
	    {
	      if (!isSimpComplex)
		{
		  j=sizeof(NDNET_IDCUMT)*((size_t)nfaces[i]+1);
		  fwrite(&j,sizeof(int),1,f);
		  fwrite(f_numVertexIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)nfaces[i]+1),f);
		  fwrite(&j,sizeof(int),1,f);

		  j=sizeof(NDNET_UINT)*((size_t)f_numVertexIndexCum[i][nfaces[i]]);
		  fwrite(&j,sizeof(int),1,f);
		  fwrite(f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)f_numVertexIndexCum[i][nfaces[i]]),f);
		  fwrite(&j,sizeof(int),1,f);
		}
	      else
		{
		  j=sizeof(NDNET_UINT)*((size_t)(i+1)*nfaces[i]);
		  fwrite(&j,sizeof(int),1,f);
		  fwrite(f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)(i+1)*nfaces[i]),f);
		  fwrite(&j,sizeof(int),1,f);
		}
	  
	    }
	}
  
      j=(1+ndims)*sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      fwrite(haveFaceFromVertex,sizeof(int),((size_t)ndims+1),f);
      fwrite(&j,sizeof(int),1,f);

      for (i=0;i<1+ndims;i++)
	{
      
	  if (haveFaceFromVertex[i])
	    {
	      j=sizeof(NDNET_IDCUMT)*((size_t)nvertex+1);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(v_numFaceIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)nvertex+1),f);
	      fwrite(&j,sizeof(int),1,f);
	  
	      j=sizeof(NDNET_UINT)*((size_t)v_numFaceIndexCum[i][nvertex]);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(v_faceIndex[i],sizeof(NDNET_UINT),((size_t)v_numFaceIndexCum[i][nvertex]),f);
	      fwrite(&j,sizeof(int),1,f);	  
	    }
	}
  
      j=(1+ndims)*(1+ndims)*sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      for (k=0;k<ndims+1;k++)
	fwrite(haveFaceFromFace[k],sizeof(int),((size_t)ndims+1),f);
      fwrite(&j,sizeof(int),1,f);

      for (i=0;i<1+ndims;i++)
	{
	  for (k=0;k<1+ndims;k++)
	    {

	      if (haveFaceFromFace[i][k])
		{
		  j=sizeof(NDNET_IDCUMT)*(1+(size_t)nfaces[i]);
		  fwrite(&j,sizeof(int),1,f);
		  fwrite(f_numFaceIndexCum[i][k],sizeof(NDNET_IDCUMT),((size_t)nfaces[i]+1),f);
		  fwrite(&j,sizeof(int),1,f);
	      
		  j=sizeof(NDNET_UINT)*((size_t)f_numFaceIndexCum[i][k][nfaces[i]]);
		  fwrite(&j,sizeof(int),1,f);
		  fwrite(f_faceIndex[i][k],sizeof(NDNET_UINT),((size_t)f_numFaceIndexCum[i][k][nfaces[i]]),f);
		  fwrite(&j,sizeof(int),1,f);	  
		}
	    }
	}

      j=sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&haveVFlags,sizeof(int),1,f);
      fwrite(&j,sizeof(int),1,f);

      if (haveVFlags)
	{
	  j=sizeof(unsigned char)*nvertex;
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(v_flag,sizeof(unsigned char),(size_t)nvertex,f);
	  fwrite(&j,sizeof(int),1,f);
	}

  
      j=sizeof(int)*(ndims+1);
      fwrite(&j,sizeof(int),1,f);
      fwrite(haveFFlags,sizeof(int),(ndims+1),f);
      fwrite(&j,sizeof(int),1,f);
  
      for (i=0;i<1+ndims;i++)
	if (haveFFlags[i])
	  {
	    j=sizeof(unsigned char)*nfaces[i];
	    fwrite(&j,sizeof(int),1,f);
	    if (j) fwrite(f_flag[i],sizeof(unsigned char),(size_t)nfaces[i],f);
	    fwrite(&j,sizeof(int),1,f);
	  }

      j=sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&ndata,sizeof(int),1,f);
      fwrite(&j,sizeof(int),1,f);
  
      for (i=0;i<ndata;i++)
	{
      
	  j=sizeof(int)+255*sizeof(char);
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(&(data[i].type),sizeof(int),1,f);
	  fwrite(data[i].name,sizeof(char)*255,1,f);
	  fwrite(&j,sizeof(int),1,f);

	  if (data[i].type==0)
	    {
	      j=sizeof(double)*nvertex;
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(data[i].data,sizeof(double),(size_t)nvertex,f);
	      fwrite(&j,sizeof(int),1,f);
	    }
	  else
	    {
	      j=sizeof(double)*nfaces[data[i].type];
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(data[i].data,sizeof(double),(size_t)nfaces[data[i].type],f);
	      fwrite(&j,sizeof(int),1,f);
	    }  
	}

      j=sizeof(int);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&nsupData,sizeof(int),1,f);
      fwrite(&j,sizeof(int),1,f);
  
      for (i=0;i<nsupData;i++)
	{
      
	  j=2*sizeof(int)+2*255*sizeof(char);
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(&(supData[i].type),sizeof(int),1,f);
	  fwrite(supData[i].name,sizeof(char)*255,1,f);
	  fwrite(&(supData[i].datasize),sizeof(int),1,f);
	  fwrite(supData[i].datatype,sizeof(char)*255,1,f);
	  fwrite(&j,sizeof(int),1,f);

	  if (supData[i].type==0)
	    {
	      j=(size_t)supData[i].datasize*nvertex;
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(supData[i].data,(size_t)supData[i].datasize,(size_t)nvertex,f);
	      fwrite(&j,sizeof(int),1,f);
	    }
	  else
	    {
	      j=(size_t)supData[i].datasize*nfaces[supData[i].type];
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(supData[i].data,(size_t)supData[i].datasize,(size_t)nfaces[supData[i].type],f);
	      fwrite(&j,sizeof(int),1,f);
	    }  
	}

 

      fclose(f);
      if (verbose>1) printf (" done.\n");
      return 0;
    }


  private:
    void setNULL()
    {
      strcpy(comment,"NO COMMENT");
      periodicity=0;
      ndims=0;ndims_net=0;
      isSimpComplex=0;
      x0=NULL;delta=NULL;
      indexSize=0;cumIndexSize=0;
      floatSize=0;
      nvertex=0;
      v_coord=NULL;
      haveVertexFromFace=NULL;
      nfaces=NULL;
      f_numVertexIndexCum=NULL;
      f_vertexIndex=NULL;
      haveFaceFromVertex=NULL;
      v_numFaceIndexCum=NULL;
      v_faceIndex=NULL;
      haveFaceFromFace=NULL;
      f_numFaceIndexCum=NULL;
      f_faceIndex=NULL;
      nsupData=0;supData=NULL;
      ndata=0;data=NULL;
      haveVFlags=0;haveFFlags=NULL;
      v_flag=NULL;f_flag=NULL;

      floatSize=sizeof(NDNET_FLOAT);
      indexSize=sizeof(NDNET_UINT);
      cumIndexSize=sizeof(NDNET_IDCUMT);
    }
    
    template <class T>
    void create(int ndims_, int ndims_net_, T nCells[])
    {   
      long k;
      isSimpComplex=1;
      ndims=ndims_;
      ndims_net=ndims_net_;
      x0=(double*)malloc(ndims*sizeof(double));
      delta=(double*)malloc(ndims*sizeof(double));

      strcpy(comment,"");
      nvertex=nCells[0];

      v_coord=NULL;

      nfaces=(NDNET_UINT*)calloc(((size_t)ndims+1),sizeof(NDNET_UINT));
      for (k=1;k<=ndims_net;k++) nfaces[k]=nCells[k];
      haveVertexFromFace=(int*)calloc(ndims+1,sizeof(int));
      f_vertexIndex = (NDNET_UINT**)calloc((1+ndims),sizeof(NDNET_UINT*));
      f_numVertexIndexCum = (NDNET_IDCUMT**)calloc((1+ndims),sizeof(NDNET_IDCUMT*));

      for (k=1;k<=ndims_net;k++) haveVertexFromFace[k]=(nCells[k]>0)?1:0;

      haveFaceFromVertex=(int*)calloc(ndims+1,sizeof(int));
      v_faceIndex = (NDNET_UINT**)calloc((1+ndims),sizeof(NDNET_UINT*));
      v_numFaceIndexCum = (NDNET_IDCUMT**)calloc((1+ndims),sizeof(NDNET_IDCUMT*));

      haveFaceFromFace =(int**)calloc(ndims+1,sizeof(int*));
      f_faceIndex = (NDNET_UINT***)calloc(ndims+1,sizeof(NDNET_UINT**));
      f_numFaceIndexCum = (NDNET_IDCUMT***)calloc(ndims+1,sizeof(NDNET_IDCUMT**));

      for (k=0;k<ndims+1;k++)
	{
	  haveFaceFromFace[k]=(int*)calloc(((size_t)ndims+1),sizeof(int));
	  f_faceIndex[k] = (NDNET_UINT**)calloc((1+ndims),sizeof(NDNET_UINT*));
	  f_numFaceIndexCum[k] = (NDNET_IDCUMT**)calloc(1+ndims,sizeof(NDNET_IDCUMT*));
	}
    
      haveVFlags=0;
      v_flag=NULL;
    
      haveFFlags=(int*)calloc(ndims+1,sizeof(int));
      f_flag=(unsigned char **)calloc(ndims+1,sizeof(unsigned char *));

      ndata=nsupData=0;
      data=NULL;
      supData=NULL;

      indexSize=sizeof(NDNET_UINT);
      cumIndexSize=sizeof(NDNET_IDCUMT);
      floatSize=sizeof(NDNET_FLOAT);
    }


    void freeData()
    {     
      int i,j;
    
      free(x0);
      free(delta);
      free(nfaces);
      free(v_coord);
      if (haveVertexFromFace !=NULL)
	for (i=0;i<=ndims;i++)
	  {
	    if (haveVertexFromFace[i])
	      {
		free(f_vertexIndex[i]);
		if (!isSimpComplex)
		  free(f_numVertexIndexCum[i]);
	      }
	  }
      free(haveVertexFromFace);
      free(f_vertexIndex);
      free(f_numVertexIndexCum);

      if (haveFaceFromVertex!=NULL)
	for (i=0;i<=ndims;i++)
	  {
	    if (haveFaceFromVertex[i])
	      {
		free(v_faceIndex[i]);
		free(v_numFaceIndexCum[i]);
	      }
	  }
      free(haveFaceFromVertex);
      free(v_faceIndex);
      free(v_numFaceIndexCum);

      if (haveFaceFromFace!=NULL)
      for (i=0;i<=ndims;i++)
	{
	  if (haveFaceFromFace[i]!=NULL)
	  for (j=0;j<=ndims;j++)
	    {
	      if (haveFaceFromFace[i][j])
		{
		  free(f_faceIndex[i][j]);
		  free(f_numFaceIndexCum[i][j]);
		}
	    }
	  free(haveFaceFromFace[i]);
	  free(f_faceIndex[i]);
	  free(f_numFaceIndexCum[i]);
	}
      free(haveFaceFromFace);
      free(f_faceIndex);
      free(f_numFaceIndexCum);

      if (haveVFlags) free(v_flag);
      if (haveFFlags!=NULL)
	for (i=0;i<=ndims;i++)
	  if (haveFFlags[i]) free(f_flag[i]);
      free(f_flag);
      free(haveFFlags);

      if (data!=NULL)
	for (i=0;i<ndata;i++)
	  {
	    free(data[i].data);
	  }
      free(data);

      if (supData!=NULL)
	for (i=0;i<nsupData;i++)
	  {
	    free(supData[i].data);
	  }
      free(supData);
      setNULL();
      fileLoaded=false;
    }
   
  };
} // IO
#include "../internal/namespace.footer"
#endif

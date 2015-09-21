#ifndef __ND_NETWORK_HXX__
#define __ND_NETWORK_HXX__

// This is ugly because it is legacy code, sorry ;(

#include <stdlib.h>
#include <stdio.h>

#include "../../tools/types/types.hxx"

#include "../../internal/namespace.header"

#define NDNET_UINT                UINT
#define NDNET_INT                 INT
#define NDNET_FLOAT               FLOAT
#define NDNET_IDCUMT              long

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

  class NDnetwork
  {
  public:

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

    template <class T1,class T2,class T3,class T4,class T5>
    NDnetwork(T1 nDims, T2 nDimsW, T3 x0_[], T4 delta_[], T5 nCells[], 
	      const char *cmt, bool periodic)
    {
      //memset(&(*this),0,sizeof(NDnetwork));

      strcpy(comment,"NO COMMENT");
      periodicity=0;
      ndims=0;ndims_net=0;
      isSimpComplex=0;
      x0=NULL;delta=NULL;
      indexSize=0;cumIndexSize=0;
      floatSize=0;
      nvertex=0;
      v_coord=NULL;
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
    }

    ~NDnetwork()
    {
      freeData();
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

  private:
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

      for (i=0;i<=ndims;i++)
	{
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
      for (i=0;i<=ndims;i++)
	if (haveFFlags[i]) free(f_flag[i]);
      free(f_flag);
      free(haveFFlags);

      for (i=0;i<ndata;i++)
	{
	  free(data[i].data);
	}
      free(data);

      for (i=0;i<nsupData;i++)
	{
	  free(supData[i].data);
	}
      free(supData);
    }
   
  };
} // IO
#include "../../internal/namespace.footer"
#endif

#ifndef __SIMPLEX_FLAGS_DEFINES_HXX__
#define __SIMPLEX_FLAGS_DEFINES_HXX__

// Flags description :
//
// A BOUNDARY simplex does not have NNEI neighbors
// A SHARED simplex is not local (i.e. the original copy is stored remotely)
// A LOCAL simplex is a simplex whose orginal copy is stored locally
// A GHOST simplex is a shared simplex with all its neighbors information correctly 
//   set locally (neighbors may be local, ghost or shadow) . Initially, a GHOST simplex 
//   shares at least one vertex with a local (i.e. non-shared) simplex.
// A SHADOW simplex does not have all the information about its neighbors locally.
//   Initially, a SHADOW simplex shares NO vertex with any shared simplex.
//   Note that shadow simplices only know their Ghost neighbors (not other shadows)
// TAGs are used to temporarily mark simplices

#define SIMPLEX_FLAG_NOTSET   (1<<0) // 1
#define SIMPLEX_FLAG_UNSAFE   (1<<1) // 2
#define SIMPLEX_FLAG_BOUNDARY (1<<2) // 4
#define SIMPLEX_FLAG_SHARED   (1<<3) // 8
#define SIMPLEX_FLAG_GHOST    (1<<4) // 16
#define SIMPLEX_FLAG_SHADOW   (1<<5) // 32

#define SIMPLEX_FLAG_TAG      (1<<7) // 128


#endif

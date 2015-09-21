#ifndef __VERTEX_FLAGS_DEFINES_HXX__
#define __VERTEX_FLAGS_DEFINES_HXX__

// Flags description :
//
// A BOUNDARY vertex is on the global boundary of the mesh (i.e. not between two processes)
// A SHARED vertex belongs to at least one shared simplex
// A LOCAL vertex is a vertex that belongs to at least one local simplex
// A GHOST vertex belongs to at least one ghost simplex, but no local simplex
// A SHADOW vertex belongs only to shadow simplices
// TAG and TAG2 are used to temporarily mark vertices

#define VERTEX_FLAG_NOTSET   (1<<0) // 1

#define VERTEX_FLAG_BOUNDARY (1<<2) // 4
#define VERTEX_FLAG_SHARED   (1<<3) // 8 
#define VERTEX_FLAG_GHOST    (1<<4) // 16
#define VERTEX_FLAG_SHADOW   (1<<5) // 32
#define VERTEX_FLAG_TAG      (1<<6) // 64
#define VERTEX_FLAG_TAG2     (1<<7) // 128

#endif

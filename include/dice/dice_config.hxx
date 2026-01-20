#ifndef __DICE_CONFIG_HXX__
#define __DICE_CONFIG_HXX__

// #define DICE_ALIGNED_MALLOC(A,N)  aligned_alloc(A,N)
#define DICE_ALIGNED_MALLOC(A, N) malloc(N)

// #define DICE_ALIGNAS(A) alignas( A )
#define DICE_ALIGNAS(A) \
    {                   \
    }

#define DICE_MESH_STATIC static
// #define DICE_MESH_STATIC {}

#endif

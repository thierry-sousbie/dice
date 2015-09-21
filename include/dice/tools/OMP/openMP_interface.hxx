#ifndef __OPENMP_INTERFACE_HXX__
#define __OPENMP_INTERFACE_HXX__

#ifdef USE_OPENMP

#include "../../tools/OMP/openmp.h"

#else

int omp_get_thread_num() {return 0;}
int omp_get_num_threads() {return 1;}
template <class T>
void omp_set_num_threads(T &n) {n=1;}
int omp_get_max_threads (void) {return 1;}
int omp_get_num_procs (void) {return 1;}

typedef char omp_lock_t;
void omp_init_lock(omp_lock_t *lock){}; 
void omp_destroy_lock(omp_lock_t *lock){};
void omp_set_lock(omp_lock_t *lock){}; 
void omp_unset_lock(omp_lock_t *lock){};

#endif //USE_OPENMP

#endif

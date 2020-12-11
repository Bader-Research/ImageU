#if (defined(SP2) || defined(CM5) || defined(MEIKO))

#include "ImageU.h"

inline
void all_gather(int* val) {
    int j;
    int n[PROCS]::;
    n[MYPROC]   = val[0];
    barrier();
    on_one 
	for (j=1 ; j<PROCS ; j++) 
	    val[j] = n[j];
    barrier();
}

inline
void all_gather_d(double* val) {
    int j;
    double n[PROCS]::;
    n[MYPROC]   = val[0];
    barrier();
    on_one 
	for (j=1 ; j<PROCS ; j++) 
	    val[j] = n[j];
    barrier();
}

inline
void all_concat(int* val) {
    int j, s;
    int n[PROCS]::;
    val[MYPROC] = val[0];
    n[MYPROC]   = val[0];
    barrier();
    for (j=1 ; j<PROCS ; j++) {
	s = (j + MYPROC) % PROCS;
	val[s] = n[s];
	barrier();
    }
}

#endif
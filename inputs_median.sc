/*                                                                    tab:8
 *
 * inputs_median.sc - Parallel Median Algorithm Data Sets
 *                    Parallel Sorting Algorithm Data Sets
 *
 * 
 * "Copyright (c) 1996 The Regents of the University of Maryland.
 * All rights reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without written agreement is
 * hereby granted, provided that the above copyright notice and the following
 * two paragraphs appear in all copies of this software.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF MARYLAND BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
 * OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 * MARYLAND HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE UNIVERSITY OF MARYLAND SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE UNIVERSITY OF MARYLAND HAS NO OBLIGATION TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS."
 *
 * Authors:             David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:             1.0
 * Creation Date:       February 6, 1996
 * Filename:            inputs_median.sc
 * History:
 */

#include "inputs_median.h"
#include <math.h>
#include "math_func.h"
#include "droute.h"
#include <limits.h>

int all_fill_array(int M, int A[PROCS]::[M], int N[PROCS]::,
		   char *filename)
{
    FILE *datafile;
    int i, j, n, tot;

    int arr[M];

    barrier();

    on_one {

	tot = 0;

	datafile = fopen(filename,"r");
	
	fscanf(datafile,"%d\n",&i);
	fscanf(datafile,"%d\n",&j);
	fprintf(outfile,"Reading input M(%d): saved with M: %d  PROCS: %d\n",
		M,i,j);
	for (i=0 ; i<PROCS ; i++) {
	    fscanf(datafile,"%d\n",&n);
	    N[i] = n;
	    tot += n;
	    for (j=0 ; j<n ; j++)
		fscanf(datafile,"%d ",&arr[j]);
	    bulk_write(&(A[i][0]), arr, n*sizeof(int));
	}

	fclose(datafile);
    }

    tot = all_bcast_i(tot);

    return(tot);
}


void fill_random(int M, int arr[M])
{
    register int
	i;

#if 1
    srandom(317*(MYPROC+17));
#endif
    
    for (i=0 ; i<M ; i++)
	arr[i] = random();
}

void all_fill_nas(int M, int A[PROCS]::[M])
{
    register int
	i, j;
    
    nas_srand();

    for (i=0 ; i<M*MYPROC ; i++) 
	j = nas_rand();

    for (i=0 ; i<M ; i++)
	A[MYPROC][i] = nas_rand();

    barrier();
}

void all_fill_cyclic(int M, int A[PROCS]::[M])
{
    register int
	i;
    
    for (i=0 ; i<M ; i++)
	A[MYPROC][i] = MYPROC + (i*PROCS);

    barrier();
}

void all_fill_pid(int M, int A[PROCS]::[M])
{
    register int
	i;
    
    for (i=0 ; i<M ; i++)
	A[MYPROC][i] = MYPROC;

    barrier();
}


void all_fill_zero(int M, int A[PROCS]::[M])
{
    register int
	i;
    
    for (i=0 ; i<M ; i++)
	A[MYPROC][i] = 0;

    barrier();
}


void all_fill_random(int M, int A[PROCS]::[M], int distrib)
/*
  distrib  entropy
     0         0
     1        31
     2        25.1
     3        16.9
     4        10.4
     5         6.2
*/
{
    register int
	i;
    
    RNG_init();

    switch (distrib) {
    case 0:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = 0;
	break;
    case 1:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = random();
	break;
    case 2:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = random() & random();
	break;
    case 3:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = random() & random() & random();
	break;
    case 4:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = random() & random() & random() & random();
	break;
    case 5:
	for (i=0 ; i<M ; i++)
	    A[MYPROC][i] = random() & random() & random() & random() & random();
	break;
    }

    barrier();
}

static inline
int hrel_lookup(int h, int M, int lastP, int *k, int idx)
{
    register int
	i,
	done;

    int addr;

    i    = 0;
    done = 0;
    while ((i<lastP)&&(!done)) 
	if (k[i] <= idx)
	    i++;
	else {
	    done = 1;
	    if (i==0)
		addr = idx;
	    else
		addr = i*h*M + idx - k[i-1];
	}
    if (!done) {
#if 0	
	addr = (PROCS-1)*h*M + MYPROC; /* <--- Might leave holes on last PE */
#else
	i = PROCS*M - k[lastP-1]; 
	i = MYPROC - PROCS + i;
	addr = (PROCS-1)*h*M + i;
#endif
    }

    return(addr);
}

void all_fill_hrel(int M, int A[PROCS]::[M], int h)
/* These are the destination addresses */
{
    register int
	i;
    int k[PROCS],
	lastP;

    barrier();

    if (h==1)
	all_fill_cyclic(M, A);
    else {
	lastP = (2 * PROCS) / h;
	for (i=0 ; i<lastP ; i++) 
	    k[i] = (int)floor( (double)(h*M) * (double)(2*PROCS - h*(i+1)) /
			       (double)(2*PROCS - h) );
	for (i=1 ; i<lastP ; i++)
	    k[i] += k[i-1];
	for (i=0 ; i<M ; i++) 
	    A[MYPROC][i] = hrel_lookup(h, M, lastP, k, PROCS*i + MYPROC +1);

    }

    barrier();

#if 0
    on_one {
	int j;
	for (j=0 ; j<M ; j++)
	    for (i=0 ; i<PROCS ; i++)
		fprintf(outfile,"Element [%3d][%5d]: %5d  (%3d, %5d)\n",
			i,j,A[i][j], A[i][j] / (h*M), A[i][j] % (h*M));
	fprintf(outfile,"\n");
	fflush(outfile);
    }
    barrier();
#endif
}


void all_fill_g_group_i(int M, int A[PROCS]::[M], int h, int g, int t)
{

  int
    i,j,k,
    p,
    offset,
    block;

  block = M/t;
  i = 0;

  barrier();
  on_one fprintf(stderr,"HERE 1xx\n");
  barrier();

  for (j=0 ; j<t ; j++) {
    p = (((PROCS/2 + j*g) % PROCS)^(MYPROC - (MYPROC % g)))*M*h + j*g*block;
    offset = (MYPROC % g) * block;
    for (k=0 ; k<(block>>1) ; k++) {
      A[MYPROC][i++] = p + offset + k;
      A[MYPROC][i++] = p + offset + block - k - 1;
    }
  }

  barrier();
  on_one fprintf(stderr,"HERE 2xx\n");
  barrier();
}

void fill_same(int M, int *arr)
{
    register int
	i;

    for (i=0 ; i<M ; i++)
	arr[i] = i;

}

void fill_linear(int M, int *arr)
{
    register int
	i;

    for (i=0 ; i<M ; i++)
	arr[i] = (MYPROC*M) + i;
}




#if (defined(T3D))
#define MAX_VAL (double)1.79769313486231e+308
#else
#define MAX_VAL (double)1.7976931348623158e+308
#endif
#define RANDOM_MAX      pow(2.0,31.0)

void fill_random_d(int M, double *arr) {
  int i,j;

  double multiplier, offset;

  srandom(23+1001*MYPROC); 

  offset     = RANDOM_MAX / 2.0;
  multiplier = MAX_VAL/offset;

  for (i=0 ; i<M ; i++)
      arr[i] = ((double)random() - offset) * multiplier;

#if 0
  barrier();
  on_one
    for (i=0 ; i<20 ; i++) {
      fprintf(outfile,"arr[%3d]: %g\n",i,arr[i]);
      fflush(outfile);
    }
  barrier();
#endif
#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
    
}


void fill_gaussian_d(int M, double *arr){
  int i,j;

  double  multiplier, offset;

  srandom(23+1001*MYPROC); 

  offset     = 2.0 * RANDOM_MAX;
  multiplier = MAX_VAL/offset;

  for (i=0 ; i<M ; i++) 
      arr[i] = ((double)random() + (double)random() +
	       (double)random() + (double)random() - offset) * multiplier;

#if 0
  barrier();
  on_one
    for (i=0 ; i<20 ; i++) {
      fprintf(outfile,"arr[%3d]: %g\n",i,arr[i]);
      fflush(outfile);
    }
  barrier();
#endif
#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
}


void fill_g_group_d(int M, double *arr, int group) {
    int
	i, j, k,
	v;

    double
	multiplier,
	t1, t2,
	offset;

    srandom(23+1001*MYPROC); 

    offset = RANDOM_MAX/2.0;
    multiplier = MAX_VAL/offset;
    t1 = RANDOM_MAX/(double)PROCS;
    v = ((MYPROC - (MYPROC % group)) + (PROCS/2)) % PROCS;

    k=0;
    for (i=0 ; i<group ; i++) {
	t2 = ((double)v * t1) - offset;
	for (j=0 ; j<(M/group) ; j++)
	    arr[k++] = (t2 + ((double)random() / (double)PROCS)) * multiplier;
	v = (++v) % PROCS; 
    } 

#if 0
  barrier();
  on_one
    for (i=0 ; i<20 ; i++) {
      fprintf(outfile,"arr[%3d]: %g\n",i,arr[i]);
      fflush(outfile);
    }
  barrier();
#endif
#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
}

void fill_bucket_d(int M, double *arr) {
  int
    i, j, k;

  double
      t1, t2,
      multiplier,
      offset;

  srandom(23+1001*MYPROC); 

  offset = RANDOM_MAX/2.0;
  multiplier = MAX_VAL/offset;
  t1 = RANDOM_MAX/(double)PROCS;

  k=0;
  for (i=0 ; i<PROCS ; i++) {
    t2 = (double)i * t1 - offset;
    for (j=0 ; j<(M/PROCS) ; j++)
      arr[k++] = (t2 + ((double)random() / (double)PROCS)) * multiplier;
  } 

#if 0
  barrier();
  on_one
    for (i=0 ; i<20 ; i++) {
      fprintf(outfile,"arr[%3d]: %g\n",i,arr[i]);
      fflush(outfile);
    }
  barrier();
#endif
#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
}

void fill_staggered_d(int M, double *arr) {
  int
    i, j, k,
    target;

  double
      multiplier,
      offset,
      t1;
  
  srandom(23+1001*MYPROC); 

  offset = RANDOM_MAX/2.0;
  multiplier = MAX_VAL/offset;

  if (MYPROC < (PROCS/2))
    target = 2*MYPROC + 1;
  else 
    target = (MYPROC - (PROCS/2))*2;

  t1 = (double)target * (RANDOM_MAX/(double)PROCS) - offset;

  k=0;
  for (i=0 ; i<M ; i++)
    arr[k++] = (t1 + ((double)random() / (double)PROCS)) * multiplier;

#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
#if 0
  barrier();
  on_one
    for (i=0 ; i<20 ; i++) {
      fprintf(outfile,"arr[%3d]: %g\n",i,arr[i]);
      fflush(outfile);
    }
  barrier();
#endif
}

void fill_zero_d(int M, double *arr){
  int i, j;

  for (i=0 ; i<M ; i++)
    arr[i] = 0.0;
#if 0
  barrier();
  for (j=0 ; j<PROCS ; j++) {
      on_proc(j) {
	  for (i=0 ; i<10 ; i++) {
	      fprintf(outfile,"PE %2d: arr[%3d]: %g\n",MYPROC,i,arr[i]);
	      fflush(outfile);
	  }
      }
      barrier();
  }
  barrier();
#endif
}


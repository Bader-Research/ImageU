/*                                                                    tab:8
 *
 * inputs_select.sc - Parallel Selection Algorithm Data Sets
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
 * Filename:            inputs_select.sc
 * History:
 */

#include "inputs_select.h"
#include <math.h>

#ifdef MEIKO
#define M_PI 3.1415926
#endif

static void all_print_n(int N[PROCS]::)
{
    register int i;

    barrier();
    
    on_one {
	for (i=0 ; i<PROCS ; i++)
	    fprintf(stdout,"N[%2d]: %12d\n",i,N[i]);
	fprintf(stdout,"\n");
    }
    barrier();
}

/*************************************************************/
int all_fill_array_control(int M, int A[PROCS]::[M], 
			   int N[PROCS]::, int total_n)
/*************************************************************/
{
    register int
	i, j;

    barrier();

    j = DIVPROCS(total_n);
    N[MYPROC] = j;

#if 0
    all_print_n(N);
#endif    

    if (j > M)
	fprintf(stderr,"ERROR: j: %d  M: %d\n",j,M);
    
    for (i=0 ; i<j ; i++)
	A[MYPROC][i] = random();

    return all_reduce_to_all_add(N[MYPROC]);
}

/*************************************************************/
int all_fill_array_linear(int M, int A[PROCS]::[M], 
			 int N[PROCS]::, int total_n)
/*************************************************************/
{
    register int
	i, j;

    barrier();

    j = (MYPROC * 2 * total_n) / (PROCS * PROCS);
    N[MYPROC] = j;

#if 0
    all_print_n(N);
#endif    

    if (j > M)
	fprintf(stderr,"ERROR: j: %d  M: %d\n",j,M);
    
    for (i=0 ; i<j ; i++)
	A[MYPROC][i] = random();

    return all_reduce_to_all_add(N[MYPROC]);
}

/*************************************************************/
int all_fill_array_exp(int M, int A[PROCS]::[M], 
			int N[PROCS]::, int total_n)
/*************************************************************/
{
    register int
	i, j;

    barrier();

    j = total_n / (1<<(MYPROC+1));
    on_proc(PROCS-1) 
	j = total_n / (1<<(PROCS-1));
	
    N[MYPROC] = j;

#if 0
    all_print_n(N);
#endif    

    if (j > M)
	fprintf(stderr,"ERROR: j: %d  M: %d\n",j,M);
    
    for (i=0 ; i<j ; i++)
	A[MYPROC][i] = random();

    return all_reduce_to_all_add(N[MYPROC]);
}

/*************************************************************/
int all_fill_array_one(int M, int A[PROCS]::[M], 
			int N[PROCS]::, int total_n)
/*************************************************************/
{
    register int
	i, j;

    barrier();

    j = 0;
#if 0
    on_proc(23 % PROCS)
#endif	
    on_one
	j = total_n;
	
    N[MYPROC] = j;

#if 0
    all_print_n(N);
#endif    

    if (j > M)
	fprintf(stderr,"ERROR: j: %d  M: %d\n",j,M);
    
    for (i=0 ; i<j ; i++)
	A[MYPROC][i] = random();

    return all_reduce_to_all_add(N[MYPROC]);
}

inline double fn_gaussian(double mean, double sd, double x) {

    return (1.0 / (sd * sqrt(2.0 * M_PI) ) *
	    exp(- pow((x-mean),2.0) / (2.0*sd*sd)));
}

#define GAUSS_END 3.0

/*************************************************************/
int all_fill_array_gauss(int M, int A[PROCS]::[M], 
			  int N[PROCS]::, int total_n)
/*************************************************************/
{
    register int
	i, j;

    double
	p,
	x_seg,
	x[PROCS],
	g[PROCS],
	tot;

    barrier();

    tot   = 0.0;
    p     = (double)PROCS;
    x_seg = 2.0 * GAUSS_END / p;

    for (i=0 ; i<PROCS ; i++) {
	x[i] = - GAUSS_END + (x_seg/2.0) + ((double)i * x_seg);
        g[i] = fn_gaussian(0.0, 1.0, x[i]);
	tot  += g[i];
    }

    j = (int)ceil((((double)total_n * g[MYPROC]) / tot));
    	
    N[MYPROC] = j;

#if 0
    all_print_n(N);
#endif    

    if (j > M)
	fprintf(stderr,"ERROR: j: %d  M: %d\n",j,M);
    
    for (i=0 ; i<j ; i++)
	A[MYPROC][i] = random();

    return all_reduce_to_all_add(N[MYPROC]);
}



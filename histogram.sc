/*									tab:8
 *
 * histogram.sc - parallel histogramming of images
 *
 * 
 * "Copyright (c) 1994,1995 The Regents of the University of Maryland.
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
 * Authors: 		David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:		1.0
 * Creation Date:	October 20, 1994
 * Filename:		histogram.sc
 * History:
 */

#include "histogram.h"

#define BITSPERINT 8*sizeof(int)

/***************************************************/
void all_histogram(int v, int w, int q, int r, 
		   int X[v][w]::[q][r], int k) 
/***************************************************/
/* Histogram the (vw x qr) = (n x n) image of "k" gray levels */
{
    int histo[PROCS]::[k];
    int A[PROCS]::[k*PROCS];

    register int
	i,
	j,
	rows,
	get_proc;

    int lastProc,
	(*locX)[r],
	*locHisto,
	*locA,
	proc_mask=PROCS-1;

    double secs;                   /* Timing marker */

    locX = tolocal(X[mapRow(v,w,MYPROC)][mapCol(v,w,MYPROC)]);
    locHisto = tolocal(histo[MYPROC]);
    locA = tolocal(A[MYPROC]);

    if (k>=PROCS) {
	rows = (int)(k/PROCS);     /* P divides k */
	on_one fprintf(outfile,"Rows per processor: %d\n",rows);
    }
    else rows = 1;                 /* k < PROCS */
  
/*****************************************************/
/* Step 0 Initialize Local Histogram Bars            */
/*****************************************************/

    for (i=0 ; i<k ; i++)
	locHisto[i] = 0;

    barrier();
    secs = get_seconds();

    all_init_timer();
    all_start_timer();
    
/*****************************************************/
/* Step 1 Histogram the local elements of Image X    */
/*****************************************************/

    for (i=0 ; i<q ; i++) 
	for (j=0 ; j<r ; j++) 
	    locHisto[locX[i][j]]++;
	
    barrier();                      /* Wait for everyone */
    /* histo is valid, p0histo is invalid */

    all_stop_timer("Step 1: local histogram");
    all_start_timer();

/***************************************************************/
/* Step 2 Transpose the P x k  histograms                      */
/***************************************************************/
    /* Transpose histo --> A */

#if (defined(CM5))
    for (i=0 ; i<PROCS ; i++) {
	if (MYPROC < k) {
	    /* Stagger prefetches using MYPROC */
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locA[rows*get_proc], &(histo[get_proc][rows*MYPROC]),
		      rows*sizeof(int));
	}
	barrier();
    }
#elif (defined(xSP2))
    if (PROCS >= k) {
	barrier();
	mpc_index(tolocal(histo), locA, rows*sizeof(int), ALLGRP);
    }
    else {
	if (MYPROC < k) {
	    for (i=0 ; i<PROCS ; i++) {
		/* Stagger prefetches using MYPROC */
		get_proc = (MYPROC + i) & proc_mask;
		bulk_get(&locA[rows*get_proc], &(histo[get_proc][rows*MYPROC]),
			 rows*sizeof(int));
	    }
	}
	sync();
	barrier();
    }
#else
    if (MYPROC < k) {
	for (i=0 ; i<PROCS ; i++) {
	    /* Stagger prefetches using MYPROC */
	    get_proc = (MYPROC + i) & proc_mask;
	    bulk_get(&locA[rows*get_proc], &(histo[get_proc][rows*MYPROC]),
		     rows*sizeof(int));
	}
    }
    sync();
    barrier();
#endif
	
    all_stop_timer("Step 2: Transpose");
    all_start_timer();

/*********************************************************/
/* Step 3 Reduce the "k/P" histograms                    */
/*********************************************************/

    if (MYPROC < k)
	for (i=1 ; i<PROCS ; i++) 
	    for (j=0 ; j<rows ; j++)
		locA[j] += locA[i*rows + j];

    barrier();
    /* A is valid */
    all_stop_timer("Step 3: Reduction");
    all_start_timer();

/*********************************************************/
/* Step 4 Collect all the histograms onto P0             */ 
/*********************************************************/

    on_proc(0) {
	lastProc = (k<PROCS) ? k : PROCS;
	for (i=1 ; i<lastProc ; i++)
	    bulk_read(&locA[rows*i], &(A[i][0]), rows*sizeof(int)); 
    }

    barrier();                     /* Wait for everyone */
    all_stop_timer("Step 4: P0 gets the entire histogram");
    
    secs = get_seconds() - secs;   /* Finish Timing */

/*****************************************************/
/* Results                                           */
/*****************************************************/

    on_proc(0) {
	fprintf(outfile,"Histogram Time: %f  for n x n = %d image on %d processors.\n",
	       secs, v*w*q*r, PROCS);
	fprintf(outfile,"  (%d colors)\n",k);
	if (RESULTS)
	    for (i=0 ; i<k ; i++)
		fprintf(outfile,"H[%3d]: %12d\n",i,locA[i]);
	fprintf(outfile,"\n");
	
    }

    all_print_timer(outfile);
    barrier(); 

}

/***************************************************/
int all_remap_histo(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int Y[v][w]::[q][r],
		    int k, int Histo_on_one[PROCS]::[k])
/***************************************************/
/* Remap the histogram of the (vw x qr) =
   (n x n) image of "k" gray levels.
   Return the new value of "k"  */
{

    register int
	i,
	j,
	done,
	bits;
    
    int (*locX)[r],
	(*locY)[r],
	histo_and,
	histo_or,
	logk,
	right_and,
	left_and,
	right_or,
	left_or,
	cut_left,
	cut_right,
	cut_tot;

    locX = tolocal(X[mapRow(v,w,MYPROC)][mapCol(v,w,MYPROC)]);
    locY = tolocal(Y[mapRow(v,w,MYPROC)][mapCol(v,w,MYPROC)]);

    logk = (int)log2(k);
    
    on_proc(0) {

	histo_and  = 0xFFFF;
	histo_or   = 0x0000;

	for (i=0 ; i<k ; i++) {
	    if (Histo_on_one[0][i] > 0) histo_and &= i;
	    if (Histo_on_one[0][i] > 0) histo_or  |= i;
	}

	i=1;
	right_and = 0;
	do {
	    bits = ~(~0<<i);
	    if ((histo_and & bits) == bits) {
		right_and = i;
		done = 0;
	    }
	    else
		done = 1;
	} while ((!done)&&(i<=logk));
	
	i=1;
	right_or = 0;
	do {
	    bits = ~(~0<<i);
	    if ((histo_or & bits) == 0) {
		right_or = i;
		done = 0;
	    }
	    else
		done = 1;
	} while ((!done)&&(i<=logk));

	cut_right = min(right_and, right_or);

/* For logk=8, what about 00000000 00000000 00000000 11111111 */
/* shift left string by BITSPERINT-logk. Then perform check. */

	i=1;
	left_and = 0;
	histo_and<<(BITSPERINT-logk);
	do {
	    bits = ~(~0>>i);
	    if ((histo_and & bits) == bits) {
		left_and = i;
		done = 0;
	    }
	    else
		done = 1;
	} while ((!done)&&(i<=logk));

/* For logk=8, what about 00000000 00000000 00000000 00011111 */

	i=1;
	left_or = 0;
	histo_or<<(BITSPERINT-logk);
	do {
	    bits = ~(~0>>i);
	    if ((histo_or & bits) == 0) {
		left_or = i;
		done = 0;
	    }
	    else
		done = 1;
	} while ((!done)&&(i<=logk));
	
	cut_left = max(left_and, left_or);

#if 1
	fprintf(outfile,"cut_left: %d  cut_right: %d \n",
		cut_left,cut_right);
#endif	

    }

    barrier();

    cut_left  = all_bcast_i(cut_left);
    cut_right = all_bcast_i(cut_right);

    cut_tot = cut_left + cut_right;

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {
	    if (cut_tot > 0) {
		/* zero out left bits */
		bits = locX[i][j] << (BITSPERINT-logk); /* leading bits */
		bits = bits & (~0>>cut_left);
		bits = bits >> (BITSPERINT-logk);       /* return back */
		/* shift right */
		locY[i][j] = bits>>cut_right;
	    }
	    else 
		locY[i][j] = locX[i][j];
	}

    barrier();

    return(k>>cut_tot);
}

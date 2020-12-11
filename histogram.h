/*									tab:8
 *
 * histogram.h - parallel histogramming of images
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
 * Filename:		histogram.h
 * History:
 */

#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

#include "ImageU.h"

void all_histogram(int v, int w, int q, int r, 
		   int X[v][w]::[q][r], int k);

int all_remap_histo(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int Y[v][w]::[q][r],
		    int k, int Histo_on_one[PROCS]::[k]);

/***************************************************/
_INLINE void all_histogram_sub(int v, int w, int q, int r, 
			       int X[v][w]::[q][r], int k,
			       int Histo_on_one[PROCS]::[k])
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
	*locA;

    locX = tolocal(X[mapRow(v,w,MYPROC)][mapCol(v,w,MYPROC)]);
    locHisto = tolocal(histo[MYPROC]);
    locA = tolocal(A[MYPROC]);

    if (k>=PROCS) {
	rows = (int)(k/PROCS);     /* P divides k */
    }
    else rows = 1;                 /* k < PROCS */
  
/*****************************************************/
/* Step 0 Initialize Local Histogram Bars            */
/*****************************************************/

    for (i=0 ; i<k ; i++)
	locHisto[i] = 0;

/*****************************************************/
/* Step 1 Histogram the local elements of Image X    */
/*****************************************************/

    for (i=0 ; i<q ; i++) 
	for (j=0 ; j<r ; j++) 
	    locHisto[locX[i][j]]++;

    barrier();                      /* Wait for everyone */
    /* histo is valid, p0histo is invalid */

/***************************************************************/
/* Step 2 Transpose the P x k  histograms                      */
/***************************************************************/
    /* Transpose histo --> A */

    for (i=0 ; i<PROCS ; i++) {
	if (MYPROC < k) {
	    /* Stagger prefetches using MYPROC */
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locA[rows*get_proc], &(histo[get_proc][rows*MYPROC]),
		      rows*sizeof(int));
	}
	barrier();
    }
	
/*********************************************************/
/* Step 3 Reduce the "k/P" histograms                    */
/*********************************************************/

    if (MYPROC < k)
	for (i=1 ; i<PROCS ; i++) 
	    for (j=0 ; j<rows ; j++)
		locA[j] += locA[i*rows + j];

    barrier();
    /* A is valid */

/*********************************************************/
/* Step 4 Collect all the histograms onto P0             */ 
/*********************************************************/

    on_proc(0) {
	lastProc = (k<PROCS) ? k : PROCS;
	for (i=1 ; i<lastProc ; i++)
	    bulk_read(&locA[rows*i], &(A[i][0]), rows*sizeof(int)); 
	for (i=0 ; i<k ; i++)
	    Histo_on_one[0][i] = locA[i];
    }

    barrier();                     /* Wait for everyone */

}


#endif

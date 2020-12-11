/*									tab:8
 *
 * columnsort.sc - parallel columnsort
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
 * Filename:		columnsort.sc
 * History:
 */

#include "columnsort.h"

void all_columnsort_int(int t, int A[PROCS]::[t], int B[PROCS]::[t]) 
{
    register int
	i,
	j,
	s,         /* s = t/p = n/(p^2) */
	elem,
	group,
	halft,
	get_proc;
	
    int tempElem,
	*locA,
	*locB;

    if (t < 2*(PROCS-1)*(PROCS-1)) {
	on_one fprintf(stderr,
		       "ERROR: ColumnSort t < 2(P-1)^2   t=%d  P=%d\n",
		       t,PROCS);
        barrier();
	exit(1);
    };

    locA = tolocal(A);
    locB = tolocal(B);

    s = (int)(t/PROCS);

    halft = (int)(t/2);

/*****************************************************/
/* Step 1 Local (QuickSort) Sort The Columns of A    */
/*****************************************************/

    fastsort(locA, t);

/*****************************************************/
/* Step 2 Unshuffle and Transpose                    */
/*****************************************************/

  /* Unshuffle A -> B */

    i = 0;
    for (elem=0 ; elem<PROCS ; elem++)
	for (group=0 ; group<s ; group++)
	    locB[i++] = locA[(PROCS*group) + elem];

  /* Transpose B -> A */

    barrier();			/* Wait for everyone */
  /* B is valid, A is invalid */

    for (i=0 ; i<PROCS ; i++) {
        get_proc = (MYPROC + i) % PROCS; /* Stagger prefetches using MYPROC */
        bulk_get(&locA[s*get_proc], &(B[get_proc][s*MYPROC]), s*sizeof(int)); 
    }

    /* Local computation can go here */

    sync();                /* Wait for all prefetches to complete */

/*****************************************************/
/* Step 3 Local (QuickSort) Sort The Columns of A    */
/*****************************************************/

    fastsort(locA, t);

/*****************************************************/
/* Step 4 Untranspose                                */
/*****************************************************/

  /* Transpose A -> B */

    barrier();			/* Wait for everyone */
  /* A is valid, B is invalid */

    for (i=0 ; i<PROCS ; i++) {
        get_proc = (MYPROC + i) % PROCS; /* Stagger prefetches using MYPROC */
        bulk_get(&locB[s*get_proc], &(A[get_proc][s*MYPROC]), s*sizeof(int)); 
    }

  /* Local computation can go here */

    sync();                /* Wait for all prefetches to complete */

    barrier();			/* Wait for everyone */
  /* B is valid, A is invalid */

/*****************************************************/
/* Step 5 Local (QuickSort) Sort The Columns of B    */
/*****************************************************/

    fastsort(locB, t);

/*****************************************************/
/* Step 6 Shift                                      */
/*****************************************************/

    barrier();			/* Wait for everyone */
  /* B is valid, A is invalid */

  /* Shift B -> A */
    bulk_get(locA, &(B[(MYPROC-1+PROCS)%PROCS][halft]), halft*sizeof(int));
  
  /* Locally shift second half */
    for (i=0 ; i < halft ; i++)
	locA[i+halft] = locB[i];

    sync();

/*********************************************************/
/* Step 7 Local (Linear Insertion) Sort The Columns of A */
/*********************************************************/

  /* Sort A */
  /* Both top and bottom halves of Data 
     on Processor 0 are already sorted */

    if (MYPROC > 0) 
	for (i=1 ; i<t ; i++) 
	    if(locA[i] < locA[i-1]) {
		tempElem = locA[i];
		for (j=i-1 ; (j>=0)&&(tempElem<locA[j]) ; j--) 
		    locA[j+1] = locA[j];
		locA[j+1] = tempElem;
	    };

/*****************************************************/
/* Step 8 UnShift                                    */
/*****************************************************/

    barrier();			/* Wait for everyone */
  /* A is valid, B is invalid */

  /* UnShift A -> B */
    bulk_get(&locB[halft], &(A[(MYPROC+1)%PROCS][0]), halft*sizeof(int));

  /* Locally unshift second half */
    for (i=0 ; i < halft ; i++)
	locB[i] = locA[i+halft];

    sync();

    barrier();			/* Wait for everyone */
  /* B is valid, A is invalid */

/*****************************************************/
/* VERIFY                                            */
/*****************************************************/

#if 0
    for (i=1 ; i<t ; i++)
	if (locB[i-1] > locB[i])
	    fprintf(outfile,"ERROR: B[%2d][%4d]: %d  B[%2d][%4d]: %d \n",
		   MYPROC,i-1,locB[i-1],MYPROC,i,locB[i]);

    on_one
	for (i=1 ; i<PROCS ; i++)
	    if (B[i-1][t-1] > B[i][0])
		fprintf(outfile,"ERROR: B[%2d][%4d]: %d  B[%2d][%4d]: %d \n",
		       i-1,t-1,B[i-1][t-1],i,0,B[i][0]);
    barrier(); 
#endif

}


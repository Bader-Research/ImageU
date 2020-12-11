/*									tab:8
 *
 * ImageU.h - Image Understanding Project
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
 * Version:		2
 * Creation Date:	October 20, 1994
 * Filename:		ImageU.h
 * History:
 *	DAB	2	Sun Dec  4 14:31:30 EST 1994
 *		Added in DO_CCGRYL
 */

#ifndef _IMAGEU_H
#define _IMAGEU_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <split-c/split-c.h>

#define VER_MAJOR 3
#define VER_MINOR 6

#define _INLINE extern inline 

#if (defined (SP1) || defined (SP2))
#define malloc(size) sc_malloc_wrapper(size)
#define free(ptr) sc_free_wrapper(ptr)
#define realloc(ptr, size) sc_realloc_wrapper(ptr, size)
#endif

#if (!defined (CM5))
#define log2(d) (log(d) / log(2.0))
#endif

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef CAT
#define CAT(a,b) a##b
#endif

#include "image_print.h"

#define DARPA_IMAGE "IUbench.im"

#define	errprnt(msg)	{ fprintf(stderr,"connComp: %s\n",msg); exit(1); }

#define IMAGE_DEFAULT 1

#define N_DEFAULT 32
#define K_DEFAULT 2

#define BARSIZE_DEFAULT 5

#define WATCH_PROC_DEFAULT 32

FILE *outfile;

int watch_proc,                 /* Processor to debug */
    OUTPUT_IMAGE,               /* Output Image Status Bits */
    UNSUPPARAM,                 /* Parameter for all_unsupported */
    DO_TEST,                    /* Test machine specifics */
    DO_UNSUP,                   /* Unsupported testing routines */
    DO_HISTO,                   /* Do a histogram? */
    DO_CCBIN,                   /* Do a binary 8-connected component? */
    DO_CCBIF,                   /* Do a binary 4-connected component? */
    DO_CCGRY,                   /* Do a grey-level 8-connected component? */
    DO_CCGRF,                   /* Do a grey-level 4-connected component? */
    DO_CCGRYL,                  /* Do a l-grey-level 8-connected component? */
    DO_SEGM,                    /* Do Image Segmentation */
    DO_GBLUR,                   /* Do Blur of image */
    DO_BENCHMARK,               /* Perform the benchmark suite? */
    RESULTS;                    /* Print verbose results? */

#define MAXLEN    80

#define MAX_TIMER 256

double timer[MAX_TIMER];
int    timer_cnt;
char   timer_msg[MAX_TIMER][MAXLEN];

_INLINE void all_init_timer() {
    barrier();
    timer_cnt = 0;
}

_INLINE void all_start_timer() {
    barrier();
    if (timer_cnt >= MAX_TIMER) {
	fprintf(stderr,"ERROR: Ran out of timers.\n");
	exit(1);
    }
    timer[timer_cnt] = get_seconds();
}

_INLINE void all_stop_timer(char* msg) {
    barrier();
    timer[timer_cnt] = get_seconds() - timer[timer_cnt];
    strcpy(timer_msg[timer_cnt], msg);
    timer_cnt++;
}

_INLINE void all_print_timer(FILE* strm) {
    register int i;
    
    on_one{
	fprintf(strm,"\n");
	fprintf(strm,"Timers:\n");
	for (i=0 ; i<timer_cnt ; i++) 
	    fprintf(strm,"%9.6f  for %s\n",timer[i],timer_msg[i]);
	fprintf(strm,"\n");
    }
}

_INLINE void all_print_timer_PS(FILE* strm) {
    register int i;
    
    on_one{
	fprintf(strm,"\n");
	fprintf(strm,"Prefix-Sums of Timers:\n");
	for (i=1 ; i<timer_cnt ; i++) 
	    timer[i] += timer[i-1];
	for (i=0 ; i<timer_cnt ; i++) 
	    fprintf(strm,"%9.6f  for %s\n",timer[i],timer_msg[i]);
	fprintf(strm,"\n");
    }
}
	    
_INLINE int mapRow(int rowsGrid, int colsGrid, int proc_id) {
/* Return the row of processor proc_id in a rowsGrid x colsGrid
   processor 2d array */
    return(proc_id / colsGrid);
}
   
_INLINE int mapCol(int rowsGrid, int colsGrid, int proc_id) {
/* Return the column of processor proc_id in a rowsGrid x colsGrid
   processor 2d array */
    return(proc_id % colsGrid);
}

_INLINE int mapProc(int rowsGrid, int colsGrid, int aRow, int aCol) {
/* Return the proc_id of the processor in a rowsGrid x colsGrid
   processor 2d array at location (aRow, aCol) */
    return(aRow * colsGrid + aCol);
}

_INLINE void assert_spread_malloc(void *spread spreadPtr)
/* Check spreadPtr to make sure it is not NULL */
{
    if (spreadPtr == NULL) {
	fprintf(stderr,"ERROR: a SpreadPtr is NULL\n");
	fflush(stderr);
	exit(1);
    }
}

_INLINE void assert_malloc(void *ptr) {
    if (ptr==NULL) {
	fprintf(stderr,"ERROR: PE%2d cannot malloc\n",MYPROC);
	fflush(stderr);
	exit(1);
    }
}

#define MODPROCS(a) ((a) & (PROCS-1))
#define DIVPROCS(a) ((a) >> PROCSLOG)

#define WHITE 0
#define BLACK 1

#define NIL   -1

typedef struct edge
{
    int a;
    int b;
} edge_t;

typedef struct elem *ELEMPTR;

typedef struct elem
{
    int color;
    int label;
    int pos;
} elem_t;

typedef struct chPair *CHPAIRPTR;

typedef struct chPair
{
    int alpha;
    int beta;
} chPair_t;

void all_get_image(int v, int w, int q, int r,
		   int X[v][w]::[q][r], int k, int);

void all_get_image_file(int v, int w, int q, int r,
			int X[v][w]::[q][r], char*);

#endif

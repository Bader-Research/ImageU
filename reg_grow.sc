/*									tab:8
 *
 * reg_grow.sc - region growing algorithm
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
 * Creation Date:	December 29, 1994
 * Filename:		reg_grow.sc
 * History:
 */

#include "reg_grow.h"
#include "math_func.h"
#include "select_b.h"
#include "columnsort.h"
#include "radixsort.h"

#define border_pix(ii,jj,v,w,q,r,i,j) \
(                               \
   ((ii==0)   && ((i)==0))   || \
   ((jj==0)   && ((j)==0))   || \
   ((jj==w-1) && ((j)==r-1)) || \
   ((ii==v-1) && ((i)==q-1))    \
) 


#if 0
/*************************************************************/
static inline void all_get_ghost_four(int ii, int jj,
			int v, int w, int q, int r, int globArr[v][w]::[q][r] , 
			int *ghostN , int *ghostS , int *ghostE , int *ghostW ,
			int nil_elem)
/*************************************************************/
{
    register int j;
    
/* Create Ghost Cells */
    
    barrier();
    
    if (ii>0) 
	bulk_get(ghostN, &(globArr[ii-1][jj][q-1][0]), r*sizeof(int));
    else
	for (j=0 ; j<r ; j++) ghostN[j] = nil_elem;
    sync();
    
    if (ii<v-1) 
	bulk_get(ghostS, &(globArr[ii+1][jj][0][0]), r*sizeof(int));
    else
	for (j=0 ; j<r ; j++) ghostS[j] = nil_elem;
    sync();
    
    if (jj>0) 
	for (j=0 ; j<q ; j++) ghostW[j] := globArr[ii][jj-1][j][r-1];
    else
	for (j=0 ; j<q ; j++) ghostW[j] = nil_elem;
    sync();
    
    if (jj<w-1) 
	for (j=0 ; j<q ; j++) ghostE[j] := globArr[ii][jj+1][j][0];
    else
	for (j=0 ; j<q ; j++) ghostE[j] = nil_elem;
    sync();
    
    barrier();
}
#endif

#define all_get_ghost_four(ii, jj, v, w, q, r, globArr, \
			   ghostN, ghostS, ghostE, ghostW, \
			   nil_elem) \
{ \
    register int j; \
     \
    barrier(); \
     \
    if (ii>0)  \
	bulk_get(ghostN, &(globArr[ii-1][jj][q-1][0]), r*sizeof(int)); \
    else \
	for (j=0 ; j<r ; j++) ghostN[j] = nil_elem; \
    sync(); \
     \
    if (ii<v-1)  \
	bulk_get(ghostS, &(globArr[ii+1][jj][0][0]), r*sizeof(int)); \
    else \
	for (j=0 ; j<r ; j++) ghostS[j] = nil_elem; \
    sync(); \
     \
    if (jj>0)  \
	for (j=0 ; j<q ; j++) ghostW[j] := globArr[ii][jj-1][j][r-1]; \
    else \
	for (j=0 ; j<q ; j++) ghostW[j] = nil_elem; \
    sync(); \
     \
    if (jj<w-1)  \
	for (j=0 ; j<q ; j++) ghostE[j] := globArr[ii][jj+1][j][0]; \
    else \
	for (j=0 ; j<q ; j++) ghostE[j] = nil_elem; \
    sync(); \
     \
    barrier(); \
} 

#if 0
/*************************************************************/
static inline void all_get_ghost_eight(int ii, int jj,
			 int v, int w, int q, int r, int globArr[v][w]::[q][r] , 
			 int *ghostN , int *ghostS , int *ghostE , int *ghostW ,
			 int *ghostNW , int *ghostNE , int *ghostSW , int *ghostSE ,
			 int nil_elem)
/*************************************************************/
{
    
/* Create Ghost Cells */
    
    all_get_ghost_four(ii,jj,v,w,q,r,globArr,
		       ghostN, ghostS, ghostE, ghostW,
		       nil_elem);
    
    *ghostNW = nil_elem;
    *ghostNE = nil_elem;
    *ghostSE = nil_elem;
    *ghostSW = nil_elem;
    
    if (ii>0) {
	if (jj>0) *ghostNW := globArr[ii-1][jj-1][q-1][r-1];
	if (jj<w-1) *ghostNE := globArr[ii-1][jj+1][q-1][0];
    };
    
    if (ii<v-1) {
	if (jj>0) *ghostSW := globArr[ii+1][jj-1][0][r-1];
	if (jj<w-1) *ghostSE := globArr[ii+1][jj+1][0][0];
    };
    
    sync();
    
    barrier();
}
#endif

#define all_get_ghost_eight(ii, jj, v, w, q, r, globArr,  \
			 ghostN, ghostS, ghostE, ghostW, \
			 ghostNW, ghostNE, ghostSW, ghostSE, \
			 nil_elem) \
{ \
     \
    all_get_ghost_four(ii,jj,v,w,q,r,globArr, \
		       ghostN, ghostS, ghostE, ghostW, \
		       nil_elem); \
     \
    *ghostNW = nil_elem; \
    *ghostNE = nil_elem; \
    *ghostSE = nil_elem; \
    *ghostSW = nil_elem; \
     \
    if (ii>0) { \
	if (jj>0) *ghostNW := globArr[ii-1][jj-1][q-1][r-1]; \
	if (jj<w-1) *ghostNE := globArr[ii-1][jj+1][q-1][0]; \
    }; \
     \
    if (ii<v-1) { \
	if (jj>0) *ghostSW := globArr[ii+1][jj-1][0][r-1]; \
	if (jj<w-1) *ghostSE := globArr[ii+1][jj+1][0][0]; \
    }; \
     \
    sync(); \
     \
    barrier(); \
} 
 

#if 0
/*************************************************************/
static inline void my_proc_nbr(int q, int r,
		 int tile[q][r], int pType[q][r],
		 int *ghostN, int *ghostS, int *ghostE, int *ghostW,
		 int ghostNW, int ghostNE, int ghostSW, int ghostSE,
		 int nil_elem,
		 register int i, register int j,
		 int *mask, int *nbrs) 
/*************************************************************/
{
/* N  */
	    switch (pType[i][j]) {
	    case 1:
	    case 2:
	    case 3:
		if (ghostN[j] == nil_elem)
		    mask[0] = 0;
		else 
		    nbrs[0] = ghostN[j];
		break;
	    default: nbrs[0] = tile[i-1][j];
	    }
/* S  */
	    switch (pType[i][j]) {
	    case 7:
	    case 8:
	    case 9:
		if (ghostS[j] == nil_elem)
		    mask[1] = 0;
		else
		    nbrs[1] = ghostS[j];
		break;
	    default: nbrs[1] = tile[i+1][j];
	    }
/* W  */
	    switch (pType[i][j]) {
	    case 1:
	    case 4:
	    case 7:
		if (ghostW[i] == nil_elem)
		    mask[2] = 0;
		else
		    nbrs[2] = ghostW[i];
		break;
	    default: nbrs[2] = tile[i][j-1];
	    }
/* E  */
	    switch (pType[i][j]) {
	    case 3:
	    case 6:
	    case 9:
		if (ghostE[i] == nil_elem)
		    mask[3] = 0;
		else
		    nbrs[3] = ghostE[i];
		break;
	    default: nbrs[3] = tile[i][j+1];
	    }
/* NW */
	    switch (pType[i][j]) {
	    case 1:
		if (ghostNW == nil_elem)
		    mask[4] = 0;
		else
		    nbrs[4] = ghostNW;
		break;
	    case 2:
	    case 3:
		if (ghostN[j-1] == nil_elem)
		    mask[4] = 0;
		else
		    nbrs[4] = ghostN[j-1];
		break;
	    case 4:
	    case 7:
		if (ghostW[i-1] == nil_elem)
		    mask[4] = 0;
		else
		    nbrs[4] = ghostW[i-1];
		break;
	    default: nbrs[4] = tile[i-1][j-1];
	    }
/* NE */
	    switch (pType[i][j]) {
	    case 3:
		if (ghostNE == nil_elem)
		    mask[5] = 0;
		else
		    nbrs[5] = ghostNE;
		break;
	    case 1:
	    case 2:
		if (ghostN[j+1] == nil_elem)
		    mask[5] = 0;
		else
		    nbrs[5] = ghostN[j+1];
		break;
	    case 6:
	    case 9:
		if (ghostE[i-1] == nil_elem)
		    mask[5] = 0;
		else
		    nbrs[5] = ghostE[i-1];
		break;
	    default: nbrs[5] = tile[i-1][j+1];
	    }
/* SW */
	    switch (pType[i][j]) {
	    case 7:
		if (ghostSW == nil_elem)
		    mask[6] = 0;
		else
		    nbrs[6] = ghostSW;
		break;
	    case 8:
	    case 9:
		if (ghostS[j-1] == nil_elem)
		    mask[6] = 0;
		else
		    nbrs[6] = ghostS[j-1];
		break;
	    case 1:
	    case 4:
		if (ghostW[i+1] == nil_elem)
		    mask[6] = 0;
		else
		    nbrs[6] = ghostW[i+1];
		break;
	    default: nbrs[6] = tile[i+1][j-1];
	    }
/* SE */
	    switch (pType[i][j]) {
	    case 9:
		if (ghostSE == nil_elem)
		    mask[7] = 0;
		else
		    nbrs[7] = ghostSE;
		break;
	    case 7:
	    case 8:
		if (ghostS[j+1] == nil_elem)
		    mask[7] = 0;
		else
		    nbrs[7] = ghostS[j+1];
		break;
	    case 3:
	    case 6:
		if (ghostE[i+1] == nil_elem)
		    mask[7] = 0;
		else
		    nbrs[7] = ghostE[i+1];
		break;
	    default: nbrs[7] = tile[i+1][j+1];
	    }
	    
}
#endif

#define my_proc_nbr(q, r, tile, pType, \
		    ghostN, ghostS, ghostE, ghostW, \
		    ghostNW, ghostNE, ghostSW, ghostSE, \
		    nil_elem, i, j, mask, nbrs)  \
\
	    switch (pType[i][j]) { \
	    case 1: \
	    case 2: \
	    case 3: \
		if (ghostN[j] == nil_elem) \
		    mask[0] = 0; \
		else \
		    nbrs[0] = ghostN[j]; \
		break; \
	    default: nbrs[0] = tile[i-1][j]; \
	    } \
	    switch (pType[i][j]) { \
	    case 7: \
	    case 8: \
	    case 9: \
		if (ghostS[j] == nil_elem) \
		    mask[1] = 0; \
		else \
		    nbrs[1] = ghostS[j]; \
		break; \
	    default: nbrs[1] = tile[i+1][j]; \
	    } \
	    switch (pType[i][j]) { \
	    case 1: \
	    case 4: \
	    case 7: \
		if (ghostW[i] == nil_elem) \
		    mask[2] = 0; \
		else \
		    nbrs[2] = ghostW[i]; \
		break; \
	    default: nbrs[2] = tile[i][j-1]; \
	    } \
	    switch (pType[i][j]) { \
	    case 3: \
	    case 6: \
	    case 9: \
		if (ghostE[i] == nil_elem) \
		    mask[3] = 0; \
		else \
		    nbrs[3] = ghostE[i]; \
		break; \
	    default: nbrs[3] = tile[i][j+1]; \
	    } \
	    switch (pType[i][j]) { \
	    case 1: \
		if (ghostNW == nil_elem) \
		    mask[4] = 0; \
		else \
		    nbrs[4] = ghostNW; \
		break; \
	    case 2: \
	    case 3: \
		if (ghostN[j-1] == nil_elem) \
		    mask[4] = 0; \
		else \
		    nbrs[4] = ghostN[j-1]; \
		break; \
	    case 4: \
	    case 7: \
		if (ghostW[i-1] == nil_elem) \
		    mask[4] = 0; \
		else \
		    nbrs[4] = ghostW[i-1]; \
		break; \
	    default: nbrs[4] = tile[i-1][j-1]; \
	    } \
	    switch (pType[i][j]) { \
	    case 3: \
		if (ghostNE == nil_elem) \
		    mask[5] = 0; \
		else \
		    nbrs[5] = ghostNE; \
		break; \
	    case 1: \
	    case 2: \
		if (ghostN[j+1] == nil_elem) \
		    mask[5] = 0; \
		else \
		    nbrs[5] = ghostN[j+1]; \
		break; \
	    case 6: \
	    case 9: \
		if (ghostE[i-1] == nil_elem) \
		    mask[5] = 0; \
		else \
		    nbrs[5] = ghostE[i-1]; \
		break; \
	    default: nbrs[5] = tile[i-1][j+1]; \
	    } \
	    switch (pType[i][j]) { \
	    case 7: \
		if (ghostSW == nil_elem) \
		    mask[6] = 0; \
		else \
		    nbrs[6] = ghostSW; \
		break; \
	    case 8: \
	    case 9: \
		if (ghostS[j-1] == nil_elem) \
		    mask[6] = 0; \
		else \
		    nbrs[6] = ghostS[j-1]; \
		break; \
	    case 1: \
	    case 4: \
		if (ghostW[i+1] == nil_elem) \
		    mask[6] = 0; \
		else \
		    nbrs[6] = ghostW[i+1]; \
		break; \
	    default: nbrs[6] = tile[i+1][j-1]; \
	    } \
	    switch (pType[i][j]) { \
	    case 9: \
		if (ghostSE == nil_elem) \
		    mask[7] = 0; \
		else \
		    nbrs[7] = ghostSE; \
		break; \
	    case 7: \
	    case 8: \
		if (ghostS[j+1] == nil_elem) \
		    mask[7] = 0; \
		else \
		    nbrs[7] = ghostS[j+1]; \
		break; \
	    case 3: \
	    case 6: \
		if (ghostE[i+1] == nil_elem) \
		    mask[7] = 0; \
		else \
		    nbrs[7] = ghostE[i+1]; \
		break; \
	    default: nbrs[7] = tile[i+1][j+1]; \
	    } 


/*************************************************************/
void all_get_pType(int v, int w, int q, int r, 
		   int pType[v][w]::[q][r]) 
/*************************************************************/
/* Fill pType with [1,9] for each pixel:

   1 | 2 | 3
   -----------
   4 | 5 | 6
   -----------
   7 | 8 | 9
*/
{
    register int
	i,
	j,
	ii,
	jj;

    int	(*locP)[r];

    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locP = tolocal(pType[ii][jj]);

    all_init_timer();

    all_start_timer();

    for (i=1 ; i<q-1 ; i++)
	for (j=1 ; j<r-1 ; j++)
	    locP[i][j] = 5;
    for (j=1 ; j<r-1 ; j++) {
	locP[0][j] = 2;
	locP[q-1][j] = 8;
    }
    for (i=1 ; i<q-1 ; i++) {
	locP[i][0] = 4;
	locP[i][r-1] = 6;
    }
    locP[0][0]     = 1;
    locP[0][r-1]   = 3;
    locP[q-1][0]   = 7;
    locP[q-1][r-1] = 9;

    all_stop_timer("Get pType");

    all_print_timer(outfile);
}

/*************************************************************/
void all_crop_border(int v, int w, int q, int r, 
		     int X[v][w]::[q][r], int Y[v][w]::[q][r],
		     int crop) 
/*************************************************************/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj,                   /* MYPROC's column */
	row,
	col;

    int	(*locX)[r],
	(*locY)[r];

    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    
    all_init_timer();

    all_start_timer();

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {
	    row = ii*q + i;
	    col = jj*r + j;
	    if ( (row < crop) || (row >= q*v - crop) ||
		 (col < crop) || (col >= r*w - crop) )
		locY[i][j] = 0;
	    else
		locY[i][j] = locX[i][j];
	}

    barrier();

    all_stop_timer("Crop Border");

    all_print_timer(outfile);
}

/*************************************************************/
void all_filter_SNF(int v, int w, int q, int r, 
		      int X[v][w]::[q][r], int k,
		      int Y[v][w]::[q][r],
		      int pType[v][w]::[q][r],
		      int steps, int epsilon)
/*************************************************************/
/* Perform SNF on interior points (do not modify 1-pixel image
   border).
*/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	d0, d1, d2,
	ii,                   /* MYPROC's row */
	jj,                   /* MYPROC's column */
	sum;

    int	(*locX)[r],
	(*locY)[r],
	(*locP)[r],
	(*locSNFin)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	tot,
	state[q][r],   /* 0=border, 1=non-fixed, 2=fixed1, 3=fixed2, 4=fixed */
	fp,
	myFP,
	tempFP,
	step,
	pfstate;
    
    int SNFin[v][w]::[q][r];
    
    double pfixed,
	limit,
	lastpf;
    
    char buf[MAXLEN];

    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    locSNFin = tolocal(SNFin[ii][jj]);
    
    all_init_timer();

    all_start_timer();

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {
	    state[i][j] = !border_pix(ii,jj,v,w,q,r,i,j);
	    locSNFin[i][j] = locX[i][j];
	}

    tot = ((q*v)-2) * ((r*w)-2);

/* constants */

    limit = 100.0;

    pfixed = 0.0;

    /* pfstate is the number of iterations that pfixed has not changed.
       After 3 iterations of no change, stop iterating. */

    pfstate = 0;

    all_stop_timer("SNF Initialization");

    for (step=1 ; ((step<=steps)&&(pfixed<limit)&&(pfstate<3)) ; step++) {

/* CHECK */
	myFP = 0;

/* Get copy of input image */

	if (step > 1)
	    for (i=0 ; i<q ; i++)
		for (j=0 ; j<r ; j++) 
		    locSNFin[i][j] = locY[i][j];


	all_start_timer();
	
/******************************************************/

/* Create Ghost Cells for SNF input image */

	all_get_ghost_eight(ii,jj,v,w,q,r,SNFin,
			    ghostN, ghostS, ghostE, ghostW,
			    &ghostNW, &ghostNE, &ghostSW, &ghostSE,
			    NIL);

/********************************************************/

	for (i=0 ; i<q ; i++)
	    for (j=0 ; j<r ; j++) {		    
		    
		d0 = locSNFin[i][j];

		if ((state[i][j]>0) && (state[i][j]<4)) {

		    for (t=0 ; t<8 ; t++)
			mask[t] = 1;
	    
		    my_proc_nbr(q, r, locSNFin, locP,
				ghostN, ghostS, ghostE, ghostW,
				ghostNW, ghostNE, ghostSW, ghostSE,
				NIL, i, j, mask, temp);
	    
		    sum = (d0 << 2) + 2;

		    d1 = temp[2];
		    d2 = temp[3];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += d1;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += d2;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += d2;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += d1;
			    }

		    d1 = temp[1];
		    d2 = temp[0];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += d1;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += d2;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += d2;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += d1;
			    }

		    d1 = temp[7];
		    d2 = temp[4];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += d1;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += d2;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += d2;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += d1;
			    }

		    d1 = temp[6];
		    d2 = temp[5];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += d1;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += d2;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += d2;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += d1;
			    }

		    locY[i][j] = sum>>2;
		    tempFP = (locY[i][j] == d0);
		    myFP += tempFP;
		    if (tempFP) {
			switch (state[i][j]) {
			case 1:
			case 2:
			case 3:
			    state[i][j]++;
			    break;
			default: fprintf(stderr,
					 "ERROR; default case reached\n");
			}
		    }
		    else {
			switch (state[i][j]) {
			case 1:
			    break;
			case 2:
			case 3:
			    state[i][j] = 1;
			    break;
			default: fprintf(stderr,
					 "ERROR; default case reached\n");
			}
		    };
		}
		else {
		    locY[i][j] = d0;
		    if (state[i][j]==4) myFP++;
		}
	    }

	sprintf(buf,"SNF iteration %3d",step);
	all_stop_timer(buf);
	
	fp = all_reduce_to_all_add(myFP);

	lastpf = pfixed;

	pfixed = 100.0 * fp / tot;

	if (lastpf==pfixed)
	    pfstate++;
	else
	    pfstate=0;

	if (RESULTS)
	    on_one 
		fprintf(outfile,"SNF3x3  epsilon= %3d   step= %2d   fixed= %8.4f\n"
			,epsilon,step,pfixed);

    }

    barrier();

    all_print_timer(outfile);
}

/*************************************************************/
void all_filter_SNF_reclass(int v, int w, int q, int r, 
			    int X[v][w]::[q][r], int k,
			    int Sel[v][w]::[q][r],
			    int Y[v][w]::[q][r],
			    int pType[v][w]::[q][r],
			    int epsilon)
/*************************************************************/
/* Perform SNF on interior points (do not modify 1-pixel image
   border).
*/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	d0, d1, d2,
	x0,
	ii,                   /* MYPROC's row */
	jj,                   /* MYPROC's column */
	sum;

    int	(*locX)[r],
	(*locY)[r],
	(*locSel)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	selN[r], selS[r], selE[q], selW[q],
	selNE, selSE, selSW, selNW,
	temp_x[8],
	temp_sel[8],
	mask[8];
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    locSel = tolocal(Sel[ii][jj]);
    
    all_init_timer();

    all_start_timer();

/******************************************************/

/* Create Ghost Cells for SNF selection image */

	all_get_ghost_eight(ii,jj,v,w,q,r,Sel,
			    selN, selS, selE, selW,
			    &selNW, &selNE, &selSW, &selSE,
			    NIL);

/********************************************************/

/* Create Ghost Cells for SNF input image */

	all_get_ghost_eight(ii,jj,v,w,q,r,X,
			    ghostN, ghostS, ghostE, ghostW,
			    &ghostNW, &ghostNE, &ghostSW, &ghostSE,
			    NIL);

/********************************************************/

	for (i=0 ; i<q ; i++)
	    for (j=0 ; j<r ; j++) {		    
		    
		if (!border_pix(ii,jj,v,w,q,r,i,j)) {

		    my_proc_nbr(q, r, locX, locP,
				ghostN, ghostS, ghostE, ghostW,
				ghostNW, ghostNE, ghostSW, ghostSE,
				NIL, i, j, mask, temp_x);
	    
		    my_proc_nbr(q, r, locSel, locP,
				selN, selS, selE, selW,
				selNW, selNE, selSW, selSE,
				NIL, i, j, mask, temp_sel);

		    x0 = locX[i][j];
		    sum = (x0 << 2) + 2;

		    d0 = locSel[i][j];
		    d1 = temp_sel[2];
		    d2 = temp_sel[3];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += temp_x[2] - x0;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += temp_x[3] - x0;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += temp_x[3] - x0;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += temp_x[2] - x0;
			    }

		    d1 = temp_sel[1];
		    d2 = temp_sel[0];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += temp_x[1] - x0;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += temp_x[0] - x0;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += temp_x[0] - x0;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += temp_x[1] - x0;
			    }

		    d1 = temp_sel[7];
		    d2 = temp_sel[4];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += temp_x[7] - x0;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += temp_x[4] - x0;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += temp_x[4] - x0;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += temp_x[7] - x0;
			    }

		    d1 = temp_sel[6];
		    d2 = temp_sel[5];
		    if ((d1-=d0) && (d2-=d0) && (d1+d2))
			if (d1+d2 > 0)
			    if (d1 < d2) {
				if (d1 <= epsilon)
				    sum += temp_x[6] - x0;
			    }
			    else {
				if (d2 <= epsilon)
				    sum += temp_x[5] - x0;
			    }
			else
			    if (d1 < d2) {
				if (d2 >= -epsilon)
				    sum += temp_x[5] - x0;
			    }
			    else {
				if (d1 >= -epsilon)
				    sum += temp_x[6] - x0;
			    }

		    locY[i][j] = sum>>2;

		}
		else
		    locY[i][j] = locX[i][j];
	    }

    all_stop_timer("SNF-reclass stage");

    barrier();

    all_print_timer(outfile);
}

/*************************************************************/
double all_find_sigma_use_sort(int v, int w, int q, int r, 
			       int X[v][w]::[q][r], int k,
			       int pType[v][w]::[q][r])
/*************************************************************/
{
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	sum,
	xsq,
	ss,
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	(*locX)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	*ssTile,
	tidx,
	bOffset;

    double sigma;

    int ssT[PROCS]::[q*r];
    int ssS[PROCS]::[q*r];

    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    ssTile = tolocal(ssT[MYPROC]);
    
    all_init_timer();

    all_start_timer();

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

    tidx = 0;
    
    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {		    
	    
	    /* Test for image border pixels */
	    if ( !border_pix(ii,jj,v,w,q,r,i,j) ) {
		
		my_proc_nbr(q, r, locX, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		sum = locX[i][j];
		xsq = sum*sum; 
		
		for (t=0 ; t<8 ; t++) {
		    sum += temp[t];
		    xsq += temp[t] * temp[t];
		}

		ss = 9 * xsq - (sum*sum);
	    }
	    else 
		ss = 0;
	    
	    ssTile[tidx++] = ss;
	}
    all_stop_timer("Find Sigma: Part 1");

    all_start_timer();
#if 0
    all_columnsort_int(tidx, ssT, ssS);
#endif    
#if 1
    all_radixsort_bitInfo(tidx, ssT, ssS, 21, 8);
#endif    
    all_stop_timer("Find Sigma: Part 2, Sort");

    all_start_timer();

    on_one {
	bOffset = 2*q*v + 2*r*w - 4;
	i = ssS[(PROCS>>1)+(bOffset/tidx)][0+(bOffset%tidx)];
	bOffset -= 1;
	i += ssS[(PROCS>>1)+(bOffset/tidx)][0+(bOffset%tidx)];
/* This is the population unbiased estimate */
/*	sigma = (double)i/(2.0*9.0*8.0); */
/* This is the standard estimate of variance  */
	sigma = (double)i/(2.0*9.0*9.0);
    }

    sigma = all_bcast_d(sigma);
    
    all_stop_timer("Find Sigma: Part 3");

    all_print_timer(outfile);

    return(sigma);
}

/*************************************************************/
double all_find_sigma_use_median(int v, int w, int q, int r, 
				 int X[v][w]::[q][r], int k,
				 int pType[v][w]::[q][r])
/*************************************************************/
{
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	sum,
	xsq,
	ss,
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	(*locX)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	*ssTile,
	tidx;

    double sigma;

    int *spread ssT;
    int *spread N;

    ssT = all_spread_malloc(PROCS,q*r*sizeof(int));
    assert_spread_malloc(ssT);
    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);
    

    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    ssTile = (int *)(ssT + MYPROC);
    
    all_init_timer();

    all_start_timer();

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

    tidx = 0;
    
    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {		    
	    
	    /* Test for image border pixels */
	    if ( !border_pix(ii,jj,v,w,q,r,i,j) ) {
		
		my_proc_nbr(q, r, locX, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		sum = locX[i][j];
		xsq = sum*sum; 
		
		for (t=0 ; t<8 ; t++) {
		    sum += temp[t];
		    xsq += temp[t] * temp[t];
		}

		ss = 9 * xsq - (sum*sum);

		ssTile[tidx++] = ss;
	    }
	}
    all_stop_timer("Find Sigma: Part 1");

    all_start_timer();

    *(N+MYPROC) = tidx;
    t = v*w*q*r - 2*q*v - 2*r*w + 4;
    i = all_select_median_unbalanced_c_i(q*r, ssT, N, t );

/* This is the population unbiased estimate */
/*	sigma = (double)i/(9.0*8.0); */
/* This is the standard estimate of variance  */

    sigma = (double)i/(9.0*9.0);
    
    all_stop_timer("Find Sigma: Part 2, Select");

    all_print_timer(outfile);

    all_spread_free(ssT);
    all_spread_free(N);
    return(sigma);
}

/*************************************************************/
void all_filter_1NN(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int k,
		    int Y[v][w]::[q][r],
		    int pType[v][w]::[q][r])
/*************************************************************/
/* Perform 1-NN filter. */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	pix,
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int T[v][w]::[q][r];

    int	(*locX)[r],
	(*locY)[r],
	(*locT)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	diff,
	num,
	step, fp, lastfp, myFP;

    char buf[MAXLEN];
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locT = tolocal(T[ii][jj]);
    locP = tolocal(pType[ii][jj]);

    all_init_timer();

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) 
	    locT[i][j] = locX[i][j];

    lastfp = 0;
    fp = 1;

    for (step=1 ; (lastfp != fp) ; step++) {

	myFP = 0;

/* Get copy of input image */

	if (step > 1)
	    for (i=0 ; i<q ; i++)
		for (j=0 ; j<r ; j++) 
		    locT[i][j] = locY[i][j];

	all_start_timer();
	
/******************************************************/

/* Create Ghost Cells for input image */

	all_get_ghost_eight(ii,jj,v,w,q,r,T,
			    ghostN, ghostS, ghostE, ghostW,
			    &ghostNW, &ghostNE, &ghostSW, &ghostSE,
			    NIL);

/********************************************************/

	for (i=0 ; i<q ; i++)
	    for (j=0 ; j<r ; j++) {		    
		    
		for (t=0 ; t<8 ; t++)
		    mask[t] = 1;
	    
		my_proc_nbr(q, r, locT, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		pix = k;

		for (t=0 ; t<8 ; t++)
		    if (mask[t]) {
			diff = abs(locT[i][j] - temp[t]);
			if (diff < pix) {
			    pix = diff;
			    num = t;
			}
		    }

		locY[i][j] = (temp[num] + locT[i][j]) / 2 ;
		if (locY[i][j] == locT[i][j])
		    myFP++;
	    }

	sprintf(buf,"1NN iteration %3d",step);
	all_stop_timer(buf);

	lastfp = fp;
	fp = all_reduce_to_all_add(myFP);

	if (RESULTS)
	    on_one 
		fprintf(outfile,"1NN   step= %2d   fixed= %8.4f\n", step,
			(double)fp/(double)(v*w*q*r));
	barrier();
    }

    all_print_timer(outfile);
}

/*************************************************************/
void all_filter_1NN_oneshot(int v, int w, int q, int r, 
			    int X[v][w]::[q][r], int k,
			    int Y[v][w]::[q][r],
			    int pType[v][w]::[q][r])
/*************************************************************/
/* Perform 1-NN filter. */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	t,
	pix,
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	(*locX)[r],
	(*locY)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	diff,
	num;
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locP = tolocal(pType[ii][jj]);

    all_init_timer();

    all_start_timer();

/******************************************************/

/* Create Ghost Cells for input image */

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

/********************************************************/

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {		    
		    
	    for (t=0 ; t<8 ; t++)
		mask[t] = 1;
	    
	    my_proc_nbr(q, r, locX, locP,
			ghostN, ghostS, ghostE, ghostW,
			ghostNW, ghostNE, ghostSW, ghostSE,
			NIL, i, j, mask, temp);

	    pix = k;

	    for (t=0 ; t<8 ; t++)
		if (mask[t]) {
		    diff = abs(locX[i][j] - temp[t]);
		    if (diff < pix) {
			pix = diff;
			num = t;
		    }
		}

	    locY[i][j] = temp[num];
	}

    barrier();

    all_stop_timer("1-NN Filter");

    all_print_timer(outfile);
}

/*************************************************************/
void all_scale_and_edge_detect(int v, int w, int q, int r,
			       int mag, 
			       int X[v][w]::[q][r],
			       int Y[v][w]::[q*mag][r*mag],
			       int pType[v][w]::[q][r],
			       int edgeColor,
			       int Ycol[v][w]::[q][r],
			       int Ymag[v][w]::[q*mag][r*mag])
/*************************************************************/
/* Perform a scale by mag, and detect edges */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj,                   /* MYPROC's column */
	x, y;

    int	(*locX)[r],
	(*locY)[r*mag],
	(*locYcol)[r],
	(*locYmag)[r*mag],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8],
	d0;
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locYcol = tolocal(Ycol[ii][jj]);
    locYmag = tolocal(Ymag[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    
    all_init_timer();

    all_start_timer();

/******************************************************/

/* Create Ghost Cells for input image */

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

/********************************************************/

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {

	    /* Test for image border pixels */
	    if ( !border_pix(ii,jj,v,w,q,r,i,j) ) {
		
		my_proc_nbr(q, r, locX, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		d0 = locX[i][j];
		
/* N */		for (y=0 ; y<mag ; y++)
		    if (d0 != temp[0])
			locY[i*mag + 0][j*mag + y] = edgeColor;
		    else 
			locY[i*mag + 0][j*mag + y] = d0;

/* S */		for (y=0 ; y<mag ; y++) 
		    if (d0 != temp[1])
			locY[i*mag + mag-1][j*mag + y] = edgeColor;
		    else 
			locY[i*mag + mag-1][j*mag + y] = d0;
/* W */		for (x=1 ; x<mag-1 ; x++) 
		    if (d0 != temp[2])
			locY[i*mag + x][j*mag + 0] = edgeColor;
		    else 
			locY[i*mag + x][j*mag + 0] = d0;
/* E */		for (x=1 ; x<mag-1 ; x++) 
		    if (d0 != temp[3])
			locY[i*mag + x][j*mag + mag-1] = edgeColor;
		    else 
			locY[i*mag + x][j*mag + mag-1] = d0;
/* Interior */
		for (x=1 ; x<mag-1 ; x++)
		    for (y=1 ; y<mag-1 ; y++)
			locY[i*mag + x][j*mag + y] = d0;
	    }
	    else
		for (x=0 ; x<mag ; x++)
		    for (y=0 ; y<mag ; y++)
			locY[i*mag + x][j*mag + y] = edgeColor;

	    for (x=0 ; x<mag ; x++)
		for (y=0 ; y<mag ; y++)
		    if (locY[i*mag + x][j*mag + y] == edgeColor)
			locYmag[i*mag + x][j*mag + y] = edgeColor;
		    else
			locYmag[i*mag + x][j*mag + y] = locYcol[i][j];
	}

    barrier();

    all_stop_timer("Scale and Edge Detection");

    all_print_timer(outfile);
}

/*************************************************************/
void all_filter_Gauss3X3(int v, int w, int q, int r, 
			 int X[v][w]::[q][r], int k,
			 int Y[v][w]::[q][r],
			 int pType[v][w]::[q][r])
/*************************************************************/
/* Perform 3x3 Gaussian convolution filter. */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	(*locX)[r],
	(*locY)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8];
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    
    all_init_timer();

    all_start_timer();

/******************************************************/

/* Create Ghost Cells for input image */

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

/********************************************************/

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {

	    /* Test for image border pixels */
	    if ( !border_pix(ii,jj,v,w,q,r,i,j) ) {
		
		my_proc_nbr(q, r, locX, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		locY[i][j] = ( (locX[i][j] << 2) +
			       ((temp[0] + temp[1] + temp[2] + temp[3])<<1) +
			       temp[4] + temp[5] + temp[6] + temp[7]) >> 4;
		    
	    }
	    else
		locY[i][j] = locX[i][j];
	}

    barrier();

    all_stop_timer("Gaussian 3x3 Convolution");

    all_print_timer(outfile);
}

/*************************************************************/
void all_filter_GaussNoise(int v, int w, int q, int r, 
			   int X[v][w]::[q][r], int k,
			   int Y[v][w]::[q][r],
			   double mean, double sigma)
/*************************************************************/
/* Add Gaussian noise N(mean, variance) to image */
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj,                   /* MYPROC's column */
	t;
    
    int	(*locX)[r],
	(*locY)[r];
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locY = tolocal(Y[ii][jj]);
    
    all_init_timer();

    all_start_timer();

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++)  {
	    t = locX[i][j] +
		(int)floor(rand_normal(mean, sigma));
	    if (t<0)   t = 0;
	    if (t>k-1) t = k-1;
	    locY[i][j] = t;
	}

    barrier();

    all_stop_timer("Add Gaussian Noise");

    all_print_timer(outfile);
}

/*************************************************************/
int all_decide_artificial(int v, int w, int q, int r, 
			  int X[v][w]::[q][r],
			  int pType[v][w]::[q][r])
/*************************************************************/
{
    
    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	pix,
	range,
	ct,
	t,
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	(*locX)[r],
	(*locP)[r],
	ghostN[r], ghostS[r], ghostE[q], ghostW[q],
	ghostNE, ghostSE, ghostSW, ghostNW,
	temp[8],
	mask[8];
    
    barrier();
    
    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locP = tolocal(pType[ii][jj]);
    
    all_init_timer();

    all_start_timer();

/******************************************************/

/* Create Ghost Cells for input image */

    all_get_ghost_eight(ii,jj,v,w,q,r,X,
			ghostN, ghostS, ghostE, ghostW,
			&ghostNW, &ghostNE, &ghostSW, &ghostSE,
			NIL);

/********************************************************/

    ct = 0;
    
    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {

	    /* Test for image border pixels */
	    if ( !border_pix(ii,jj,v,w,q,r,i,j) ) {
		
		my_proc_nbr(q, r, locX, locP,
			    ghostN, ghostS, ghostE, ghostW,
			    ghostNW, ghostNE, ghostSW, ghostSE,
			    NIL, i, j, mask, temp);

		pix = locX[i][j];
		range = 0;

		for (t=0 ; t<8 ; t++)
		    if (temp[t] != pix)
			range = 1;

		if (!range)
		    ct++;
		    
	    }
	}

    range = all_reduce_to_all_add(ct);
    ct = (q*v - 2)*(r*w - 2);

    all_stop_timer("Decide Artificial");

    all_print_timer(outfile);

    if (((double)range / (double)ct) >= 0.50)
	return(1);

    return(0);
}


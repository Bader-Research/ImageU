/*									tab:8
 *
 * sorting.h - fast sequential sorting routines
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
 * Filename:		sorting.h
 * History:
 */


#ifndef _SORTING_H
#define _SORTING_H

#include "ImageU.h"
#include "intpair.h"
#include "mach_specifics.h"

int intcompare(int *, int *);
int elemcompare(ELEMPTR, ELEMPTR);
int chpaircompare(CHPAIRPTR, CHPAIRPTR);
int edgecompare(edge_t *, edge_t *);
int intpaircompare(intpair_t *, intpair_t *);
int binsearch(chPair_t *, int , int , int );
void radixsort(int *, int);
void radixsort_19(int *, int);

_INLINE unsigned bits(unsigned x, int k, int j) {
/* Return the j bits which appear k bits from the right in x */
    return (x>>k) & ~(~0<<j);
}

_INLINE void elemcopy(ELEMPTR a, ELEMPTR b) {
/* Deep copy of b into a */

    if ((a==NULL)||(b==NULL)) return;
    a->color = b->color;
    a->label = b->label;
    a->pos   = b->pos;
}

_INLINE void chPairCopy(CHPAIRPTR a, CHPAIRPTR b) {
/* Deep copy of b into a */

    a->alpha = b->alpha;
    a->beta  = b->beta;
}

_INLINE void radixsort_elem(ELEMPTR a, int n) {
/* Radix sort a list of n elem_pointers, a[], with labels between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	pass;

    int	nextpass,
	*count,
	*bitArr;

    ELEMPTR b;
    
    if ((b = (ELEMPTR)malloc(n*sizeof(elem_t)))==NULL)
	fprintf(stderr,"ERROR: radix temp array could not be malloc'ed\n");

    if ((bitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: bitArr could not be malloc'ed\n");

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */

    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i].label,pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) elemcopy(&b[--count[bitArr[i]]], &a[i]); 

	nextpass = pass+1;

	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i].label,nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) elemcopy(&a[--count[bitArr[i]]], &b[i]);
    }
    free(bitArr);
    free(count);
    free(b);
}

_INLINE void radixsort_chPair(CHPAIRPTR a, int n) {
/* Radix sort a list of n chPair_pointers, a[], with labels between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	pass;

    int	nextpass,
	*count,
	*bitArr;

    CHPAIRPTR b;
    
    if ((b = (CHPAIRPTR)malloc(n*sizeof(chPair_t)))==NULL)
	fprintf(stderr,"ERROR: radix temp array could not be malloc'ed\n");

    if ((bitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: bitArr could not be malloc'ed\n");

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */

    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i].alpha,pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) chPairCopy(&b[--count[bitArr[i]]], &a[i]); 

	nextpass = pass+1;

	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i].alpha,nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) chPairCopy(&a[--count[bitArr[i]]], &b[i]);
    }
    free(bitArr);
    free(count);
    free(b);
}

_INLINE int chLookup(chPair_t * chList, int chSize, int oldLabel)
/* Search chList.Alpha for oldLabel, and return beta if it exists,
   or oldLabel if it does not */
/* since chListAlpha is sorted, we can use a binary search */
{
    register int found_idx;

    found_idx = binsearch(chList, 0, chSize-1, oldLabel);
    if (found_idx < 0)
	return(oldLabel);
    else
	return(chList[found_idx].beta);

}

_INLINE int fill_chList(chPair_t * chList, CHPAIRPTR changeList, int len) 
/* Copy the array changeList into two integer arrays. Copy ONLY
   the unique pairs (changeList is ALREADY sorted by alpha's.
   Return the number of unique pairs. */
{
    register int i;
    register int count = 0;

    if (len > 0) {
	chList[0].alpha = changeList[0].alpha;
	chList[0].beta  = changeList[0].beta;
	count++;
    }

    for (i=1 ; i < len ; i++) {
	if (changeList[i].alpha > chList[count-1].alpha) {
	    chList[count].alpha = changeList[i].alpha;
	    chList[count].beta  = changeList[i].beta;
	    count++;
	}
    }

    return(count);
}

_INLINE void radixsort_edge(edge_t * aL, int n) {
/* Radix sort a list of n edges, a[], with labels between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	temp,
	pass;

    int	nextpass,
	*count,
	*bitArr;

    edge_t *bL;
    
    if ((bL = (edge_t *)malloc(n*sizeof(edge_t)))==NULL)
	fprintf(stderr,"ERROR: radix temp array could not be malloc'ed\n");

    if ((bitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: bitArr could not be malloc'ed\n");

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */

    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   aL[] and bL[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(aL[i].a,pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) {
	    temp = --count[bitArr[i]];
	    bL[temp].a = aL[i].a; 
	    bL[temp].b = aL[i].b;
	}

	nextpass = pass+1;

	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(bL[i].a,nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) {
	    temp = --count[bitArr[i]];
	    aL[temp].a = bL[i].a;
	    aL[temp].b = bL[i].b;
	}
    }
    free(bitArr);
    free(count);
    free(bL);
}

_INLINE void radixsort_intpair(intpair_t * aL, int n) {
/* Radix sort a list of n intpairs, a[], with labels between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	temp,
	pass;

    int	nextpass,
	*count,
	*bitArr;

    intpair_t *bL;
    
    if ((bL = (intpair_t *)malloc(n*sizeof(intpair_t)))==NULL)
	fprintf(stderr,"ERROR: radix temp array could not be malloc'ed\n");

    if ((bitArr = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: bitArr could not be malloc'ed\n");

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */

    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   aL[] and bL[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(aL[i].a,pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) {
	    temp = --count[bitArr[i]];
	    bL[temp].a = aL[i].a; 
	    bL[temp].b = aL[i].b;
	}

	nextpass = pass+1;

	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(bL[i].a,nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) {
	    temp = --count[bitArr[i]];
	    aL[temp].a = bL[i].a;
	    aL[temp].b = bL[i].b;
	}
    }
    free(bitArr);
    free(count);
    free(bL);
}

#define insertsort(a,b) insertsort_i(a,b)

_INLINE void insertsort_i(int *A, int n) {

#define DATA_TYPE int
    
    register DATA_TYPE item;
    register int i,j;
    
    for (i=1 ; i<n ; i++) {
	item = A[i];
	j = i-1;
	while ((j>=0)&&(item < A[j])) {
	    A[j+1] = A[j];
	    j--;
	}
	A[j+1] = item;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_d(double *A, int n) {

#define DATA_TYPE double
    
    register DATA_TYPE item;
    register int i,j;
    
    for (i=1 ; i<n ; i++) {
	item = A[i];
	j = i-1;
	while ((j>=0)&&(item < A[j])) {
	    A[j+1] = A[j];
	    j--;
	}
	A[j+1] = item;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_elem(elem_t *A, int n) {

#define DATA_TYPE elem_t
    
    register DATA_TYPE item;
    register int i,j,k;
    
    for (i=1 ; i<n ; i++) {
	item.color = A[i].color;
	item.label = A[i].label;
	item.pos   = A[i].pos;
	j = i-1;
	while ((j>=0)&&(item.label < A[j].label)) {
	    A[k = j+1].color = A[j].color;
	    A[k].label       = A[j].label;
	    A[k].pos         = A[j].pos;
	    j--;
	}
	A[k = j+1].color = item.color;
	A[k = j+1].label = item.label;
	A[k = j+1].pos   = item.pos;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_chPair(chPair_t *A, int n) {

#define DATA_TYPE chPair_t
    
    register DATA_TYPE item;
    register int i,j, k;
    
    for (i=1 ; i<n ; i++) {
	item.alpha = A[i].alpha;
	item.beta  = A[i].beta;
	j = i-1;
	while ((j>=0)&&(item.alpha < A[j].alpha)) {
	    A[k = j+1].alpha = A[j].alpha;
	    A[k].beta        = A[j].beta;
	    j--;
	}
	A[k = j+1].alpha = item.alpha;
	A[k].beta        = item.beta;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_edge(edge_t *A, int n) {

#define DATA_TYPE edge_t
    
    register DATA_TYPE item;
    register int i,j, k;
    
    for (i=1 ; i<n ; i++) {
	item.a = A[i].a;
	item.b = A[i].b;
	j = i-1;
	while ((j>=0)&&
	       ((item.a < A[j].a) || ((item.a==A[j].a)&&(item.b < A[j].b)))) {
	    A[k = j+1].a = A[j].a;
	    A[k].b       = A[j].b;
	    j--;
	}
	A[k = j+1].a = item.a;
	A[k].b       = item.b;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_intpair(intpair_t *A, int n) {

#define DATA_TYPE intpair_t
    
    register DATA_TYPE item;
    register int i,j, k;
    
    for (i=1 ; i<n ; i++) {
	item.a = A[i].a;
	item.b = A[i].b;
	j = i-1;
	while ((j>=0)&&
	       ((item.a < A[j].a) || ((item.a==A[j].a)&&(item.b < A[j].b)))) {
	    A[k = j+1].a = A[j].a;
	    A[k].b       = A[j].b;
	    j--;
	}
	A[k = j+1].a = item.a;
	A[k].b       = item.b;
    }

#undef DATA_TYPE    
}

_INLINE void fastsort(int* arr, int nel) {

    if (nel>=RADIXSORT_INT_BREAKPT)
	radixsort(arr,nel); 
    else
	insertsort(arr,nel);
}

_INLINE void fastsort_19(int* arr, int nel) {

    if (nel>=RADIXSORT_INT_BREAKPT)
	radixsort_19(arr,nel); 
    else
	insertsort_i(arr,nel);
}

_INLINE void fastsort_elem(ELEMPTR arr, int nel) {

    if (nel>=RADIXSORT_ELEM_BREAKPT)
	radixsort_elem(arr,nel);
    else
	insertsort_elem(arr,nel);
}
	      

_INLINE void fastsort_chPair(CHPAIRPTR arr, int nel) {

    if (nel>=RADIXSORT_CHPAIR_BREAKPT)
	radixsort_chPair(arr,nel);
    else
        insertsort_chPair(arr,nel);
}
	      
_INLINE void fastsort_edge(edge_t *arr, int nel) {

    if (nel>=RADIXSORT_EDGE_BREAKPT)
	radixsort_edge(arr,nel);
    else
	insertsort_edge(arr,nel);
}

_INLINE void fastsort_intpair(intpair_t *arr, int nel) {

    if (nel>=RADIXSORT_EDGE_BREAKPT)
	radixsort_intpair(arr,nel);
    else
	insertsort_intpair(arr,nel);
}


#endif



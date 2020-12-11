/*									tab:8
 *
 * sorting.sc - fast sequential sorting routines
 *
 * 
 * "Copyright (c) 1994,1995,1996 The Regents of the University of Maryland.
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
 * Filename:		sorting.sc
 * History:
 */

#include "sorting.h"

int intcompare(int *i, int *j)
{
    return(*i - *j);
}

int elemcompare(ELEMPTR i, ELEMPTR j)
{
    return(i->label - j->label);
}

int chpaircompare(CHPAIRPTR i, CHPAIRPTR j)
{
    return(i->alpha - j->alpha);
}

int edgecompare(edge_t *i, edge_t *j)
{
    int temp;
    temp = i->a - j->a;
    return ((temp!=0) ? temp : i->b - j->b);
}

int intpaircompare(intpair_t *i, intpair_t *j)
{
    int temp;
    temp = i->a - j->a;
    return ((temp!=0) ? temp : i->b - j->b);
}

int binsearch(chPair_t * searchList, int min_idx, int max_idx, int target) 
/* Search searchList recursively using binary search to find target.
   Return -1 if it does not exist, or target's index if it does. */
{
    register int mid_idx;

    mid_idx = min_idx + (max_idx - min_idx) / 2;
    if (max_idx < min_idx)
	return(-1);
    else if (target < searchList[mid_idx].alpha)
	return binsearch(searchList, min_idx,     mid_idx - 1, target);
    else if (target > searchList[mid_idx].alpha)
	return binsearch(searchList, mid_idx + 1, max_idx    , target);
    else
	return mid_idx;
}



void radixsort(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
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
	bitArr[n],
	b[n];

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

	nextpass = pass+1;
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];
    }
    free(count);
}

void radixsort_19(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
/* Assume m = 19 */

    register int
	i,
	j;

    int	*count,
	bitArr[n],
	b[n];

    count = (int*)malloc(1024*sizeof(int));
    assert_malloc(count);

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (j=0 ; j<1024 ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],0,10)]++;
    for (j=1 ; j<1024 ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

    for (j=0 ; j<512 ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],10,9)]++;
    for (j=1 ; j<512 ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];

    free(count);
}


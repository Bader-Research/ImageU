/*                                                                    tab:8
 *
 * select_seq.sc - Sequential Selection Algorithm
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
 * Filename:            select_seq.sc
 * History:
 */

#include "select_seq.h"
#include "sorting.h"
#include <math.h>

#define DEFAULT_R       7
#define SELECT_CUTOFF  22

int select_mom_i(int* list, int n, int k) {

#define DATA_TYPE int
    
    register DATA_TYPE
	*s1, *s2;

    register
	blocks,
	left;

    DATA_TYPE
	result,
	*sublist,
	*listptr,
	*listptr_med;
    
    if (n < SELECT_CUTOFF) {
	insertsort_i(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    sublist = (DATA_TYPE *)malloc(blocks*sizeof(DATA_TYPE));
    assert_malloc(sublist);
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_i(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_i(sublist,blocks,blocks/2);
    
    free(sublist);
    
    left = partition_i(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_i(list,left,k);
	else
	    if (k > left) 
		result = select_mom_i(list+left,n-left,k-left);
    }
    else {
	insertsort_i(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}

double select_mom_d(double* list, int n, int k) {

#define DATA_TYPE double
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*sublist,
	*listptr,
	*listptr_med;
    
    if (n < SELECT_CUTOFF) {
	insertsort_d(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;

    sublist = (DATA_TYPE *)malloc(blocks*sizeof(DATA_TYPE));
    assert_malloc(sublist);

    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_d(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_d(sublist,blocks,blocks/2);
    
    free(sublist);
    
    left = partition_d(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_d(list,left,k);
	else
	    if (k > left) 
		result = select_mom_d(list+left,n-left,k-left);
    }
    else {
	insertsort_d(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}


int select_mom_alloc_i(int* list, int n, int k, int* sublist) {

#define DATA_TYPE int
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*listptr,
	*listptr_med;
    
    if (n < SELECT_CUTOFF) {
	insertsort_i(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_i(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_alloc_i(sublist,blocks,blocks/2,s2);
    
    left = partition_i(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_alloc_i(list,left,k,sublist);
	else
	    if (k > left) 
		result = select_mom_alloc_i(list+left,n-left,k-left,sublist);
    }
    else {
	insertsort_i(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}

double select_mom_alloc_d(double* list, int n, int k, double* sublist) {

#define DATA_TYPE double
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*listptr,
	*listptr_med;
    
    if (n < SELECT_CUTOFF) {
	insertsort_d(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_d(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_alloc_d(sublist,blocks,blocks/2,s2);
    
    left = partition_d(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_alloc_d(list,left,k,sublist);
	else
	    if (k > left) 
		result = select_mom_alloc_d(list+left,n-left,k-left,sublist);
    }
    else {
	insertsort_d(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}


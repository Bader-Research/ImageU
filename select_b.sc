/*                                                                    tab:8
 *
 * select_b.sc - Parallel Selection Algorithm
 *
 * 
 * "Copyright (c) 1995 The Regents of the University of Maryland.
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
 * Creation Date:       April 17, 1995
 * Filename:            select_b.sc
 * History:
 */

#include "select_b.h"

#define THRESH_SEL  (PROCS*PROCS)
#define PROFILE_DETAILED 0

static int binsearch_gt(int * list, int min_idx, int max_idx, int target) 
/* return the highest index of the array list with element <= target */
/* return min_idx-1 if all elements of list > target */
/* "THE PRICE IS RIGHT" */    
{
    register int mid_idx;

#if PROFILE_DETAILED
    fprintf(outfile,"PE %2d: in binsearch_gt (%5d, %5d): %5d\n",
	    MYPROC, min_idx, max_idx, target);
#endif
    
    if (max_idx < min_idx)
	return (min_idx-1);
    if (min_idx == max_idx)
	if (list[min_idx] <= target)
	    return (min_idx);
	else
	    return (min_idx-1);
    mid_idx = min_idx + (max_idx - min_idx) / 2;
    if (list[mid_idx] <= target)                 /* This person bid under */
	return binsearch_gt(list, mid_idx+1, max_idx,   target);
    else                                         /* This person bid over  */
	return binsearch_gt(list, min_idx,   mid_idx-1, target);
}


void all_LoadBalance_b_i(int M, int *spread A,
			 int *spread N, int total_n) {
/* Load Balance in place */

#define DATA_TYPE int

    register int
/*	get_proc, */
	j, t;
    
    int locN[PROCS],
	q[PROCS],
	d[PROCS],
	src[PROCS],
	snk[PROCS],
	src_rnk[PROCS],
	snk_rnk[PROCS],
	src_cnt, snk_cnt,
	my_extra_elems,
	my_extra_elem_idx,
	l, r,
	my_rank_first, my_rank_last,
	l_rank_first, l_rank_last,
	l_size, l_offset,
	done;

    DATA_TYPE *locA;

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_LoadBalance_b: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    barrier();
#endif

    locA = (DATA_TYPE *)(A+MYPROC);

#if 1
    locN[0] = *(N+MYPROC);
    all_concat(locN);
#else
    for (j=0 ; j<PROCS ; j++) {
	get_proc = MODPROCS(MYPROC + j);
	locN[get_proc] = *(N+get_proc);
	barrier();
    }
#endif

    t  = DIVPROCS(total_n);

    if (t*PROCS < total_n)
	t++;

    for (j=0 ; j<PROCS-1 ; j++) 
	q[j] = t;

    q[PROCS-1] = total_n - (t * (PROCS-1));

    for (j=0 ; j<PROCS ; j++) {
	d[j]   = locN[j] - q[j];
	src[j] = (d[j]>0);
	snk[j] = (d[j]<0);
    }

    src_cnt=0;	
    snk_cnt=0;
    for (j=0 ; j<PROCS ; j++) {
	if (src[j]) {
	    src_cnt += d[j];
	    src_rnk[j] = src_cnt;
	}
	else
	    src_rnk[j] = 0;

	if (snk[j]) {
	    snk_cnt -= d[j];
	    snk_rnk[j] = snk_cnt;
	}
	else
	    snk_rnk[j] = 0;
    }

#if 1
    if (src[MYPROC]) {
	my_rank_first = src_rnk[MYPROC] - d[MYPROC] + 1;
	my_rank_last  = src_rnk[MYPROC];
	my_extra_elems = d[MYPROC];
	my_extra_elem_idx = q[MYPROC];

	/* Calculate left processor index */
	l = -1;
	done = 0;
	while (!done) {
	    l++;
	    if (snk[l] && (my_rank_first <= snk_rnk[l]))
		    done = 1;
	}

	l_rank_first = snk_rnk[l] + d[l] + 1;
	l_rank_last  = snk_rnk[l];

	l_offset = locN[l] + (my_rank_first - l_rank_first);
	l_size   = min(l_rank_last, my_rank_last) - my_rank_first + 1;

	bulk_store(((DATA_TYPE *global)(A+l))+l_offset,
		   locA+my_extra_elem_idx,
		   l_size*sizeof(DATA_TYPE));

	if (l_size < my_extra_elems) {

	    /* Calculate right processor index */
	    r = l;
	    done = 0;
	    while (!done) {
		r++;
		if (snk[r] && (my_rank_last <= snk_rnk[r]))
		    done = 1;
	    }
	
	    /* if (l<r) */
	    my_extra_elem_idx += l_size;
	    my_extra_elems    -= l_size;

	    while (l+1 < r) {
		l++;
		if (snk[l]) {
		    l_offset = locN[l];
		    l_size   = -d[l];

		    bulk_store(((DATA_TYPE *global)(A+l))+l_offset,
			       locA + my_extra_elem_idx,
			       l_size*sizeof(DATA_TYPE));
		    my_extra_elem_idx += l_size;
		    my_extra_elems    -= l_size;
		}
	    }

	    l_offset = locN[r];
	    l_size   = my_extra_elems;

	    bulk_store(((DATA_TYPE *global)(A+r))+l_offset,
		       locA + my_extra_elem_idx,
		       l_size*sizeof(DATA_TYPE));
	}

    }

    all_store_sync(); 
    barrier();

#endif
#if 0
    if (src[MYPROC]) {
	my_rank_first = src_rnk[MYPROC] - d[MYPROC] + 1;
	my_rank_last  = src_rnk[MYPROC];
	my_extra_elems = d[MYPROC];
	my_extra_elem_idx = q[MYPROC];

	/* Calculate left processor index */
	l = -1;
	done = 0;
	while (!done) {
	    l++;
	    if (snk[l] && (my_rank_first <= snk_rnk[l]))
		    done = 1;
	}

	l_rank_first = snk_rnk[l] + d[l] + 1;
	l_rank_last  = snk_rnk[l];

	l_offset = locN[l] + (my_rank_first - l_rank_first);
	l_size   = min(l_rank_last, my_rank_last) - my_rank_first + 1;

	bulk_write(((DATA_TYPE *global)(A+l))+l_offset,
		   locA +my_extra_elem_idx,
		   l_size*sizeof(DATA_TYPE));

	if (l_size < my_extra_elems) {

	    /* Calculate right processor index */
	    r = l;
	    done = 0;
	    while (!done) {
		r++;
		if (snk[r] && (my_rank_last <= snk_rnk[r]))
		    done = 1;
	    }
	
	    /* if (l<r) */
	    my_extra_elem_idx += l_size;
	    my_extra_elems    -= l_size;

	    while (l+1 < r) {
		l++;
		if (snk[l]) {
		    l_offset = locN[l];
		    l_size   = -d[l];

		    bulk_write(((DATA_TYPE *global)(A+l))+l_offset,
			       locA + my_extra_elem_idx,
			       l_size*sizeof(DATA_TYPE));
		    my_extra_elem_idx += l_size;
		    my_extra_elems    -= l_size;
		}
	    }

	    l_offset = locN[r];
	    l_size   = my_extra_elems;

	    bulk_write(((DATA_TYPE *global)(A+r))+l_offset,
		       locA + my_extra_elem_idx,
		       l_size*sizeof(DATA_TYPE));
	}
    }
    
    barrier();

#endif
#if 0
    if (src[MYPROC]) {
	my_rank_first = src_rnk[MYPROC] - d[MYPROC] + 1;
	my_rank_last  = src_rnk[MYPROC];
	my_extra_elems = d[MYPROC];
	my_extra_elem_idx = q[MYPROC];

	/* Calculate left processor index */
	l = -1;
	done = 0;
	while (!done) {
	    l++;
	    if (snk[l] && (my_rank_first <= snk_rnk[l]))
		    done = 1;
	}

	l_rank_first = snk_rnk[l] + d[l] + 1;
	l_rank_last  = snk_rnk[l];

	l_offset = locN[l] + (my_rank_first - l_rank_first);
	l_size   = min(l_rank_last, my_rank_last) - my_rank_first + 1;

	bulk_put(((DATA_TYPE *global)(A+l))+l_offset, locA + my_extra_elem_idx,
		   l_size*sizeof(DATA_TYPE));

	if (l_size < my_extra_elems) {

	    /* Calculate right processor index */
	    r = l;
	    done = 0;
	    while (!done) {
		r++;
		if (snk[r] && (my_rank_last <= snk_rnk[r]))
		    done = 1;
	    }
	
	    /* if (l<r) */
	    my_extra_elem_idx += l_size;
	    my_extra_elems    -= l_size;

	    while (l+1 < r) {
		l++;
		if (snk[l]) {
		    l_offset = locN[l];
		    l_size   = -d[l];

		    bulk_put(((DATA_TYPE *global)(A+l))+l_offset,
			     locA + my_extra_elem_idx,
			     l_size*sizeof(DATA_TYPE));
		    my_extra_elem_idx += l_size;
		    my_extra_elems    -= l_size;
		}
	    }

	    l_offset = locN[r];
	    l_size   = my_extra_elems;

	    bulk_put(((DATA_TYPE *global)(A+r))+l_offset, locA + my_extra_elem_idx,
		       l_size*sizeof(DATA_TYPE));
	}

    }

    sync();
    barrier();

#endif

    *(N+MYPROC) = q[MYPROC];

    barrier();

#undef DATA_TYPE
}

void all_LoadBalance_b_d(int M, double *spread A,
			 int *spread N, int total_n) {
/* Load Balance in place */

#define DATA_TYPE double

    register int
/*	get_proc, */
	j, t;
    
    int locN[PROCS],
	q[PROCS],
	d[PROCS],
	src[PROCS],
	snk[PROCS],
	src_rnk[PROCS],
	snk_rnk[PROCS],
	src_cnt, snk_cnt,
	my_extra_elems,
	my_extra_elem_idx,
	l, r,
	my_rank_first, my_rank_last,
	l_rank_first, l_rank_last,
	l_size, l_offset,
	done;

    DATA_TYPE *locA;

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_LoadBalance_b: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    barrier();
#endif

    locA = (DATA_TYPE *)(A+MYPROC);

#if 1
    locN[0] = *(N+MYPROC);
    all_concat(locN);
#else
    for (j=0 ; j<PROCS ; j++) {
	get_proc = MODPROCS(MYPROC + j);
	locN[get_proc] = *(N+get_proc);
	barrier();
    }
#endif

    t  = DIVPROCS(total_n);

    if (t*PROCS < total_n)
	t++;

    for (j=0 ; j<PROCS-1 ; j++) 
	q[j] = t;

    q[PROCS-1] = total_n - (t * (PROCS-1));

    for (j=0 ; j<PROCS ; j++) {
	d[j]   = locN[j] - q[j];
	src[j] = (d[j]>0);
	snk[j] = (d[j]<0);
    }

    src_cnt=0;	
    snk_cnt=0;
    for (j=0 ; j<PROCS ; j++) {
	if (src[j]) {
	    src_cnt += d[j];
	    src_rnk[j] = src_cnt;
	}
	else
	    src_rnk[j] = 0;

	if (snk[j]) {
	    snk_cnt -= d[j];
	    snk_rnk[j] = snk_cnt;
	}
	else
	    snk_rnk[j] = 0;
    }

    if (src[MYPROC]) {
	my_rank_first = src_rnk[MYPROC] - d[MYPROC] + 1;
	my_rank_last  = src_rnk[MYPROC];
	my_extra_elems = d[MYPROC];
	my_extra_elem_idx = q[MYPROC];

	/* Calculate left processor index */
	l = -1;
	done = 0;
	while (!done) {
	    l++;
	    if (snk[l] && (my_rank_first <= snk_rnk[l]))
		    done = 1;
	}

	l_rank_first = snk_rnk[l] + d[l] + 1;
	l_rank_last  = snk_rnk[l];

	l_offset = locN[l] + (my_rank_first - l_rank_first);
	l_size   = min(l_rank_last, my_rank_last) - my_rank_first + 1;

	bulk_store(((DATA_TYPE *global)(A+l))+l_offset,
		   locA+my_extra_elem_idx,
		   l_size*sizeof(DATA_TYPE));

	if (l_size < my_extra_elems) {

	    /* Calculate right processor index */
	    r = l;
	    done = 0;
	    while (!done) {
		r++;
		if (snk[r] && (my_rank_last <= snk_rnk[r]))
		    done = 1;
	    }
	
	    /* if (l<r) */
	    my_extra_elem_idx += l_size;
	    my_extra_elems    -= l_size;

	    while (l+1 < r) {
		l++;
		if (snk[l]) {
		    l_offset = locN[l];
		    l_size   = -d[l];

		    bulk_store(((DATA_TYPE *global)(A+l))+l_offset,
			       locA + my_extra_elem_idx,
			       l_size*sizeof(DATA_TYPE));
		    my_extra_elem_idx += l_size;
		    my_extra_elems    -= l_size;
		}
	    }

	    l_offset = locN[r];
	    l_size   = my_extra_elems;

	    bulk_store(((DATA_TYPE *global)(A+r))+l_offset,
		       locA + my_extra_elem_idx,
		       l_size*sizeof(DATA_TYPE));
	}

    }

    all_store_sync(); 
    barrier();

    *(N+MYPROC) = q[MYPROC];

    barrier();

#undef DATA_TYPE
}

/*************************************************************/
int all_select_b(int M, int *spread A,
		 int *spread N, int total_n, int i)
/*************************************************************/
{

    int
	meds[PROCS],
	mom,
	j, k, l,
	t, v,
	result,
	*seqval;

    int *locA,
	*locN;

#if PROFILE_DETAILED
    double secs;
#endif
    
    barrier();
#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_b: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    barrier();
#endif

    if (total_n < THRESH_SEL) {
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
        assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (int*)malloc(total_n*sizeof(int));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (int *global)(A+j), locN[j]*sizeof(int));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    fastsort(seqval,t);
	    result = seqval[i];
	    free(seqval);
	}
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_bcast(result);
	free(locN);
      }
    else {

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
	all_LoadBalance_b_i(M, A, N, total_n);
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	barrier();
#endif    

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (int *)(A+MYPROC);

	fastsort(locA, *locN); 
	*meds = locA[(*locN + 1) >> 1];

	all_gather(meds);
	on_one fastsort(meds, PROCS);

	mom     = all_bcast(meds[PROCS>>1]);
	k       = binsearch_gt(locA, 0, *locN - 1, mom);
	l       = k+1;
	t       = all_reduce_to_all_add(l);
	v       = (i <= t);
	j       = (v ? t : (total_n - t) );

	total_n = j;

	if (v)
	    *locN = l ;
	else {
	    *locN -= l;
#if 1
	    bcopy(locA + l, locA, (*locN)*sizeof(int));
#else
	    for (j=0 ; j< *locN ; j++)
		locA[j] = locA[j+l];
#endif
	}

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Sel2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_select_b(M, A, N, total_n, (v ? i : (i-t)));
    }

    return (result);
}

int all_select_median_b(int M, int *spread A)
{
    int *spread N;
    int t;

    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);
    
    *(N+MYPROC) = M;
    t = PROCS*M;
    
#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_median_b (%d)\n",MYPROC,t);
	fflush(outfile);
    }
    barrier();
#endif

    return all_select_b(M, A, N, t, (t+1)>>1);
}

int all_select_median_unbalanced_b(int M, int *spread A,
				   int *spread N, int total_n)
{
    return all_select_b(M, A, N, total_n, (total_n+1)>>1);
}


/*                                                                    tab:8
 *
 * select_b_nas.sc - Parallel Selection Algorithm for NAS IS
 *
 * 
 * "Copyright (c) 1995,1996 The Regents of the University of Maryland.
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
 * Filename:            select_b_nas.sc
 * History:
 */

#include "select_b_nas.h"

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


/*************************************************************/
static int all_select_b_nas(int M, int *spread A,
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
	fprintf(outfile,"PE %2d: all_select_b_nas: total_n: %d  i: %d\n",
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
	    fastsort_19(seqval,t);
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
	all_LoadBalance_b(M, A, N, total_n);
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

	fastsort_19(locA, *locN); 
	*meds = locA[(*locN + 1) >> 1];

	all_gather(meds);
	on_one fastsort_19(meds, PROCS);

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
	    bcopy(locA+l , locA, (*locN) * sizeof(int));
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
	result = all_select_b_nas(M, A, N, total_n, (v ? i : (i-t)));
    }

    return (result);
}

int all_select_median_b_nas(int M, int *spread A)
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
	fprintf(outfile,"PE %2d: all_select_median_b_nas (%d)\n",MYPROC,t);
	fflush(outfile);
    }
    barrier();
#endif

    return all_select_b_nas(M, A, N, t, (t+1)>>1);
}



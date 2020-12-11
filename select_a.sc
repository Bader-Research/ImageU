/*                                                                    tab:8
 *
 * select_a.sc - Parallel Selection Algorithm
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
 * Filename:            select_a.sc
 * History:
 */

#include "select_a.h"

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
void all_LoadBalance_a(int M, int *spread A, int *spread LB,
		     int *spread N, int total_n)
/*************************************************************/
{
    register int
/*	get_proc, */
	j;
    
    int locN[PROCS],
	ps[PROCS],
	q,
	l, r,
	lc, rc,
	l_n, r_n,
	l_off, t,
	*locA,
	*locLB;

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_LoadBalance: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    barrier();
#endif

    locA = (int *)(A + MYPROC);
    locLB = (int *)(LB + MYPROC);

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

    ps[0] = locN[0];

    for (j=1 ; j<PROCS ; j++)
	ps[j] = ps[j-1] + locN[j];

    q  = DIVPROCS(total_n);

    if (q*PROCS < total_n)
	q++;

    lc = MYPROC * q;
    rc = (MYPROC+1) * q;

    on_proc(PROCS-1) {
	t = (q*PROCS) - total_n;
	q = q-t;             /* since total_n > P^2, this should be +ve */
	rc = total_n;
    }
	
    l = 0;
    while (lc >= ps[l])
	l++;

    r = 0;
    while (rc >  ps[r])
	r++;

    l_n   = min(q , ps[l] - lc);
    l_off = lc - (ps[l] - locN[l]);

    bulk_read(locLB, (int *global)(A+l)+l_off, l_n*sizeof(int));
    t = l_n;
    
    if (l<r) {
	while (l+1 < r) {
	    l++;
	    bulk_read(locLB+t, (int *global)(A+l), locN[l]*sizeof(int));
	    t += locN[l];
	}

	r_n   = rc - ps[r-1];
	bulk_read(locLB + t, (int *global)(A+r), r_n*sizeof(int));
	t += r_n;

    }

    if (t!=q)
	fprintf(stderr,"ERROR: t != q\n");

    *(N+MYPROC) = t;

    barrier();
}


/*************************************************************/
int all_select_a(int M, int *spread A, int *spread LB,
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
	*locLB,
	*locN;

#if PROFILE_DETAILED
    double secs;
#endif
    
    barrier();
#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select:     total_n: %d  i: %d\n",
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
	  sync();
	  fastsort(seqval,t);
	  result = seqval[i];
	  free(seqval);
	}
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Seq1  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq1:  %9.6f\n",secs);
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
	all_LoadBalance_a(M, A, LB, N, total_n);
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"LB1   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB1:   %9.6f\n",secs);
	}
	barrier();
#endif    

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (int *)(A+MYPROC);
	locLB = (int *)(LB+MYPROC);

	fastsort(locLB, *locN); 
	*meds = locLB[(*locN + 1) >> 1];

	all_gather(meds);
	on_one fastsort(meds, PROCS);

	mom     = all_bcast(meds[PROCS>>1]);
	k       = binsearch_gt(locLB, 0, *locN - 1, mom);
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
	    bcopy(locLB + l, locLB, (*locN)*sizeof(int));
#else
	    for (j=0 ; j< *locN ; j++)
		locLB[j] = locLB[j+l];
#endif	    
	}

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Sel1  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Sel1:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_select_a(M, LB, A, N, total_n, (v ? i : (i-t)));

    }

    return (result);
}

int all_select_median_a(int M, int *spread A)
{
    int *spread N;
    int *spread LB; 
    int t;

    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);
    
    *(N+MYPROC) = M;
    t = PROCS*M;

    LB = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(LB);
    
    return all_select_a(M, A, LB, N, t, (t+1)>>1);
}

int all_select_median_unbalanced_a(int M, int *spread A,
				   int *spread N, int total_n)
{
    int *spread LB;

    LB = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(LB);

    return all_select_a(M, A, LB, N, total_n, (total_n+1)>>1);
}

int all_select_median_unbalanced_make_data(int M, int *spread A,
					   int *spread N, int total_n)
{

    FILE *datafile;
    int i,j,n;

    int arr[M];

    barrier();

    on_one {

	datafile = fopen("dfile.dat","w+");

	fprintf(datafile,"%d\n",M);
	fprintf(datafile,"%d\n",PROCS);
	for (i=0 ; i<PROCS ; i++) {
	    n = *(N+i);
	    fprintf(datafile,"%d\n",n);
	    bulk_read(arr, (int *global)(A+i), n*sizeof(int));
	    for (j=0 ; j<n ; j++)
		fprintf(datafile,"%d ",arr[j]);
	    fprintf(datafile,"\n");
	}

	fclose(datafile);
    }

    barrier();
    
}

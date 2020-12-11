/*                                                                    tab:8
 *
 * select_c.sc - Parallel Selection Algorithm
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
 * Filename:            select_c.sc
 * History:
 */

#include "select_c.h"
#include "select_b.h"
#include "select_seq.h"

static int THRESH_SEL;
#define THRESH_SEL_INIT (PROCS*PROCS)
static int SELECT_MIN_MAX;
#define SELECT_MIN_MAX_INIT (PROCS)

/* NOTE:
 *  You must call all_select_c_init() before a call to all_select_c...()
 */

#define PROFILE_DETAILED 0
#define PROFILE_GENERAL  1

#if PROFILE_GENERAL
static double time_select, time_lb;
#endif

int all_select_c_i(int M, int *spread A,
		     int *spread N, int total_n, int i) {

#define DATA_TYPE int

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_i(seqval, t, i);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_i(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (DATA_TYPE *)(A+MYPROC);

	*meds = select_mom_i(locA, *locN, *locN >> 1);

	all_gather(meds);

	on_one mom = select_mom_i(meds, PROCS, PROCS>>1);

	mom     = all_bcast(mom);
	k       = partition_unk_piv_i(locA, *locN, mom);
	t       = all_reduce_to_all_add(k);

	if (t==total_n) {
	    result = mom;
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		*locN    = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
		*locN   -= k;
		
		bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
	    }

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	    result = all_select_c_i(M, A, N, total_n, i);
	}
    }

#if PROFILE_GENERAL
#endif
    return (result);

#undef DATA_TYPE
}

double all_select_c_d(int M, double *spread A,
		      int *spread N, int total_n, int i) {

#define DATA_TYPE double

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;
    
#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    barrier();
#endif

    if (total_n <= THRESH_SEL) {
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_d(seqval, t, i);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_bcast_d(result);
	free(locN);
      }
    else {

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_d(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (DATA_TYPE *)(A+MYPROC);

	*meds = select_mom_d(locA, *locN, *locN >> 1);

	all_gather_d(meds);

	on_one mom = select_mom_d(meds, PROCS, PROCS>>1);

	mom     = all_bcast_d(mom);
	k       = partition_unk_piv_d(locA, *locN, mom);
	t       = all_reduce_to_all_add(k);

	if (t==total_n) {
	    result = mom;
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		*locN    = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
		*locN   -= k;

		bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
	    }

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	    result = all_select_c_d(M, A, N, total_n, i);
	}
    }

#if PROFILE_GENERAL
#endif

    return (result);

#undef DATA_TYPE
}

int all_select_c_alloc_i_orig(int M, int *spread A,
			      int *spread N, int total_n, int i,
			      int *sublist) {

#define DATA_TYPE int

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_alloc_i(seqval, t, i,sublist);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_i(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (DATA_TYPE *)(A+MYPROC);

	*meds = select_mom_alloc_i(locA, *locN, *locN >> 1, sublist);

	all_gather(meds);

	on_one mom = select_mom_alloc_i(meds, PROCS, PROCS>>1,sublist);

	mom     = all_bcast(mom);
	k       = partition_unk_piv_i(locA, *locN, mom);
	t       = all_reduce_to_all_add(k);

	if (t==total_n) {
	    result = mom;
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		*locN    = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
		*locN   -= k;

		bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
	    }

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	    result = all_select_c_alloc_i(M, A, N, total_n, i,sublist);
	}
    }

    return (result);

#undef DATA_TYPE
}

int all_select_min_alloc_i(int M, int *spread A,
			   int *spread N, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    int
	j, t;

    DATA_TYPE
	piv,
	result,
	*locA,
	*seqval;
    
    locA = (DATA_TYPE *)(A+MYPROC);
    t = *(N+MYPROC);

    if (i==1) {
	result = locA[0];
	for (j=1 ; j<t ; j++)
	    result = min(result,locA[j]);
	result = all_reduce_to_all_min(result);
    }
    else{
	piv = select_mom_alloc_i(locA, t, i, sublist);
	j   = partition_i(locA, t, piv);

	barrier();

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 X*sizeof(DATA_TYPE));
		t += X;
	    }
	    sync();
	    result = select_mom_alloc_i(seqval, t, i,sublist);
	    free(seqval);
	}

	result = all_bcast(result);
    }

    return (result);
}

int all_select_max_alloc_i(int M, int *spread A,
			   int *spread N, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    int
	j, t,
	offset,
	last;

    DATA_TYPE
	piv,
	result,
	*locA,
	*seqval;
    

    locA = (DATA_TYPE *)(A+MYPROC);
    t = *(N+MYPROC);

    if (i==total_n) {
	result = locA[0];
	for (j=1 ; j<t ; j++)
	    result = max(result,locA[j]);
	result = all_reduce_to_all_max(result);
    }
    else {
	piv = select_mom_alloc_i(locA, t, i, sublist);
	j   = partition_i(locA, t, piv);

	barrier();

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    offset = M - X;

	    t = 0;
	    for (j=0; j<PROCS-1 ; j++) {
		bulk_get(seqval+t,
			 (DATA_TYPE *global)(A+j) + offset,
			 X*sizeof(DATA_TYPE));
		t += X;
	    }

	    offset = *(N+(PROCS-1)) - X;
	    bulk_get(seqval+t,
		     (DATA_TYPE *global)(A+(PROCS-1)) + offset,
		     X*sizeof(DATA_TYPE));
	    t += X;

	    sync();
	    j = i - (total_n - X*PROCS);
	    result = select_mom_alloc_i(seqval, t, j,sublist);
	    free(seqval);
	}

	result = all_bcast(result);
    }
    return (result);

#undef DATA_TYPE
}

double all_select_min_alloc_d(int M, double *spread A,
			      int *spread N, int total_n, int i,
			      double *sublist, int X) {
#define DATA_TYPE double

    int
	j, t;

    DATA_TYPE
	piv,
	result,
	*locA,
	*seqval;
    

    locA = (DATA_TYPE *)(A+MYPROC);
    t = *(N+MYPROC);
    
    if (i==1) {
	result = locA[0];
	for (j=1 ; j<t ; j++)
	    result = min(result,locA[j]);
	result = all_reduce_to_all_dmin(result);
    }
    else {
	piv = select_mom_alloc_d(locA, t, i, sublist);
	j   = partition_d(locA, t, piv);

	barrier();

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 X*sizeof(DATA_TYPE));
		t += X;
	    }
	    sync();
	    result = select_mom_alloc_d(seqval, t, i,sublist);
	    free(seqval);
	}

	result = all_bcast_d(result);
    }
    return (result);

#undef DATA_TYPE
}

double all_select_max_alloc_d(int M, double *spread A,
			   int *spread N, int total_n, int i,
			   double *sublist, int X) {
#define DATA_TYPE double

    int
	j, t,
	offset,
	last;

    DATA_TYPE
	piv,
	result,
	*locA,
	*seqval;
    

    locA = (DATA_TYPE *)(A+MYPROC);
    t = *(N+MYPROC);

    if (i==total_n) {
	result = locA[0];
	for (j=1 ; j<t ; j++)
	    result = max(result,locA[j]);
	result = all_reduce_to_all_dmax(result);
    }
    else {
	piv = select_mom_alloc_d(locA, t, i, sublist);
	j   = partition_d(locA, t, piv);

	barrier();

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    offset = M - X;

	    t = 0;
	    for (j=0; j<PROCS-1 ; j++) {
		bulk_get(seqval+t,
			 (DATA_TYPE *global)(A+j) + offset,
			 X*sizeof(DATA_TYPE));
		t += X;
	    }

	    offset = *(N+(PROCS-1)) - X;
	    bulk_get(seqval+t,
		     (DATA_TYPE *global)(A+(PROCS-1)) + offset,
		     X*sizeof(DATA_TYPE));
	    t += X;

	    sync();
	    j = i - (total_n - X*PROCS);
	    result = select_mom_alloc_d(seqval, t, j,sublist);
	    free(seqval);
	}

	result = all_bcast_d(result);
    }
    return (result);

#undef DATA_TYPE
}

int all_select_c_alloc_i(int M, int *spread A,
			 int *spread N, int total_n, int i,
			 int *sublist) {

#define DATA_TYPE int

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_alloc_i(seqval, t, i,sublist);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_i(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	barrier();
#endif    

	if (i <= SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    result = all_select_min_alloc_i(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"MIN   n: %12d  ",total_n);
		fprintf(outfile,"Time for MIN:   %9.6f\n",secs);
	    }
	    barrier();
#endif
	}
	else if (i >= total_n - SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    result = all_select_max_alloc_i(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"MAX   n: %12d  ",total_n);
		fprintf(outfile,"Time for MAX:   %9.6f\n",secs);
	    }
	    barrier();
#endif
	}
	else {
	
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    locN = (int *)(N+MYPROC);
	    locA = (DATA_TYPE *)(A+MYPROC);

	    *meds = select_mom_alloc_i(locA, *locN, *locN >> 1, sublist);

	    all_gather(meds);

	    on_one mom = select_mom_alloc_i(meds, PROCS, PROCS>>1,sublist);

	    mom     = all_bcast(mom);
	    k       = partition_unk_piv_i(locA, *locN, mom);
	    t       = all_reduce_to_all_add(k);

	    if (t==total_n) {
		result = mom;
#if PROFILE_GENERAL
		barrier();
		secs_gen = get_seconds() - secs_gen;
		time_select += secs_gen;
		barrier();
#endif    
#if PROFILE_DETAILED
		barrier();
		secs = get_seconds() - secs;
		on_one {
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		barrier();
#endif    
	    }
	    else {
		if (i <= t) {
		    total_n  = t;
		    *locN    = k;
		}
		else {
		    total_n -= t;
		    i       -= t;
		    *locN   -= k;

		    bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
		}

#if PROFILE_GENERAL
		barrier();
		secs_gen = get_seconds() - secs_gen;
		time_select += secs_gen;
		barrier();
#endif    
#if PROFILE_DETAILED
		barrier();
		secs = get_seconds() - secs;
		on_one {
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		barrier();
#endif    
		result = all_select_c_alloc_i(M, A, N, total_n, i,sublist);
	    }
	}
    }

    return (result);

#undef DATA_TYPE
}

double all_select_c_alloc_d_orig(int M, double *spread A,
				 int *spread N, int total_n, int i,
				 double* sublist) {

#define DATA_TYPE double

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;
    
#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    barrier();
#endif

    if (total_n <= THRESH_SEL) {
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_alloc_d(seqval, t, i,sublist);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_bcast_d(result);
	free(locN);
      }
    else {

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_d(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
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
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	locN = (int *)(N+MYPROC);
	locA = (DATA_TYPE *)(A+MYPROC);

	*meds = select_mom_alloc_d(locA, *locN, *locN >> 1, sublist);

	all_gather_d(meds);

	on_one mom = select_mom_alloc_d(meds, PROCS, PROCS>>1, sublist);

	mom     = all_bcast_d(mom);
	k       = partition_unk_piv_d(locA, *locN, mom);
	t       = all_reduce_to_all_add(k);

#if 0
	barrier();
	for (j=0 ; j<PROCS ; j++) {
	  barrier();
	  on_proc(j) {
	    fprintf(outfile,"total: %d  k: %d  t:%d  i: %d\n",
		    total_n,k,t,i);
	    fflush(outfile);
	  }
	  barrier();
	}
	barrier();
#endif

	if (t==total_n) {
	    result = mom;
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		*locN    = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
		*locN   -= k;
		
		bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
	    }
	    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    barrier();
#endif    
	    result = all_select_c_alloc_d(M, A, N, total_n, i,sublist);
	}
    }

    return (result);

#undef DATA_TYPE
}


double all_select_c_alloc_d(int M, double *spread A,
			    int *spread N, int total_n, int i,
			    double* sublist) {

#define DATA_TYPE double

    int	j, k, 
	t,
	*locN;

    DATA_TYPE
	meds[PROCS],
	mom,
	result,
	*seqval,
	*locA;
    
#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

    barrier();

#if PROFILE_DETAILED
    barrier();
    on_one {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    barrier();
#endif

    if (total_n <= THRESH_SEL) {
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    
	locN = (int*) malloc(PROCS * sizeof(int));
	assert_malloc(locN);

	locN[0] = *(N+MYPROC);
	all_gather(locN);

	on_one {
	    seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	    assert_malloc(seqval);

	    t = 0;
	    for (j=0; j<PROCS ; j++) {
		bulk_get(seqval+t, (DATA_TYPE *global)(A+j),
			 locN[j]*sizeof(DATA_TYPE));
		t += locN[j];
	    }
#if 1
	    if (t != total_n)
		fprintf(stderr,"ERROR: t != total_n\n");
#endif	
	    sync();
	    result = select_mom_alloc_d(seqval, t, i,sublist);
	    free(seqval);
	}
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_select += secs_gen;
	barrier();
#endif    
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
	}
	barrier();
#endif    
	result = all_bcast_d(result);
	free(locN);
      }
    else {

#if PROFILE_DETAILED
	barrier();
	secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds();
#endif    

	all_LoadBalance_b_d(M, A, N, total_n);

#if PROFILE_GENERAL
	barrier();
	secs_gen = get_seconds() - secs_gen;
	time_lb += secs_gen;
	barrier();
#endif    
#if PROFILE_DETAILED
	barrier();
	secs = get_seconds() - secs;
	on_one {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	barrier();
#endif    

	if (i <= SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    result = all_select_min_alloc_d(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"MIN   n: %12d  ",total_n);
		fprintf(outfile,"Time for MIN:   %9.6f\n",secs);
	    }
	    barrier();
#endif
	}
	else if (i >= total_n - SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    result = all_select_max_alloc_d(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds() - secs_gen;
	    time_select += secs_gen;
	    barrier();
#endif    
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds() - secs;
	    on_one {
		fprintf(outfile,"MAX   n: %12d  ",total_n);
		fprintf(outfile,"Time for MAX:   %9.6f\n",secs);
	    }
	    barrier();
#endif
	}
	else {
	
#if PROFILE_DETAILED
	    barrier();
	    secs = get_seconds();
#endif    
#if PROFILE_GENERAL
	    barrier();
	    secs_gen = get_seconds();
#endif    

	    locN = (int *)(N+MYPROC);
	    locA = (DATA_TYPE *)(A+MYPROC);

	    *meds = select_mom_alloc_d(locA, *locN, *locN >> 1, sublist);

	    all_gather_d(meds);

	    on_one mom = select_mom_alloc_d(meds, PROCS, PROCS>>1, sublist);
	    
	    mom     = all_bcast_d(mom);
	    k       = partition_unk_piv_d(locA, *locN, mom);
	    t       = all_reduce_to_all_add(k);

#if 0
	    barrier();
	    for (j=0 ; j<PROCS ; j++) {
		barrier();
		on_proc(j) {
		    fprintf(outfile,"total: %d  k: %d  t:%d  i: %d\n",
			    total_n,k,t,i);
		    fflush(outfile);
		}
		barrier();
	    }
	    barrier();
#endif

	    if (t==total_n) {
		result = mom;
#if PROFILE_GENERAL
		barrier();
		secs_gen = get_seconds() - secs_gen;
		time_select += secs_gen;
		barrier();
#endif    
#if PROFILE_DETAILED
		barrier();
		secs = get_seconds() - secs;
		on_one {
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		barrier();
#endif    
	    }
	    else {
		if (i <= t) {
		    total_n  = t;
		    *locN    = k;
		}
		else {
		    total_n -= t;
		    i       -= t;
		    *locN   -= k;

		    bcopy(locA+k, locA, (*locN)*sizeof(DATA_TYPE));
		}

#if PROFILE_GENERAL
		barrier();
		secs_gen = get_seconds() - secs_gen;
		time_select += secs_gen;
		barrier();
#endif    
#if PROFILE_DETAILED
		barrier();
		secs = get_seconds() - secs;
		on_one {
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		barrier();
#endif    
		result = all_select_c_alloc_d(M, A, N, total_n, i,sublist);
	    }
	}
    }

    return (result);

#undef DATA_TYPE
}

void all_select_c_init() {
    THRESH_SEL = THRESH_SEL_INIT;
    SELECT_MIN_MAX = SELECT_MIN_MAX_INIT;
}

int all_select_median_c_i(int M, int *spread A) {

#define DATA_TYPE int

    int *spread N;
    DATA_TYPE *sublist;
    DATA_TYPE result;
    int t;


    all_select_c_init();

    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);

    sublist = (DATA_TYPE *)malloc(M*sizeof(DATA_TYPE));
    assert_malloc(sublist);
    
    *(N+MYPROC) = M;
    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_lb = 0.0;
#endif    
    result = all_select_c_alloc_i(M, A, N, t, (t>>1)+1,sublist);
#if PROFILE_GENERAL
    barrier();
    on_one 
	fprintf(outfile,"Select_c general profile: Select: %9.6f  LB: %9.6f\n",
		time_select, time_lb);
    barrier();
#endif
    free(sublist);
    return (result);

#undef DATA_TYPE
}

double all_select_median_c_d(int M, double *spread A)
{

#define DATA_TYPE double

    int *spread N;
    DATA_TYPE *sublist;
    DATA_TYPE result;
    int t;

    all_select_c_init();

    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);
    
    sublist = (DATA_TYPE *)malloc(M*sizeof(DATA_TYPE));
    assert_malloc(sublist);
    
    *(N+MYPROC) = M;
    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_lb = 0.0;
#endif    
    result = all_select_c_alloc_d(M, A, N, t, (t>>1) + 1,sublist);
#if PROFILE_GENERAL
    barrier();
    on_one 
	fprintf(outfile,"Select_c general profile: Select: %9.6f  LB: %9.6f\n",
		time_select, time_lb);
    barrier();
#endif
    free(sublist);
    return (result);

#undef DATA_TYPE
}

int all_select_median_unbalanced_c_i(int M, int *spread A,
				       int *spread N, int total_n) {
    all_select_c_init();
    return all_select_c_i(M, A, N, total_n, (total_n>>1) + 1);

}

double all_select_median_unbalanced_c_d(int M, double *spread A,
					     int *spread N, int total_n) {
    all_select_c_init();
    return all_select_c_d(M, A, N, total_n, (total_n>>1) + 1);
}


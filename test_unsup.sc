/*                                                                    tab:8
 *
 * test_unsup.sc - Unsupported testing routines
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
 * Creation Date:       May 2, 1995
 * Filename:            test_unsup.sc
 * History:
 */

#include "test_unsup.h"
#include "select_a.h"
#include "select_b.h"
#include "select_c.h"
#include "inputs_select.h"
#include "inputs_median.h"
#include "math_func.h"
#include "droute.h"
#include "radixsort.h"

#if 1
#define RADIX1
#endif

#if (defined(CM5) || defined(MEIKO) || defined(SP2) || defined(T3D))
#define RADIX3
#endif

#if (defined(SP2))
#define RADIX4
#endif

#define RADIX5
#define RADIX6

#define INPUT_HREL 1
#define INPUT_GGRP 2
#define GGRP_H     1
#define GGRP_G     8
#define GGRP_T     8

#define BENCH_I_N   1
#define BENCH_I_R   1
#define BENCH_I_S   1
#define BENCH_I_Z   1
#define BENCH_I_C   1
#define BENCH_I_N20 1

#define TEST_MACH_LIMIT_LOG 20
#define TEST_MACH_LIMIT (1<<TEST_MACH_LIMIT_LOG)
#define NAS_REPEAT    3
#define WARM_UP_COUNT 4
#define WARM_UP       1
#define DO_NAS        0
#define DO_DOUBLES    0
#define DO_HREL       1
#define DO_SELECT     0

void all_unsupported_select_d()
{
    register int
	i, j;

    double result;

    double *locA;
    double *spread A;

    double secs;

    A = all_spread_malloc(PROCS,TEST_MACH_LIMIT*sizeof(double));
    assert_spread_malloc(A);
    locA = (double *)(A + MYPROC);

    RNG_init();

#if WARM_UP
    for (i=0 ; i<WARM_UP_COUNT ; i++) {
      barrier();
      j = TEST_MACH_LIMIT;
      fill_random_d(j, locA);
      barrier();
      secs = get_seconds();
      result = all_select_median_c_d(j, A);  
      barrier();
      secs = get_seconds() - secs;
      on_one {
	fprintf(outfile,"dis%d: %9.6f %g\n",i,secs,result);
	fflush(outfile);
      }
      barrier();
    }
#endif
    
#if 1
    for (i = (1<<12) ; i <= TEST_MACH_LIMIT ; i = (i<<1) ) {

	barrier();

/* Test median random double */

	fill_random_d(i, locA);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [R] doubles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();

/* Test median gaussian double */

	fill_gaussian_d(i, locA);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [G] doubles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();

#if 0
/* Test median 2-g double */

	fill_g_group_d(i, locA, 2);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [2-G] dbles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();
#endif
	
#if 0
/* Test median bucket double */

	fill_bucket_d(i, locA);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [B] doubles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();
#endif
	
/* Test median staggered double */

	fill_staggered_d(i, locA);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [S] doubles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();
	
#if 0
/* Test median zero entropy double */

	fill_zero_d(i, locA);

	barrier();
	secs = get_seconds();
	result = all_select_median_c_d(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,
		    "Median_c [%d]::[%d] [Z] doubles is: %g  Time: %9.6f\n",
		    PROCS, i, result, secs);
	    fflush(outfile);
	}
	barrier();
#endif	
    }
#endif
    
    all_spread_free(A);
}


void all_unsupported_select()
{
    register int
	i, j;

    int
	n,
	*locA;
    int *spread N;
    int *spread A;

    double secs;

#if DO_DOUBLES
    all_unsupported_select_d();
#endif

    N = all_spread_malloc(PROCS,sizeof(int));
    assert_spread_malloc(N);

    A = all_spread_malloc(PROCS,TEST_MACH_LIMIT*sizeof(int));
    assert_spread_malloc(A);
    locA = (int *)(A + MYPROC);

    RNG_init();

#if WARM_UP
    for (j=0 ; j<WARM_UP_COUNT ; j++) {
      barrier();

      i = TEST_MACH_LIMIT;
      n = (i>>2);
      n = all_fill_array_linear(i, A, N, n);
      secs = get_seconds();
      all_LoadBalance_b(i,A,N,n);
      secs = get_seconds() - secs;
      on_one {
	fprintf(outfile,"dis%d n: %12d LB_B: %9.6f\n",j,n,secs);
	fflush(outfile);
      }

      barrier();
    }
#endif

#if DO_NAS
    barrier();

    if (TEST_MACH_LIMIT_LOG + PROCSLOG >= 23) {
	int *spread B;

	i = DIVPROCS(1<<23);

	B = all_spread_malloc(PROCS,i*sizeof(int));
	assert_spread_malloc(B);

	all_fill_nas(i, B);

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  barrier();
	  bcopy((int*)(B+MYPROC),(int*)(A+MYPROC),i*sizeof(int));

	  barrier();
	  secs = get_seconds();
	  j = all_select_median_b(i, A);
	  barrier();
	  secs = get_seconds() - secs;

	  on_one {
	    fprintf(outfile,"Median_b of [%d]::[%d] NAS is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  barrier();
	  bcopy((int*)(B+MYPROC),(int*)(A+MYPROC),i*sizeof(int));

	  barrier();
	  secs = get_seconds();
	  j = all_select_median_b_nas(i, A);
	  barrier();
	  secs = get_seconds() - secs;

	  on_one {
	    fprintf(outfile,"Median_b_nas of [%d]::[%d] NAS is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  barrier();
	  bcopy((int*)(B+MYPROC),(int*)(A+MYPROC),i*sizeof(int));
	
	  barrier();
	  secs = get_seconds();
	  j = all_select_median_c(i, A);
	  barrier();
	  secs = get_seconds() - secs;
	  
	  on_one {
	    fprintf(outfile,"Median_c of [%d]::[%d] NAS is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}
	all_spread_free(B);
      }
    barrier();
#endif    

#if 1
    for (i = (1<<12) ; i <= TEST_MACH_LIMIT ; i = (i<<1) ) {

	barrier();

/* Test median same */

	fill_same(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_a(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_a of [%d]::[%d] same ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();
	
/* Test median same */

	fill_same(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_b(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_b of [%d]::[%d] same ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median same */

	fill_same(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_c(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_c of [%d]::[%d] same ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();
	
/* Test median linear */

	fill_linear(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_a(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_a of [%d]::[%d] linear ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median linear */

	fill_linear(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_b(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_b of [%d]::[%d] linear ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median linear */

	fill_linear(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_c(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_c of [%d]::[%d] linear ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median random */

	fill_random(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_a(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_a of [%d]::[%d] random ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median random */

	fill_random(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_b(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_b of [%d]::[%d] random ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

/* Test median random */

	fill_random(i, locA);

	barrier();
	secs = get_seconds();
	j = all_select_median_c(i, A);
	barrier();
	secs = get_seconds() - secs;

	on_one {
	    fprintf(outfile,"Median_c of [%d]::[%d] random ints is: %9d  Time: %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	}
	barrier();

    }
#endif
    
#if 0
/* Test median band5-512 */

    i = (1<<13);
    n = all_fill_array(i, A, N, "im/band5.512.dat");

    barrier();
    secs = get_seconds();
    j = all_select_median_unbalanced_a(i, A, N, n);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,"Median_a of [%d]::[%d] band 5 (512x512) problem: %9d  Time: %9.6f\n",
		PROCS, i, j,secs);
	fflush(outfile);
    }
    barrier();

/* Test median band5-512 */

    i = (1<<13);
    n = all_fill_array(i, A, N, "im/band5.512.dat");

    barrier();
    secs = get_seconds();
    j = all_select_median_unbalanced_b(i, A, N, n);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,"Median_b of [%d]::[%d] band 5 (512x512) problem: %9d  Time: %9.6f\n",
		PROCS, i, j,secs);
	fflush(outfile);
    }
    barrier();

/* Test median b5-1024 */

    i = (1<<15);
    n = all_fill_array(i, A, N, "im/b5.1024.dat");

    barrier();
    secs = get_seconds();
    j = all_select_median_unbalanced_a(i, A, N, n);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,"Median_a of [%d]::[%d] b5 (1024x1024) problem: %9d  Time: %9.6f\n",
		PROCS, i, j,secs);
	fflush(outfile);
    }
    barrier();

/* Test median band5-512 */

    i = (1<<15);
    n = all_fill_array(i, A, N, "im/b5.1024.dat");

    barrier();
    secs = get_seconds();
    j = all_select_median_unbalanced_b(i, A, N, n);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,"Median_b of [%d]::[%d] b5 (1024x1024) problem: %9d  Time: %9.6f\n",
		PROCS, i, j,secs);
	fflush(outfile);
    }
    barrier();

#endif

    all_spread_free(N);
    all_spread_free(A);
}


void check_sort(int q, int A[PROCS]::[q]) {
    int i;
    for (i=1; i<q ; i++) 
	if (A[MYPROC][i-1] > A[MYPROC][i]) {
	    fprintf(stderr,
		    "ERROR: PE%2d: A[%2d][%6d] > A[%2d][%6d] (%6d,%6d)\n",
		    MYPROC,MYPROC,i-1,MYPROC,i,A[MYPROC][i-1],A[MYPROC][i]);
	    exit(1);
	}
    if (MYPROC<PROCS-1)
	if (A[MYPROC][q-1]>A[MYPROC+1][0]) {
	    fprintf(stderr,
		    "ERROR: PE%2d: A[%2d][%6d] > A[%2d][%6d] (%6d,%6d)\n",
		    MYPROC,MYPROC,q-1,MYPROC+1,0,A[MYPROC][q-1],A[MYPROC+1][0]);
	    exit(1);
	}
}

double calc_bandwidth_hrel(int M, int h, double secs) {

    double bw;
    int b1, b2;

    b1 = M/PROCS + PROCS + 1;
    b2 = h*(M/PROCS) + PROCS;
    bw = (1.0 - (1.0/(double)PROCS));    /* only (p-1)/p is sent via network */
    bw *= (double)(PROCS*(b1 + b2));     /* elements of data sent            */
    bw *= 2.0;                           /*              and received        */
    bw *= (double)sizeof(intpair_t);     /* bytes per element                */
    bw /= secs;                          /* bytes per second                 */
    bw /= 1000000.0;                     /* MB/s                             */
    return (bw);    
}

void all_unsupported_hrel1(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel1  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel1  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel1  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel1  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_hrel2(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det2(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel2  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det2(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel2  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det2(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel2  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det2(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel2  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_hrel3(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det3(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel3  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det3(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel3  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det3(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel3  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det3(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel3  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_hrel4(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det4(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel4  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det4(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel4  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det4(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel4  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_h_rel_det4(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"hrel4  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
      	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_linperm(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_linear_perm(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linperm  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_linear_perm(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linperm  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_linear_perm(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linperm  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_linear_perm(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linperm  h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_linperm_H(int M, int test_input) {
    int *spread Addr;
    int *spread Keys;
    int *spread Routed;

    int h;

    double secs;

    Addr = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Addr);
    Keys = all_spread_malloc(PROCS,M*sizeof(int));
    assert_spread_malloc(Keys);

/*****************************************************************/
/*  Test hrel 1                                                  */
/*****************************************************************/

    barrier();
    h = 1;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_single_phase(M, h, Keys, Addr, Routed); 
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linpermH h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();
    
/*****************************************************************/
/*  Test hrel 2                                                  */
/*****************************************************************/

    barrier();
    h = 2;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_single_phase(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linpermH h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 4                                                  */
/*****************************************************************/

    barrier();
    h = 4;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_single_phase(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linpermH h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();

/*****************************************************************/
/*  Test hrel 8                                                  */
/*****************************************************************/

    barrier();
    h = 8;

    Routed = all_spread_malloc(PROCS,h*M*sizeof(int));
    assert_spread_malloc(Routed);
    barrier();

    switch (test_input) {
      case INPUT_HREL : all_fill_hrel(M, Addr, h); break;
      case INPUT_GGRP : all_fill_g_group_i(M, Addr,
					   GGRP_H, GGRP_G, GGRP_T); break;
    }
    all_fill_pid(M, Keys);
    all_fill_zero(h*M, Routed);

    barrier();
    secs = get_seconds();
    all_route_single_phase(M, h, Keys, Addr, Routed);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"linpermH h,procs,size,time,bw(MB/s) %2d %3d %12d %7.3f %7.3f \n",
		h, PROCS, M, secs, calc_bandwidth_hrel(M,h,secs));
	fflush(outfile);
    }
    barrier();
    all_spread_free(Routed);
    barrier();


/*****************************************************************/
/* Free mem                                                      */
/*****************************************************************/
    all_spread_free(Keys);
    all_spread_free(Addr);

}

void all_unsupported_hrel(int intparam, int test_input) {
    register int
	i;

    int n;

    int b[4];
    int      *locA;

    double secs;

    int       *spread A;


    if (intparam==999) {

	/* Junk runs */
	all_unsupported_hrel3(PROCS*256, test_input);
	all_unsupported_hrel3(PROCS*256, test_input);

	for (i=(1<<10) ; i < TEST_MACH_LIMIT ; i = (i<<1)) {
#if (defined(RADIX1))
	    all_unsupported_hrel1(i, test_input);
#endif
#if (defined(RADIX3))
	    all_unsupported_hrel3(i, test_input);
#endif
#if (defined(RADIX4))
	    all_unsupported_hrel4(i, test_input);
#endif	    
	}
	return;
    }

    if ((intparam>=1010)&&(intparam<=1020)) {

	/* Junk runs */
	all_unsupported_hrel1(PROCS*256, test_input);
	all_unsupported_hrel1(PROCS*256, test_input);

	i = (1<<(intparam-1000));
#if (defined(RADIX1))
	all_unsupported_hrel1(i, test_input);
#endif
#if (defined(RADIX3))
	all_unsupported_hrel3(i, test_input);
#endif
#if (defined(RADIX4))
	all_unsupported_hrel4(i, test_input);
#endif	    
	return;
    }

    if (intparam==9999) {

	/* Junk runs */
      all_unsupported_linperm(PROCS*256, test_input);
      all_unsupported_linperm(PROCS*256, test_input);

      for (i=(1<<10) ; i < TEST_MACH_LIMIT ; i = (i<<1)) {
	all_unsupported_linperm  (i, test_input);
	all_unsupported_linperm_H(i, test_input);
	all_unsupported_linperm  (i, test_input);
	all_unsupported_linperm_H(i, test_input);
      }
      
      return;
    }

    if (intparam > 0)
	i = (1<<intparam);
    else
	i = DIVPROCS(1<<23);

    A = all_spread_malloc(PROCS,i*sizeof(int));
    assert_spread_malloc(A);

    locA  = tolocal(A[MYPROC]);

#if 0
/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 8;
    b[1] = 8;
    b[2] = 8;
    b[3] = 8;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 4, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 9;
    b[1] = 9;
    b[2] = 9;
    b[3] = 5;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 4, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 10;
    b[1] = 10;
    b[2] = 10;
    b[3] = 2;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 4, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 11;
    b[1] = 11;
    b[2] = 10;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 3, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 12;
    b[1] = 12;
    b[2] = 8;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 3, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 13;
    b[1] = 13;
    b[2] = 6;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 3, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 14;
    b[1] = 14;
    b[2] = 4;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 3, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 15;
    b[1] = 15;
    b[2] = 2;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 3, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();

/*****************************************************************/
/*  Test Radixsort with NAS input                                */
/*****************************************************************/

    barrier();

    all_fill_random(i, A, 1);
    b[0] = 16;
    b[1] = 16;
    b[2] = 0;
    b[3] = 0;

    barrier();
    secs = get_seconds();
    all_radixsort_passInfo(i, A, A, 2, b);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"[%2d,%2d,%2d,%2d], (%d) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		b[0],b[1],b[2],b[3],i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if BENCH_I_N
/*****************************************************************/
/*  Radixsort with NAS input                                     */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif
    
#if (defined(RADIX2))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif
    
#if (defined(RADIX5))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX6))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) NAS            on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#endif

#if BENCH_I_R
/*****************************************************************/
/*  Radixsort with random entropy=31                             */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX2))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX5))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX6))
    all_fill_random(i, A, 1);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) random (e=31)  on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#endif

#if BENCH_I_S
/*****************************************************************/
/*  Radixsort with random entropy=6.2                            */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX2))
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))    
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif


#if (defined(RADIX5))    
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif


#if (defined(RADIX6))    
    all_fill_random(i, A, 5);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) random (e=6.2) on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#endif

#if BENCH_I_Z
/*****************************************************************/
/*  Radixsort with random entropy=0                              */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX2))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#if (defined(RADIX5))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#if (defined(RADIX6))
    all_fill_random(i, A, 0);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) random (e=0)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#endif

#if BENCH_I_C
/*****************************************************************/
/*  Radixsort with cyclic                                        */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX2))
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))    
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX5))    
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX6))    
    all_fill_cyclic(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) cyclic         on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#endif
    
#if BENCH_I_N20
/*****************************************************************/
/*  Radixsort (20 bit) with NAS input                            */
/*****************************************************************/

    barrier();

#if (defined(RADIX1))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 1);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"1 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX2))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 2);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"2 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX3))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 3);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"3 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif

#if (defined(RADIX4))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 4);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"4 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#if (defined(RADIX5))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 5);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"5 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#if (defined(RADIX6))
    all_fill_nas(i, A);

    barrier();
    secs = get_seconds();
    all_radixsort20(i, A, A, 6);
    barrier();
    secs = get_seconds() - secs;

    on_one {
	fprintf(outfile,
		"6 Radixsort of (%7d) NAS (20-bit)   on %d procs, Time: %7.3f (us/ky/P = %6.2f)\n",
		i, PROCS, secs, (secs*1000000.0)/(double)i);
	fflush(outfile);
    }
    barrier();
    check_sort(i,A);
    barrier();
#endif    

#endif
     
}

void all_unsupported(int intparam) {

#if DO_HREL
    all_unsupported_hrel(intparam, INPUT_HREL);
/*    all_unsupported_linperm(intparam); */
#endif

#if DO_SELECT
    all_unsupported_select();
#endif
}

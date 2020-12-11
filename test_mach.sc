/*                                                                    tab:8
 *
 * test_mach.sc - Test Parallel Machine Specifics
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
 * Creation Date:       April 18, 1995
 * Filename:            test_mach.sc
 * History:
 */

#include "test_mach.h"
#include "math_func.h"

#define TEST_MACH_LIMIT 256
#define TEST_MACH_STEP  1
#define TEST_MACH_TRIES 10
#define TEST_XPOSE_SIZE (PROCS*PROCS*PROCS)
#define TEST_SLOT_SIZE  (PROCS*PROCS)
#define LARGE_DOUBLE    1000000000.0
#define TEST_ALG        2
#define TEST_TIMES      (TEST_MACH_LIMIT / TEST_MACH_STEP) + 1

#define TEST_TYPES      4
/************************/
#define TYPE_INT  0
#define TYPE_ELEM 1
#define TYPE_CHPR 2
#define TYPE_EDGE 3

static void test_fill_int(int M, int arr[M])
{
    register int
	i;

    for (i=0 ; i<M ; i++)
	arr[i] = random();

}

static void test_fill_elem(int M, elem_t arr[M])
{
    register int
	i;

    for (i=0 ; i<M ; i++) {
	arr[i].color = random();
	arr[i].label = random();
	arr[i].pos   = random();
    }
}

static void test_fill_chPair(int M, chPair_t arr[M])
{
    register int
	i;

    for (i=0 ; i<M ; i++) {
	arr[i].alpha = random();
	arr[i].beta  = random();
    }
}

static void test_fill_edge(int M, edge_t arr[M])
{
    register int
	i;

    for (i=0 ; i<M ; i++) {
	arr[i].a = random();
	arr[i].b = random();
    }
}

inline static
double test_radixsort_i(int *originalList, int elems) {

#define DATA_TYPE int    

    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (int *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	radixsort(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE 
}

inline static
double test_radixsort_elem(elem_t *originalList, int elems) {

#define DATA_TYPE elem_t

    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	radixsort_elem(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE 
}

inline static
double test_radixsort_chPair(chPair_t *originalList, int elems) {

#define DATA_TYPE chPair_t    

    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	radixsort_chPair(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE 
}

inline static
double test_radixsort_edge(edge_t *originalList, int elems) {

#define DATA_TYPE edge_t

    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	radixsort_edge(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE 
}

inline static
double test_insertsort_i(int *originalList, int elems) {

#define DATA_TYPE int
    
    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	insertsort_i(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE    
}

inline static
double test_insertsort_elem(elem_t *originalList, int elems) {

#define DATA_TYPE elem_t
    
    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	insertsort_elem(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE    
}

inline static
double test_insertsort_chPair(chPair_t *originalList, int elems) {

#define DATA_TYPE chPair_t
    
    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	insertsort_chPair(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE    
}

inline static
double test_insertsort_edge(edge_t *originalList, int elems) {

#define DATA_TYPE edge_t
    
    double
	total_time, over_time;
    int loop;
    DATA_TYPE *list;

    list = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
    assert_malloc(list);
    
    total_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
	insertsort_edge(list,elems);
    }
    total_time = get_seconds() - total_time;

    over_time = get_seconds();
    for (loop=0 ; loop<TEST_MACH_TRIES ; loop++) {
	bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    }
    over_time = get_seconds() - over_time;

    total_time -= over_time;
    total_time /= (double)TEST_MACH_TRIES;
    
    free(list);
    return (total_time);

#undef DATA_TYPE    
}

void all_test_mach()
{
    register int
	i, j,
	s, get_proc;

    int	rIntBkpt    = 0,
	rElemBkpt   = 0,
	rChPairBkpt = 0,
	rEdgeBkpt   = 0,
	mXposeB     = 0;

    double
	tm[TEST_TYPES][TEST_TIMES][TEST_ALG];

    int cnt[TEST_TYPES][TEST_TIMES],
	tmax[TEST_TYPES];
   
    int      *locA;
    elem_t   *locB;
    chPair_t *locC;
    edge_t   *locD;
    int      *locX1,
             *locX2;

    int       A[PROCS]::[TEST_MACH_LIMIT];
    elem_t    B[PROCS]::[TEST_MACH_LIMIT];
    chPair_t  C[PROCS]::[TEST_MACH_LIMIT];
    edge_t    D[PROCS]::[TEST_MACH_LIMIT];
    int      X1[PROCS]::[TEST_XPOSE_SIZE];
    int      X2[PROCS]::[TEST_XPOSE_SIZE];

    double secs, best_time[2];

    FILE* machfile;

    locA  = tolocal( A[MYPROC]);
    locB  = tolocal( B[MYPROC]);
    locC  = tolocal( C[MYPROC]);
    locD  = tolocal( D[MYPROC]);
    locX1 = tolocal(X1[MYPROC]);
    locX2 = tolocal(X2[MYPROC]);

    RNG_init();

    barrier();

    on_one {
/******* WARMUP *********************************************/
	j = 0;
	for (i=10; i<TEST_MACH_LIMIT ; i+=TEST_MACH_STEP) {
	    test_fill_int(i, locA);
	    tm[TYPE_INT][j][0] = test_insertsort_i(locA, i);
	    tm[TYPE_INT][j][1] = test_radixsort_i(locA, i);
	    cnt[TYPE_INT][j] = i;
	    j++;
	}
	tmax[TYPE_INT] = j;
/******* WARMUP *********************************************/

	j = 0;
	for (i=10; i<TEST_MACH_LIMIT ; i+=TEST_MACH_STEP) {
	    test_fill_int(i, locA);
	    tm[TYPE_INT][j][0] = test_insertsort_i(locA, i);
	    tm[TYPE_INT][j][1] = test_radixsort_i(locA, i);
	    cnt[TYPE_INT][j] = i;
	    fprintf(outfile," int[%3d]: ins: %9.6f  rs: %9.6f\n",
		    cnt[TYPE_INT][j],
		    tm[TYPE_INT][j][0], tm[TYPE_INT][j][1]);
	    j++;
	}
	tmax[TYPE_INT] = j;

	j = 0;
	for (i=10; i<TEST_MACH_LIMIT ; i+=TEST_MACH_STEP) {
	    test_fill_elem(i, locB);
	    tm[TYPE_ELEM][j][0] = test_insertsort_elem(locB, i);
	    tm[TYPE_ELEM][j][1] = test_radixsort_elem(locB, i);
	    cnt[TYPE_ELEM][j] = i;
	    fprintf(outfile,"elem[%3d]: ins: %9.6f  rs: %9.6f\n",
		    cnt[TYPE_ELEM][j],
		    tm[TYPE_ELEM][j][0], tm[TYPE_ELEM][j][1]);
	    j++;
	}
	tmax[TYPE_ELEM] = j;

	j = 0;
	for (i=10; i<TEST_MACH_LIMIT ; i+=TEST_MACH_STEP) {
	    test_fill_chPair(i, locC);
	    tm[TYPE_CHPR][j][0] = test_insertsort_chPair(locC, i);
	    tm[TYPE_CHPR][j][1] = test_radixsort_chPair(locC, i);
	    cnt[TYPE_CHPR][j] = i;
	    fprintf(outfile,"chpr[%3d]: ins: %9.6f  rs: %9.6f\n",
		    cnt[TYPE_CHPR][j],
		    tm[TYPE_CHPR][j][0], tm[TYPE_CHPR][j][1]);
	    j++;
	}
	tmax[TYPE_CHPR] = j;

	j = 0;
	for (i=10; i<TEST_MACH_LIMIT ; i+=TEST_MACH_STEP) {
	    test_fill_edge(i, locD);
	    tm[TYPE_EDGE][j][0] = test_insertsort_edge(locD, i);
	    tm[TYPE_EDGE][j][1] = test_radixsort_edge(locD, i);
	    cnt[TYPE_EDGE][j] = i;
	    fprintf(outfile,"edge[%3d]: ins: %9.6f  rs: %9.6f\n",
		    cnt[TYPE_EDGE][j],
		    tm[TYPE_EDGE][j][0], tm[TYPE_EDGE][j][1]);
	    j++;
	}
	tmax[TYPE_EDGE] = j;

/*******************************************************************/

	i = 0;
	while ((i<tmax[TYPE_INT]) &&
	       (tm[TYPE_INT][i][0] < tm[TYPE_INT][i][1]))
	    i++;
	rIntBkpt = cnt[TYPE_INT][i];

	i = 0;
	while ((i<tmax[TYPE_ELEM]) &&
	       (tm[TYPE_ELEM][i][0] < tm[TYPE_ELEM][i][1]))
	    i++;
	rElemBkpt = cnt[TYPE_ELEM][i];

	i = 0;
	while ((i<tmax[TYPE_CHPR]) &&
	       (tm[TYPE_CHPR][i][0] < tm[TYPE_CHPR][i][1]))
	    i++;
	rChPairBkpt = cnt[TYPE_CHPR][i];

	i = 0;
	while ((i<tmax[TYPE_EDGE]) &&
	       (tm[TYPE_EDGE][i][0] < tm[TYPE_EDGE][i][1]))
	    i++;
	rEdgeBkpt = cnt[TYPE_EDGE][i];
	
/*******************************************************************/

    }

    barrier();

    s = TEST_SLOT_SIZE;
    
    best_time[0] = LARGE_DOUBLE;
    for (j=0 ; j<TEST_MACH_TRIES ; j++) {
	test_fill_int(TEST_XPOSE_SIZE, locX1);
	barrier();
	secs = get_seconds();

	for (i=0 ; i<PROCS ; i++) {
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locX2[s*get_proc],&(X1[get_proc][s*MYPROC]),
		      s*sizeof(int));
	    barrier();
	}

	
	secs = get_seconds() - secs;
	best_time[0] = min(secs, best_time[0]);
    }

    best_time[1] = LARGE_DOUBLE;
    for (j=0 ; j<TEST_MACH_TRIES ; j++) {
	test_fill_int(TEST_XPOSE_SIZE, locX1);
	barrier();
	secs = get_seconds();

	for (i=0 ; i<PROCS ; i++) {
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locX2[s*get_proc],&(X1[get_proc][s*MYPROC]),
		      s*sizeof(int));
	}

	secs = get_seconds() - secs;
	best_time[1] = min(secs, best_time[1]);
    }


    on_one{
	if (best_time[0] <= best_time[1])
	    mXposeB = 1;
	else
	    mXposeB = 0;
    }
	
    barrier();

    on_one {
	machfile = fopen(MACH_FILENAME,"w+");
	fprintf(machfile,"#ifndef _MACH_SPECIFICS_H\n");
	fprintf(machfile,"#define _MACH_SPECIFICS_H\n\n");
	fprintf(machfile,"#define RADIXSORT_INT_BREAKPT\t\t%4d\n",rIntBkpt);
	fprintf(machfile,"#define RADIXSORT_ELEM_BREAKPT\t\t%4d\n",rElemBkpt);
	fprintf(machfile,"#define RADIXSORT_CHPAIR_BREAKPT\t%4d\n",rChPairBkpt);
	fprintf(machfile,"#define RADIXSORT_EDGE_BREAKPT\t\t%4d\n",rEdgeBkpt);
	fprintf(machfile,"#define MATRIX_XPOSE_BARRIER\t\t%4d\n",mXposeB);
	fprintf(machfile,"\n#endif\n\n");
	fflush(machfile);
	fclose(machfile);
    }

    barrier();

}


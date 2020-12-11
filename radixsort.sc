/*									tab:8
 *
 * radixsort.sc - Integer sort by radix
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
 * Creation Date:	March 11, 1995
 * Filename:		radixsort.sc
 * History:
 */

/* Radix Sort in Split-C */

#include "radixsort.h"
#include "droute.h"

/****************************************************/
static inline
void all_countsort_safe(int q,
			int Key[PROCS]::[q],
			int Sorted[PROCS]::[q],
			int R,
			int bitOff, int m,
			int route)
/****************************************************/
/* R (range)      must be a multiple of PROCS */
/* q (elems/proc) must be a multiple of PROCS */
{

    register int
	j,
	k,
	s,
	get_proc,
	offset;
    
    int *lKey,
	*lAddr,
	*lIndex,
	*lScanTran;

    intpair_t
	*lIntScanTran,
	*lScans;

    int Addr[PROCS]::[q],
	Index[PROCS]::[R],
	ScanTran[PROCS]::[R],
	addrBits[q];

    intpair_t
	IntScanTran[PROCS]::[R],
	Scans[PROCS]::[R];
    
    lKey   = tolocal(Key[MYPROC]);
    lAddr  = tolocal(Addr[MYPROC]);
    lIndex = tolocal(Index[MYPROC]);
    lScanTran = tolocal(ScanTran[MYPROC]);

    lIntScanTran = tolocal(IntScanTran[MYPROC]);
    lScans = tolocal(Scans[MYPROC]);

    barrier();

    if (R < PROCS)
	R = PROCS;

    s = R/PROCS;

    for (k=0 ; k<R ; k++)
	lIndex[k] = 0;

    for (k=0 ; k<q ; k++) {
	addrBits[k] = bits(lKey[k],bitOff,m);
	lIndex[addrBits[k]]++;
    }

/* Perform Matrix Transpose of Index -> ScanTran */    

    barrier();			/* Wait for everyone */

    for (k=0 ; k<PROCS ; k++) {
        get_proc = (MYPROC + k) % PROCS; /* Stagger prefetches using MYPROC */
        bulk_read(&lScanTran[s*get_proc],
		  &(Index[get_proc][s*MYPROC]),
		  s*sizeof(int));
	barrier();
    }

/* Calculate local Prefix Sums */    

    for (j=0 ; j<s ; j++)
	for (k=1 ; k<PROCS ; k++)
	    lScanTran[k*s + j] += lScanTran[(k-1)*s + j];

/* Create IntLeaveScan by interleaving ScanTran[.][*] with ScanTran[P-1][*]*/

    for (j=0 ; j<s ; j++)
	for (k=0 ; k<PROCS ; k++) {
	    lIntScanTran[k*s + j].a = lScanTran[k*s + j];
	    lIntScanTran[k*s + j].b = lScanTran[(PROCS-1)*s + j];
	}

    barrier();
    
/* Perform Inverse Matrix Transpose of IntScanTran -> Scans */
    
    for (k=0 ; k<PROCS ; k++) {
        get_proc = (MYPROC + k) % PROCS; /* Stagger prefetches using MYPROC */
        bulk_read(&lScans[s*get_proc],
		  &(IntScanTran[get_proc][s*MYPROC]),
		  s*sizeof(intpair_t));
	barrier();
    }

    offset = 0;

    for (k=0 ; k<R ; k++) {
	lIndex[k] = (lScans[k].a - lIndex[k]) + offset;
	offset  += lScans[k].b;
    }

    for (k=0 ; k<q ; k++) {
	lAddr[k] = lIndex[addrBits[k]];
	lIndex[addrBits[k]]++;
    }

    switch (route) {
    case 1 : all_route_h_rel_det(q, 1, Key, Addr, Sorted); break;
    case 2 : all_route_h_rel_det2(q, 1, Key, Addr, Sorted); break;
    case 3 : all_route_h_rel_det3(q, 1, Key, Addr, Sorted); break;
    case 4 : all_route_h_rel_det4(q, 1, Key, Addr, Sorted); break;
    case 5 : all_route_linear_perm(q, 1, Key, Addr, Sorted); break;
    case 6 : all_route_single_phase(q, 1, Key, Addr, Sorted); break;
    }
}


/****************************************************/
void all_countsort(int q,
		   int Key[PROCS]::[q],
		   int Sorted[PROCS]::[q],
		   int R,
		   int bitOff, int m,
		   int route)
/****************************************************/
/* R (range)      must be a multiple of PROCS */
/* q (elems/proc) must be a multiple of PROCS */
{

    register int
	j,
	k,
	s,
	get_proc,
	offset;
    
    int *lKey,
	*lAddr,
	*lIndex,
	*lScanTran,
	*t1ptr,
	*t2ptr,
	*t3ptr;

    intpair_t
	*lIntScanTran,
	*lScans,
	*p1ptr,
	*p2ptr;

    int *spread Addr;
    int *spread Index;
    int *spread ScanTran;
    int *addrBits;

    intpair_t *spread IntScanTran;
    intpair_t *spread Scans;

    Addr = all_spread_malloc(PROCS,q*sizeof(int));
    assert_spread_malloc(Addr);
    Index = all_spread_malloc(PROCS,R*sizeof(int));
    assert_spread_malloc(Index);
    ScanTran = all_spread_malloc(PROCS,R*sizeof(int));
    assert_spread_malloc(ScanTran);
    addrBits = (int*)malloc(q*sizeof(int));
    assert_malloc(addrBits);

    IntScanTran = all_spread_malloc(PROCS,R*sizeof(intpair_t));
    assert_spread_malloc(IntScanTran);
    Scans = all_spread_malloc(PROCS,R*sizeof(intpair_t));
    assert_spread_malloc(Scans);
    
    lKey   = (int *)(Key + MYPROC);
    lAddr  = (int *)(Addr + MYPROC);
    lIndex = (int *)(Index + MYPROC);
    lScanTran = (int *)(ScanTran + MYPROC);

    lIntScanTran = (intpair_t *)(IntScanTran + MYPROC);
    lScans = (intpair_t *)(Scans + MYPROC);

    barrier();

    if (R < PROCS)
	R = PROCS;

    s = DIVPROCS(R);

    t1ptr = lIndex-1;
    t2ptr = lIndex+R;
    while (++t1ptr<t2ptr)
	*t1ptr = 0;

    for (k=0 ; k<q ; k++) {
	*(t1ptr = addrBits+k) = bits(*(lKey+k),bitOff,m);
	(*(lIndex + *t1ptr))++;
    }

/* Perform Matrix Transpose of Index -> ScanTran */    

    barrier();			/* Wait for everyone */

#if (defined(CM5))
    for (k=0 ; k<PROCS ; k++) {
        get_proc = MODPROCS(MYPROC + k); /* Stagger prefetches using MYPROC */
        bulk_read(lScanTran + (s*get_proc),
		 (int *global)(Index + get_proc) + (s*MYPROC),
		  s*sizeof(int));
	barrier();
    }
#else
    for (k=0 ; k<PROCS ; k++) {
	get_proc = MODPROCS(MYPROC + k); /* Stagger prefetches using MYPROC */
        bulk_get(lScanTran + (s*get_proc),
		 (int *global)(Index + get_proc) + (s*MYPROC),
		 s*sizeof(int));
    }
    sync();
    barrier();
#endif
    
/* Calculate local Prefix Sums */    

    for (j=0 ; j<s ; j++) {
	t2ptr = t1ptr = lScanTran+j;
	t2ptr += s;
	t3ptr = t1ptr + PROCS*s;
	while (t2ptr < t3ptr) {
	    *t2ptr += *t1ptr;
	    t1ptr = t2ptr;
	    t2ptr += s;
	}
    }

/* Create IntLeaveScan by interleaving ScanTran[.][*] with ScanTran[P-1][*]*/

    t2ptr = lScanTran + (PROCS-1)*s;
    for (j=0 ; j<s ; j++) {
	t1ptr = lScanTran + j;
	p1ptr = lIntScanTran + j;
	p2ptr = p1ptr + s*PROCS;
	while (p1ptr < p2ptr) {
	    p1ptr->a = *t1ptr;
	    p1ptr->b = *t2ptr;
	    t1ptr += s;
	    p1ptr += s;
	}
	t2ptr++;
    }

    barrier();
    
/* Perform Inverse Matrix Transpose of IntScanTran -> Scans */
    
#if (defined(CM5))
    for (k=0 ; k<PROCS ; k++) {
        get_proc = MODPROCS(MYPROC + k); /* Stagger prefetches using MYPROC */
        bulk_read(lScans + (s*get_proc),
		 (intpair_t *global)(IntScanTran + get_proc) + (s*MYPROC),
		  s*sizeof(intpair_t));
	barrier();
    }
#else
    for (k=0 ; k<PROCS ; k++) {
        get_proc = MODPROCS(MYPROC + k); /* Stagger prefetches using MYPROC */
        bulk_get(lScans + (s*get_proc),
		 (intpair_t *global)(IntScanTran + get_proc) + (s*MYPROC),
		 s*sizeof(intpair_t));
    }
    sync();
    barrier();
#endif
    
    offset = 0;

    for (k=0 ; k<R ; k++) {
	*(lIndex + k) = ((p1ptr=lScans+k)->a - *(lIndex+k)) + offset;
	offset  += p1ptr->b;
    }

    for (k=0 ; k<q ; k++) {
	*(lAddr + k) = *(t1ptr = (lIndex + *(addrBits+k)));
	(*t1ptr)++;
    }

    all_spread_free(Scans);
    all_spread_free(IntScanTran);

    free(addrBits);
    all_spread_free(ScanTran);
    all_spread_free(Index);

    switch (route) {
    case 1 : all_route_h_rel_det(q, 1, Key, Addr, Sorted); break;
    case 2 : all_route_h_rel_det2(q, 1, Key, Addr, Sorted); break;
    case 3 : all_route_h_rel_det3(q, 1, Key, Addr, Sorted); break;
    case 4 : all_route_h_rel_det4(q, 1, Key, Addr, Sorted); break;
    case 5 : all_route_linear_perm(q, 1, Key, Addr, Sorted); break;
    case 6 : all_route_single_phase(q, 1, Key, Addr, Sorted); break;
    }

    all_spread_free(Addr);
}


/****************************************************/
inline void all_radixsort20(int q,
			    int Keys[PROCS]::[q],
			    int Sorted[PROCS]::[q],
			    int route)
/****************************************************/
{
    all_countsort(q, Keys,   Sorted, (1<<10),  0, 10, route);
    all_countsort(q, Sorted, Sorted, (1<<10), 10, 10, route);
}

/****************************************************/
inline void all_radixsort_passInfo(int q,
				   int Keys[PROCS]::[q],
				   int Sorted[PROCS]::[q],
				   int passes,
				   int *b)
/****************************************************/
{
    register int
	i,
	j;

    j=0;
    all_countsort(q, Keys, Sorted, (1<<b[0]), 0, b[0], 1);

    for (i=1 ; i<passes ; i++) 
	all_countsort(q, Sorted, Sorted, max((1<<b[i]),PROCS),
		      (j+=b[i-1]), b[i], 1);

}










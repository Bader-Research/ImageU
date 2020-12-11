/*									tab:8
 *
 * radixsort.h - Integer sort by radix
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
 * Filename:		radixsort.h
 * History:
 */

#ifndef _RADIXSORT_H
#define _RADIXSORT_H

#include "ImageU.h"
#include "sorting.h"
#include "intpair.h"

void all_countsort(int q,
		   int Key[PROCS]::[q],
		   int Sorted[PROCS]::[q],
		   int R,
		   int bitOff, int m,
		   int route);

/****************************************************/
_INLINE void all_radixsort(int q,
			   int Keys[PROCS]::[q],
			   int Sorted[PROCS]::[q],
			   int route)
/****************************************************/
{
/*  all_radixsort_bits(q, Keys, Sorted, sizeof(int)<<3, route); */
    
    all_countsort(q, Keys,   Sorted, (1<<11),  0, 11, route);
    all_countsort(q, Sorted, Sorted, (1<<11), 11, 11, route);
    all_countsort(q, Sorted, Sorted, (1<<10), 22, 10, route);
}

/****************************************************/
_INLINE void all_radixsort_bits(int q,
				int Keys[PROCS]::[q],
				int Sorted[PROCS]::[q],
				int numBits,
				int route)
/****************************************************/
{
    register int
	i,
	w, R, m,
	rng,
	passes;

    w = numBits;            /* The number of bits in the key */
    m = w>>2;               /* radix of 10 bits */
    R = 1<<m;               /* R=1024 */
    passes = w / m;         /* Total number of passes needed */

    all_countsort(q, Keys, Sorted, R, 0, m, route);  /* Unroll loop */

    for (i=1 ; i<passes ; i++) 
	all_countsort(q, Sorted, Sorted, R, i*m, m, route);

/* Range must be a multiple of processors */
    i = m*passes;
    if (i != w) /* Finish up sorting */  {
	rng = max( 1<<(w-i) , PROCS);
	all_countsort(q, Sorted, Sorted, rng, i, w-i, route);
    }    
}

/****************************************************/
_INLINE void all_radixsort_bitInfo(int q,
				   int Keys[PROCS]::[q],
				   int Sorted[PROCS]::[q],
				   int w, int m)
/****************************************************/
{
    register int
	i,
	R,
	rng,
	passes;

    R = 1 << m;
    passes = w / m;         /* Total number of passes needed */

    all_countsort(q, Keys, Sorted, R, 0, m, 1);  /* Unroll loop */

    for (i=1 ; i<passes ; i++) 
	all_countsort(q, Sorted, Sorted, R, i*m, m, 1);

/* Range must be a multiple of processors */
    i = m*passes;
    if (i != w) /* Finish up sorting */  {
	rng = max( 1<<(w-i) , PROCS);
	all_countsort(q, Sorted, Sorted, rng, i, w-i, 1);
    }    
}


#endif

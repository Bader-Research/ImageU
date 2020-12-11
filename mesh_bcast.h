/*									tab:8
 *
 * mesh_bcast.h - support for connected components
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
 * Authors: 		David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:		1.0
 * Creation Date:	April 13, 1995
 * Filename:		mesh_bcast.h
 * History:
 */

#ifndef _MESH_BCAST_H
#define _MESH_BCAST_H

#include "ImageU.h"

_INLINE int greyCode(int x) {
    return ((x>>1) ^ x );
}

_INLINE int greyInvCode(int x)
{
    register int
	z = x;
    
    while (x) {
	x = x>>1;
	z ^= x;
    }
	
    return(z);
}

_INLINE int greyEncoder(int vIt, int wIt, int i, int j)
{
    return (greyCode(i)<<wIt) | greyCode(j);
}

_INLINE void greyDecoder(int vIt, int wIt, int* i, int* j, int x)
{
    *i = greyInvCode( x>>wIt );
    *j = greyInvCode( x & (~(~0 << wIt)) );
}

void all_mesh_bcast(int v, int w, int l,
		    int vIt, int wIt,		       
		    int ii, int jj, int grpManV, int grpManW,
		    chPair_t chList[v][w]::[l],
		    int changeSize[v][w]::);

#endif


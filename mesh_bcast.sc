/*									tab:8
 *
 * mesh_bcast.sc - support for connected components
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
 * Filename:		mesh_bcast.sc
 * History:
 */

#include "mesh_bcast.h"

/*************************************************************/
void all_mesh_bcast(int v, int w, int l,
		    int vIt, int wIt,		       
		    int ii, int jj, int grpManV, int grpManW,
		    chPair_t chList[v][w]::[l],
		    int changeSize[v][w]::)
/*************************************************************/
{
    int
	orgV, orgW,       /* original of my group in real coords */
	g_grey,           /* my grey code address inside my group */
	g_GM,             /* the Group Manager's grey code offset */
	g_label,          /* my final "hypercube" label inside my group */
	i,
	send_i, send_j,
	d;

    orgV    = (ii>>vIt) << vIt;
    orgW    = (jj>>wIt) << wIt;
    g_grey  = greyEncoder(vIt, wIt, ii-orgV, jj-orgW);
    g_GM    = greyEncoder(vIt, wIt, grpManV-orgV, grpManW-orgW);
    g_label = g_grey ^ g_GM;
    d       = vIt + wIt;

    for (i=0 ; i<d ; i++) {
	if (g_label < (1<<i)) {
	    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
	    changeSize[send_i+orgV][send_j+orgW] = changeSize[ii][jj];
	}
	barrier();
    }

    for (i=0 ; i<d ; i++) {
	if ((g_label < (1<<i)) && (changeSize[ii][jj] > 0)) {
	    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
	    bulk_write(&(chList[send_i+orgV][send_j+orgW][0]),
		       tolocal(chList[ii][jj]),
		       changeSize[ii][jj]*sizeof(chPair_t));
	}
	barrier();
    }
}



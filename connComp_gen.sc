/*									tab:8
 *
 * connComp_gen.sc - auxiliary routines for connected components
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
 * Version:		2
 * Creation Date:	October 20, 1994
 * Filename:		connComp_gen.sc
 * History:
 *	DAB	2	Wed Dec 14 12:54:08 EST 1994
 *		Modified BFS_auxGraph - move thru auxIndex in reverse
 */

#include "connComp_gen.h"

int print_tile(int* tile, int* labels, int q, int r)
/* Print my tile which is q x r */
{
    register int
	i, j;

    fprintf(outfile,"PE %2d: Printing tile\n",MYPROC);
    fprintf(outfile,"      IMAGE                LABELS\n");

    for (i=0 ; i<q ; i++) {
	for (j=0 ; j<r ; j++)
	    fprintf(outfile,"%3d ",tile[i*r + j]);
	fprintf(outfile,"     ");
	for (j=0 ; j<r ; j++)
	    fprintf(outfile,"%3d ",labels[i*r + j]);
	fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n\n");
}




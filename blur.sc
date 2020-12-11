/*									tab:8
 *
 * blur.sc - blur an image
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
 * Creation Date:	February 15, 1995
 * Filename:		blur.sc
 * History:
 */

#include "blur.h"
#include "math_func.h"

/*************************************************************/
void all_blur_image(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int k,
		    int Y[v][w]::[q][r],
		    int sigma)
/*************************************************************/
{
    double blur_sigma;

    int Z[v][w]::[q][r];          /* pType */

    char buf[MAXLEN];
    char buf2[MAXLEN];

    RNG_init();
    
    all_get_pType(v,w,q,r,Z); 

    blur_sigma = (double)sigma / 0.375;
    
    all_filter_GaussNoise(v,w,q,r,X,k,Y,0,blur_sigma); 

    all_filter_Gauss3X3(v,w,q,r,Y,k,Y,Z); 

    on_one {
	if (OUTPUT_IMAGE) {
	    sprintf(buf,"Noise and blur of grey image, sigma = %d",sigma);
	    sprintf(buf2,"imageBlur.%d",sigma);
	    fprintf(outfile,"Writing Output Image... ");
	    fflush(outfile);	    
	    image_print_im(buf, buf2,
			   v,w,q,r,k,Y);
	    fprintf(outfile,"Done.\n");
	    fflush(outfile);
	}
    }

}

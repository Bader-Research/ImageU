/*									tab:8
 *
 * image_print.sc - output images
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
 * Author: 		David A. Bader <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:		1
 * Creation Date:	November 2, 1995
 * Filename:		image_print.sc
 * History:
 */

#include "image_print.h"
#include "ImageU.h"

#define is_bit_set(x, b)  ((x) & (1<<(b)))

char* image_add_tag(char* outfilename, int image_t)
{
    char* buf;

    buf = (char*)malloc(MAXLEN*sizeof(char));
    assert_malloc(buf);
    
    switch (image_t) {
    case VISTA_IMAGE : sprintf(buf,"%s.v",outfilename);
	break;
    case PS_IMAGE : sprintf(buf,"%s.eps",outfilename);
	break;
    case GIF_IMAGE : sprintf(buf,"%s.gif",outfilename);
	break;
    case PGM_IMAGE : sprintf(buf,"%s.pgm",outfilename);
	break;
#if X11
    case X11_IMAGE : sprintf(buf,"%s",outfilename);
	break;
#endif	
    default : fprintf(stderr,"ERROR: Bad image type (%d)\n",image_t);
	exit(1);
    }

    return(buf);
}

void image_print_im(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r])
{
    if (is_bit_set(OUTPUT_IMAGE,VISTA_IMAGE))
	print_vista_im(desc,image_add_tag(outfilename,VISTA_IMAGE),
		       v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PS_IMAGE))
	print_ps_im(desc,image_add_tag(outfilename,PS_IMAGE),
		    v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,GIF_IMAGE))
	print_gif_im(desc,image_add_tag(outfilename,GIF_IMAGE),
		     v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PGM_IMAGE))
	print_pgm_im(desc,image_add_tag(outfilename,PGM_IMAGE),
		     v,w,q,r,k,image);
#if X11
    if (is_bit_set(OUTPUT_IMAGE,X11_IMAGE))
	print_x11_im(desc,image_add_tag(outfilename,X11_IMAGE),
		     v,w,q,r,k,image);
#endif    
}

void image_print_cc_orig(char *desc, char *outfilename,
			 int v, int w, int q, int r, int k,
			 int image[v][w]::[q][r])
{
    if (is_bit_set(OUTPUT_IMAGE,VISTA_IMAGE))
	print_vista_cc_orig(desc,image_add_tag(outfilename,VISTA_IMAGE),
			    v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PS_IMAGE))
	print_ps_cc_orig(desc,image_add_tag(outfilename,PS_IMAGE),
			 v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,GIF_IMAGE))
	print_gif_cc_orig(desc,image_add_tag(outfilename,GIF_IMAGE),
			 v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PGM_IMAGE))
	print_pgm_cc_orig(desc,image_add_tag(outfilename,PGM_IMAGE),
			 v,w,q,r,k,image);
#if X11
    if (is_bit_set(OUTPUT_IMAGE,X11_IMAGE))
	print_x11_cc_orig(desc,image_add_tag(outfilename,X11_IMAGE),
			 v,w,q,r,k,image);
#endif    
}

void image_print_cc(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r])
{
    if (is_bit_set(OUTPUT_IMAGE,VISTA_IMAGE))
	print_vista_cc(desc,image_add_tag(outfilename,VISTA_IMAGE),
		       v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PS_IMAGE))
	print_ps_cc(desc,image_add_tag(outfilename,PS_IMAGE),
		    v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,GIF_IMAGE))
	print_gif_cc(desc,image_add_tag(outfilename,GIF_IMAGE),
		     v,w,q,r,k,image);
    if (is_bit_set(OUTPUT_IMAGE,PGM_IMAGE))
	print_pgm_cc(desc,image_add_tag(outfilename,PGM_IMAGE),
		     v,w,q,r,k,image);
#if X11
    if (is_bit_set(OUTPUT_IMAGE,X11_IMAGE))
	print_x11_cc(desc,image_add_tag(outfilename,X11_IMAGE),
		     v,w,q,r,k,image);
#endif    
}

int image_print_real_color(char *desc, char *outfilename,
			   int v, int w, int q, int r, int k,
			   int maxLabel,
			   int Labels[v][w]::[q][r], 
			   int X[v][w]::[q][r])
{
    int i=0;
    
    if (is_bit_set(OUTPUT_IMAGE,VISTA_IMAGE))
	i = print_vista_real_color(desc,image_add_tag(outfilename,VISTA_IMAGE),
				   v,w,q,r,k,
				   maxLabel, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,PS_IMAGE))
	i = print_ps_real_color(desc,image_add_tag(outfilename,PS_IMAGE),
				v,w,q,r,k,
				maxLabel, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,GIF_IMAGE))
	i = print_gif_real_color(desc,image_add_tag(outfilename,GIF_IMAGE),
				 v,w,q,r,k,
				 maxLabel, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,PGM_IMAGE))
	i = print_pgm_real_color(desc,image_add_tag(outfilename,PGM_IMAGE),
				 v,w,q,r,k,
				 maxLabel, Labels, X);
#if X11
    if (is_bit_set(OUTPUT_IMAGE,X11_IMAGE))
	i = print_x11_real_color(desc,image_add_tag(outfilename,X11_IMAGE),
				 v,w,q,r,k,
				 maxLabel, Labels, X);
#endif
    
    return (i);
}

int image_print_real_color_on_one(char *desc, char *outfilename,
				  int v, int w, int q, int r, int k,
				  int Labels[v][w]::[q][r], 
				  int X[v][w]::[q][r])
{
    int i=0;
    
    if (is_bit_set(OUTPUT_IMAGE,VISTA_IMAGE))
	i = print_vista_real_color_on_one(desc,
					  image_add_tag(outfilename,VISTA_IMAGE),
					  v,w,q,r,k, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,PS_IMAGE))
	i = print_ps_real_color_on_one(desc,
				       image_add_tag(outfilename,PS_IMAGE),
				       v,w,q,r,k, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,GIF_IMAGE))
	i = print_gif_real_color_on_one(desc,
				       image_add_tag(outfilename,GIF_IMAGE),
				       v,w,q,r,k, Labels, X);
    if (is_bit_set(OUTPUT_IMAGE,PGM_IMAGE))
	i = print_pgm_real_color_on_one(desc,
				       image_add_tag(outfilename,PGM_IMAGE),
				       v,w,q,r,k, Labels, X);
#if X11
    if (is_bit_set(OUTPUT_IMAGE,X11_IMAGE))
	i = print_x11_real_color_on_one(desc,
				       image_add_tag(outfilename,X11_IMAGE),
				       v,w,q,r,k, Labels, X);
#endif
    
    return (i);
}


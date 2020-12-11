/*									tab:8
 *
 * image_print.h - output images 
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
 * Creation Date:	November 2, 1995
 * Filename:		image_print.h
 * History:
 */

#ifndef _IMAGE_PRINT_H
#define _IMAGE_PRINT_H

#include "print_vista.h"
#include "print_ps.h"
#include "print_gif.h"
#include "print_pgm.h"
#if X11
#include "print_X11.h"
#endif

#define  VISTA_IMAGE 1
#define  PS_IMAGE    2
#define  GIF_IMAGE   3
#define  PGM_IMAGE   4
#if X11
#define  X11_IMAGE   5
#endif

void image_print_im(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r]);

void image_print_cc(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r]);

int image_print_real_color(char *desc, char *outfilename,
			   int v, int w, int q, int r, int k,
			   int maxLabel,
			   int Labels[v][w]::[q][r],
			   int X[v][w]::[q][r]);

int image_print_real_color_on_one(char *desc, char *outfilename,
				  int v, int w, int q, int r, int k,
				  int Labels[v][w]::[q][r],
				  int X[v][w]::[q][r]);

/* image_print_real_color returns the number of labels */


#endif




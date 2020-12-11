/*									tab:8
 *
 * print_pgm.h - output images in Portable Graymap (PGM) format
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
 * Creation Date:	October 20, 1994
 * Filename:		print_pgm.h
 * History:
 */

#ifndef _PRINT_PGM_H
#define _PRINT_PGM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <split-c/split-c.h>

void print_pgm_im(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r]);

void print_pgm_cc(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r]);

int print_pgm_real_color(char *desc, char *outfilename,
			   int v, int w, int q, int r, int k,
			   int maxLabel,
			   int Labels[v][w]::[q][r],
			   int X[v][w]::[q][r]);

int print_pgm_real_color_on_one(char *desc, char *outfilename,
				  int v, int w, int q, int r, int k,
				  int Labels[v][w]::[q][r],
				  int X[v][w]::[q][r]);

/* print_pgm_real_color returns the number of labels */

#endif




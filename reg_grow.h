/*									tab:8
 *
 * reg_grow.h - region growing algorithm
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
 * Creation Date:	December 29, 1994
 * Filename:		reg_grow.h
 * History:
 */

#ifndef _REGGROW_H
#define _REGGROW_H

#include "ImageU.h"
#include "queue.h"
#include "connComp_gen.h"

void all_filter_SNF(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int k,
		    int Y[v][w]::[q][r],
		    int pType[v][w]::[q][r],
		    int steps, int epsilon);

void all_filter_SNF_reclass(int v, int w, int q, int r, 
			    int X[v][w]::[q][r], int k,
			    int Sel[v][w]::[q][r],
			    int Y[v][w]::[q][r],
			    int pType[v][w]::[q][r],
			    int epsilon);

void all_crop_border(int v, int w, int q, int r, 
		     int X[v][w]::[q][r], int Y[v][w]::[q][r],
		     int crop);

void all_get_pType(int v, int w, int q, int r, 
		   int pType[v][w]::[q][r]);

double all_find_sigma_use_sort(int v, int w, int q, int r, 
		   int X[v][w]::[q][r], int k,
		   int pType[v][w]::[q][r]);

double all_find_sigma_use_median(int v, int w, int q, int r, 
			     int X[v][w]::[q][r], int k,
			     int pType[v][w]::[q][r]);

void all_filter_1NN(int v, int w, int q, int r, 
		    int X[v][w]::[q][r], int k,
		    int Y[v][w]::[q][r],
		    int pType[v][w]::[q][r]);

void all_filter_1NN_oneshot(int v, int w, int q, int r, 
			    int X[v][w]::[q][r], int k,
			    int Y[v][w]::[q][r],
			    int pType[v][w]::[q][r]);

void all_scale_and_edge_detect(int v, int w, int q, int r,
			       int mag, 
			       int X[v][w]::[q][r],
			       int Y[v][w]::[q*mag][r*mag],
			       int pType[v][w]::[q][r],
			       int edgeColor,
			       int Ycol[v][w]::[q][r],
			       int Ymag[v][w]::[q*mag][r*mag]);

void all_filter_Gauss3X3(int v, int w, int q, int r, 
			 int X[v][w]::[q][r], int k,
			 int Y[v][w]::[q][r],
			 int pType[v][w]::[q][r]);

void all_filter_GaussNoise(int v, int w, int q, int r, 
			   int X[v][w]::[q][r], int k,
			   int Y[v][w]::[q][r],
			   double mean, double sigma);

int all_decide_artificial(int v, int w, int q, int r, 
			  int X[v][w]::[q][r],
			  int pType[v][w]::[q][r]);

#endif






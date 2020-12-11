/*                                                                    tab:8
 *
 * inputs_median.h - Parallel Median Algorithm Data Sets
 *                   Parallel Sorting Algorithm Data Sets
 *
 * 
 * "Copyright (c) 1996 The Regents of the University of Maryland.
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
 * Authors:             David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:             1.0
 * Creation Date:       February 6, 1996
 * Filename:            inputs_median.h
 * History:
 */

#ifndef _INPUTS_MEDIAN_H
#define _INPUTS_MEDIAN_H

#include "ImageU.h"

int all_fill_array(int M, int A[PROCS]::[M], int N[PROCS]::,
		   char *filename);

void all_fill_nas(int M, int A[PROCS]::[M]);
void all_fill_cyclic(int M, int A[PROCS]::[M]);
void all_fill_pid(int M, int A[PROCS]::[M]);
void all_fill_random(int M, int A[PROCS]::[M], int distrib);
void all_fill_hrel(int M, int A[PROCS]::[M], int h);
void all_fill_g_group_i(int M, int A[PROCS]::[M], int h, int g, int t);

void fill_same(int , int *);
void fill_linear(int , int *);
void fill_random(int , int *);

void fill_random_d(int , double *);
void fill_gaussian_d(int , double *);
void fill_g_group_d(int , double *, int);
void fill_bucket_d(int , double *);
void fill_staggered_d(int , double *);
void fill_zero_d(int , double *);

#endif



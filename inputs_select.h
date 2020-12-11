/*                                                                    tab:8
 *
 * inputs_select.h - Parallel Selection Algorithm Data Sets
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
 * Filename:            inputs_select.h
 * History:
 */

#ifndef _INPUTS_SELECT_H
#define _INPUTS_SELECT_H

#include "ImageU.h"

int all_fill_array_control(int M, int A[PROCS]::[M], 
			   int N[PROCS]::, int total_n);

int all_fill_array_linear(int M, int A[PROCS]::[M], 
			  int N[PROCS]::, int total_n);

int all_fill_array_exp(int M, int A[PROCS]::[M], 
		       int N[PROCS]::, int total_n);

int all_fill_array_one(int M, int A[PROCS]::[M], 
		       int N[PROCS]::, int total_n);

int all_fill_array_gauss(int M, int A[PROCS]::[M], 
			 int N[PROCS]::, int total_n);

#endif

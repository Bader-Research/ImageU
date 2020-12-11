/*                                                                    tab:8
 *
 * droute.h - Deterministic Routing of H-relations
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
 * Authors:             David A. Bader   <dbader@umiacs.umd.edu>
 *                      David Helman     <helman@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:             1.0
 * Creation Date:       August 1, 1995
 * Filename:            droute.h
 * History:
 */

#ifndef _DROUTE_H
#define _DROUTE_H

#include "ImageU.h"
#include "intpair.h"

int all_route_h_rel_det(int q, int h,
			int Keys[PROCS]::[q],
			int Addr[PROCS]::[q],
			int Routed[PROCS]::[h*q]);

int all_route_h_rel_det2(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed);

int all_route_h_rel_det3(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed);

#if (defined(SP2))
int all_route_h_rel_det4(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed);
#endif

int all_route_linear_perm(int q, int h,
			  int *spread Keys,
			  int *spread Addr,
			  int *spread Routed);

int all_route_single_phase(int q, int h,
			   int *spread Keys,
			   int *spread Addr,
			   int *spread Routed);
    
#endif

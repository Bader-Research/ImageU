/*									tab:8
 *
 * connComp_gen.h - auxiliary routines for connected components
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
 * Filename:		connComp_gen.h
 * History:
 */

#ifndef _CONNCOMPGEN_H
#define _CONNCOMPGEN_H

#include "ImageU.h"
#include "mesh_bcast.h"
#include "sorting.h"

typedef struct tileHook *TILEHOOKPTR;

struct tileHook
{
    int label;
    int size;
    int* addr;
};

_INLINE int fill_tileHook(TILEHOOKPTR *tileHookArr,
			 CHPAIRPTR tempHookArr, int pairs) {
/* Sort the tempHookArr (label, addr) by label. */
/* Copy the data into structure tileHook Arr */    
/* return the number of hooks */

    register int
	i,
	j,
	hooks;

    TILEHOOKPTR hookPtr;

    if (pairs == 0) return(0);
    
    fastsort_chPair(tempHookArr, pairs);

    hooks = 1;

    for (i=1 ; i<pairs ; i++)
	if (tempHookArr[i].alpha > tempHookArr[i-1].alpha)
	    hooks++;

    if ((*tileHookArr = (TILEHOOKPTR)malloc(hooks*sizeof(struct tileHook)))==NULL)
	fprintf(stderr,"ERROR: tileHookArr could not be malloc'ed\n");

    hookPtr = *tileHookArr;  /* Alias for tileHookArr */
    
    /* Count sizes of each hook */

    for (i=1 ; i<hooks ; i++)
	hookPtr[i].size = 0;

    hookPtr[0].size = 1;
    j = 0;
    
    for (i=1 ; i<pairs ; i++) {
	if (tempHookArr[i].alpha > tempHookArr[i-1].alpha)
	    j++;
	hookPtr[j].size++;
    }

    /* Allocate memory for addresses on each hook */
    
    for (i=0 ; i<hooks ; i++) 
	if ((hookPtr[i].addr = (int*)malloc((hookPtr[i].size)*
						sizeof(int))) == NULL)
	    fprintf(stderr,"ERROR: tileHook(sub)Arr could not be malloc'ed\n");

    /* Fill Out Hook Structure's addresses */

    for (i=1 ; i<hooks ; i++)
	hookPtr[i].size = 0;

    j = 0;
    hookPtr[j].size    = 1;
    hookPtr[j].label   = tempHookArr[0].alpha;
    hookPtr[j].addr[0] = tempHookArr[0].beta;
    
    for (i=1 ; i<pairs ; i++) {
	if (tempHookArr[i].alpha > tempHookArr[i-1].alpha) {
	    j++;
	    hookPtr[j].label   = tempHookArr[i].alpha;
	    hookPtr[j].addr[0] = tempHookArr[i].beta;
	}
	else 
	    hookPtr[j].addr[hookPtr[j].size] = tempHookArr[i].beta;

	hookPtr[j].size++;
    }

    return(hooks);
}

_INLINE void destroy_tileHook(TILEHOOKPTR tileHookArr, int size) {
    register int i;

    for (i=0 ; i<size ; i++)
	free( tileHookArr[i].addr );
}

_INLINE void create_hooks(CHPAIRPTR tileHookOrig, TILEHOOKPTR tileHookArr, int hooks) {
/* Create a record of the original labels and a hook address into the tile */
    register int i;

    for (i=0 ; i<hooks ; i++) {
	tileHookOrig[i].alpha = tileHookArr[i].label;
	tileHookOrig[i].beta  = (tileHookArr[i].addr)[0];
    }
}

int print_tile(int*, int*, int, int);

_INLINE void addAuxEdge(edge_t *auxG, int *auxC, int aVal, int bVal)
/* Add the undirected edges <a,b> and <b,a>, into the auxiliary graph */
{

/* If a == b, do nothing */

    if (aVal==bVal) return;

    auxG[*auxC].a = aVal;
    auxG[*auxC].b = bVal;
    (*auxC)++;
    auxG[*auxC].a = bVal;
    auxG[*auxC].b = aVal;
    (*auxC)++;
}

typedef struct equeue_tag {
    int count;
    int front;
    int rear;
    int MAXQ;
    int *entry;
} eQueue_type;

_INLINE void eQInitialize(eQueue_type *equeue_ptr, int n) {
    equeue_ptr->count = 0;
    equeue_ptr->front = 0;
    equeue_ptr->rear = -1;
    equeue_ptr->MAXQ = n;
    if ((equeue_ptr->entry = (int*)(malloc(equeue_ptr->MAXQ*sizeof(int))))
	==NULL)
	fprintf(stderr,"ERROR: eQueue entries could not be malloc'ed\n");
}

_INLINE void eQDestruct(eQueue_type *equeue_ptr) {
    free(equeue_ptr->entry);
}

_INLINE void eQAddQueue(eQueue_type *equeue_ptr, int i) {
    if (equeue_ptr->count >= equeue_ptr->MAXQ) {
	fprintf(stderr,"Queue is FULL\n");
	exit(1);
    }
    else {
	equeue_ptr->count++;
	equeue_ptr->rear = (equeue_ptr->rear + 1) % equeue_ptr->MAXQ;
	equeue_ptr->entry[equeue_ptr->rear] = i;
    }
}

_INLINE void eQDeleteQueue(eQueue_type *equeue_ptr, int *i) {
    if (equeue_ptr->count <= 0) {
	fprintf(stderr,"Queue is EMPTY\n");
	exit(1);
    }
    else {
	equeue_ptr->count--;
	*i = equeue_ptr->entry[equeue_ptr->front];
	equeue_ptr->front = (equeue_ptr->front + 1) % equeue_ptr->MAXQ;
    }
}

_INLINE int eQSize(eQueue_type *equeue_ptr) {
    return (equeue_ptr->count);
}

_INLINE int eQEmpty(eQueue_type *equeue_ptr) {
    return (equeue_ptr->count <= 0);
}

_INLINE int eQFull(eQueue_type *equeue_ptr) {
    return (equeue_ptr->count >= equeue_ptr->MAXQ);
}


_INLINE void BFS_auxGraph(edge_t *auxGraph, int edgeCount, int n, int *auxIndex,
		  int *cc, int *inLabels) {

    eQueue_type Q;
    int i;
    register int
	j, x,
	newLabel;

    int	*mark;
    
    if ((mark = (int*)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: mark arr could not be malloc'ed\n");
    for (i=0 ; i<n ; i++)
	mark[i] = ((auxIndex[i]==NIL) ? 1 : 0);

    eQInitialize(&Q,n);

    for (x=0 ; x<n ; x++)
	if (!mark[x]) {
	    eQAddQueue(&Q, x);
	    mark[x] = 1;
	    newLabel = inLabels[x];
	    do {
		eQDeleteQueue(&Q, &i);
		
		cc[i] = newLabel;

		if ((j = auxIndex[i]) != NIL) 
		    while ((j<edgeCount) && (auxGraph[j].a == i)) {
			if (!mark[auxGraph[j].b]) {
			    eQAddQueue(&Q, auxGraph[j].b);
			    mark[auxGraph[j].b] = 1;
			}
			j++;
		    }
	    } while (!eQEmpty(&Q));
	}

    eQDestruct(&Q);
    free(mark);

}

#endif

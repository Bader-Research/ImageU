/*									tab:8
 *
 * queue.h - queue data structures and algorithms
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
 * Filename:		queue.h
 * History:
 */

#ifndef _QUEUE_H
#define _QUEUE_H

#include "ImageU.h"

typedef struct queue_tag {
    int count;
    int front;
    int rear;
    int MAXQ;
    int *entryI;
    int *entryJ;
} Queue_type;

_INLINE void QInitialize(Queue_type *queue_ptr, int q, int r) {
    queue_ptr->count = 0;
    queue_ptr->front = 0;
    queue_ptr->rear = -1;
    queue_ptr->MAXQ = q*r;
    if ((queue_ptr->entryI = (int*)(malloc(queue_ptr->MAXQ*sizeof(int))))
	==NULL)
	fprintf(stderr,"ERROR: entryI could not be malloc'ed\n");
    if ((queue_ptr->entryJ = (int*)(malloc(queue_ptr->MAXQ*sizeof(int))))
	==NULL)
	fprintf(stderr,"ERROR: entryJ could not be malloc'ed\n");
}

_INLINE void QDestruct(Queue_type *queue_ptr) {
    free(queue_ptr->entryI);
    free(queue_ptr->entryJ);
}

_INLINE void QAddQueue(int i, int j, Queue_type *queue_ptr) {
    if (queue_ptr->count >= queue_ptr->MAXQ) {
	fprintf(stderr,"Queue is FULL\n");
	exit(1);
    }
    else {
	queue_ptr->count++;
	queue_ptr->rear = (queue_ptr->rear + 1) % queue_ptr->MAXQ;
	queue_ptr->entryI[queue_ptr->rear] = i;
	queue_ptr->entryJ[queue_ptr->rear] = j;
    }
}

_INLINE void QDeleteQueue(int *i, int *j, Queue_type *queue_ptr) {
    if (queue_ptr->count <= 0) {
	fprintf(stderr,"Queue is EMPTY\n");
	exit(1);
    }
    else {
	queue_ptr->count--;
	*i = queue_ptr->entryI[queue_ptr->front];
	*j = queue_ptr->entryJ[queue_ptr->front];
	queue_ptr->front = (queue_ptr->front + 1) % queue_ptr->MAXQ;
    }
}

_INLINE int QSize(Queue_type *queue_ptr) {
    return queue_ptr->count;
}

_INLINE int QEmpty(Queue_type *queue_ptr) {
    return (queue_ptr->count <= 0);
}

_INLINE int QFull(Queue_type *queue_ptr) {
    return (queue_ptr->count >= queue_ptr->MAXQ);
}

#endif

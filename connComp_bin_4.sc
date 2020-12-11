/*									tab:8
 *
 * connComp_bin_4.sc - 4-connected components of binary images
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
 * Filename:		connComp_bin_4.sc
 * History:
 */

#include "connComp_bin_4.h"

/*******************************************************/
/* Sequential BFS Connected Components of a Graph in C */
/*******************************************************/

static inline void BFS_binary_4(int v, int w, int q, int r, int ii, int jj,
	 int (*image)[r], int (*mark)[r], int (*lab)[r]) {
    Queue_type Q;
    register int
	i, j, x, y,
	newLabel;
    
    Q.count = 0;
    Q.front = 0;
    Q.rear = -1;
    Q.MAXQ = q*r;
    if ((Q.entryI = (int*)(malloc(Q.MAXQ*sizeof(int))))==NULL)
	fprintf(stderr,"ERROR: entryI could not be malloc'ed\n");
    if ((Q.entryJ = (int*)(malloc(Q.MAXQ*sizeof(int))))==NULL)
	fprintf(stderr,"ERROR: entryJ could not be malloc'ed\n");

    for (x=0 ; x<q ; x++)
	for (y=0 ; y<r ; y++) {
	    if (!mark[x][y]) {
		if (Q.count >= Q.MAXQ) {
		    fprintf(stderr,"Queue is FULL\n");
		    exit(1);
		}
		Q.count++;
		Q.rear = (Q.rear + 1) % Q.MAXQ;
		Q.entryI[Q.rear] = x;
		Q.entryJ[Q.rear] = y;

		mark[x][y] = BLACK;
		newLabel = (ii*q + x)*q*v + (jj*r + y) + 1;
		do {
		    if (Q.count <= 0) {
			fprintf(stderr,"Queue is EMPTY\n");
			exit(1);
		    }
		    Q.count--;
		    i = Q.entryI[Q.front];
		    j = Q.entryJ[Q.front];
		    Q.front = (Q.front + 1) % Q.MAXQ;

		    lab[i][j] = newLabel; /* Visit(i,j) */

		    if ((j<r-1) && !mark[i][j+1]) {
			if (Q.count >= Q.MAXQ) {
			    fprintf(stderr,"Queue is FULL\n");
			    exit(1);
			}
			Q.count++;
			Q.rear = (Q.rear + 1) % Q.MAXQ;
			Q.entryI[Q.rear] = i;
			Q.entryJ[Q.rear] = j+1;
			mark[i][j+1] = BLACK;
		    }
		    if ((i<q-1) && !mark[i+1][j]) {
			if (Q.count >= Q.MAXQ) {
			    fprintf(stderr,"Queue is FULL\n");
			    exit(1);
			}
			Q.count++;
			Q.rear = (Q.rear + 1) % Q.MAXQ;
			Q.entryI[Q.rear] = i+1;
			Q.entryJ[Q.rear] = j;
			mark[i+1][j] = BLACK;
		    }
		    if ((j>0) && !mark[i][j-1]) {
			if (Q.count >= Q.MAXQ) {
			    fprintf(stderr,"Queue is FULL\n");
			    exit(1);
			}
			Q.count++;
			Q.rear = (Q.rear + 1) % Q.MAXQ;
			Q.entryI[Q.rear] = i;
			Q.entryJ[Q.rear] = j-1;
			mark[i][j-1] = BLACK;
		    }
		    if ((i>0) && !mark[i-1][j]) {
			if (Q.count >= Q.MAXQ) {
			    fprintf(stderr,"Queue is FULL\n");
			    exit(1);
			}
			Q.count++;
			Q.rear = (Q.rear + 1) % Q.MAXQ;
			Q.entryI[Q.rear] = i-1;
			Q.entryJ[Q.rear] = j;
			mark[i-1][j] = BLACK;
		    }
		} while (!(Q.count <= 0));
	    }
	}
    free(Q.entryI);
    free(Q.entryJ);
}

static inline void BFS_binary_4_recolor(int q, int r, int x, int y, int (*lab)[r],
		    int (*image)[r], int (*mark)[r], int newLabel) {

    Queue_type Q;
    register int i, j;
    
    if (!mark[x][y]) {
	mark[x][y] = BLACK;

	Q.count = 1;
	Q.front = 0;
	Q.rear = 0;
	Q.MAXQ = q*r;
        if ((Q.entryI = (int*)(malloc(Q.MAXQ*sizeof(int))))==NULL)
    	    fprintf(stderr,"ERROR: entryI could not be malloc'ed\n");
        if ((Q.entryJ = (int*)(malloc(Q.MAXQ*sizeof(int))))==NULL)
            fprintf(stderr,"ERROR: entryJ could not be malloc'ed\n");

	Q.entryI[Q.rear] = x;
	Q.entryJ[Q.rear] = y;
	
	do {
	    if (Q.count <= 0) {
		fprintf(stderr,"Queue is EMPTY\n");
		exit(1);
	    }
	    Q.count--;
	    i = Q.entryI[Q.front];
	    j = Q.entryJ[Q.front];
	    Q.front = (Q.front + 1) % Q.MAXQ;
	    
	    lab[i][j] = newLabel;

	    if ((j<r-1) && !mark[i][j+1]) {
		if (Q.count >= Q.MAXQ) {
		    fprintf(stderr,"Queue is FULL\n");
		    exit(1);
		}
		Q.count++;
		Q.rear = (Q.rear + 1) % Q.MAXQ;
		Q.entryI[Q.rear] = i;
		Q.entryJ[Q.rear] = j+1;
		mark[i][j+1] = BLACK;
	    }
	    if ((i<q-1) && !mark[i+1][j]) {
		if (Q.count >= Q.MAXQ) {
		    fprintf(stderr,"Queue is FULL\n");
		    exit(1);
		}
		Q.count++;
		Q.rear = (Q.rear + 1) % Q.MAXQ;
		Q.entryI[Q.rear] = i+1;
		Q.entryJ[Q.rear] = j;
		mark[i+1][j] = BLACK;
	    }
	    if ((j>0) && !mark[i][j-1]) {
		if (Q.count >= Q.MAXQ) {
		    fprintf(stderr,"Queue is FULL\n");
		    exit(1);
		}
		Q.count++;
		Q.rear = (Q.rear + 1) % Q.MAXQ;
		Q.entryI[Q.rear] = i;
		Q.entryJ[Q.rear] = j-1;
		mark[i][j-1] = BLACK;
	    }
	    if ((i>0) && !mark[i-1][j]) {
		if (Q.count >= Q.MAXQ) {
		    fprintf(stderr,"Queue is FULL\n");
		    exit(1);
		}
		Q.count++;
		Q.rear = (Q.rear + 1) % Q.MAXQ;
		Q.entryI[Q.rear] = i-1;
		Q.entryJ[Q.rear] = j;
		mark[i-1][j] = BLACK;
	    }
	} while (!(Q.count <= 0));
	free(Q.entryI);
	free(Q.entryJ);
    }
}

/*************************************************************/
void all_connected_components_binary_four_sub(int v, int w, int q, int r, 
			      int X[v][w]::[q][r], int k,
			      int Label[v][w]::[q][r]) 
/*************************************************************/
/* Find the connected components of a (vw x qr) = (n x n) binary image
   on a parallel machine with PROCS processors */
#ifndef _COMPILE_LITE
{

    ELEMPTR border_right,     /* The right border to be converted into graph problem */
	    border_left;      /* The  left border to be converted into graph problem */
elem_t
	ShadowBorder[v][w]::[q*v],    /* Shadow Manager's border */
	*locShadowBorder;             /* ptr to shadow manager's border */


    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	logv,                 /* log2(v) */
	logw,                 /* log2(w) */
	vIt,                  /* The iteration num of combining up/down blocks together */
	wIt,                  /*      "        "             left/right blocks together */
	group_manager,        /* Boolean = at this iteration, I do the sequential work */
	shadow_manager,       /* Boolean = at this iteration, I do the shadow work */
	grpManV,              /* if not a group manager, this is the row of my manager */
	grpManW,              /* if not a group manager, this is the col of my manager */
	length,               /* length of the border being merged */
	borderLen,            /* size of the total border merge graph problem */
	mark[q][r],           /* Markers for initial connected components on tile */
	hooks,                /* Number of cc's on our tile */
	beta,                 /* temp variable used in merge update */
	*ccOut,               /* results of initial connected components on a tile */
	*old_label_right,     /* Saved value for original right label */
	*new_label_right,     /* New   value for original right label */
	*old_label_left,      /* Saved value for original left  label */
	*new_label_left,      /* New   value for original left  label */
	*color_left,          /* The pixel value of the left border */
	*color_right,         /* The pixel value of the right border */
	*input_labelling,     /* The labelling passed into connected components */
	changeSize[v][w]::,   /* Proc. (v,w)'s # of (alpha,beta) changes from merge */
	chSize,               /* Local copy of our group manager's changeSize */
	(*locX)[r],           /* local pointer to our subimage */
	(*locLabel)[r],       /* local pointer to our subimage's labels */
	*ccLabel;             /* local pointer to our subimage's labels (linear) */
    
#if 1
    int
	orgV, orgW,       /* original of my group in real coords */
	g_grey,           /* my grey code address inside my group */
	g_GM,             /* the Group Manager's grey code offset */
	g_label,          /* my final "hypercube" label inside my group */
	send_i,	send_j;
#endif

    CHPAIRPTR changeList,     /* Array of changes created by group manager */
	tileHookOrig,         /* Array of original (label,addr) hooks */
	tempHookArr;          /* Temp array of border labels and their addresses */

    TILEHOOKPTR tileHookArr;  /* Array of tile border hooks */

chPair_t
	chList[v][w]::[2*q*v],    /* Group Manager's change list */
	*locChList;               /* ptr to our local list */

#if 0
    char buf[MAXLEN];
#endif

#if 0
    double secs;              /* Timing marker */
#endif

    edge_t *auxGraph;
    int *auxIndex,
	edgeCount;
		
#if 0
    on_one fprintf(outfile,"4-Connected Components of a Binary Image:\n");
#endif

    if (k != 2)
	on_one 
	    fprintf(stderr,"ERROR: A binary image is needed for this algorithm\n");

    barrier();

#if 0
    all_init_timer();

    all_start_timer();
    secs = get_seconds();     /* Begin timing */
#endif

    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locLabel = tolocal(Label[ii][jj]);
    locChList = tolocal(chList[ii][jj]);
    locShadowBorder = tolocal(ShadowBorder[ii][jj]);

    ccLabel = &(locLabel[0][0]);
    
    logv = (int)log2((double)v);
    logw = (int)log2((double)w);
  
/* New tile connected components */

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {
	    locLabel[i][j] = 0;
	    mark[i][j] = (locX[i][j] ? WHITE : BLACK);
	}

    BFS_binary_4(v,w,q,r,ii,jj, locX, mark, locLabel);
    
/* Create tileHookArr with size >= 2q + 2r - 4 (border of tile) */    

    borderLen = 2*(q+r) - 4;

    if ((tempHookArr = (CHPAIRPTR)malloc(borderLen*sizeof(chPair_t)))==NULL)
	fprintf(stderr,"ERROR: tempHookArr could not be malloc'ed\n");

/* For each pixel in the border of the tile, add it to the list if it
   is a unique label > 0 */

    chSize = 0;
    
    for (j=0 ; j<r ; j++) {
	if (locLabel  [0][j] > 0) {
	    tempHookArr[chSize].alpha = locLabel  [0][j];
	    tempHookArr[chSize].beta  = j;
	    chSize++;
	}
	if (locLabel[q-1][j] > 0) {
	    tempHookArr[chSize].alpha = locLabel[q-1][j];
	    tempHookArr[chSize].beta  = (q-1)*r +j;
	    chSize++;
	}
    }

    for (i=1 ; i<q-1 ; i++) {
	if (locLabel[i][0]   > 0) {
	    tempHookArr[chSize].alpha = locLabel[i][0];
	    tempHookArr[chSize].beta  = i*r;
	    chSize++;
	}
	if (locLabel[i][r-1] > 0) {
	    tempHookArr[chSize].alpha = locLabel[i][r-1];
	    tempHookArr[chSize].beta  = i*r + r-1;
	    chSize++;
	}
    }

    hooks = fill_tileHook(&tileHookArr, tempHookArr, chSize);

    if ((tileHookOrig = (CHPAIRPTR)malloc(hooks*sizeof(chPair_t)))==NULL)
	fprintf(stderr,"ERROR: tileHookOrig could not be malloc'ed\n");

    create_hooks(tileHookOrig, tileHookArr, hooks);

#if 0
    all_stop_timer("Initialization");
#endif
    
    barrier();

    length = q*v; /* n = q*v */

    if ((border_right = (ELEMPTR)malloc(length*sizeof(elem_t)))==NULL)
	fprintf(stderr,"ERROR: border_right could not be malloc'ed\n");
    if ((border_left  = (ELEMPTR)malloc(length*sizeof(elem_t)))==NULL)
	fprintf(stderr,"ERROR: border_left could not be malloc'ed\n");
    if ((old_label_right  = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: old_label_right could not be malloc'ed\n");
    if ((new_label_right  = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: new_label_right could not be malloc'ed\n");
    if ((old_label_left   = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: old_label_left could not be malloc'ed\n");
    if ((new_label_left   = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: new_label_left could not be malloc'ed\n");
    if ((color_left       = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: color_left could not be malloc'ed\n");
    if ((color_right      = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: color_right could not be malloc'ed\n");

    if ((ccOut = (int*)malloc(2*length * sizeof(int)))==NULL)  
	fprintf(stderr, "ERROR: ccOut list could not be malloc'ed\n");

    if ((input_labelling = (int*)malloc(2*length * sizeof(int)))==NULL)  
	fprintf(stderr, "ERROR: input_labelling list could not be malloc'ed\n");

    if ((changeList = (CHPAIRPTR)malloc(2*length * sizeof(chPair_t)))==NULL)  
	fprintf(stderr, "ERROR: changeList could not be malloc'ed\n");
    
    if ((auxGraph = (edge_t *)malloc(6*2*length*sizeof(edge_t)))==NULL)
	fprintf(stderr, "ERROR: auxGraph could not be malloc'ed\n");
    if ((auxIndex = (int *)malloc(2*length*sizeof(int)))==NULL)
	fprintf(stderr, "ERROR: auxIndex could not be malloc'ed\n");

/* Combine tiles across processors */

    vIt = 0;
    wIt = 0;

    while ((vIt<logv)||(wIt<logw)) {

	/* Combine left/right */

	if (wIt < logw) {

#if 0
	    all_start_timer();
#endif

	    wIt++;

	    length = q * (1 << vIt);

/* if I am a manager this round: */
/* Then - prefetch the border from w_index ++ */
/*      - find the connected components */
/*      - broadcast the label changes to my clients */
/*      - update my pixels */

/* Group Manager: Last wIt bits of jj = 0111.1 (a zero and wIt-1 1's)
		  Last vIt bits of ii each = 0 */

	    group_manager = ( ((jj & ~(~0<<wIt)) == (1<<(wIt-1))-1)
			      && ((ii & ((1<<vIt)-1) ) == 0) );
	    
/* Shadow Manager: The processor one to the right of the group manager.
                  Last wIt bits of jj = 1000.0 (a one and wIt-1 0's)
		  Last vIt bits of ii each = 0 */

	    shadow_manager = ( ((jj & ~(~0<<wIt)) == (1<<(wIt-1)))
			      && ((ii & ((1<<vIt)-1) ) == 0) );
	    
	    if (group_manager) {

/* Get border from right processor, and set up label arrays */
		
		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			border_left   [(j*q) + i].color :=     X[ii+j][jj][i][r-1];
			border_left   [(j*q) + i].label := Label[ii+j][jj][i][r-1];
		    }

		    for (i=0 ; i<q ; i++) {
			border_left   [(j*q) + i].pos   = (j*q) + i;
			border_right  [(j*q) + i].pos   = (j*q) + i;
		    }

		    sync();
		    for (i=0 ; i<q ; i++) {
			color_left     [(j*q) + i] = border_left [(j*q) + i].color;
			old_label_left [(j*q) + i] = border_left [(j*q) + i].label;
		    }
		}
		
/* Prefetch Right Border */

		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			border_right  [(j*q) + i].color :=     X[ii+j][jj+1][i][0];
			border_right  [(j*q) + i].label := Label[ii+j][jj+1][i][0];
		    }
		}
		
/* Sort each border list by label, and then move through the list linearly, 
   attaching each pair that are of the same label */
		
		fastsort_elem(border_left, length);

		sync();

		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			color_right    [(j*q) + i] = border_right[(j*q) + i].color;
			old_label_right[(j*q) + i] = border_right[(j*q) + i].label;
		    }
		}
	    }
	    
	    if (shadow_manager) {

/* Get border from processors, and set up label arrays */
		
		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			locShadowBorder[(j*q) + i].color :=     X[ii+j][jj][i][0];
			locShadowBorder[(j*q) + i].label := Label[ii+j][jj][i][0];
		    }

		    for (i=0 ; i<q ; i++) 
			locShadowBorder[(j*q) + i].pos   = (j*q) + i;

		    sync();
		}
		
		fastsort_elem(locShadowBorder,length);
	    }
	    
	    barrier();

/* Group Manager should retrieve results from Shadow manager */

	    if (group_manager) {

		bulk_get(border_right, &(ShadowBorder[ii][jj+1][0]), 
			     length*sizeof(elem_t)); 

		borderLen = 2*length;

		/* Each vertex can have at most 5 edges incident to it */
		/* There are at most 2*(5*borderLen - 4) edges in the graph */
		edgeCount = 0;
		
/* Add edges */
/* First, add edges stringing along pixels of the same label down each of 
   the borders */
/* Then add edges representing links between the two borders */
		
		for (i=1 ; i<length ; i++) {
		    if ((border_left[i-1].label == border_left[i].label)
			&& (border_left[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				border_left[i-1].pos,
				border_left[i].pos);
		}

		sync();
		
		for (i=1 ; i<length ; i++) {
		    if ((border_right[i-1].label == border_right[i].label)
			&& (border_right[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				length + border_right[i-1].pos,
				length + border_right[i].pos);
		}
		
/* Scan down the left border and add edges wherever it meets the right bord */
		
		for (i=0 ; i<length ; i++) 
		    if (color_left[i]==color_right[i]) 
			addAuxEdge(auxGraph, &edgeCount, i, length + i);

		if (edgeCount > 0) {
		/* Sort the edges, then create an auxIndex */
		    fastsort_edge(auxGraph, edgeCount);

		    for (i=0 ; i<borderLen ; i++)
			auxIndex[i] = NIL;
		
		    auxIndex[auxGraph[0].a] = 0;

		    if (edgeCount > 1)
			for (i=1 ; i<edgeCount ; i++) 
			    if (auxGraph[i].a > auxGraph[i-1].a)
				auxIndex[auxGraph[i].a] = i;

		    for (i=0 ; i<length ; i++) {
			input_labelling[i]        = old_label_left[i];
			input_labelling[i+length] = old_label_right[i];
			ccOut[i]        = input_labelling[i];
			ccOut[i+length] = input_labelling[i+length];
		    }
		
		    BFS_auxGraph(auxGraph, edgeCount, borderLen, auxIndex,
				 ccOut, input_labelling);
		}

		chSize = 0;
		
		if (edgeCount > 0) {

		    for (i=0 ; i<length ; i++) {
			new_label_left[i]  = ccOut[i];
			new_label_right[i] = ccOut[length+i];
		    }

/* Compare old_label's with new_label's and propogate the differences */

		    for (i=0 ; i<length ; i++) {
			if (old_label_left[i] != new_label_left[i]) {
			    changeList[chSize].alpha = old_label_left[i];
			    changeList[chSize].beta  = new_label_left[i];
			    chSize++;
			}
			if (old_label_right[i] != new_label_right[i]) {
			    changeList[chSize].alpha = old_label_right[i];
			    changeList[chSize].beta  = new_label_right[i];
			    chSize++;
			}
		    }

		    fastsort_chPair(changeList, chSize);

		    chSize = fill_chList(locChList, changeList, chSize);
		}

		changeSize [ii][jj] = chSize;
		grpManV = ii;
		grpManW = jj;

	    }
	    else {  /* Figure out who my group_manager is */
		/* ii with the last vIt bits set to 0 */
		grpManV = (ii >> vIt) << vIt;
		/* jj with the last wIt bits set to 0111.1 (a zero and wIt-1 1's) */
		grpManW = (((jj >> wIt) << wIt) | ((1<<(wIt-1))-1));
	    } 
	    
	    barrier();

#if 0
/* Prefetch Change List from Group Manager */
	    
	    if (!group_manager) {
		chSize = changeSize[grpManV][grpManW];
		
/* Fix this prefetch later to be a better broadcast algorithm */

		if (chSize > 0) {
		    bulk_get(locChList, &(chList[grpManV][grpManW][0]), 
			     chSize*sizeof(chPair_t)); 
		    sync();
		}

	    }
#endif
#if 0
	    all_mesh_bcast(v, w, (q*v) << 1,
			   vIt, wIt,		       
			   ii, jj, grpManV, grpManW,
			   chList, changeSize);
	    chSize = changeSize[ii][jj];
#endif
#if 1
	    orgV    = (ii>>vIt) << vIt;
	    orgW    = (jj>>wIt) << wIt;
	    g_grey  = greyEncoder(vIt, wIt, ii-orgV, jj-orgW);
	    g_GM    = greyEncoder(vIt, wIt, grpManV-orgV, grpManW-orgW);
	    g_label = g_grey ^ g_GM;
	    j       = vIt + wIt;

	    for (i=0 ; i<j ; i++) {
		if (g_label < (1<<i)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
		    changeSize[send_i+orgV][send_j+orgW]
			= changeSize[ii][jj];
		}
		barrier();
	    }

	    chSize = changeSize[ii][jj];

	    for (i=0 ; i<j ; i++) {
		if ((g_label < (1<<i)) && (chSize > 0)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) ); 
		    bulk_write(&(chList[send_i+orgV][send_j+orgW][0]),
			       tolocal(chList[ii][jj]),
			       chSize*sizeof(chPair_t));
		}
		barrier();
	    }	    
#endif

/* Make changes */

/* We have a sorted array of changes of size chSize,
           a data structure of tile border hooks.
   Recolor only the border of each tile by binary searching the
   label change for each hook's label, and then updating the label
   at each of the addresses in the addr. list. */
	    
	    for (i=0 ; i<hooks ; i++) {          /* For Each Hook */
		beta = chLookup(locChList, chSize, tileHookArr[i].label);
		if (tileHookArr[i].label != beta) {
		    for (j=0 ; j < tileHookArr[i].size ; j++) 
			ccLabel[(tileHookArr[i].addr)[j]] = beta;
		    tileHookArr[i].label = beta;
		}
	    }

#if 0
	    sprintf(buf,"Merge Phase %d",vIt+wIt);
	    all_stop_timer(buf);
#endif

	    barrier();

	}
	
/* Combine up/down */
/* left == up , right == down */
	
	if (vIt < logv) {

#if 0
	    all_start_timer();
#endif
	    
	    vIt++;

	    length = r * (1 << wIt);

/* if I am a manager this round: */
/* Then - prefetch the border from v_index ++ */
/*      - find the connected components */
/*      - broadcast the label changes to my clients */
/*      - update my pixels */

/* Group Manager: Last vIt bits of ii = 0111.1 (a zero and vIt-1 1's)
		  Last wIt bits of jj each = 0 */

	    group_manager = ( ((ii & ~(~0<<vIt)) == (1<<(vIt-1))-1)
			      && ((jj & ((1<<wIt)-1) ) == 0) );

/* Shadow Manager: The processor one to the bottom of the group manager.
                  Last vIt bits of ii = 1000.0 (a one and vIt-1 0's)
		  Last wIt bits of jj each = 0 */

	    shadow_manager = ( ((ii & ~(~0<<vIt)) == (1<<(vIt-1)))
			      && ((jj & ((1<<wIt)-1) ) == 0) );
	    
	    if (group_manager) {

/* Get border from right processor, and set up label arrays */
		
		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			border_left   [(j*r) + i].color :=     X[ii ][jj+j][q-1][i];
			border_left   [(j*r) + i].label := Label[ii ][jj+j][q-1][i];
		    }

		    for (i=0 ; i<r ; i++) {
			border_left   [(j*r) + i].pos   = (j*r) + i;
			border_right  [(j*r) + i].pos   = (j*r) + i;
		    }
		    sync();
		    for (i=0 ; i<r ; i++) {
			color_left     [(j*r) + i] = border_left [(j*r) + i].color;
			old_label_left [(j*r) + i] = border_left [(j*r) + i].label;
		    }
		}
		
/* Prefetch Right Border */

		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			border_right  [(j*r) + i].color :=     X[ii+1][jj+j][0][i];
			border_right  [(j*r) + i].label := Label[ii+1][jj+j][0][i];
		    }
		}
		
/* Sort each border list by label, and then move through the list linearly, 
   attaching each pair that are of the same label */
		
		fastsort_elem(border_left, length);

		sync();

		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			color_right    [(j*r) + i] = border_right[(j*r) + i].color;
			old_label_right[(j*r) + i] = border_right[(j*r) + i].label;
		    }
		}
	    }	    

	    if (shadow_manager) {

/* Get border from processors, and set up label arrays */
		
		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			locShadowBorder[(j*r) + i].color :=     X[ii][jj+j][0][i];
			locShadowBorder[(j*r) + i].label := Label[ii][jj+j][0][i];
		    }

		    for (i=0 ; i<r ; i++) 
			locShadowBorder[(j*r) + i].pos   = (j*r) + i;

		    sync();
		}
		
		fastsort_elem(locShadowBorder,length);
	    }
	    
	    barrier();

/* Group Manager should retrieve results from Shadow manager */

	    if (group_manager) {

		bulk_get(border_right, &(ShadowBorder[ii+1][jj][0]), 
			     length*sizeof(elem_t)); 

/* Create Adjacency List problem */

		borderLen = 2*length;

		/* Each vertex can have at most 5 edges incident to it */
		/* There are at most 2*(5*borderLen - 4) edges in the graph */
		edgeCount = 0;
		
/* Add edges */
/* First, add edges stringing along pixels of the same label down each of 
   the borders */
/* Then add edges representing links between the two borders */
		
		for (i=1 ; i<length ; i++) {
		    if ((border_left[i-1].label == border_left[i].label)
			&& (border_left[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				border_left[i-1].pos,
				border_left[i].pos);
		}

		sync();
		
		for (i=1 ; i<length ; i++) {
		    if ((border_right[i-1].label == border_right[i].label)
			&& (border_right[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				length + border_right[i-1].pos,
				length + border_right[i].pos);
		}
		
/* Scan down the left border and add edges wherever it meets the right bord */
		
		for (i=0 ; i<length ; i++)
		    if (color_left[i]==color_right[i]) 
			    addAuxEdge(auxGraph, &edgeCount, i, length + i );
		
		if (edgeCount > 0) {
		/* Sort the edges, then create an auxIndex */
		    fastsort_edge(auxGraph, edgeCount);

		    for (i=0 ; i<borderLen ; i++)
			auxIndex[i] = NIL;
		
		    auxIndex[auxGraph[0].a] = 0;

		    if (edgeCount > 1)
			for (i=1 ; i<edgeCount ; i++) 
			    if (auxGraph[i].a > auxGraph[i-1].a)
				auxIndex[auxGraph[i].a] = i;

		    for (i=0 ; i<length ; i++) {
			input_labelling[i]        = old_label_left[i];
			input_labelling[i+length] = old_label_right[i];
			ccOut[i]        = input_labelling[i];
			ccOut[i+length] = input_labelling[i+length];
		    }
		
		    BFS_auxGraph(auxGraph, edgeCount, borderLen, auxIndex,
				 ccOut, input_labelling);
		}

/* Compare old_label's with new_label's and propogate the differences */
			        
		chSize = 0;
		
		if (edgeCount > 0) {

		    for (i=0 ; i<length ; i++) {
			new_label_left[i]  = ccOut[i];
			new_label_right[i] = ccOut[length+i];
		    }

		    for (i=0 ; i<length ; i++) {
			if (old_label_left[i] != new_label_left[i]) {
			    changeList[chSize].alpha = old_label_left[i];
			    changeList[chSize].beta  = new_label_left[i];
			    chSize++;
			}
			if (old_label_right[i] != new_label_right[i]) {
			    changeList[chSize].alpha = old_label_right[i];
			    changeList[chSize].beta  = new_label_right[i];
			    chSize++;
			}
		    }
		
		    fastsort_chPair(changeList, chSize);

		    chSize = fill_chList(locChList, changeList, chSize);

		}

		changeSize [ii][jj] = chSize;
		grpManV = ii;
		grpManW = jj;

	    }
	    else {  /* Figure out who my group_manager is */
		/* ii with the last vIt bits set to 0111.1 (a zero and vIt-1 1's) */
		grpManV = (((ii >> vIt) << vIt) | ((1<<(vIt-1))-1));
		/* jj with the last wIt bits set to 0 */
		grpManW = (jj >> wIt) << wIt;
	    }

	    barrier();

#if 0
/* Prefetch Change List from Group Manager */
	    
	    if (!group_manager) {

		chSize = changeSize[grpManV][grpManW];

/* Fix this prefetch later to be a better broadcast algorithm */

		if (chSize > 0) {
		    bulk_get(locChList, &(chList[grpManV][grpManW][0]), 
			     chSize*sizeof(chPair_t)); 
		    sync();
		}
	    }
#endif
#if 0
	    all_mesh_bcast(v, w, (q*v) << 1,
			   vIt, wIt,		       
			   ii, jj, grpManV, grpManW,
			   chList, changeSize);
	    chSize = changeSize[ii][jj];
#endif
#if 1
	    orgV    = (ii>>vIt) << vIt;
	    orgW    = (jj>>wIt) << wIt;
	    g_grey  = greyEncoder(vIt, wIt, ii-orgV, jj-orgW);
	    g_GM    = greyEncoder(vIt, wIt, grpManV-orgV, grpManW-orgW);
	    g_label = g_grey ^ g_GM;
	    j       = vIt + wIt;

	    for (i=0 ; i<j ; i++) {
		if (g_label < (1<<i)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
		    changeSize[send_i+orgV][send_j+orgW]
			= changeSize[ii][jj];
		}
		barrier();
	    }

	    chSize = changeSize[ii][jj];

	    for (i=0 ; i<j ; i++) {
		if ((g_label < (1<<i)) && (chSize > 0)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) ); 
		    bulk_write(&(chList[send_i+orgV][send_j+orgW][0]),
			       tolocal(chList[ii][jj]),
			       chSize*sizeof(chPair_t));
		}
		barrier();
	    }	    
#endif

/* Make changes */

/* We have a sorted array of changes of size chSize,
           a data structure of tile border hooks.
   Recolor only the border of each tile by binary searching the
   label change for each hook's label, and then updating the label
   at each of the addresses in the addr. list. */
	    
	    for (i=0 ; i<hooks ; i++) {          /* For Each Hook */
		beta = chLookup(locChList, chSize, tileHookArr[i].label);
		if (tileHookArr[i].label != beta) {
		    for (j=0 ; j < tileHookArr[i].size ; j++)
			ccLabel[(tileHookArr[i].addr)[j]] = beta;
		    tileHookArr[i].label = beta;
		}
	    }

#if 0
	    sprintf(buf,"Merge Phase %d",vIt+wIt);
	    all_stop_timer(buf);
#endif

	    barrier();

	}
    }
    
/* We must relabel the interior of our tile by searching through
   tileHookOrig, and if any labels have changed, perform a depth-first search
   and recolor connected pixels */

#if 0
    all_start_timer();
#endif
          
    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) 
	    mark[i][j] = (locX[i][j] ? WHITE : BLACK);
	    
    for (i=0 ; i<hooks ; i++)
	if (tileHookOrig[i].alpha != tileHookArr[i].label)
	    /* Update component */
	    /* This will check all adjacent pixels */

            BFS_binary_4_recolor(q,r,
			       tileHookOrig[i].beta / r,
			       tileHookOrig[i].beta % r,
			       locLabel, locX, mark, tileHookArr[i].label);

#if 0
    all_stop_timer("Update Interior Labels");
#endif

    destroy_tileHook(tileHookArr,hooks);    
    free(tileHookArr);
    free(tileHookOrig);
    free(tempHookArr);

    free(auxIndex);
    free(auxGraph);
    free(input_labelling);
    free(ccOut);
    free(border_right);
    free(border_left);
    free(color_right);
    free(color_left);
    free(old_label_right);
    free(new_label_right);
    free(old_label_left);
    free(new_label_left);

    barrier();

#if 0
    secs = get_seconds() - secs;  /* Finish Timing */
#endif

/*****************************************************/
/* Results                                           */
/*****************************************************/

/* DONE - output results of connected components */

 
#if 0
    on_one {
	fprintf(outfile,"4-Connected Components Time: %f  ",secs);
	fprintf(outfile,"for n x n = %d image on %d processors.\n", v*w*q*r, PROCS);

	if (RESULTS) {
	    fprintf(outfile,"Original Binary Image:\n");
	    for (i=0 ; i<v ; i++) 
		for (ii=0 ; ii < q ; ii++) {
		    for (j=0 ; j<w ; j++)
			for (jj=0 ; jj < r ; jj++)
			    fprintf(outfile,"%1d ",X[i][j][ii][jj]);
		    fprintf(outfile,"\n");
		}
	    fprintf(outfile,"\n");

	    fprintf(outfile,"4-Connected Components:\n");
	    for (i=0 ; i<v ; i++) 
		for (ii=0 ; ii < q ; ii++) {
		    for (j=0 ; j<w ; j++)
			for (jj=0 ; jj < r ; jj++)
			    fprintf(outfile,"%3d ",Label[i][j][ii][jj]);
		    fprintf(outfile,"\n");
		}
	    fprintf(outfile,"\n");
	}

	if (OUTPUT_IMAGE) {
	    fprintf(outfile,"Writing Output Image... ");
	    fflush(outfile);
	    image_print_cc("4-connected components of binary image",
			   "image4ccbin",
			   v,w,q,r,k,Label);
	    fprintf(outfile,"Done.\n");
	    fflush(outfile);
	}

	all_print_timer(outfile);
	    
    }
#endif
    
}
#else
{}
#endif

/*************************************************************/
void all_connected_components_binary_four(int v, int w, int q, int r, 
			      int X[v][w]::[q][r], int k)
/*************************************************************/
/* Find the connected components of a (vw x qr) = (n x n) binary image
   on a parallel machine with PROCS processors */
#ifndef _COMPILE_LITE
{

    ELEMPTR border_right,     /* The right border to be converted into graph problem */
	    border_left;      /* The  left border to be converted into graph problem */
elem_t
	ShadowBorder[v][w]::[q*v],    /* Shadow Manager's border */
	*locShadowBorder;             /* ptr to shadow manager's border */


    register int
	i,                    /* Loop variable */
	j,                    /* Loop variable */
	ii,                   /* MYPROC's row */
	jj;                   /* MYPROC's column */

    int	logv,                 /* log2(v) */
	logw,                 /* log2(w) */
	vIt,                  /* The iteration num of combining up/down blocks together */
	wIt,                  /*      "        "             left/right blocks together */
	group_manager,        /* Boolean = at this iteration, I do the sequential work */
	shadow_manager,       /* Boolean = at this iteration, I do the shadow work */
	grpManV,              /* if not a group manager, this is the row of my manager */
	grpManW,              /* if not a group manager, this is the col of my manager */
	length,               /* length of the border being merged */
	borderLen,            /* size of the total border merge graph problem */
	mark[q][r],           /* Markers for initial connected components on tile */
	hooks,                /* Number of cc's on our tile */
	beta,                 /* temp variable used in merge update */
	*ccOut,               /* results of initial connected components on a tile */
	*old_label_right,     /* Saved value for original right label */
	*new_label_right,     /* New   value for original right label */
	*old_label_left,      /* Saved value for original left  label */
	*new_label_left,      /* New   value for original left  label */
	*color_left,          /* The pixel value of the left border */
	*color_right,         /* The pixel value of the right border */
	*input_labelling,     /* The labelling passed into connected components */
	changeSize[v][w]::,   /* Proc. (v,w)'s # of (alpha,beta) changes from merge */
	chSize,               /* Local copy of our group manager's changeSize */
	(*locX)[r],           /* local pointer to our subimage */
	(*locLabel)[r],       /* local pointer to our subimage's labels */
	*ccLabel;             /* local pointer to our subimage's labels (linear) */
    
#if 1
    int
	orgV, orgW,       /* original of my group in real coords */
	g_grey,           /* my grey code address inside my group */
	g_GM,             /* the Group Manager's grey code offset */
	g_label,          /* my final "hypercube" label inside my group */
	send_i,	send_j;
#endif

    CHPAIRPTR changeList,     /* Array of changes created by group manager */
	tileHookOrig,         /* Array of original (label,addr) hooks */
	tempHookArr;          /* Temp array of border labels and their addresses */

    TILEHOOKPTR tileHookArr;  /* Array of tile border hooks */

chPair_t
	chList[v][w]::[2*q*v],    /* Group Manager's change list */
	*locChList;               /* ptr to our local list */

    char buf[MAXLEN];

    double secs;              /* Timing marker */

    edge_t *auxGraph;
    int *auxIndex,
	edgeCount;
		
    int Label[v][w]::[q][r];

    on_one fprintf(outfile,"4-Connected Components of a Binary Image:\n");

    if (k != 2)
	on_one 
	    fprintf(stderr,"ERROR: A binary image is needed for this algorithm\n");

    barrier();

    all_init_timer();

    all_start_timer();
    secs = get_seconds();     /* Begin timing */

    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    locX = tolocal(X[ii][jj]);
    locLabel = tolocal(Label[ii][jj]);
    locChList = tolocal(chList[ii][jj]);
    locShadowBorder = tolocal(ShadowBorder[ii][jj]);

    ccLabel = &(locLabel[0][0]);
    
    logv = (int)log2((double)v);
    logw = (int)log2((double)w);
  
/* New tile connected components */

    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) {
	    locLabel[i][j] = 0;
	    mark[i][j] = (locX[i][j] ? WHITE : BLACK);
	}

    BFS_binary_4(v,w,q,r,ii,jj, locX, mark, locLabel);
    
/* Create tileHookArr with size >= 2q + 2r - 4 (border of tile) */    

    borderLen = 2*(q+r) - 4;

    if ((tempHookArr = (CHPAIRPTR)malloc(borderLen*sizeof(chPair_t)))==NULL)
	fprintf(stderr,"ERROR: tempHookArr could not be malloc'ed\n");

/* For each pixel in the border of the tile, add it to the list if it
   is a unique label > 0 */

    chSize = 0;
    
    for (j=0 ; j<r ; j++) {
	if (locLabel  [0][j] > 0) {
	    tempHookArr[chSize].alpha = locLabel  [0][j];
	    tempHookArr[chSize].beta  = j;
	    chSize++;
	}
	if (locLabel[q-1][j] > 0) {
	    tempHookArr[chSize].alpha = locLabel[q-1][j];
	    tempHookArr[chSize].beta  = (q-1)*r +j;
	    chSize++;
	}
    }

    for (i=1 ; i<q-1 ; i++) {
	if (locLabel[i][0]   > 0) {
	    tempHookArr[chSize].alpha = locLabel[i][0];
	    tempHookArr[chSize].beta  = i*r;
	    chSize++;
	}
	if (locLabel[i][r-1] > 0) {
	    tempHookArr[chSize].alpha = locLabel[i][r-1];
	    tempHookArr[chSize].beta  = i*r + r-1;
	    chSize++;
	}
    }

    hooks = fill_tileHook(&tileHookArr, tempHookArr, chSize);

    if ((tileHookOrig = (CHPAIRPTR)malloc(hooks*sizeof(chPair_t)))==NULL)
	fprintf(stderr,"ERROR: tileHookOrig could not be malloc'ed\n");

    create_hooks(tileHookOrig, tileHookArr, hooks);

    all_stop_timer("Initialization");
    
    barrier();

    length = q*v; /* n = q*v */

    if ((border_right = (ELEMPTR)malloc(length*sizeof(elem_t)))==NULL)
	fprintf(stderr,"ERROR: border_right could not be malloc'ed\n");
    if ((border_left  = (ELEMPTR)malloc(length*sizeof(elem_t)))==NULL)
	fprintf(stderr,"ERROR: border_left could not be malloc'ed\n");
    if ((old_label_right  = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: old_label_right could not be malloc'ed\n");
    if ((new_label_right  = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: new_label_right could not be malloc'ed\n");
    if ((old_label_left   = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: old_label_left could not be malloc'ed\n");
    if ((new_label_left   = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: new_label_left could not be malloc'ed\n");
    if ((color_left       = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: color_left could not be malloc'ed\n");
    if ((color_right      = (int*)malloc(length*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: color_right could not be malloc'ed\n");

    if ((ccOut = (int*)malloc(2*length * sizeof(int)))==NULL)  
	fprintf(stderr, "ERROR: ccOut list could not be malloc'ed\n");

    if ((input_labelling = (int*)malloc(2*length * sizeof(int)))==NULL)  
	fprintf(stderr, "ERROR: input_labelling list could not be malloc'ed\n");

    if ((changeList = (CHPAIRPTR)malloc(2*length * sizeof(chPair_t)))==NULL)  
	fprintf(stderr, "ERROR: changeList could not be malloc'ed\n");
    
    if ((auxGraph = (edge_t *)malloc(6*2*length*sizeof(edge_t)))==NULL)
	fprintf(stderr, "ERROR: auxGraph could not be malloc'ed\n");
    if ((auxIndex = (int *)malloc(2*length*sizeof(int)))==NULL)
	fprintf(stderr, "ERROR: auxIndex could not be malloc'ed\n");

/* Combine tiles across processors */

    vIt = 0;
    wIt = 0;

    while ((vIt<logv)||(wIt<logw)) {

	/* Combine left/right */

	if (wIt < logw) {

	    all_start_timer();

	    wIt++;

	    length = q * (1 << vIt);

/* if I am a manager this round: */
/* Then - prefetch the border from w_index ++ */
/*      - find the connected components */
/*      - broadcast the label changes to my clients */
/*      - update my pixels */

/* Group Manager: Last wIt bits of jj = 0111.1 (a zero and wIt-1 1's)
		  Last vIt bits of ii each = 0 */

	    group_manager = ( ((jj & ~(~0<<wIt)) == (1<<(wIt-1))-1)
			      && ((ii & ((1<<vIt)-1) ) == 0) );
	    
/* Shadow Manager: The processor one to the right of the group manager.
                  Last wIt bits of jj = 1000.0 (a one and wIt-1 0's)
		  Last vIt bits of ii each = 0 */

	    shadow_manager = ( ((jj & ~(~0<<wIt)) == (1<<(wIt-1)))
			      && ((ii & ((1<<vIt)-1) ) == 0) );
	    
	    if (group_manager) {

/* Get border from right processor, and set up label arrays */
		
		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			border_left   [(j*q) + i].color :=     X[ii+j][jj][i][r-1];
			border_left   [(j*q) + i].label := Label[ii+j][jj][i][r-1];
		    }

		    for (i=0 ; i<q ; i++) {
			border_left   [(j*q) + i].pos   = (j*q) + i;
			border_right  [(j*q) + i].pos   = (j*q) + i;
		    }

		    sync();
		    for (i=0 ; i<q ; i++) {
			color_left     [(j*q) + i] = border_left [(j*q) + i].color;
			old_label_left [(j*q) + i] = border_left [(j*q) + i].label;
		    }
		}
		
/* Prefetch Right Border */

		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			border_right  [(j*q) + i].color :=     X[ii+j][jj+1][i][0];
			border_right  [(j*q) + i].label := Label[ii+j][jj+1][i][0];
		    }
		}
		
/* Sort each border list by label, and then move through the list linearly, 
   attaching each pair that are of the same label */
		
		fastsort_elem(border_left, length);

		sync();

		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			color_right    [(j*q) + i] = border_right[(j*q) + i].color;
			old_label_right[(j*q) + i] = border_right[(j*q) + i].label;
		    }
		}
	    }
	    
	    if (shadow_manager) {

/* Get border from processors, and set up label arrays */
		
		for (j=0 ; j < (1<<vIt) ; j++) {
		    for (i=0 ; i<q ; i++) {
			locShadowBorder[(j*q) + i].color :=     X[ii+j][jj][i][0];
			locShadowBorder[(j*q) + i].label := Label[ii+j][jj][i][0];
		    }

		    for (i=0 ; i<q ; i++) 
			locShadowBorder[(j*q) + i].pos   = (j*q) + i;

		    sync();
		}
		
		fastsort_elem(locShadowBorder,length);
	    }
	    
	    barrier();

/* Group Manager should retrieve results from Shadow manager */

	    if (group_manager) {

		bulk_get(border_right, &(ShadowBorder[ii][jj+1][0]), 
			     length*sizeof(elem_t)); 

		borderLen = 2*length;

		/* Each vertex can have at most 5 edges incident to it */
		/* There are at most 2*(5*borderLen - 4) edges in the graph */
		edgeCount = 0;
		
/* Add edges */
/* First, add edges stringing along pixels of the same label down each of 
   the borders */
/* Then add edges representing links between the two borders */
		
		for (i=1 ; i<length ; i++) {
		    if ((border_left[i-1].label == border_left[i].label)
			&& (border_left[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				border_left[i-1].pos,
				border_left[i].pos);
		}

		sync();
		
		for (i=1 ; i<length ; i++) {
		    if ((border_right[i-1].label == border_right[i].label)
			&& (border_right[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				length + border_right[i-1].pos,
				length + border_right[i].pos);
		}
		
/* Scan down the left border and add edges wherever it meets the right bord */
		
		for (i=0 ; i<length ; i++) 
		    if (color_left[i]==color_right[i]) 
			addAuxEdge(auxGraph, &edgeCount, i, length + i);

		if (edgeCount > 0) {
		/* Sort the edges, then create an auxIndex */
		    fastsort_edge(auxGraph, edgeCount);

		    for (i=0 ; i<borderLen ; i++)
			auxIndex[i] = NIL;
		
		    auxIndex[auxGraph[0].a] = 0;

		    if (edgeCount > 1)
			for (i=1 ; i<edgeCount ; i++) 
			    if (auxGraph[i].a > auxGraph[i-1].a)
				auxIndex[auxGraph[i].a] = i;

		    for (i=0 ; i<length ; i++) {
			input_labelling[i]        = old_label_left[i];
			input_labelling[i+length] = old_label_right[i];
			ccOut[i]        = input_labelling[i];
			ccOut[i+length] = input_labelling[i+length];
		    }
		
		    BFS_auxGraph(auxGraph, edgeCount, borderLen, auxIndex,
				 ccOut, input_labelling);
		}

		chSize = 0;
		
		if (edgeCount > 0) {

		    for (i=0 ; i<length ; i++) {
			new_label_left[i]  = ccOut[i];
			new_label_right[i] = ccOut[length+i];
		    }

/* Compare old_label's with new_label's and propogate the differences */

		    for (i=0 ; i<length ; i++) {
			if (old_label_left[i] != new_label_left[i]) {
			    changeList[chSize].alpha = old_label_left[i];
			    changeList[chSize].beta  = new_label_left[i];
			    chSize++;
			}
			if (old_label_right[i] != new_label_right[i]) {
			    changeList[chSize].alpha = old_label_right[i];
			    changeList[chSize].beta  = new_label_right[i];
			    chSize++;
			}
		    }

		    fastsort_chPair(changeList, chSize);

		    chSize = fill_chList(locChList, changeList, chSize);
		}

		changeSize [ii][jj] = chSize;
		grpManV = ii;
		grpManW = jj;

	    }
	    else {  /* Figure out who my group_manager is */
		/* ii with the last vIt bits set to 0 */
		grpManV = (ii >> vIt) << vIt;
		/* jj with the last wIt bits set to 0111.1 (a zero and wIt-1 1's) */
		grpManW = (((jj >> wIt) << wIt) | ((1<<(wIt-1))-1));
	    } 
	    
	    barrier();

#if 0
/* Prefetch Change List from Group Manager */
	    
	    if (!group_manager) {
		chSize = changeSize[grpManV][grpManW];
		
/* Fix this prefetch later to be a better broadcast algorithm */

		if (chSize > 0) {
		    bulk_get(locChList, &(chList[grpManV][grpManW][0]), 
			     chSize*sizeof(chPair_t)); 
		    sync();
		}

	    }
#endif
#if 0
	    all_mesh_bcast(v, w, (q*v) << 1,
			   vIt, wIt,		       
			   ii, jj, grpManV, grpManW,
			   chList, changeSize);
	    chSize = changeSize[ii][jj];
#endif
#if 1
	    orgV    = (ii>>vIt) << vIt;
	    orgW    = (jj>>wIt) << wIt;
	    g_grey  = greyEncoder(vIt, wIt, ii-orgV, jj-orgW);
	    g_GM    = greyEncoder(vIt, wIt, grpManV-orgV, grpManW-orgW);
	    g_label = g_grey ^ g_GM;
	    j       = vIt + wIt;

	    for (i=0 ; i<j ; i++) {
		if (g_label < (1<<i)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
		    changeSize[send_i+orgV][send_j+orgW]
			= changeSize[ii][jj];
		}
		barrier();
	    }

	    chSize = changeSize[ii][jj];

	    for (i=0 ; i<j ; i++) {
		if ((g_label < (1<<i)) && (chSize > 0)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) ); 
		    bulk_write(&(chList[send_i+orgV][send_j+orgW][0]),
			       tolocal(chList[ii][jj]),
			       chSize*sizeof(chPair_t));
		}
		barrier();
	    }	    
#endif

/* Make changes */

/* We have a sorted array of changes of size chSize,
           a data structure of tile border hooks.
   Recolor only the border of each tile by binary searching the
   label change for each hook's label, and then updating the label
   at each of the addresses in the addr. list. */
	    
	    for (i=0 ; i<hooks ; i++) {          /* For Each Hook */
		beta = chLookup(locChList, chSize, tileHookArr[i].label);
		if (tileHookArr[i].label != beta) {
		    for (j=0 ; j < tileHookArr[i].size ; j++) 
			ccLabel[(tileHookArr[i].addr)[j]] = beta;
		    tileHookArr[i].label = beta;
		}
	    }

	    sprintf(buf,"Merge Phase %d",vIt+wIt);
	    all_stop_timer(buf);

	    barrier();

	}
	
/* Combine up/down */
/* left == up , right == down */
	
	if (vIt < logv) {

	    all_start_timer();
	    
	    vIt++;

	    length = r * (1 << wIt);

/* if I am a manager this round: */
/* Then - prefetch the border from v_index ++ */
/*      - find the connected components */
/*      - broadcast the label changes to my clients */
/*      - update my pixels */

/* Group Manager: Last vIt bits of ii = 0111.1 (a zero and vIt-1 1's)
		  Last wIt bits of jj each = 0 */

	    group_manager = ( ((ii & ~(~0<<vIt)) == (1<<(vIt-1))-1)
			      && ((jj & ((1<<wIt)-1) ) == 0) );

/* Shadow Manager: The processor one to the bottom of the group manager.
                  Last vIt bits of ii = 1000.0 (a one and vIt-1 0's)
		  Last wIt bits of jj each = 0 */

	    shadow_manager = ( ((ii & ~(~0<<vIt)) == (1<<(vIt-1)))
			      && ((jj & ((1<<wIt)-1) ) == 0) );
	    
	    if (group_manager) {

/* Get border from right processor, and set up label arrays */
		
		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			border_left   [(j*r) + i].color :=     X[ii ][jj+j][q-1][i];
			border_left   [(j*r) + i].label := Label[ii ][jj+j][q-1][i];
		    }

		    for (i=0 ; i<r ; i++) {
			border_left   [(j*r) + i].pos   = (j*r) + i;
			border_right  [(j*r) + i].pos   = (j*r) + i;
		    }
		    sync();
		    for (i=0 ; i<r ; i++) {
			color_left     [(j*r) + i] = border_left [(j*r) + i].color;
			old_label_left [(j*r) + i] = border_left [(j*r) + i].label;
		    }
		}
		
/* Prefetch Right Border */

		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			border_right  [(j*r) + i].color :=     X[ii+1][jj+j][0][i];
			border_right  [(j*r) + i].label := Label[ii+1][jj+j][0][i];
		    }
		}
		
/* Sort each border list by label, and then move through the list linearly, 
   attaching each pair that are of the same label */
		
		fastsort_elem(border_left, length);

		sync();

		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			color_right    [(j*r) + i] = border_right[(j*r) + i].color;
			old_label_right[(j*r) + i] = border_right[(j*r) + i].label;
		    }
		}
	    }	    

	    if (shadow_manager) {

/* Get border from processors, and set up label arrays */
		
		for (j=0 ; j < (1<<wIt) ; j++) {
		    for (i=0 ; i<r ; i++) {
			locShadowBorder[(j*r) + i].color :=     X[ii][jj+j][0][i];
			locShadowBorder[(j*r) + i].label := Label[ii][jj+j][0][i];
		    }

		    for (i=0 ; i<r ; i++) 
			locShadowBorder[(j*r) + i].pos   = (j*r) + i;

		    sync();
		}
		
		fastsort_elem(locShadowBorder,length);
	    }
	    
	    barrier();

/* Group Manager should retrieve results from Shadow manager */

	    if (group_manager) {

		bulk_get(border_right, &(ShadowBorder[ii+1][jj][0]), 
			     length*sizeof(elem_t)); 

/* Create Adjacency List problem */

		borderLen = 2*length;

		/* Each vertex can have at most 5 edges incident to it */
		/* There are at most 2*(5*borderLen - 4) edges in the graph */
		edgeCount = 0;
		
/* Add edges */
/* First, add edges stringing along pixels of the same label down each of 
   the borders */
/* Then add edges representing links between the two borders */
		
		for (i=1 ; i<length ; i++) {
		    if ((border_left[i-1].label == border_left[i].label)
			&& (border_left[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				border_left[i-1].pos,
				border_left[i].pos);
		}

		sync();
		
		for (i=1 ; i<length ; i++) {
		    if ((border_right[i-1].label == border_right[i].label)
			&& (border_right[i].label > 0))
			addAuxEdge(auxGraph, &edgeCount,
				length + border_right[i-1].pos,
				length + border_right[i].pos);
		}
		
/* Scan down the left border and add edges wherever it meets the right bord */
		
		for (i=0 ; i<length ; i++)
		    if (color_left[i]==color_right[i]) 
			    addAuxEdge(auxGraph, &edgeCount, i, length + i );
		
		if (edgeCount > 0) {
		/* Sort the edges, then create an auxIndex */
		    fastsort_edge(auxGraph, edgeCount);

		    for (i=0 ; i<borderLen ; i++)
			auxIndex[i] = NIL;
		
		    auxIndex[auxGraph[0].a] = 0;

		    if (edgeCount > 1)
			for (i=1 ; i<edgeCount ; i++) 
			    if (auxGraph[i].a > auxGraph[i-1].a)
				auxIndex[auxGraph[i].a] = i;

		    for (i=0 ; i<length ; i++) {
			input_labelling[i]        = old_label_left[i];
			input_labelling[i+length] = old_label_right[i];
			ccOut[i]        = input_labelling[i];
			ccOut[i+length] = input_labelling[i+length];
		    }
		
		    BFS_auxGraph(auxGraph, edgeCount, borderLen, auxIndex,
				 ccOut, input_labelling);
		}

/* Compare old_label's with new_label's and propogate the differences */
			        
		chSize = 0;
		
		if (edgeCount > 0) {

		    for (i=0 ; i<length ; i++) {
			new_label_left[i]  = ccOut[i];
			new_label_right[i] = ccOut[length+i];
		    }

		    for (i=0 ; i<length ; i++) {
			if (old_label_left[i] != new_label_left[i]) {
			    changeList[chSize].alpha = old_label_left[i];
			    changeList[chSize].beta  = new_label_left[i];
			    chSize++;
			}
			if (old_label_right[i] != new_label_right[i]) {
			    changeList[chSize].alpha = old_label_right[i];
			    changeList[chSize].beta  = new_label_right[i];
			    chSize++;
			}
		    }
		
		    fastsort_chPair(changeList, chSize);

		    chSize = fill_chList(locChList, changeList, chSize);

		}

		changeSize [ii][jj] = chSize;
		grpManV = ii;
		grpManW = jj;

	    }
	    else {  /* Figure out who my group_manager is */
		/* ii with the last vIt bits set to 0111.1 (a zero and vIt-1 1's) */
		grpManV = (((ii >> vIt) << vIt) | ((1<<(vIt-1))-1));
		/* jj with the last wIt bits set to 0 */
		grpManW = (jj >> wIt) << wIt;
	    }

	    barrier();

#if 0
/* Prefetch Change List from Group Manager */
	    
	    if (!group_manager) {

		chSize = changeSize[grpManV][grpManW];

/* Fix this prefetch later to be a better broadcast algorithm */

		if (chSize > 0) {
		    bulk_get(locChList, &(chList[grpManV][grpManW][0]), 
			     chSize*sizeof(chPair_t)); 
		    sync();
		}
	    }
#endif
#if 0
	    all_mesh_bcast(v, w, (q*v) << 1,
			   vIt, wIt,		       
			   ii, jj, grpManV, grpManW,
			   chList, changeSize);
	    chSize = changeSize[ii][jj];
#endif
#if 1
	    orgV    = (ii>>vIt) << vIt;
	    orgW    = (jj>>wIt) << wIt;
	    g_grey  = greyEncoder(vIt, wIt, ii-orgV, jj-orgW);
	    g_GM    = greyEncoder(vIt, wIt, grpManV-orgV, grpManW-orgW);
	    g_label = g_grey ^ g_GM;
	    j       = vIt + wIt;

	    for (i=0 ; i<j ; i++) {
		if (g_label < (1<<i)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) );
		    changeSize[send_i+orgV][send_j+orgW]
			= changeSize[ii][jj];
		}
		barrier();
	    }

	    chSize = changeSize[ii][jj];

	    for (i=0 ; i<j ; i++) {
		if ((g_label < (1<<i)) && (chSize > 0)) {
		    greyDecoder(vIt, wIt, &send_i, &send_j, g_grey^(1<<i) ); 
		    bulk_write(&(chList[send_i+orgV][send_j+orgW][0]),
			       tolocal(chList[ii][jj]),
			       chSize*sizeof(chPair_t));
		}
		barrier();
	    }	    
#endif

/* Make changes */

/* We have a sorted array of changes of size chSize,
           a data structure of tile border hooks.
   Recolor only the border of each tile by binary searching the
   label change for each hook's label, and then updating the label
   at each of the addresses in the addr. list. */
	    
	    for (i=0 ; i<hooks ; i++) {          /* For Each Hook */
		beta = chLookup(locChList, chSize, tileHookArr[i].label);
		if (tileHookArr[i].label != beta) {
		    for (j=0 ; j < tileHookArr[i].size ; j++)
			ccLabel[(tileHookArr[i].addr)[j]] = beta;
		    tileHookArr[i].label = beta;
		}
	    }

	    sprintf(buf,"Merge Phase %d",vIt+wIt);
	    all_stop_timer(buf);

	    barrier();

	}
    }
    
/* We must relabel the interior of our tile by searching through
   tileHookOrig, and if any labels have changed, perform a depth-first search
   and recolor connected pixels */

    all_start_timer();
          
    for (i=0 ; i<q ; i++)
	for (j=0 ; j<r ; j++) 
	    mark[i][j] = (locX[i][j] ? WHITE : BLACK);
	    
    for (i=0 ; i<hooks ; i++)
	if (tileHookOrig[i].alpha != tileHookArr[i].label)
	    /* Update component */
	    /* This will check all adjacent pixels */

            BFS_binary_4_recolor(q,r,
			       tileHookOrig[i].beta / r,
			       tileHookOrig[i].beta % r,
			       locLabel, locX, mark, tileHookArr[i].label);

    all_stop_timer("Update Interior Labels");

    destroy_tileHook(tileHookArr,hooks);    
    free(tileHookArr);
    free(tileHookOrig);
    free(tempHookArr);

    free(auxIndex);
    free(auxGraph);
    free(input_labelling);
    free(ccOut);
    free(border_right);
    free(border_left);
    free(color_right);
    free(color_left);
    free(old_label_right);
    free(new_label_right);
    free(old_label_left);
    free(new_label_left);

    barrier();

    secs = get_seconds() - secs;  /* Finish Timing */

/*****************************************************/
/* Results                                           */
/*****************************************************/

/* DONE - output results of connected components */

 
    on_one {
	fprintf(outfile,"4-Connected Components Time: %f  ",secs);
	fprintf(outfile,"for n x n = %d image on %d processors.\n", v*w*q*r, PROCS);

	if (RESULTS) {
	    fprintf(outfile,"Original Binary Image:\n");
	    for (i=0 ; i<v ; i++) 
		for (ii=0 ; ii < q ; ii++) {
		    for (j=0 ; j<w ; j++)
			for (jj=0 ; jj < r ; jj++)
			    fprintf(outfile,"%1d ",X[i][j][ii][jj]);
		    fprintf(outfile,"\n");
		}
	    fprintf(outfile,"\n");

	    fprintf(outfile,"4-Connected Components:\n");
	    for (i=0 ; i<v ; i++) 
		for (ii=0 ; ii < q ; ii++) {
		    for (j=0 ; j<w ; j++)
			for (jj=0 ; jj < r ; jj++)
			    fprintf(outfile,"%3d ",Label[i][j][ii][jj]);
		    fprintf(outfile,"\n");
		}
	    fprintf(outfile,"\n");
	}

	if (OUTPUT_IMAGE) {
	    fprintf(outfile,"Writing Output Image... ");
	    fflush(outfile);
	    image_print_cc("4-connected components of binary image",
			   "image4ccbin",
			   v,w,q,r,k,Label);
	    fprintf(outfile,"Done.\n");
	    fflush(outfile);
	}

	all_print_timer(outfile);
	    
    }
    
}
#else
{}
#endif

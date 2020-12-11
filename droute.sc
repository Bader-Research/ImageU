/*                                                                    tab:8
 *
 * droute.sc - Deterministic Routing of H-relations
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
 * Filename:            droute.sc
 * History:
 */

#include "droute.h"
#include "sorting.h"

#define CHECK_HREL  0
#define TDEBUG      0

typedef struct record {
    struct record *next; 
    intpair_t *bin;
} record_t;


void all_check_hrel(int q, int* rout, int M, int h) 
{
    int i, k;

    int *spread cnt;
    int *mycnt;

    barrier();
	
    cnt = all_spread_malloc(PROCS, PROCS*sizeof(int));
    assert_spread_malloc(cnt);

    mycnt = (int*)(cnt + MYPROC);

    for (i=0 ; i<PROCS ; i++)
	mycnt[i] = 0;

    for (i=0 ; i<q ; i++) 
	if ((rout[i] < 0) || (rout[i] > PROCS))
	    fprintf(outfile,
		    "ERROR: (%2d), i:%4d  rout[i]:%4d M:%d  h:%1d\n",
		    MYPROC,i,rout[i],M,h);
	else 
	    mycnt[rout[i]]++;

#if 0
    if ((MYPROC < PROCS-1)||(MYPROC==PROCS-1)) {
	...
    }
    else {
      on_proc(PROCS-1) {
	for (i=0 ; i<h*M ; i++) 
	  if (rout[i] > PROCS)
	    fprintf(outfile,
		    "ERROR: (%2d), i:%4d  rout[i]:%4d M:%d  h:%1d\n",
		    MYPROC,i,rout[i],M,h);
	  else
	    if (rout[i] > 0)
	      mycnt[rout[i]]++;
	    
      }
    }
#endif

    barrier();

    on_one {
	fprintf(outfile,"\t");
	for (k=0 ; k<PROCS ; k++)
	    fprintf(outfile,"%3d ",k);
	fprintf(outfile,"\n\n");

	for (i=0 ; i<PROCS ; i++) {
	    fprintf(outfile,"%2d : \t",i);
	    for (k=0 ; k<PROCS ; k++)
		fprintf(outfile,"%3d ",*((int *global)(cnt+i) + k));
	    fprintf(outfile,"\n");
	}
	fflush(outfile);
    }
    barrier();

    all_spread_free(cnt);

}


/****************************************************/
inline int all_route_h_rel_det_safe(int q, int h,
				    int Keys[PROCS]::[q],
				    int Addr[PROCS]::[q],
				    int Routed[PROCS]::[h*q])
/****************************************************/
/* q (elems/proc) must be a multiple of PROCS */
{
    register int
	i, j, k,
	l, d, m,
	get_proc;

    int b1, b2;

    int *lKeys,
	*lAddr,
	*lRouted,
	tag[PROCS];

    intpair_t
	*lC, *lD;

    intpair_t C[PROCS]::[h*q + PROCS*(PROCS+1)];
    intpair_t D[PROCS]::[h*q + PROCS*(PROCS+1)];

    lKeys = tolocal(Keys[MYPROC]);
    lAddr = tolocal(Addr[MYPROC]);
    lRouted = tolocal(Routed[MYPROC]);

    lC = tolocal(C[MYPROC]);
    lD = tolocal(D[MYPROC]);
    
    barrier();

/* Init First Phase */

    b1 = q/PROCS + PROCS + 1;

    for (i=0 ; i<PROCS ; i++) {
	lC[i*b1].a = 0;         /* zero slot counters */
	tag[i]     = i;         /* stagger address tags */
    }

/* Place data into arrays */

    for (i=0 ; i<q ; i++) {
	l = lAddr[i] / (h*q); /* get key -> proc address */
	d = tag[l];           /* get current slot for address */
	j = d*b1;             /* slot base address */
	m = ++lC[j].a;        /* get slot offset + 1 */
	if (m >= b1)
	    fprintf(stderr,"ERROR: Slot filled (1)\n");

	lC[j + m].a = lKeys[i]; /* place Key data */
	lC[j + m].b = lAddr[i]; /* place Key addr */
	tag[l] = (d + 1) % PROCS;  /* next slot */
    }

/* Transpose C -> D */    
	
    barrier();		

    for (i=0 ; i<PROCS ; i++) {
        get_proc = (MYPROC + i) % PROCS; 
        bulk_read(&lD[b1*get_proc],
		  &(C[get_proc][b1*MYPROC]),
		  b1*sizeof(intpair_t));
	barrier();
    }

/* Init Second Phase */

    b2 = h*(q/PROCS) + PROCS;

    for (i=0 ; i<PROCS ; i++)
	lC[i*b2].a = 0;         /* zero slot counters */

/* Place data into arrays */

    for (i=0 ; i<PROCS ; i++) {
	j = i*b1;
	l = 1;
	while (l <= lD[j].a) {
	    d = (lD[j+l].b) / (h*q);
	    k = d*b2;
	    m = ++lC[k].a;
	    if (m >= b2)
		fprintf(stderr,"ERROR: Slot filled (2)\n");

	    lC[k+m].a = lD[j+l].a;
	    lC[k+m].b = lD[j+l].b;
	    l++;
	}
    }

/* Transpose C -> D */    
	
    barrier();		

    for (i=0 ; i<PROCS ; i++) {
        get_proc = (MYPROC + i) % PROCS; 
        bulk_read(&lD[b2*get_proc],
		  &(C[get_proc][b2*MYPROC]),
		  b2*sizeof(intpair_t));
	barrier();
    }

/* Unpack D */

    k = 0;
    for (i=0 ; i<PROCS ; i++) {
	j = i*b2;           
	l = 1;
	while (l <= lD[j].a) {
	    d = lD[j + l].b % (h*q);
	    lRouted[d] = lD[j + l].a;
	    k++;
	    l++;
	}
    }

    return(k);
}

/****************************************************/
int all_route_h_rel_det(int q, int h,
			int Keys[PROCS]::[q],
			int Addr[PROCS]::[q],
			int Routed[PROCS]::[h*q])
/****************************************************/
/* q (elems/proc) must be a multiple of PROCS */
{
    register int
	i, j, k,
	l, d, m,
	get_proc;

    int b1, b2;
    int LOGQ;

    int *lKeys,
	*lAddr,
	*lRouted,
	tag[PROCS],
	*tptr;

    intpair_t
	*lC, *lD,
	*t1ptr, *t2ptr;

    intpair_t *spread C;
    intpair_t *spread D;

    C = all_spread_malloc(PROCS,(h*q + PROCS*(PROCS+1))*sizeof(intpair_t));
    assert_spread_malloc(C);
    D = all_spread_malloc(PROCS,(h*q + PROCS*(PROCS+1))*sizeof(intpair_t));
    assert_spread_malloc(D);
	
    lKeys = tolocal(Keys[MYPROC]);
    lAddr = tolocal(Addr[MYPROC]);
    lRouted = tolocal(Routed[MYPROC]);

    lC = (intpair_t *)(C+MYPROC);
    lD = (intpair_t *)(D+MYPROC);
    
    barrier();

/* Init First Phase */

    LOGQ = (int)log2((double)(h*q));
    
    b1 = DIVPROCS(q) + (PROCS/2) + 1;

    j = -b1;
    for (i=0 ; i<PROCS ; i++) {
	(lC + (j+=b1))->a = 0;   /* zero slot counters */
#if 0
	*(tag+i) = i;            /* stagger address tags */
#endif	
	*(tag+i) = (i+MYPROC) % PROCS;    /* stagger address tags */
    }

/* Place data into arrays */

    for (i=0 ; i<q ; i++) {
	d = *(tptr = tag + (*(lAddr+i) >> LOGQ) );         
	j = d*b1;           
	m = ++((lC+j)->a);      
	if (m >= b1)
	    fprintf(stderr,"ERROR: Slot filled (1)\n");
	k = j+m;
	(t1ptr = lC + k)->a = *(lKeys+i); 
	t1ptr->b = *(lAddr+i); 
	*tptr = MODPROCS(d + 1);
    }

/* Transpose C -> D */    
	
    barrier();		

#if (defined(CM5))
    for (i=0 ; i<PROCS ; i++) {
        get_proc = MODPROCS(MYPROC + i);
        bulk_read(lD + (b1*get_proc),
		 (intpair_t *global)(C+get_proc) + (b1*MYPROC),
		  b1*sizeof(intpair_t));
	barrier();
    }
#elif (defined(SP2))
    barrier();
    mpc_index(lC, lD, b1*sizeof(intpair_t), ALLGRP);
#else
    for (i=0 ; i<PROCS ; i++) {
        get_proc = MODPROCS(MYPROC + i);
        bulk_get(lD + (b1*get_proc),
		 (intpair_t *global)(C+get_proc) + (b1*MYPROC),
		 b1*sizeof(intpair_t));
    }
    sync();
    barrier();
#endif

/* Init Second Phase */

#if 0
    b2 = DIVPROCS(h*q) + PROCS;
#endif    
    b2 = DIVPROCS(h*q) + PROCS/2 ;

    j = -b2;
    for (i=0 ; i<PROCS ; i++)
	(lC+(j+=b2))->a = 0;      

/* Place data into arrays */

    for (i=0 ; i<PROCS ; i++) {
	j = i*b1;
	l = 1;
	while (l <= (lD+j)->a) {
	    k = (((lD+j+l)->b) >> LOGQ)*b2;
	    m = ++((t1ptr = lC+k)->a);
	    if (m >= b2)
		fprintf(stderr,"ERROR: Slot filled (2)\n");

	    (t1ptr += m)->a = (t2ptr = lD+j+l)->a;
	    t1ptr->b = t2ptr->b;
	    l++;
	}
    }

/* Transpose C -> D */    

    barrier();		

#if (defined(CM5))
    for (i=0 ; i<PROCS ; i++) {
        get_proc = MODPROCS(MYPROC + i);
        bulk_read(lD + b2*get_proc,
		 (intpair_t *global)(C+get_proc) + (b2*MYPROC),
		  b2*sizeof(intpair_t));
	barrier();
    }
#elif (defined(SP2))
    barrier();
    mpc_index(lC, lD, b2*sizeof(intpair_t), ALLGRP);
#else
    for (i=0 ; i<PROCS ; i++) {
        get_proc = MODPROCS(MYPROC + i);
        bulk_get(lD + b2*get_proc,
		 (intpair_t *global)(C+get_proc) + (b2*MYPROC),
		 b2*sizeof(intpair_t));
    }
    sync();
    barrier();
#endif

/* Unpack D */

    k = 0;
    for (i=0 ; i<PROCS ; i++) {
	j = i*b2;           
	l = 0;
	while (++l <= (lD+j)->a) {
	    d = ((t1ptr = lD+j+l)->b) & ((h*q)-1);
	    *(lRouted+d) = t1ptr->a;
	    k++;
	}
    }
    all_spread_free(D);
    all_spread_free(C);

#if 0
    barrier();
    on_one {
	fprintf(outfile,"hrel1 Output on PE %3d:\n",MYPROC);
	fprintf(outfile,"p = %3d  h = %2d  q = %12d  (h*q: %d j: %d)\n",
		PROCS, h, q, h*q, k);
	for (i=0 ; i<k ; i++)
	    fprintf(outfile,"elem[%12d]: %12d\n",i,
		    *((int *)(Routed+MYPROC)+i));
	fflush(outfile);
    }
    barrier();
#endif

#if CHECK_HREL
    all_check_hrel(k, (int*)(Routed+MYPROC), q, h);
#endif

    return(k);
}

int all_route_h_rel_det2(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed) {
    register int
	i, 
	LOGQ,
	entries,
	values,
	dest;

    intpair_t
	*A_ptr,
	*buck_ptr,
	*buff_ptr,
	*B_ptr,
	*bin_ptr,
	*out_ptr,
	*t_ptr;
    int *in_ptr,
	*rank_ptr,
	*t1_ptr,
	*myout_ptr;

    int *L_out_counts;

    intpair_t *global dest_ptr;
    intpair_t *global G_out_ptr;
    int *global G_out_counts;

    intpair_t *spread temp_ptr;
    intpair_t *spread rec_ptr;
    intpair_t *spread unload_ptr;

    intpair_t **B_bin_ptr;

    record_t **loc_ptr; 
    record_t *record_ptr;

    intpair_t **B_bin;

    record_t *Info;
    record_t **Loc;

    intpair_t *spread A;
    intpair_t *spread B;
    intpair_t *spread C;
    intpair_t *spread Buffer1;
    intpair_t *spread Buffer2;

    int *spread Out_Counts;

    int buck_size = (int)(q/PROCS) + 2*PROCS + 1;
    int bin_size = h*buck_size + 1;

    B_bin = (intpair_t **)malloc(PROCS*sizeof(intpair_t *));
    assert_malloc(B_bin);
    Info = (record_t *)malloc(PROCS*sizeof(record_t));
    assert_malloc(Info);
    Loc = (record_t **)malloc(PROCS*sizeof(record_t *));
    assert_malloc(Loc);
    A = all_spread_malloc(PROCS,PROCS*buck_size*sizeof(intpair_t));
    assert_spread_malloc(A);
    B = all_spread_malloc(PROCS,PROCS*bin_size*sizeof(intpair_t));
    assert_spread_malloc(B);
    Buffer1 = all_spread_malloc(PROCS,bin_size*sizeof(intpair_t));
    assert_spread_malloc(Buffer1);
    Buffer2 = all_spread_malloc(PROCS,bin_size*sizeof(intpair_t));
    assert_spread_malloc(Buffer2);
    Out_Counts = all_spread_malloc(PROCS,PROCS*sizeof(int));
    assert_spread_malloc(Out_Counts);

/* Distribute the contents of In into the bins of A */

    /* First initialize things: */

    LOGQ = (int)log2((double)(h*q));
    in_ptr = (int *) (Keys+MYPROC) - 1;
    rank_ptr = (int *) (Addr+MYPROC) - 1;

    A_ptr = (intpair_t *)(A+MYPROC) - buck_size;
    
    for (i=0,record_ptr=Info;i<PROCS;i++) {
	*(Loc + i) = record_ptr;
	record_ptr->next = record_ptr + 1;
	(record_ptr++)->bin = (A_ptr+=buck_size);
    }

    (Info +PROCS - 1)->next = Info;

    /*Now, distribute the contents of Keys: */

    t1_ptr = rank_ptr + q;
    while (++rank_ptr <= t1_ptr) {
	(++((*(loc_ptr = Loc + (*rank_ptr >> LOGQ)))->bin))->b =
	    *(++in_ptr);
	((*loc_ptr)->bin)->a = *rank_ptr;
	*loc_ptr = (*loc_ptr)->next;
    }

    barrier();

/* Transpose the slices of A, sorting the incoming data into the appropriate
   bins of B: */


    record_ptr = Info;
    rec_ptr = Buffer1;     /*The buffer currently being broadcast to */
    unload_ptr = Buffer2;  /*The buffer currently being unloaded */
    A_ptr = (intpair_t *)(A+MYPROC);
    dest = MODPROCS(MYPROC+1);
    buck_ptr = A_ptr + dest*buck_size;
    buck_ptr->a = entries = (record_ptr+dest)->bin - buck_ptr;
    dest_ptr = (intpair_t *global)(rec_ptr+dest); 

    bulk_store(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 

    barrier();

    B_ptr = (intpair_t *)(B+MYPROC) - bin_size;
    B_bin_ptr = B_bin;

    for (i=0;i<PROCS;i++)
	*(B_bin_ptr++) = (B_ptr+=bin_size);

    buck_ptr = A_ptr + MYPROC*buck_size; 
    entries = (record_ptr+MYPROC)->bin - buck_ptr;
    t_ptr = buck_ptr + entries;
    while (++buck_ptr <= t_ptr)
	*((*(B_bin + ((buck_ptr->a) >> LOGQ)))++) = *buck_ptr;

    /* Swap buffers: */

    temp_ptr = unload_ptr;
    unload_ptr = rec_ptr;
    rec_ptr = temp_ptr;

    barrier();
    
    all_store_sync(); 

    for (i=2 ; i<PROCS ; i++) {
	dest = MODPROCS(MYPROC+i);
	buck_ptr = A_ptr + dest*buck_size;
	buck_ptr->a = (entries = (record_ptr+dest)->bin - buck_ptr);  
	dest_ptr = (intpair_t *global)(rec_ptr+dest); 

	bulk_store(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 
	barrier();
	buff_ptr = (intpair_t *)(unload_ptr + MYPROC);
	entries = buff_ptr->a;
	t_ptr = buff_ptr + entries;
	while (++buff_ptr <= t_ptr)
	    *((*(B_bin + ((buff_ptr->a) >> LOGQ)))++) = *buff_ptr;

        /* Swap buffers: */

	temp_ptr = unload_ptr;
	unload_ptr = rec_ptr;
	rec_ptr = temp_ptr;
	all_store_sync(); 
    }

    buff_ptr = (intpair_t *)(unload_ptr + MYPROC);
    entries = buff_ptr->a;
    t_ptr = buff_ptr + entries;
    while ((++buff_ptr) <= t_ptr)
	*((*(B_bin + ((buff_ptr->a) >> LOGQ)))++) = *buff_ptr;

    barrier();

/* Transpose the slices of B into C */

    C = A;

    L_out_counts = (int *)(Out_Counts + MYPROC);
    G_out_counts = (int *global)(Out_Counts + 
				 (MODPROCS(MYPROC + PROCS - 1)));

    B_ptr = (intpair_t *)(B + MYPROC);
    out_ptr = (intpair_t *)(C + MYPROC);
    bin_ptr = B_ptr + MYPROC*bin_size;
    values = *(B_bin + MYPROC) - bin_ptr;
    *(G_out_counts++) :- values;    

    for (i=0;i<values;i++)
	*(out_ptr + i) = *(bin_ptr + i);

    all_store_sync();

/* Now route the remaining elements */

    for (i=1;i<PROCS;i++) {
	values = *(L_out_counts++); 
	dest = MODPROCS(MYPROC + i);
	bin_ptr = B_ptr + dest*bin_size;
	entries = *(B_bin + dest) - bin_ptr;
	*(G_out_counts++) :- values + entries;   
	G_out_ptr = (intpair_t *global)(C + dest) + values;
	bulk_store(G_out_ptr,bin_ptr,entries*sizeof(intpair_t));
	all_store_sync(); 
    }

    values = *L_out_counts; 

    myout_ptr = (int*)(Routed+MYPROC);
    out_ptr   = (intpair_t *)(C+MYPROC) - 1;

    t_ptr = out_ptr + values;
    while (++out_ptr <= t_ptr)
	*(myout_ptr + ((out_ptr->a) & ((h*q)-1))) = out_ptr->b;

    barrier();

    free(B_bin);
    free(Loc);
    free(Info);

    all_spread_free(A);
    all_spread_free(B);
    all_spread_free(Buffer1);
    all_spread_free(Buffer2);
    all_spread_free(Out_Counts);

#if CHECK_HREL
    all_check_hrel(values, (int*)(Routed+MYPROC), q, h);
#endif

    return values;
}


int all_route_h_rel_det3(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed) {

    int buck_size = DIVPROCS(q)   + PROCS/2 + 2;
    int hq        = h*q;
    int bin_size  = DIVPROCS(hq) + PROCS/2 + 2;

    register int
	i,
	t,
	LOGQ,
	entries,
	dest,
	hqmask;

    int pbins,
	values;
    
    intpair_t
	*A_ptr,
	*buck_ptr,
	*buff_ptr,
	*B_ptr,
	*bin_ptr,
	*s_ptr,
	*t_ptr,
	*C_ptr,
	*D_ptr,
	*p1_ptr;
    
    int	*in_ptr,
	*rank_ptr,
	*t1_ptr,
	*myout_ptr;


    intpair_t *global dest_ptr;
    
    intpair_t *spread A;
    intpair_t *spread B; 
    intpair_t *spread C; 
    intpair_t *spread D; 
    
    intpair_t **C_bin_ptr;
    
    record_t **loc_ptr; 
    record_t *record_ptr;
    
    intpair_t **C_bin;
    record_t *Info;
    record_t **Loc;
    
    intpair_t *spread Buffer1;
    intpair_t *spread Buffer2;

#if TDEBUG
    double secs;
#endif

#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    pbins = PROCS * bin_size;

    C_bin =  (intpair_t **)malloc(PROCS*sizeof(intpair_t *));
    assert_malloc(C_bin);
    Info = (record_t *)malloc(PROCS*sizeof(record_t));
    assert_malloc(Info);
    Loc = (record_t **)malloc(PROCS*sizeof(record_t *));
    assert_malloc(Loc);
    
    Buffer1 = all_spread_malloc(PROCS,pbins*sizeof(intpair_t));
    assert_spread_malloc(Buffer1);
    Buffer2 = all_spread_malloc(PROCS,pbins*sizeof(intpair_t));
    assert_spread_malloc(Buffer2);
    
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) malloc \t\t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

/* Distribute the contents of In into the bins of A */
    
    /* First initialize things: */

#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    LOGQ = (int)log2((double)(h*q));
    A = Buffer1;
    in_ptr = (int *)(Keys + MYPROC) - 1;
    rank_ptr = (int *)(Addr + MYPROC) - 1;
    A_ptr =  (intpair_t *)(A + MYPROC) - buck_size;
    
    for (i=0,record_ptr=Info ; i<PROCS ; i++) {
#if 0
	*(Loc + i) = record_ptr;
#endif	
	*(Loc + MODPROCS(i+MYPROC)) = record_ptr;
	record_ptr->next = record_ptr + 1;
	(record_ptr++)->bin = (A_ptr+=buck_size);
    }
    
    (Info + PROCS - 1)->next = Info;
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) init_ptrs \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

    /*Now, distribute the contents of Keys: */
    
#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    t1_ptr = rank_ptr + q;
    while (++rank_ptr <= t1_ptr) {
	(++((*(loc_ptr = Loc + (*rank_ptr >> LOGQ)))->bin))->b =
	    *(++in_ptr);
	((*loc_ptr)->bin)->a = *rank_ptr;
	*loc_ptr = (*loc_ptr)->next;
    }
    
    barrier();
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) distibute \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

/* Transpose the slices of A into  B: */
  
#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    B = Buffer2;
    record_ptr = Info;    
    A_ptr = (intpair_t *)(A + MYPROC); 
    B_ptr = (intpair_t *)(B + MYPROC); 
    buck_ptr = A_ptr + (t = MYPROC*buck_size);
    buff_ptr = B_ptr + t;
    buff_ptr->a = entries = (record_ptr+MYPROC)->bin - buck_ptr;
    bcopy(buck_ptr+1, buff_ptr+1, entries*sizeof(intpair_t));
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) init_trans A \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

#if (defined(T3D))
    for (i=1;i<PROCS;i++) {
	dest = MYPROC^i; 
	buck_ptr = A_ptr + dest*buck_size;
	buck_ptr->a = entries = (record_ptr+dest)->bin - buck_ptr;  
	dest_ptr = (intpair_t *global)(B+dest) + t;   
	bulk_write(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 
	barrier();
    }
#elif (defined(MEIKO))
    for (i=1;i<PROCS;i++) {
	dest = MYPROC^i; 
	buck_ptr = A_ptr + dest*buck_size;
	buck_ptr->a = entries = (record_ptr+dest)->bin - buck_ptr;  
	dest_ptr = (intpair_t *global)(B+dest) + t;   
	bulk_write(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 
    }
#elif (defined(SP2))
    for (i=1;i<PROCS;i++) {
	dest = MODPROCS(MYPROC+i); 
	buck_ptr = A_ptr + dest*buck_size;
	buck_ptr->a = entries = (record_ptr+dest)->bin - buck_ptr;  
	dest_ptr = (intpair_t *global)(B+dest) + t;   
	bulk_write(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 
    }
    barrier();
#else
    for (i=1;i<PROCS;i++) {
	dest = MODPROCS(MYPROC+i);
	buck_ptr = A_ptr + dest*buck_size;
	buck_ptr->a = entries = (record_ptr+dest)->bin - buck_ptr;  
	dest_ptr = (intpair_t *global)(B+dest) + t;   
	bulk_store(dest_ptr,buck_ptr,(entries+1)*sizeof(intpair_t)); 
    }
    
    all_store_sync();
    barrier();
#endif
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) trans A \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

/* Sort the Contents of B into the appropriate bins of C: */
    
#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    C = Buffer1;
    s_ptr = (C_ptr = (intpair_t *)(C+MYPROC)) - bin_size ;
    C_bin_ptr = C_bin;
    
    for (i=0;i<PROCS;i++)
	*(C_bin_ptr + i) = (s_ptr+=bin_size);
    
    t_ptr = B_ptr + (PROCS*buck_size);
    s_ptr = B_ptr - buck_size;
    while ((buff_ptr = (s_ptr += buck_size)) < t_ptr) {
	p1_ptr = buff_ptr + (buff_ptr->a);
	while (++buff_ptr <= p1_ptr)
	    *(++(*(C_bin + ((buff_ptr->a)>>LOGQ)))) = *buff_ptr;
    }
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) copy B->C \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

/* Transpose the slices of C into D */ 
    
#if TDEBUG
    barrier();
    secs = get_seconds();
#endif
    D = Buffer2;
    D_ptr = (intpair_t *)(D + MYPROC); 
    bin_ptr = C_ptr + (t = MYPROC*bin_size);
    buff_ptr = D_ptr + t;
    buff_ptr->a = entries = *(C_bin + MYPROC) - bin_ptr;
    bcopy(bin_ptr+1, buff_ptr+1, entries*sizeof(intpair_t));

#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) init_trans B \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

#if (defined(T3D))
    for (i=1;i<PROCS;i++) {
	dest = MYPROC^i;
	bin_ptr = C_ptr + dest*bin_size;
	bin_ptr->a = (entries = *(C_bin + dest) - bin_ptr);  
	dest_ptr = (intpair_t *global)(D + dest) + t;   
	bulk_write(dest_ptr,bin_ptr,(entries+1)*sizeof(intpair_t)); 
	barrier();
    }
#elif (defined(MEIKO))
    for (i=1;i<PROCS;i++) {
	dest = MYPROC^i;
	bin_ptr = C_ptr + dest*bin_size;
	bin_ptr->a = (entries = *(C_bin + dest) - bin_ptr);  
	dest_ptr = (intpair_t *global)(D + dest) + t;   
	bulk_write(dest_ptr,bin_ptr,(entries+1)*sizeof(intpair_t)); 
    }
#elif (defined(SP2))
    for (i=1;i<PROCS;i++) {
	dest = MODPROCS(MYPROC+i); 
	bin_ptr = C_ptr + dest*bin_size;
	bin_ptr->a = (entries = *(C_bin + dest) - bin_ptr);  
	dest_ptr = (intpair_t *global)(D + dest) + t;   
	bulk_write(dest_ptr,bin_ptr,(entries+1)*sizeof(intpair_t)); 
    }
    barrier();
#else
    for (i=1;i<PROCS;i++) {
	dest = MODPROCS(MYPROC+i);
	bin_ptr = C_ptr + dest*bin_size;
	bin_ptr->a = (entries = *(C_bin + dest) - bin_ptr);  
	dest_ptr = (intpair_t *global)(D + dest) + t;   
	bulk_store(dest_ptr,bin_ptr,(entries+1)*sizeof(intpair_t)); 
    }
    
    all_store_sync();
    barrier();
#endif
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) trans B \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

/* Move the Contents of D into Routed */
    
#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    values = 0;
    hqmask = hq - 1;
    myout_ptr = (int *)(Routed + MYPROC);  
    t_ptr = D_ptr + pbins;
    s_ptr = D_ptr - bin_size;
    while ((buff_ptr = (s_ptr += bin_size)) < t_ptr) {
	values += (entries = (buff_ptr->a));
	p1_ptr = buff_ptr + entries; 
	while (++buff_ptr <= p1_ptr)
#if (defined(SP2))
	    {
		i = (buff_ptr->a) & hqmask;
		if (i < hq) 
		    *(myout_ptr + i) = buff_ptr->b;
		else {
		    fprintf(stderr,"ERROR: (%3d)  h:%d q:%d i:%d\n",
			    MYPROC, h, q, i);
		    fflush(stderr);
		}
	    }
#else
	*(myout_ptr + ((buff_ptr->a) & hqmask)) = buff_ptr->b;
#endif
    }
    
    
    barrier();
    
#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) D->Routed \t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif

#if TDEBUG
    barrier();
    secs = get_seconds();
#endif

    free(C_bin);
    free(Loc);
    free(Info);

    all_spread_free(Buffer1);
    all_spread_free(Buffer2);

#if TDEBUG
    barrier();
    secs = get_seconds() - secs;
    on_one {
	fprintf(outfile,"P(%3d) q(%7d) h(%2d) free \t\t %9.6f \n",PROCS,q,h,secs);
	fflush(outfile);
    }
    barrier();
#endif
    
#if CHECK_HREL
    all_check_hrel(values, (int*)(Routed+MYPROC), q, h);
#endif

    return values;
    
}

#if (defined(SP2))
int all_route_h_rel_det4(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed) {

    int buck_size = DIVPROCS(q) + PROCS/2 + 1;
    int bin_size = DIVPROCS(h*q) + PROCS/2 ;
    
    register int
	i,
	LOGQ,
	entries,
	values;
    
    intpair_t
	*A_ptr,
	*buck_ptr,
	*buff_ptr,
	*B_ptr,
	*bin_ptr,
	*s_ptr,
	*t_ptr,
	*C_ptr,
	*D_ptr,
	*p1_ptr;
    
    int	*in_ptr,
	*rank_ptr,
	*t1_ptr,
	*myout_ptr;
    
    intpair_t *spread A;
    intpair_t *spread B; 
    intpair_t *spread C; 
    intpair_t *spread D; 
    
    intpair_t **C_bin_ptr;
    
    record_t **loc_ptr; 
    record_t *record_ptr;
    
    intpair_t **C_bin;
    record_t *Info;
    record_t **Loc;
    
    intpair_t *spread Buffer1;
    intpair_t *spread Buffer2;

    C_bin = (intpair_t **)malloc(PROCS*sizeof(intpair_t *));
    assert_malloc(C_bin);
    Info = (record_t *)malloc(PROCS*sizeof(record_t));
    assert_malloc(Info);
    Loc = (record_t **)malloc(PROCS*sizeof(record_t *));
    assert_malloc(Loc);
    
    Buffer1 = all_spread_malloc(PROCS,PROCS*bin_size*sizeof(intpair_t));
    assert_spread_malloc(Buffer1);
    Buffer2 = all_spread_malloc(PROCS,PROCS*bin_size*sizeof(intpair_t));
    assert_spread_malloc(Buffer2);
    
/* Distribute the contents of In into the bins of A */
    
    /* First initialize things: */
    
    LOGQ = (int)log2((double)(h*q));
    A = Buffer1;
    in_ptr = (int *)(Keys + MYPROC) - 1;
    rank_ptr = (int *)(Addr + MYPROC) - 1;
    A_ptr = (intpair_t *)(A + MYPROC) - buck_size;
    
    for (i=0,record_ptr=Info;i<PROCS;i++) {
	*(Loc + MODPROCS(i + MYPROC)) = record_ptr;
	record_ptr->next = record_ptr + 1;
	(record_ptr++)->bin = (A_ptr+=buck_size);
    }
    
    (Info + PROCS - 1)->next = Info;
    
    /*Now, distribute the contents of Keys: */
    
    t1_ptr = rank_ptr + q;
    while (++rank_ptr <= t1_ptr) {
	(++((*(loc_ptr = Loc + (*rank_ptr >> LOGQ)))->bin))->b =
	    *(++in_ptr);
	((*loc_ptr)->bin)->a = *rank_ptr;
	*loc_ptr = (*loc_ptr)->next;
    }
    
    barrier();
    
/* Transpose the slices of A into  B: */
    
    B = Buffer2;
    record_ptr = Info;    
    A_ptr = (intpair_t *)(A + MYPROC); 
    B_ptr = (intpair_t *)(B + MYPROC); 
    buck_ptr = A_ptr;
    for (i=0 ; i<PROCS ; i++,(buck_ptr += buck_size))
	buck_ptr->a = (record_ptr + i)->bin - buck_ptr;
    
#if 1
    barrier();
#endif
    mpc_index(A_ptr,B_ptr,buck_size*(sizeof(intpair_t)),ALLGRP);
    
    
/* Sort the Contents of B into the appropriate bins of C: */
    
    C = Buffer1;
    s_ptr = (C_ptr = (intpair_t *)(C + MYPROC)) - bin_size ;
    C_bin_ptr = C_bin;
    
    for (i=0;i<PROCS;i++)
	*(C_bin_ptr + i) = (s_ptr+=bin_size);
    
    t_ptr = B_ptr + PROCS*buck_size;
    s_ptr = B_ptr - buck_size;
    while ((buff_ptr = (s_ptr += buck_size)) < t_ptr) {
	p1_ptr = buff_ptr + (buff_ptr->a);
	while (++buff_ptr <= p1_ptr)
	    *(++(*(C_bin + ((buff_ptr->a)>>LOGQ)))) = *buff_ptr;
	}
    
    barrier();
    
/* Transpose the slices of C into D */ 
    
    D = Buffer2;
    D_ptr = (intpair_t *)(D + MYPROC); 
    bin_ptr = C_ptr;
    for (i=0;i<PROCS;i++,(bin_ptr += bin_size))
	bin_ptr->a = *(C_bin + i) - bin_ptr;
    
#if 1
    barrier();
#endif
    mpc_index(C_ptr,D_ptr,bin_size*(sizeof(intpair_t)),ALLGRP);
    
    /* Move the Contents of D into Routed */
    
    values = 0;
    myout_ptr = (int *)(Routed + MYPROC);  
    t_ptr = D_ptr + PROCS*bin_size;
    s_ptr = D_ptr - bin_size;
    while ((buff_ptr = (s_ptr += bin_size)) < t_ptr) {
	values += (entries = (buff_ptr->a));
	p1_ptr = buff_ptr + entries; 
	while (++buff_ptr <= p1_ptr)
	    *(myout_ptr + ((buff_ptr->a) & ((h*q)-1))) = buff_ptr->b;
    }
    
    barrier();
    
    free(C_bin);
    free(Loc);
    free(Info);
    
    all_spread_free(Buffer1);
    all_spread_free(Buffer2);
    
#if 0
    barrier();
    on_one {
	fprintf(outfile,"hrel4 Output on PE %3d:\n",MYPROC);
	fprintf(outfile,"p = %3d  h = %2d  q = %12d  (h*q: %d j: %d)\n",
		PROCS, h, q, h*q, values);
	for (i=0 ; i<values ; i++)
	    fprintf(outfile,"elem[%12d]: %12d\n",i,
		    *((int *)(Routed+MYPROC)+i));
	fflush(outfile);
    }
    barrier();
#endif

    return values;
    
}
#else
int all_route_h_rel_det4(int q, int h,
			 int *spread Keys,
			 int *spread Addr,
			 int *spread Routed) {}
#endif

int all_route_linear_perm(int q, int h,
			  int *spread Keys,
			  int *spread Addr,
			  int *spread Routed)
{
    int j, hq;
    register int i, hq_mask;

    int *lAddr,
	*lKeys,
	*lRout,
	**commMatrix,
	*lcommS,
	LOGQ,
	*offset,
	*send_offset;
    
    int *spread commMatrixS;

    intpair_t *spread sendArr;
    intpair_t *spread recvArr;
    intpair_t *global dest;
    intpair_t *loc,
	*lsend,
	*lsend2,
	*lrecv,
	**bptr;

    hq = h*q;
    
    commMatrixS = all_spread_malloc(PROCS, PROCS*sizeof(int));
    assert_spread_malloc(commMatrixS);

    sendArr = all_spread_malloc(PROCS, q*sizeof(intpair_t));
    assert_spread_malloc(sendArr);

    recvArr = all_spread_malloc(PROCS, hq*sizeof(intpair_t));
    assert_spread_malloc(recvArr);

    commMatrix = (int **)malloc(PROCS*sizeof(int*));
    assert_malloc(commMatrix);

    for (i=0 ; i<PROCS ; i++) {
	commMatrix[i] = (int*)malloc(PROCS*sizeof(int));
	assert_malloc(commMatrix[i]);
    }
	
    offset = (int*)malloc(PROCS*sizeof(int));
    assert_malloc(offset);

    send_offset = (int*)malloc(PROCS*sizeof(int));
    assert_malloc(send_offset);

    lsend2 = (intpair_t *)malloc(q*sizeof(intpair_t));
    assert_malloc(lsend2);

    bptr = (intpair_t **)malloc(PROCS*sizeof(intpair_t *));
    assert_malloc(bptr);

    lAddr  = (int *)(Addr + MYPROC);
    lKeys  = (int *)(Keys + MYPROC);
    lRout  = (int *)(Routed + MYPROC);
    lcommS = (int *)(commMatrixS + MYPROC);
    lsend  = (intpair_t *)(sendArr + MYPROC);
    lrecv  = (intpair_t *)(recvArr + MYPROC);

    barrier();

    LOGQ = (int)log2((double)(hq));

    for (i=0 ; i<PROCS ; i++) 
	lcommS[i] = 0;

    for (i=0 ; i<q ; i++) {
	lsend[i].a = lKeys[i];
	j = lsend[i].b = lAddr[i];

	lcommS[j >> LOGQ]++;
    }

    barrier();

    for (i=0 ; i<PROCS ; i++) {
	j = MYPROC ^ i;
	bulk_read(commMatrix[j], commMatrixS +j,
		  PROCS*sizeof(int));
    }

    bptr[0] = lsend2;
    for (i=1 ; i<PROCS ; i++) 
	bptr[i] = bptr[i-1] + lcommS[i-1];
    for (i=0 ; i<q ; i++) {
	j = (lsend[i].b>>LOGQ);
	bptr[j]->a = lsend[i].a;
	bptr[j]->b = lsend[i].b;
	bptr[j]++;
    }
    
    
/*    sync(); */
    barrier();

    for (i=0 ; i<PROCS ; i++) {
	offset[i] = 0;
	for (j=0 ; j<MYPROC ; j++)
	    offset[i] += commMatrix[j][i];
    }

    send_offset[0] = 0;
    for (i=1 ; i<PROCS ; i++)
	send_offset[i] = send_offset[i-1] + lcommS[i-1];

#if (defined(CM5))
    for (i=0 ; i<PROCS ; i++) {
	j = MYPROC^i;
	dest = (intpair_t *global)(recvArr + j);
	dest += offset[j];
	bulk_store(dest, lsend2+send_offset[j],
		   commMatrix[MYPROC][j]*sizeof(intpair_t));
	barrier();
    }
    all_store_sync();
    barrier();
#elif (defined(SP2))
    for (i=0 ; i<PROCS ; i++) {
	j = MODPROCS(MYPROC+i);
	dest = (intpair_t *global)(recvArr + j);
	dest += offset[j];
	bulk_write(dest, lsend2+send_offset[j],
		   commMatrix[MYPROC][j]*sizeof(intpair_t));
    }
    barrier();
#elif (defined(T3D))
    for (i=0 ; i<PROCS ; i++) {
	j = MYPROC^i;
	dest = (intpair_t *global)(recvArr + j);
	dest += offset[j];
	bulk_write(dest, lsend2+send_offset[j],
		   commMatrix[MYPROC][j]*sizeof(intpair_t));
	barrier();
    }
#elif (defined(MEIKO))
    for (i=0 ; i<PROCS ; i++) {
	j = MYPROC^i;
	dest = (intpair_t *global)(recvArr + j);
	dest += offset[j];
	bulk_write(dest, lsend2+send_offset[j],
		   commMatrix[MYPROC][j]*sizeof(intpair_t));
    }
#else
    for (i=0 ; i<PROCS ; i++) {
	j = MYPROC^i;
	dest = (intpair_t *global)(recvArr + j);
	dest += offset[j];
	bulk_write(dest, lsend2+send_offset[j],
		   commMatrix[MYPROC][j]*sizeof(intpair_t));
    }
    barrier();
#endif

    j = 0;
    for (i=0 ; i<PROCS ; i++)
	j += commMatrix[i][MYPROC];

    hq_mask = hq - 1;
    for (i=0 ; i<j ; i++) 
	*(lRout + (lrecv[i].b & hq_mask)) = lrecv[i].a;

    barrier();

    free(bptr);
    free(lsend2);
    free(send_offset);
    free(offset);
    for (i=0 ; i<PROCS ; i++) 
	free(commMatrix[i]);
    free(commMatrix);

    all_spread_free(recvArr);
    all_spread_free(sendArr);
    all_spread_free(commMatrixS);

#if CHECK_HREL
    all_check_hrel(j, lRout, q, h);
#endif
	
    return(j);
}


int all_route_single_phase(int q, int h, int *spread Keys,
			   int *spread Addr, int *spread Routed) 
     
{ 
  register int
      i, t,
      LOGQ,
      entries,
      values,
      dest;

  int PROCS_m1  = PROCS - 1;
  int bin_size  = ((h*q) >> PROCSLOG) + PROCS/2 + 2;
  int size_m1   = q*h - 1;    
  
  double time1, time2, time3, time4, time5;
  double tot_time;
  
  intpair_t
      *A_ptr,
      *buff_ptr,
      *B_ptr,
      *t_ptr,
      *temp_ptr;
  
  int
      *in_ptr,
      *rank_ptr,
      *ti_ptr,
      *out_ptr,
      *Count,
      *Prefix,
      *hist_ptr;
  
  intpair_t *global dest_ptr;
  
  intpair_t *spread A;
  intpair_t *spread B; 
  
  intpair_t **Loc;
  
  intpair_t *spread Buffer1;
  intpair_t *spread Buffer2;
  
  int *spread Hist1;
  int *spread Hist2;
  
  Loc = (intpair_t **) malloc(PROCS*sizeof(intpair_t *));
  assert_malloc(Loc);
  
  Count = (int *) malloc(PROCS*sizeof(int));
  assert_malloc(Count);

  Prefix = (int *) malloc(PROCS*sizeof(int));
  assert_malloc(Prefix);
  
  Buffer1 = all_spread_malloc(PROCS,PROCS*bin_size*sizeof(intpair_t));
  assert_spread_malloc(Buffer1);
  
  Buffer2 = all_spread_malloc(PROCS,PROCS*bin_size*sizeof(intpair_t));
  assert_spread_malloc(Buffer2);

  Hist1 = all_spread_malloc(PROCS,PROCS*sizeof(int));
  assert_spread_malloc(Hist1);

  Hist2 = all_spread_malloc(PROCS,PROCS*sizeof(int));  
  assert_spread_malloc(Hist2);
  
  /* Sort the keys by destination: */
  
  barrier();
  
  time1 = get_seconds();
  
  LOGQ = (int) ((log((double)(h*q)))/(log(2.0)));
  A = Buffer1;
  in_ptr = (int *)(Keys + MYPROC) - 1;
  rank_ptr = (int *)(Addr + MYPROC) - 1;
  
  for (i=0 ; i<PROCS ; i++)
      *(Count + i) = 0;
  
  ti_ptr = rank_ptr + q;
  while (++rank_ptr <= ti_ptr) 
      ++(*(Count +  (*rank_ptr >> LOGQ)));
  
  *Loc = A_ptr = (intpair_t *)(A + MYPROC);
  
  hist_ptr = (int *)(Hist1 + MYPROC); 
  
  *Prefix = 0;
  
  for (i=1 ; i<PROCS ; i++) {
      *(Loc + i) = A_ptr + (*(Prefix + i) = *(Prefix + i - 1) + 
			    (*(hist_ptr + i - 1) = *(Count + i - 1)));
  }
  
  *(hist_ptr + PROCS - 1) = *(Count + PROCS - 1);
  
  rank_ptr = (int *)(Addr + MYPROC) - 1;
  
  ti_ptr = rank_ptr + q;
  while (++rank_ptr <= ti_ptr) {
      (temp_ptr = ((*(Loc + (*rank_ptr >> LOGQ)))++))->a = *rank_ptr;
      temp_ptr->b = *(++in_ptr);
   }
  
  barrier();
  time1 = get_seconds() - time1;    
  
  
  /* Transpose the Histogram: */
  
  time2 = get_seconds();
  
  for (i=0 ; i<PROCS ; i++) {
      dest = (MYPROC+i) & PROCS_m1;
      bulk_store((int *global)(Hist2+dest) + MYPROC,
		 hist_ptr + dest,
		 sizeof(int)); 
  }
  
  all_store_sync(); 
  
  hist_ptr = (int *)(Hist2 + MYPROC); 
  
  for (i=1 ; i<PROCS ; i++)
      *(hist_ptr + i) += *(hist_ptr + i - 1);
  
  values = *(hist_ptr + PROCS_m1);
  
  for (i=PROCS_m1;i>0;i--) 
    *(hist_ptr + i) = *(hist_ptr + i - 1); 
  
  *(hist_ptr) = 0; 
  
  for (i=0;i<PROCS;i++) {
      dest = ((MYPROC+i) & PROCS_m1);
      bulk_store((int *global)(Hist1+dest) + MYPROC,
		 hist_ptr + dest,
		 sizeof(int)); 
  }
  
  all_store_sync(); 
  
  barrier();
  time2 = get_seconds() - time2;    
  
  
  /* Route the Contents of A into B */ 
  
  time3 = get_seconds();
  
  B = Buffer2;
  A_ptr = (intpair_t *)(A + MYPROC); 
  hist_ptr = (int *)(Hist1 + MYPROC);     
  
  for (i=0 ; i<PROCS ; i++) {
      dest = MYPROC^i;
      buff_ptr = A_ptr + *(Prefix + dest);
      dest_ptr = (intpair_t *global)(B + dest) + *(hist_ptr + dest);   
      bulk_store(dest_ptr,
		 buff_ptr,
		 (*(Count + dest))*sizeof(intpair_t)); 
  }
  
  all_store_sync(); 
  
  barrier();
  time3 = get_seconds() - time3;    
  
  /* Unpack the Received Data: */
  
  time4 = get_seconds();
  
  out_ptr = (int *)(Routed + MYPROC);  
  buff_ptr = (B_ptr = (intpair_t *)(B + MYPROC)) - 1;
  t_ptr = B_ptr + values;
  
  while (++buff_ptr < t_ptr)
      *(out_ptr + ((buff_ptr->a) & size_m1)) = buff_ptr->b;
  
  
  barrier();
  time4 = get_seconds() - time4;    
  
  
  tot_time = time1+time2+time3+time4;
  
  free(Loc);
  free(Count);
  free(Prefix);
  
  all_spread_free(Buffer1);
  all_spread_free(Buffer2);
  all_spread_free(Hist1);   
  all_spread_free(Hist2);    

#if CHECK_HREL
  all_check_hrel(values, out_ptr, q, h);
#endif

  return(values);

}


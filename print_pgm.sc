/*									tab:8
 *
 * print_pgm.sc - output images to Portable Graymap (PGM) format
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
 * Author: 		David A. Bader <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:		1
 * Creation Date:	November 2, 1995
 * Filename:		print_pgm.sc
 * History:
 */

#include "print_pgm.h"
#include "connComp_gen.h"
#include "radixsort.h"

static void
print_pgm_preamble(FILE *outfile, char* desc, int n, int m, int k)
{
    fprintf(outfile,"P5\n");
    fprintf(outfile,"# %s\n",desc);
    fprintf(outfile,"%d %d \n",m, n);
    fprintf(outfile,"%d\n",k-1);
}

static void
print_pgm_preamble_with_labels(FILE *outfile,
				 char* desc, int n, int m, int k, int l)
{
    fprintf(outfile,"P5\n");
    fprintf(outfile,"# %s, %d labels\n",desc,l);
    fprintf(outfile,"%d %d \n",m, n);
    fprintf(outfile,"%d\n",k-1);
}

void print_pgm_im(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r])
{
    register int
	i, ii, j, jj;

    FILE *outfile;

    outfile = fopen(outfilename,"w+");

    print_pgm_preamble(outfile,desc,v*q,r*w,k);
    
    for (i=0 ; i<v ; i++) 
	for (ii=0 ; ii < q ; ii++) 
	    for (j=0 ; j<w ; j++)
		for (jj=0 ; jj < r ; jj++)
		    fprintf(outfile,"%c",(unsigned char)(image[i][j][ii][jj]));

    fclose(outfile);

}

void print_pgm_cc_orig(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r])
{
    register int
	i, ii, j, jj;

    register double val;

    FILE *outfile;

    outfile = fopen(outfilename,"w+");
    
    print_pgm_preamble(outfile,desc,v*q,r*w,k);

    for (i=0 ; i<v ; i++) 
	for (ii=0 ; ii < q ; ii++) 
	    for (j=0 ; j<w ; j++)
		for (jj=0 ; jj < r ; jj++) {
		    val = ((double)image[i][j][ii][jj]) /
			((double)(v*w*q*r) + 1.0);
		    val *= 255.0;
		    val = floor(val);
		    fprintf(outfile,"%c",(unsigned char)((int)val));
		}

    fclose(outfile);
}

void print_pgm_cc(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int image[v][w]::[q][r])
{
    register int
	i, ii, j, jj;

    register double val;

    FILE *outfile;

    outfile = fopen(outfilename,"w+");
    
    print_pgm_preamble(outfile,desc,v*q,r*w,k);

/* Instead of linear mapping of labels to [0..255], use
   255 (sin (1E+7 * x) + 1.0) / 2 */

    for (i=0 ; i<v ; i++) 
	for (ii=0 ; ii < q ; ii++) 
	    for (j=0 ; j<w ; j++)
		for (jj=0 ; jj < r ; jj++) {
		    val = sin(12345678900.0*(double)image[i][j][ii][jj]) + 1.0;
		    val *= (255.0 / 2.0);
		    val = floor(val);
		    fprintf(outfile,"%c",(unsigned char)((int)val));
		}

    fclose(outfile);
}

static int binArrSearch(int * searchList, int min_idx, int max_idx, int target) 
{
    register int mid_idx;

    mid_idx = min_idx + (max_idx - min_idx) / 2;
    if (max_idx < min_idx)
	return(-1);
    else if (target < searchList[mid_idx])
	return binArrSearch(searchList, min_idx,     mid_idx - 1, target);
    else if (target > searchList[mid_idx])
	return binArrSearch(searchList, mid_idx + 1, max_idx    , target);
    else
	return mid_idx;
}

static inline int lookup_label(int * intArr, int arrSize, int label)
{
    return ( binArrSearch(intArr, 0, arrSize-1, label) );
}

#if 0
int lookup_label(int *labels, int n, int lab) {
    register int i=0;

    while ((i<n) && (lab > labels[i]))
	i++;

    if (lab == labels[i])
	return(i);
    else return(0);
}
#endif

int print_pgm_real_color(char *desc, char *outfilename,
			   int v, int w, int q, int r, int k,
			   int maxLabel,
			   int Labels[v][w]::[q][r], 
			   int X[v][w]::[q][r])
{

    register int
	i, j, ii, jj,
	idx, get_proc;

    int n, minL, maxL,
	l, t, rows, steps;

    int RVAL = 1<<10;    /* Number of buckets per loop */
    
    int (*locLabel)[r],
	(*locX)[r],
	(*locRealIm)[r];

    int RealIm[v][w]::[q][r];
    
    intpair_t
	RColor[PROCS]::[RVAL],
	Xpose[PROCS]::[RVAL],
	*locRColor,
	*locXpose;
    
    FILE *outfile;

    on_one
	outfile = fopen(outfilename,"w+");

    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);
    
    locX = tolocal(X[ii][jj]);
    locLabel = tolocal(Labels[ii][jj]);
    locRealIm = tolocal(RealIm[ii][jj]);
    locRColor = tolocal(RColor[MYPROC]);
    locXpose = tolocal(Xpose[MYPROC]);

    l=0;
    n = v*w*q*r;
    rows = DIVPROCS(RVAL);
    steps = maxLabel / RVAL;
    if (steps*RVAL != maxLabel)
	steps++;

    barrier();

    for (t=0 ; t<steps ; t++) {
	for (i=0 ; i<RVAL ; i++) {
	    locRColor[i].a = 0;
	    locRColor[i].b = 0;
	}
	minL = t*RVAL;
	maxL = (t+1)*RVAL;

	for (i=0 ; i < q ; i++) 
	    for (j=0 ; j < r ; j++) 
		if ((locLabel[i][j] >= minL)&&(locLabel[i][j] < maxL)) {
		    idx = locLabel[i][j] % RVAL;
		    locRColor[idx].a++;
		    locRColor[idx].b += locX[i][j];
		}

	barrier();
	
	for (i=0 ; i<PROCS ; i++) {
	    /* Stagger prefetches using MYPROC */
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locXpose[rows*get_proc],
		      &(RColor[get_proc][rows*MYPROC]),
		      rows*sizeof(intpair_t));
	    barrier();
	}

	for (i=0 ; i<rows ; i++)
	    for (j=1 ; j<PROCS ; j++) {
		locXpose[j*rows + i].a += locXpose[(j-1)*rows + i].a;
		locXpose[j*rows + i].b += locXpose[(j-1)*rows + i].b;
	    }

	barrier();

	for (i=0 ; i<rows ; i++)
	    if (locXpose[(PROCS-1)*rows + i].a > 0)
		locRColor[(MYPROC*rows) + i].a =
		    locXpose[(PROCS-1)*rows + i].b /
		    locXpose[(PROCS-1)*rows + i].a;
	    else
		locRColor[(MYPROC*rows) + i].a = 0;

	barrier();

	for (i=1 ; i<PROCS ; i++) {
	    /* Stagger prefetches using MYPROC */
	    get_proc = (MYPROC + i) % PROCS;
	    bulk_read(&locRColor[rows*get_proc],
		      &(RColor[get_proc][rows*get_proc]),
		      rows*sizeof(intpair_t));
	    barrier();
	}

	on_one
	    for (i=0; i<RVAL ; i++)
		if (locRColor[i].a > 0)
		    l++;

	/* Color output */
	for (i=0 ; i < q ; i++) 
	    for (j=0 ; j < r ; j++) 
		if ((locLabel[i][j] >= minL)&&(locLabel[i][j] < maxL)) {
		    idx = locLabel[i][j] % RVAL;
		    locRealIm[i][j] = locRColor[idx].a;
		}
    }

    barrier();

    on_one{

	print_pgm_preamble_with_labels(outfile,desc,v*q,r*w,k,l);

	for (i=0 ; i<v ; i++) 
	    for (ii=0 ; ii < q ; ii++) 
		for (j=0 ; j<w ; j++)
		    for (jj=0 ; jj < r ; jj++) {
			fprintf(outfile,"%c",
				(unsigned char)(RealIm[i][j][ii][jj]));
		    }

	fclose(outfile);
    }

    barrier();
    return(1);
}

int print_pgm_real_color_on_one(char *desc, char *outfilename,
		    int v, int w, int q, int r, int k,
		    int Labels[v][w]::[q][r], 
		    int X[v][w]::[q][r])
{
    register int
	i, ii, j, jj,
	idx, count;
    int n, l, num;
    int* mapLabel;
    int* mapColor;
    
    edge_t *colorArr;

    register double sum;

    FILE *outfile;

    outfile = fopen(outfilename,"w+");
    n = v*w*q*r;

    if ((colorArr = (edge_t *)malloc(n*sizeof(edge_t)))==NULL)
	fprintf(stderr,"ERROR: cannot malloc colorArr\n");
    if ((mapColor = (int *)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: cannot malloc mapColor\n");
    if ((mapLabel = (int *)malloc(n*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: cannot malloc mapLabel\n");

    idx = 0;

    for (i=0 ; i<v ; i++) 
	for (ii=0 ; ii < q ; ii++) 
	    for (j=0 ; j<w ; j++)
		for (jj=0 ; jj < r ; jj++) {
		    colorArr[idx].a := Labels[i][j][ii][jj];
		    colorArr[idx].b := X[i][j][ii][jj];
		    idx++;
		}

    sync();

    fastsort_edge(colorArr, idx);

    l = 0;
    count = 0;
    sum = 0.0;

    for (i=0 ; i<idx-1 ; i++) {
	count++;
	sum += (double)colorArr[i].b;
	if (colorArr[i].a != colorArr[i+1].a) {
	    mapLabel[l] = colorArr[i].a;
	    mapColor[l] = (int)floor(sum / (double)count);
#if 0
	    fprintf(stdout,"(%5d) Mapping Label %4d to color %3d\n",
		    l,mapLabel[l], mapColor[l]);
#endif
	    l++;
	    count = 0;
	    sum = 0.0;
	}
    }

    count++;
    sum += colorArr[idx-1].b;
    mapLabel[l] = colorArr[idx-1].a;
    mapColor[l] = (int)floor(sum / (double)count);
    l++;

    print_pgm_preamble_with_labels(outfile,desc,v*q,r*w,k,l);

    for (i=0 ; i<v ; i++) 
	for (ii=0 ; ii < q ; ii++) 
	    for (j=0 ; j<w ; j++)
		for (jj=0 ; jj < r ; jj++) {
		    num = lookup_label(mapLabel,l,Labels[i][j][ii][jj]);
		    fprintf(outfile,"%c",(unsigned char)(mapColor[num]));
		}

    fclose(outfile);

    free(mapLabel);
    free(mapColor);
    free(colorArr);

    return(l);
}


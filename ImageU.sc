/*									tab:8
 *
 * ImageU.sc - Image Understanding Project
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
 * Version:		2
 * Creation Date:	October 20, 1994
 * Filename:		ImageU.sc
 * History:
 *	DAB	2	Sun Dec  4 14:03:17 EST 1994
 *		Added in all_connected_components_grey_l() call
 */

#include "ImageU.h"
#include "test_mach.h"
#include "math_func.h"

/***************************************************/
void all_get_image(int v, int w, int q, int r,
		   int X[v][w]::[q][r], int k, int image_option) 
/***************************************************/
{
    register int
	i, j,          /* each local pixel */
	ii, jj,        /* my logical grid */
	I, J,          /* a pixel's position */
	barSize;

    int	(*locX)[r];

    double n, m,
	m1[6], m2[6];

    locX = tolocal(X[mapRow(v,w,MYPROC)][mapCol(v,w,MYPROC)]);

    barSize = BARSIZE_DEFAULT;

    n = (double) (v * q);
    m = (double) (w * r);
    
    m1[0] = 0.125 * n;
    m1[1] = 0.375 * n;
    m1[2] = 0.625 * n;
    m1[3] = 0.875 * n;
    m1[4] = 0.5   * n;
    m1[5] = 0.75  * n;

    m2[0] = 0.125 * m;
    m2[1] = 0.375 * m;
    m2[2] = 0.625 * m;
    m2[3] = 0.875 * m;
    m2[4] = 0.5   * m;
    m2[5] = 0.75  * m;

    ii   = mapRow(v,w,MYPROC);
    jj   = mapCol(v,w,MYPROC);

/**********************************************/
/* Initialize Image Array X                   */
/**********************************************/

    RNG_init();

    for (i=0 ; i<q ; i++) {
	I = ii*q + i;
	for (j=0 ; j<r ; j++) {
	    J = jj*r + j;
	    switch (image_option) {
	    case 1 : locX[i][j] =    /* Rows */
		((I % (2*barSize)) < barSize); break;
	    case 2 : locX[i][j] =    /* Columns */
		((J % (2*barSize)) < barSize); break;
	    case 3 : locX[i][j] =    /* Diagonals / */
		(((I+J) % (2*barSize)) < barSize); break;
	    case 4 : locX[i][j] =    /* Diagonals \ */
		(((I>J ? I-J : J-I) % (2*barSize)) < barSize); break;
	    case 5 : locX[i][j] =    /* Cross */
		((I>=m1[0])&&(I<=m1[3])&&(J>=m2[1])&&(J<=m2[2]) ||
		 (I>=m1[1])&&(I<=m1[2])&&(J>=m2[0])&&(J<=m2[3]));
		break;
	    case 6 : locX[i][j] =    /* Disc */
		sqrt((I-m1[4])*(I-m1[4]) + (J-m2[4])*(J-m2[4]))
		<= (n>m? m2[5] : m1[5])*0.4;
		break;
	    case 7 : locX[i][j] =    /* Concentric Circles */
		(int)floor(sqrt((I-m1[4])*(I-m1[4]) + (J-m2[4])*(J-m2[4])))
		% (2*barSize) < barSize;
		break;
	    case 8 : locX[i][j] =    /* Four Boxes */
		((I>=m1[0])&&(I<=m1[3])&&(J>=m2[0])&&(J<=m2[3])) &&
		(!((I>=m1[0])&&(I<=m1[3])&&(J>=m2[1])&&(J<=m2[2]) ||
		 (I>=m1[1])&&(I<=m1[2])&&(J>=m2[0])&&(J<=m2[3])));
		break;
	    case 9 : locX[i][j] =    /* Spiral */
		(int)floor(sqrt(( J<=m2[4] ?  /* left half? */
		   ((I-m1[4]-barSize)*(I-m1[4]-barSize)) :
		   ((I-m1[4]+barSize)*(I-m1[4]+barSize)) ) +
		   (J-m2[4])*(J-m2[4]))) % (barSize) < barSize/2;
		break;
	    case 101 : locX[i][j] = random() & (k-1); break;   /* random */
	    case 102 : locX[i][j] = (i&(k-1)); break;          /* rows */
	    case 103 : locX[i][j] = (j&(k-1)); break;          /* columns */
	    case 104 : locX[i][j] = (((i+j)&3) == 2); break ;  /* diagonals / */
	    case 105 : locX[i][j] = (((i>j ? i-j : j-i)&3)==2);/* diagonals \ */
		break;
	    case 106 : locX[i][j] = 1; break;                  /* all one */
	    case 107 : locX[i][j] = 0; break;                  /* all zero */
	    case 240 : locX[i][j] = (random()%100)<40?1:0; break;
	    case 250 : locX[i][j] = (random()%100)<50?1:0; break;
	    case 260 : locX[i][j] = (random()%100)<60?1:0; break;
		
	    default : locX[i][j] = 1;
	    }
	}
    }

    barrier();
}

/***************************************************/
void all_get_image_file(int v, int w, int q, int r,
			int X[v][w]::[q][r], char* infilename)
/***************************************************/
{
    register int
	i, j,          /* tile position */
	ii, jj;        /* logical grid position */

    FILE* infile;
    unsigned char anum;

    barrier();

    on_one {

	infile = fopen(infilename,"r");

	if (infile == NULL) errprnt("input file does not exist.");

	rewind(infile);

/***************************/
/* Read Image Array X      */
/***************************/

	for (ii = 0 ; ii<v ; ii++)
	    for (i=0 ; i<q ; i++) 
		for (jj = 0 ; jj<w ; jj++) 
		    for (j=0 ; j<r ; j++) {
			fscanf(infile,"%c",&anum);
			X[mapRow(v,w,mapProc(v,w,ii,jj))]
			    [mapCol(v,w,mapProc(v,w,ii,jj))]
			    [i][j] = (int)anum;
		    }

	fclose(infile);
    }

    barrier();
}

/***************************************************/
splitc_main(int argc, char* argv[])
/***************************************************/
{
    int n,                          /* Image size  */
	v,                          /* processor array length */
	w,                          /* processor array width */
	q,                          /* tile length */
	r,                          /* tile width  */
	k,                          /* number of gray levels */
	l,                          /* local relaxation value, if used */
	zVal,                       /* variance input */
	image_option,               /* The image number to use */
	i;

    int (*spread X)[q][r];          /* the image   */

    int (*spread Y)[q][r];          /* intermediate image */

    char *outfilename = NULL,
	*infilename = NULL,
	*s,**argvv = argv;
   
    image_option = IMAGE_DEFAULT;
    n            = N_DEFAULT;
    k            = K_DEFAULT;
    watch_proc   = WATCH_PROC_DEFAULT;
    RESULTS      = 0;
    OUTPUT_IMAGE = 0;
    UNSUPPARAM   = 0;
    DO_TEST      = 0;
    DO_UNSUP     = 0;
    DO_HISTO     = 0;
    DO_CCBIN     = 0;
    DO_CCBIF     = 0;
    DO_CCGRY     = 0;
    DO_CCGRF     = 0;
    DO_CCGRYL    = 0;
    DO_SEGM      = 0;
    DO_GBLUR     = 0;
    DO_BENCHMARK = 0;  /* This overrides all other arguments */

    while (--argc > 0 && (*++argvv)[0] == '-')
	for (s=argvv[0]+1; *s != '\0'; s++)
	    switch (*s) {
	    case 'Z': DO_BENCHMARK = 1; break;
	    case 'V': OUTPUT_IMAGE = OUTPUT_IMAGE^(1<<VISTA_IMAGE); break;
	    case 'P': OUTPUT_IMAGE = OUTPUT_IMAGE^(1<<PS_IMAGE); break;
	    case 'F': OUTPUT_IMAGE = OUTPUT_IMAGE^(1<<GIF_IMAGE); break;
	    case 'M': OUTPUT_IMAGE = OUTPUT_IMAGE^(1<<PGM_IMAGE); break;
#if X11		
	    case 'X': OUTPUT_IMAGE = OUTPUT_IMAGE^(1<<X11_IMAGE); break;
#endif		
	    case 'v': RESULTS      = 1; break;
	    case 'h': DO_HISTO     = 1; break;
	    case 'b': DO_CCBIN     = 1; break;
	    case 'B': DO_CCBIF     = 1; break;
	    case 'g': DO_CCGRY     = 1; break;
	    case 'G': DO_CCGRF     = 1; break;
	    case 's': DO_SEGM      = 1; break;
	    case 't': DO_TEST      = 1; break;
	    case 'U': DO_UNSUP     = 1; /* Optional parameter */
		if (argc > 1)
		    if ((*++argvv)[0] != '-') {
			argc--;
			UNSUPPARAM = atoi(*argvv);
		    }
		    else
			argvv--;
		break;
	    case 'l': if (argc <= 1) 
		errprnt("local relaxation (l) expected after -l");
		argc--; l = atoi(*++argvv); DO_CCGRYL = 1; break;
	    case 'n': if (argc <= 1) 
		errprnt("size n of nxn image expected after -n");
		argc--; n = atoi(*++argvv); break;
	    case 'k': if (argc <= 1) 
		errprnt("number of colors expected after -k");
		argc--; k = atoi(*++argvv); break;
	    case 'i': if (argc <= 1) 
		errprnt("image number expected after -i");
		argc--; image_option = atoi(*++argvv); break;
	    case 'f': if (argc <= 1) 
		errprnt("input filename expected after -f (e.g. -f filename)");
		argc--;
		infilename = (char *)malloc(MAXLEN*sizeof(char));
		strcpy(infilename, *++argvv); break;
	    case 'w': if (argc <= 1) 
		errprnt("processor id to watch expected after -w");
		argc--; watch_proc = atoi(*++argvv); break;
	    case 'z': if (argc <= 1) 
		errprnt("variance expected after -z");
		argc--; zVal = atoi(*++argvv); DO_GBLUR = 1; break;
	    case 'o': if (argc <= 1) 
		errprnt("output filename expected after -o (e.g. -o filename)");
		argc--;
		outfilename = (char *)malloc(MAXLEN*sizeof(char));
		strcpy(outfilename, *++argvv); break;
	    default: errprnt("illegal option");
	    }
    if (argc != 0)
	errprnt(": [-h] [-n image_size] [-k colors] [-i image_number | -f filename] [-o outfile] [-b] [-B] [-g] [-G] [-Z] [-VPFM]");
    
    if (outfilename == NULL)
	outfile = stdout;
    else
	outfile = fopen(outfilename,"a+");

    on_one fprintf(outfile,"Running [v %d.%d]\n",VER_MAJOR,VER_MINOR);
    
    if (DO_BENCHMARK) {
/* BENCHMARK RUN */
	OUTPUT_IMAGE = 0;
	
	for (i=1 ; i<=9 ; i++)
	    for (n=128 ; n<=1024 ; n*= 2) {
		if (n < PROCS) n = PROCS;       /* Make sure that PROCS <= n */

		v = (int)pow(2.0 , floor( (double)PROCSLOG / 2.0 ));
		w = (int)pow(2.0 ,  ceil( (double)PROCSLOG / 2.0 ));
		q = n/v;
		r = n/w;
		k = 2;

		on_one {
		    fprintf(outfile,"Test Image Number: %d\n",i);
		    fprintf(outfile,"Image is %d x %d, with processor grid:",n,n);
		    fprintf(outfile," [%d,%d] and tiles: [%d,%d]\n",v,w,q,r);
		    fprintf(outfile,"k colors: %d \n",k);
		}

		X = all_spread_malloc(v*w, q*r*sizeof(int));
		assert_spread_malloc(X);

		on_one fprintf(outfile,"Creating Image... ");
		on_one fflush(outfile);
		all_get_image(v,w,q,r,X,k,i);
		on_one fprintf(outfile,"Done.\n");
		on_one fflush(outfile);

		all_connected_components_binary(v,w,q,r,X,k); 
		all_connected_components_binary_four(v,w,q,r,X,k); 
		all_spread_free(X);
		on_one fflush(outfile);
	    }

	/* Run DARPA Benchmark Image */
	n = 512;
	k = 256;
	infilename = DARPA_IMAGE;

	if (n < PROCS) n = PROCS;       /* Make sure that PROCS <= n */

	v = (int)pow(2.0 , floor( (double)PROCSLOG / 2.0 ));
	w = (int)pow(2.0 ,  ceil( (double)PROCSLOG / 2.0 ));
	q = n/v;
	r = n/w;

	on_one {
	    fprintf(outfile,"Image Filename: %s\n",infilename);	    
	    fprintf(outfile,"Image is %d x %d, with processor grid:",n,n);
	    fprintf(outfile," [%d,%d] and tiles: [%d,%d]\n",v,w,q,r);
	    fprintf(outfile,"k colors: %d \n",k);
	}

	X = all_spread_malloc(v*w, q*r*sizeof(int));
	assert_spread_malloc(X);

	if (infilename != NULL) {
	    on_one fprintf(outfile,"Loading Image... ");
	    on_one fflush(outfile);
	    all_get_image_file(v,w,q,r,X,infilename);
	    on_one fprintf(outfile,"Done.\n");
	    on_one fflush(outfile);
	}

	all_connected_components_grey(v,w,q,r,X,k);
	all_connected_components_grey_four(v,w,q,r,X,k);

	all_spread_free(X);
    }
    else if (DO_TEST) 
	all_test_mach();
    else if (DO_UNSUP) 
	all_unsupported(UNSUPPARAM);
    else { /* Use USER parameters */

	if (n < PROCS) n = PROCS;       /* Make sure that PROCS <= n */

	v = (int)pow(2.0 , floor( (double)PROCSLOG / 2.0 ));
	w = (int)pow(2.0 ,  ceil( (double)PROCSLOG / 2.0 ));
	q = n/v;
	r = n/w;

	on_one {
	    if (infilename != NULL)
		fprintf(outfile,"Image Filename: %s\n",infilename);	    
	    else
		fprintf(outfile,"Test Image Number: %d\n",image_option);
	    fprintf(outfile,"Image is %d x %d, with processor grid:",n,n);
	    fprintf(outfile," [%d,%d] and tiles: [%d,%d]\n",v,w,q,r);
	    fprintf(outfile,"k colors: %d \n",k);
	}

/* X is the image */
	X = all_spread_malloc(v*w, q*r*sizeof(int));
	assert_spread_malloc(X);

	if (infilename != NULL) {
	    on_one fprintf(outfile,"Loading Image... ");
	    on_one fflush(outfile);
	    all_get_image_file(v,w,q,r,X,infilename);
	    on_one fprintf(outfile,"Done.\n");
	    on_one fflush(outfile);
	}
	else {
	    on_one fprintf(outfile,"Creating Image... ");
	    on_one fflush(outfile);
	    all_get_image(v,w,q,r,X,k,image_option);
	    on_one fprintf(outfile,"Done.\n");
	    on_one fflush(outfile);
	}

	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Input Image... ");
		fflush(outfile);
		image_print_im("Original image",
			       "image-orig",
			       v,w,q,r,k,X);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
	
	if (DO_HISTO) all_histogram(v,w,q,r,X,k);
	if (DO_CCBIN) all_connected_components_binary(v,w,q,r,X,k); 
	if (DO_CCBIF) all_connected_components_binary_four(v,w,q,r,X,k); 
	if (DO_CCGRY) all_connected_components_grey(v,w,q,r,X,k);
	if (DO_CCGRF) all_connected_components_grey_four(v,w,q,r,X,k);
	if (DO_CCGRYL) all_connected_components_grey_l(v,w,q,r,X,k,l);
	if (DO_SEGM)   all_segment_image(v,w,q,r,X,k);
	if (DO_GBLUR) {
	    Y = all_spread_malloc(v*w, q*r*sizeof(int));
	    assert_spread_malloc(Y);
	    
	    all_blur_image(v,w,q,r,X,k,Y,zVal);
	    all_blur_image(v,w,q,r,X,k,Y,2);
	    all_blur_image(v,w,q,r,X,k,Y,3);
	    all_blur_image(v,w,q,r,X,k,Y,4);
	    all_blur_image(v,w,q,r,X,k,Y,5);
	    all_blur_image(v,w,q,r,X,k,Y,6);
	    all_blur_image(v,w,q,r,X,k,Y,7);
	    
	    all_spread_free(Y);
	}

/*********************/
/* Destroy image     */
/*********************/

	all_spread_free(X);
    }

	
    if (outfilename != NULL) {
	fflush(outfile);
	on_one fclose(outfile);
    }
#if X11
    if (OUTPUT_IMAGE & (1<<X11_IMAGE)) {
	on_one {
	    fprintf(stdout,"Press return to exit. \n");
	    fflush(stdout);
	    i = getchar();
	}
    }
#endif    
	

}

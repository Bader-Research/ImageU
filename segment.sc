/*									tab:8
 *
 * segment.sc - segment an image
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
 * Creation Date:	February 9, 1995
 * Filename:		segment.sc
 * History:
 */

#include "segment.h"

#define SNFLIMIT  200
#define SDTIMES   2
#define MAG       2

/*************************************************************/
void all_segment_image_conservative(int v, int w, int q, int r, 
		       int X[v][w]::[q][r], int k)
/*************************************************************/
{
    int epsilon,
	artif;

    double sigma,
	sd;
    
#if 0
    int Histo[1]::[k];            /* Histogram */
#endif
    
    int Z[v][w]::[q][r];          /* pType */
    int Y[v][w]::[q][r];          /* Intermediate Image */
    int Y2[v][w]::[q][r];          /* Intermediate Image */
    int outImage[v][w]::[q*MAG][r*MAG]; /* Output segmented image */
    int Ymag[v][w]::[q*MAG][r*MAG]; /* Output segmented image colors */
   
#if 0
    all_histogram_sub(v,w,q,r,X,k,Histo);
    k = all_remap_histo(v,w,q,r,X,Y,k,Histo); 
#endif

    all_get_pType(v,w,q,r,Z);

    artif = all_decide_artificial(v,w,q,r,X,Z);

    on_one {
	fprintf(outfile,"Image is ");
	if (artif)
	    fprintf(outfile,"artificial\n");
	else
	    fprintf(outfile,"real\n"); 
    }

    if (artif) {
	all_init_timer();
	all_start_timer();
	
	if (k==2)
	    all_connected_components_binary_four_sub(v,w,q,r,X,k,Y); 
	else
	    all_connected_components_grey_four_sub(v,w,q,r,X,k,Y);

	all_stop_timer("Connected Components");
	all_print_timer(outfile);

	on_one {
	  if (OUTPUT_IMAGE) {
	    fprintf(outfile,"Writing Segmented Image... ");
	    fflush(outfile);
	    image_print_cc("4-connected components segmentation",
			   "imageSegmented",
			   v,w,q,r,k,Y);
	    fprintf(outfile,"Done.\n");
	    fflush(outfile);
	  }
	}

    }
    else {

	sigma = all_find_sigma_use_median(v,w,q,r,X,k,Z);
	sd = sqrt(sigma);
	epsilon = (int)ceil(SDTIMES * sd); 

	if (RESULTS)
	    on_one
		fprintf(outfile,"sd: %9.2f  epsilon: %d\n",sd,epsilon);

#if 0
	all_filter_SNF_reclass(v,w,q,r,X,k,X,Y,Z,epsilon); 
	all_filter_SNF_reclass(v,w,q,r,X,k,Y,Y,Z,epsilon); 
	all_filter_SNF_reclass(v,w,q,r,X,k,Y,Y,Z,epsilon); 
	all_filter_SNF_reclass(v,w,q,r,X,k,Y,Y,Z,epsilon);
#endif	
    
#if 0
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image 1... ");
		fflush(outfile);
		image_print_im("after reclass", "image1",
			       v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif

	all_filter_SNF(v,w,q,r,X,k,Y,Z, 4, 0);

#if 0
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image 2... ");
		fflush(outfile);
		image_print_im("after SNF1", "image2",
			       v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif

	all_filter_SNF(v,w,q,r,Y,k,Y,Z, SNFLIMIT, epsilon); 

#if 0
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image 3... ");
		fflush(outfile);
		image_print_im("after SNF2", "image3",
			       v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif

	all_filter_SNF(v,w,q,r,Y,k,Y,Z, SNFLIMIT, 0);

#if 0
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image 4... ");
		fflush(outfile);
		image_print_im("after SNF3", "image4",
			       v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif

	all_filter_1NN(v,w,q,r,Y,k,Y,Z);

	all_crop_border(v,w,q,r,Y,Y,3);

#if 1
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image (After enhance)... ");
		fflush(outfile);
		image_print_im("after enhancement process",
			       "imageSNF", v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif

	all_init_timer();
	all_start_timer();

	all_connected_components_grey_l_sub(v,w,q,r,Y,k,epsilon,Y2); 

	all_stop_timer("Connected Components");
	all_print_timer(outfile);

	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Segmented Output Image... ");
		fflush(outfile);
		image_print_real_color_on_one(
		    "SNF Segmented Image",
		    "imageSegmented",
		    v,w,q,r,k,Y2,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
	
#if 0
	all_scale_and_edge_detect(v,w,q,r,MAG,Y2,outImage,Z,0,Y,Ymag); 

	if (OUTPUT_IMAGE) {
	    barrier();
	    on_one {
		fprintf(outfile,"Writing Edge Colored Output Image... ");
		fflush(outfile);
	    }
	    image_print_real_color(
		"Edge Detected SNF/Segmented Image",
		"imageZedge",
		v,w,q*MAG,r*MAG,k,v*w*q*r, outImage,Ymag);
	    on_one {
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	    barrier();
	}
#endif
	
    }
}

/*************************************************************/
void all_segment_image(int v, int w, int q, int r, 
		       int X[v][w]::[q][r], int k)
/*************************************************************/
{
    int epsilon,
	artif;

    double sigma,
	sd;
    
    int Z[v][w]::[q][r];          /* pType */
    int Y[v][w]::[q][r];          /* Intermediate Image */
   
    all_get_pType(v,w,q,r,Z);

    artif = all_decide_artificial(v,w,q,r,X,Z);

    on_one {
	fprintf(outfile,"Image is ");
	if (artif)
	    fprintf(outfile,"artificial\n");
	else
	    fprintf(outfile,"real\n"); 
    }

    if (artif) {
	all_init_timer();
	all_start_timer();
	
	if (k==2)
	    all_connected_components_binary_four_sub(v,w,q,r,X,k,Y); 
	else
	    all_connected_components_grey_four_sub(v,w,q,r,X,k,Y);

	all_stop_timer("Connected Components");
	all_print_timer(outfile);

	on_one  {
	  if (OUTPUT_IMAGE) {
	    fprintf(outfile,"Writing Segmented Image... ");
	    fflush(outfile);
	    image_print_cc("4-connected components segmentation",
			   "imageSegmented",
			   v,w,q,r,k,Y);
	    fprintf(outfile,"Done.\n");
	    fflush(outfile);
	  }
	}
    }
    else {

	sigma = all_find_sigma_use_median(v,w,q,r,X,k,Z);
	sd = sqrt(sigma);
	epsilon = (int)ceil(SDTIMES * sd); 

	if (RESULTS)
	    on_one
		fprintf(outfile,"sd: %9.2f  epsilon: %d\n",sd,epsilon);

	all_filter_SNF(v,w,q,r,X,k,Y,Z, 4, 0);
	all_filter_SNF(v,w,q,r,Y,k,Y,Z, SNFLIMIT, epsilon); 
	all_filter_SNF(v,w,q,r,Y,k,Y,Z, SNFLIMIT, 0);
	all_filter_1NN(v,w,q,r,Y,k,Y,Z);
	all_crop_border(v,w,q,r,Y,Y,3);

#if 1
	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Output Image (After enhance)... ");
		fflush(outfile);
		image_print_im("after enhancement process",
			       "imageSNF", v,w,q,r,k,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
#endif
#if 1
	on_one {
	  fprintf(outfile,"Performing delta-Connected Components: d = %d \n",
		  epsilon);
	  fflush(outfile);
	}
#endif

	all_init_timer();
	all_start_timer();

	all_connected_components_grey_l_sub(v,w,q,r,Y,k,epsilon,Z); 

	all_stop_timer("delta-Connected Components");
	all_print_timer(outfile);

	on_one {
	    if (OUTPUT_IMAGE) {
		fprintf(outfile,"Writing Segmented Output Image... ");
		fflush(outfile);
		image_print_real_color_on_one(
		    "SNF Segmented Image",
		    "imageSegmented",
		    v,w,q,r,k,Z,Y);
		fprintf(outfile,"Done.\n");
		fflush(outfile);
	    }
	}
    }
}


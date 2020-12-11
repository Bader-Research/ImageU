#ifndef _MACH_SPECIFICS_H
#define _MACH_SPECIFICS_H

#if (defined (CM5))
#define RADIXSORT_INT_BREAKPT		  99
#define RADIXSORT_ELEM_BREAKPT		  69
#define RADIXSORT_CHPAIR_BREAKPT	  81
#define RADIXSORT_EDGE_BREAKPT		  70
#define MATRIX_XPOSE_BARRIER		   1

#elif (defined (SP2))
#define RADIXSORT_INT_BREAKPT             91
#define RADIXSORT_ELEM_BREAKPT            67
#define RADIXSORT_CHPAIR_BREAKPT          72
#define RADIXSORT_EDGE_BREAKPT            72
#define MATRIX_XPOSE_BARRIER               1

#elif (defined(MEIKO))
#define RADIXSORT_INT_BREAKPT             71
#define RADIXSORT_ELEM_BREAKPT            57
#define RADIXSORT_CHPAIR_BREAKPT          67
#define RADIXSORT_EDGE_BREAKPT            56
#define MATRIX_XPOSE_BARRIER               0

#elif (defined (T3D)
#define RADIXSORT_INT_BREAKPT           1760
#define RADIXSORT_ELEM_BREAKPT           960
#define RADIXSORT_CHPAIR_BREAKPT        1110
#define RADIXSORT_EDGE_BREAKPT          1010
#define MATRIX_XPOSE_BARRIER               0

#else

#define RADIXSORT_INT_BREAKPT		  99
#define RADIXSORT_ELEM_BREAKPT		  69
#define RADIXSORT_CHPAIR_BREAKPT	  81
#define RADIXSORT_EDGE_BREAKPT		  70
#define MATRIX_XPOSE_BARRIER		   0

#endif

#endif


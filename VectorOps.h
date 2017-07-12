/*	Matrix and Vector subroutines plus some NSPCG extras		*/

/*	extras - John B. Pormann, Duke University			*/
/*		 NSF/ERC for Emerging Cardiovascular Technologies	*/
/*		 see http://www.ee.duke.edu/people/jpormann.html	*/

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */

/* RCSID: $Id: VectorOps.h 14 2007-05-11 14:57:55Z jbp $ */


#ifndef _VECTOR_H
#define _VECTOR_H

/* vector "package" */
typedef struct {
	int   size;
	real* data;
	domain_t dtype;
} vector;
typedef struct {
	int  size;
	int* data;
	domain_t dtype;
} ivector;
typedef struct {
	int   size;
	byte* data;
	domain_t dtype;
} bvector;
/* use _INIT_VECTOR to initialize ***ALL*** vector structures!! */
#define _INIT_VECTOR {0,NULL,-1}

/* sparse matrix structure */
typedef struct {
	int   type;
	int   rows,cols;
	int   maxnz;
	int   csep;
	domain_t dtype;
	int   msgtag;
	int*  jcoef;
	real* coef;
} sparse;
/* use _INIT_SPARSE to initialize ***ALL*** sparse structures!! */
#define _INIT_SPARSE {0,0,0,0,0,-1,MPI_ANY_TAG,NULL,NULL}

/* the sparse matrix type field, supported by NSPCG */
#define _PRIMARY   1
#define _PRIMARY_U 1
#define _BAND_S    2
#define _BAND_U    3
#define _BAND      3
#define _BANDED    3
#define _COORD_S   4
#define _COORD_U   5
#define _COORD     5
/* not exactly supported by NSPCG */
#define _PRIMARY_S 6
/* definitely not supported by NSPCG */
#define _NULL      0
#define _STENCIL   10
#define _FULL      11
#define _ZERO      12

/* boundary condition flags for STN matrices */
#define _BC_NULL   0
#define _BC_ROW    1
#define _BC_COL    2
#define _BC_FULL   4

/* Function Prototypes: */
/* - memory allocation and freeing */
int    vecalloc( vector* v1, domain_t dt );
int    vecfree( vector* v1 );
int    ivecalloc( ivector* v1, domain_t dt );
int    ivecfree( ivector* v1 );
int    bvecalloc( bvector* v1, domain_t dt );
int    bvecfree( bvector* v1 );
int    spralloc( sparse* M, domain_t dtp, int type, int r, int c, int maxnz );
int    sprfree( sparse* M );

/* - basic memory-type operations */
int    vcopy( vector x, vector y );
int    vcopyi( vector x, int i, vector y, int j, int l );
int    vswap( vector* xp, vector* yp );
int    vzero( vector v );
int    vfill( real alpha, vector v );
int    ivcopy( ivector x, ivector y );
int    ivcopyi( ivector x, int i, ivector y, int j, int l );
int    ivswap( ivector* xp, ivector* yp );
int    ivzero( ivector v );
int    ivfill( int alpha, ivector v );
int    bvcopy( bvector x, bvector y );
int    bvcopyi( bvector x, int i, bvector y, int j, int l );
int    bvswap( bvector* xp, bvector* yp );
int    bvzero( bvector v );
int    bvfill( byte alpha, bvector v );
int    sswap( sparse* M1, sparse* M2 );
int    scopy( sparse M1, sparse M2 );
int    szero( sparse M1 );

/* - basic mathematical operations */
real   vmax( vector x );
real   vmin( vector x );
real   innerprod( vector x, vector y );
#define dotprod(X,Y) innerprod(X,Y)
real   norm1( vector x );
real   norm2( vector x );
real   normI( vector x );
int    vscale( real beta, vector x );
int    vrecip( real beta, vector x );
int    vvscale( real beta, vector x, vector y );
int    vvrecip( real beta, vector x, vector y );
int    voffset( real alpha, vector x );
int    vvmpy( real alpha, vector x, vector y, vector z );
int    vvmpyb( real alpha, vector x, vector y, real beta, vector z );
int    vvdiv( real alpha, vector x, vector y, vector z );
int    vvdivb( real alpha, vector x, vector y, real beta, vector z );
int    vvadd( real alpha, vector x, vector y );
int    vvaddb( real alpha, vector x, real beta, vector y );
int    sscale( real alpha, sparse M1 );
int    ssadd( real alpha, sparse M1, sparse M2 );
int    ssadddiag( real alpha, vector D1, sparse M2 );
/* **** Level 2 (sparse) BLAS */
int    sprvec( real alpha, sparse M, vector x, real beta, vector y );

/* **** helpers */
real   sprget( sparse A, int i, int j );
int    sprput( sparse A, int i, int j, real v );
void   sprnonz( sparse A, int i, int* n, int* lst );

/* - file I/O procedures */
#define tokensep " \t\n"
int    vecreadinfo( char* filename, int* sz );
int    vecread( char* filename, vector v1 );
int    vecwrite( char* filename, vector v1 );
int    ivecreadinfo( char* filename, int* sz );
int    ivecread( char* filename, ivector v1 );
int    ivecwrite( char* filename, ivector v1 );
int    bvecreadinfo( char* filename, int* sz );
int    bvecread( char* filename, bvector v1 );
int    bvecwrite( char* filename, bvector v1 );

/* - some simple-minded stuff */
extern char* prnt_hdr;
extern char* prnt_ftr;
extern char* prnt_sep;
extern char* prnt_fmt;
extern char* prnt_ifmt;
void vecprint( vector V );
void ivecprint( ivector V );
void bvecprint( bvector V );

/* this is really an internal function */
void swapendian( void* ptr, int dsize, int num );

#endif


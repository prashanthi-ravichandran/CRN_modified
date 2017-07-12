/*	Basic Vector functions for C				*/

/*	John Pormann, NSF/ERC-CECT at Duke University		*/
/*		jbp@erc-sparc.mc.duke.edu			*/

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */



#include "CardioWave.h"
#include "VectorOps.h"

#define TAG_SEQ 1234

/* for the print functions */
char*  prnt_hdr  = "[ ";
char*  prnt_ftr  = "]\n";
char*  prnt_sep  = ";\n ";
char*  prnt_fmt  = "%6.2lf ";
char*  prnt_ifmt = "%6d ";
char*  prnt_bfmt = "%3d ";

char* VectorOps_RCSID = "$Id: VectorOps.c 14 2007-05-11 14:57:55Z jbp $";

static logical SequentialFirsttime = true;

int vecalloc( vector* v1, domain_t dt ) {
	v1->size  = DomainLocalSize( dt );
	v1->dtype = dt;
	if( v1->size == 0 ) {
		v1->data = NULL;
	} else {
		v1->data = (real*)malloc( v1->size*sizeof(real) );
		if( v1->data == NULL ) {
			return( -1 );
		}
		memset( (char*)v1->data, 0, v1->size*sizeof(real) );

	}

	return( 0 );
}

int ivecalloc( ivector* v1, domain_t dt ) {
	v1->size  = DomainLocalSize( dt );
	v1->dtype = dt;
	if( v1->size == 0 ) {
		v1->data = NULL;
	} else {
		v1->data = (int*)malloc( v1->size*sizeof(int) );
		if( v1->data == NULL ) {
			return( -1 );
		}
		memset( (char*)v1->data, 0, v1->size*sizeof(int) );

	}

	return( 0 );
}

int bvecalloc( bvector* v1, domain_t dt ) {
	v1->size  = DomainLocalSize( dt );
	v1->dtype = dt;
	if( v1->size == 0 ) {
		v1->data = NULL;
	} else {
		v1->data = (byte*)malloc( v1->size*sizeof(byte) );
		if( v1->data == NULL ) {
			return( -1 );
		}
		memset( (char*)v1->data, 0, v1->size*sizeof(byte) );

	}

	return( 0 );
}

int vecfree( vector* v1 ) {
	if( v1->data != NULL ) {
		free( v1->data );
	}

	v1->size  = 0;
	v1->data  = NULL;
	v1->dtype = -1;

	return( 0 );
}

int ivecfree( ivector* v1 ) {
	if( v1->data != NULL ) {
		free( v1->data );
	}

	v1->size  = 0;
	v1->data  = NULL;
	v1->dtype = -1;

	return( 0 );
}

int bvecfree( bvector* v1 ) {
	if( v1->data != NULL ) {
		free( v1->data );
	}

	v1->size  = 0;
	v1->data  = NULL;
	v1->dtype = -1;

	return( 0 );
}

/* copy x to y */
int vcopy( vector x, vector y ) {
	int size;

	size = x.size;
	if( y.size != size ) {
		return( -1 );
	}

	memcpy( (char*)y.data, (char*)x.data, size*sizeof(real) );

	return( 0 );
}

int ivcopy( ivector x, ivector y ) {
	int size;

	size = x.size;
	if( y.size != size ) {
		return( -1 );
	}

	memcpy( (char*)y.data, (char*)x.data, size*sizeof(int) );

	return( 0 );
}

int bvcopy( bvector x, bvector y ) {
	int size;

	size = x.size;
	if( y.size != size ) {
		return( -1 );
	}

	memcpy( (char*)y.data, (char*)x.data, size*sizeof(byte) );

	return( 0 );
}

/* copy x from index i to y at index j, length l  */
int vcopyi( vector x, int i, vector y, int j, int l ) {
	int s;

	s = x.size;
	if( (i+l) > s ) {
		return( -1 );
	}
	s = y.size;
	if( (j+l) > s ) {
		return( -1 );
	}

	memcpy( (char*)(y.data+j), (char*)(x.data+i), l*sizeof(real) );

	return( 0 );
}

/* swap x and y */
int vswap( vector* xp, vector* yp ) {
	int size;
	real* temp;
	domain_t itemp;

	size     = yp->size;
	yp->size = xp->size;
	xp->size = size;

	temp     = xp->data;
	xp->data = yp->data;
	yp->data = temp;

	itemp     = xp->dtype;
	xp->dtype = yp->dtype;
	yp->dtype = itemp;

	return( 0 );
}

int ivswap( ivector* xp, ivector* yp ) {
	int size;
	int* temp;
	domain_t itemp;

	size     = yp->size;
	yp->size = xp->size;
	xp->size = size;

	temp     = xp->data;
	xp->data = yp->data;
	yp->data = temp;

	itemp     = xp->dtype;
	xp->dtype = yp->dtype;
	yp->dtype = itemp;

	return( 0 );
}

int bvswap( bvector* xp, bvector* yp ) {
	int size;
	byte* temp;
	domain_t itemp;

	size     = yp->size;
	yp->size = xp->size;
	xp->size = size;

	temp     = xp->data;
	xp->data = yp->data;
	yp->data = temp;

	itemp     = xp->dtype;
	xp->dtype = yp->dtype;
	yp->dtype = itemp;

	return( 0 );
}

/* zero out vector */
int vzero( vector v ) {
	memset( (char*)v.data, 0, v.size*sizeof(real) );

	return( 0 );
}

int ivzero( ivector v ) {
	memset( (char*)v.data, 0, v.size*sizeof(int) );

	return( 0 );
}

int bvzero( bvector v ) {
	memset( (char*)v.data, 0, v.size*sizeof(byte) );

	return( 0 );
}

/* fill in the vector */
int vfill( real alpha, vector v ) {
	int i,m;
	int size = v.size;
	real* vv = v.data;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		vv[i]   = alpha;
		vv[i+1] = alpha;
		vv[i+2] = alpha;
		vv[i+3] = alpha;
	}
	for(;i<size;i++) {
		vv[i] = alpha;
	}

	return( 0 );
}

int ivfill( int alpha, ivector v ) {
	int i,m;
	int size = v.size;
	int* vv = v.data;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		vv[i]   = alpha;
		vv[i+1] = alpha;
		vv[i+2] = alpha;
		vv[i+3] = alpha;
	}
	for(;i<size;i++) {
		vv[i] = alpha;
	}

	return( 0 );
}

int bvfill( byte alpha, bvector v ) {
	int i,m;
	int size = v.size;
	byte* vv = v.data;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		vv[i]   = alpha;
		vv[i+1] = alpha;
		vv[i+2] = alpha;
		vv[i+3] = alpha;
	}
	for(;i<size;i++) {
		vv[i] = alpha;
	}

	return( 0 );
}

/* inner or dot product */
real innerprod( vector x, vector y ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real* yv = y.data;
	real sum,finalsum;

	if( size != y.size ) {
		return( -1 );
	}

	m = size - size%4;

	/* compute local piece */
	sum = 0.0;
	for(i=0;i<m;i+=4) {
		sum += xv[i]*yv[i] + xv[i+1]*yv[i+1] + xv[i+2]*yv[i+2] 
		     + xv[i+3]*yv[i+3];
	}
	for(;i<size;i++) {
		sum += xv[i]*yv[i];
	}

	MPI_Allreduce( &sum, &finalsum, 1, MPI_SSREAL, MPI_SUM,
		MPI_COMM_WORLD );

	return( finalsum );
}

/* returns the maximum value in x */
real vmax( vector x ) {
	int i;
	int size = x.size;
	real* xv = x.data;
	real mx,finalmx;

	mx = xv[0];
	for(i=1;i<size;i++) {
		if( xv[i] > mx ) {
			mx = xv[i];
		}
	}

	MPI_Allreduce( &mx, &finalmx, 1, MPI_SSREAL, MPI_MAX,
		MPI_COMM_WORLD );

	return( finalmx );
}

/* returns the minimum value in x */
real vmin( vector x ) {
	int i;
	int size = x.size;
	real* xv = x.data;
	real mn,finalmn;

	mn = xv[0];
	for(i=1;i<size;i++) {
		if( xv[i] < mn ) {
			mn = xv[i];
		}
	}

	MPI_Allreduce( &mn, &finalmn, 1, MPI_SSREAL, MPI_MIN,
		MPI_COMM_WORLD );

	return( finalmn );
}

/* 1-norm = sum(abs(x[i])) */
/* - should go from smallest components to largest to avoid rounding */
/*   errors (ie. 1+0.005+0.005 may go to 1.00 instead of 1.01)       */
real norm1( vector x ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real sum = 0.0, finalsum;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		sum += fabs(xv[i]) + fabs(xv[i+1]) + fabs(xv[i+2]) 
		     + fabs(xv[i+3]);
	}
	for(;i<size;i++) {
		sum += fabs(xv[i]);
	}

	MPI_Allreduce( &sum, &finalsum, 1, MPI_SSREAL, MPI_SUM,
		MPI_COMM_WORLD );

	return( finalsum );
}
int bnorm1( bvector x ) {
	int i,m;
	int size = x.size;
	byte* xv = x.data;
	int sum = 0, finalsum;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		sum += abs(xv[i]) + abs(xv[i+1]) + abs(xv[i+2]) 
		     + abs(xv[i+3]);
	}
	for(;i<size;i++) {
		sum += abs(xv[i]);
	}

	MPI_Allreduce( &sum, &finalsum, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD );

	return( finalsum );
}
int inorm1( ivector x ) {
	int i,m;
	int size = x.size;
	int* xv = x.data;
	int sum = 0, finalsum;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		sum += abs(xv[i]) + abs(xv[i+1]) + abs(xv[i+2]) 
		     + abs(xv[i+3]);
	}
	for(;i<size;i++) {
		sum += abs(xv[i]);
	}

	MPI_Allreduce( &sum, &finalsum, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD );

	return( finalsum );
}

/* 2-norm = sqrt(sum(x[i]^2)) */
/* - should go from smallest components to largest to avoid rounding */
/*   errors (ie. 1+0.005+0.005 may go to 1.00 instead of 1.01)       */
real norm2( vector x ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real sum = 0.0, finalsum;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		sum += xv[i]*xv[i] + xv[i+1]*xv[i+1] + xv[i+2]*xv[i+2] 
		     + xv[i+3]*xv[i+3];
	}
	for(;i<size;i++) {
		sum += xv[i]*xv[i];
	}

	MPI_Allreduce( &sum, &finalsum, 1, MPI_SSREAL, MPI_SUM,
		MPI_COMM_WORLD );

	return( sqrt(finalsum) );
}

/* Inf-norm = max(abs(x[i])) */
real normI( vector x ) {
	int i;
	int size = x.size;
	real* xv = x.data;
	real mx  = -9.9e9, finalmx;
	real av;

	mx = fabs( xv[0] );
	for(i=1;i<size;i++) {
		av = fabs( xv[i] );
		if( av > mx ) {
			mx = av;
		}
	}

	MPI_Allreduce( &mx, &finalmx, 1, MPI_SSREAL, MPI_MAX,
		MPI_COMM_WORLD );

	return( finalmx );
}

/* x = beta*x */
int vscale( real beta, vector x ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;

	if( beta == 1.0 ) {
		return( 0 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		xv[i]   *= beta;
		xv[i+1] *= beta;
		xv[i+2] *= beta;
		xv[i+3] *= beta;
	}
	for(;i<size;i++) {
		xv[i] *= beta;
	}

	return( 0 );
}

/* y = alpha*x*y */
int vvscale( real alpha, vector x, vector y ) {
	int i,m;
	int size = y.size;
	real* xv = x.data;
	real* yv = y.data;
	real  x0,x1,x2,x3;

	if( x.size == 1 ) {
		/* assume that this is a scalar */
		if( alpha*xv[0] == 1.0 ) {
			return( 0 );
		}
		vscale( alpha*xv[0], y );
		return( 0 );
	}
	if( size != x.size ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		yv[i]   *= x0;
		yv[i+1] *= x1;
		yv[i+2] *= x2;
		yv[i+3] *= x3;
	}
	for(;i<size;i++) {
		yv[i] *= alpha*xv[i];
	}

	return( 0 );
}

/* y = alpha*(1/y) */
int vrecip( real alpha, vector x ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real  x0,x1,x2,x3;

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = 1.0/xv[i];
		x1 = 1.0/xv[i+1];
		x2 = 1.0/xv[i+2];
		x3 = 1.0/xv[i+3];
		xv[i]   = alpha*x0;
		xv[i+1] = alpha*x1;
		xv[i+2] = alpha*x2;
		xv[i+3] = alpha*x3;
	}
	for(;i<size;i++) {
		xv[i] = alpha/xv[i];
	}

	return( 0 );
}

/* y = alpha*(1/x)*y */
int vvrecip( real alpha, vector x, vector y ) {
	int i,m;
	int size = y.size;
	real* xv = x.data;
	real* yv = y.data;
	real  x0,x1,x2,x3;

	if( x.size == 1 ) {
		/* assume that this is a scalar */
		vscale( alpha/xv[0], y );
		return( 0 );
	}
	if( size != x.size ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha/xv[i];
		x1 = alpha/xv[i+1];
		x2 = alpha/xv[i+2];
		x3 = alpha/xv[i+3];
		yv[i]   *= x0;
		yv[i+1] *= x1;
		yv[i+2] *= x2;
		yv[i+3] *= x3;
	}
	for(;i<size;i++) {
		yv[i] *= alpha/xv[i];
	}

	return( 0 );
}

/* x = beta+x */
int voffset( real beta, vector x ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;

	if( beta == 0.0 ) {
		return( 0 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		xv[i]   += beta;
		xv[i+1] += beta;
		xv[i+2] += beta;
		xv[i+3] += beta;
	}
	for(;i<size;i++) {
		xv[i] += beta;
	}

	return( 0 );
}

/* z = alpha*x.*y + z */
int vvmpy( real alpha, vector x, vector y, vector z ) {
	int i,m;
	int size = z.size;
	real* xv = x.data;
	real* yv = y.data;
	real* zv = z.data;
	real  x0,x1,x2,x3;

	if( (x.size==1) and (size==y.size) ) {
		vvadd( alpha*x.data[0], y, z );
		return( 0 );
	}
	if( size != y.size ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		zv[i]   += x0*yv[i];
		zv[i+1] += x1*yv[i+1];
		zv[i+2] += x2*yv[i+2];
		zv[i+3] += x3*yv[i+3];
	}
	for(;i<size;i++) {
		zv[i] += alpha*xv[i]*yv[i];
	}

	return( 0 );
}

/* z = alpha*x.*y + beta*z */
int vvmpyb( real alpha, vector x, vector y, real beta, vector z ) {
	int i,m;
	int size = z.size;
	real* xv = x.data;
	real* yv = y.data;
	real* zv = z.data;
	real  x0,x1,x2,x3;

	if( (x.size==1) and (size==y.size) ) {
		vvaddb( alpha*x.data[0], y, beta, z );
		return( 0 );
	}
	if( size != y.size ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		zv[i]   = x0*yv[i]   + beta*zv[i];
		zv[i+1] = x1*yv[i+1] + beta*zv[i+1];
		zv[i+2] = x2*yv[i+2] + beta*zv[i+2];
		zv[i+3] = x3*yv[i+3] + beta*zv[i+3];
	}
	for(;i<size;i++) {
		zv[i] = alpha*xv[i]*yv[i] + beta*zv[i];
	}

	return( 0 );
}

/* z = alpha*x./y + z */
int vvdiv( real alpha, vector x, vector y, vector z ) {
	int i,m;
	int size = z.size;
	real* xv = x.data;
	real* yv = y.data;
	real* zv = z.data;
	real  x0,x1,x2,x3;

	if( (size!=x.size) or (size!=y.size) ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		zv[i]   += x0/yv[i];
		zv[i+1] += x1/yv[i+1];
		zv[i+2] += x2/yv[i+2];
		zv[i+3] += x3/yv[i+3];
	}
	for(;i<size;i++) {
		zv[i] += alpha*xv[i]/yv[i];
	}

	return( 0 );
}

/* z = alpha*x./y + beta*z */
int vvdivb( real alpha, vector x, vector y, real beta, vector z ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real* yv = y.data;
	real* zv = z.data;
	real  x0,x1,x2,x3;

	if( (size!=y.size) or (size!=z.size) ) {
		return( -1 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		zv[i]   = x0/yv[i]   + beta*zv[i];
		zv[i+1] = x1/yv[i+1] + beta*zv[i+1];
		zv[i+2] = x2/yv[i+2] + beta*zv[i+2];
		zv[i+3] = x3/yv[i+3] + beta*zv[i+3];
	}
	for(;i<size;i++) {
		zv[i] = alpha*xv[i]/yv[i] + beta*zv[i];
	}

	return( 0 );
}

/* y = alpha*x + y */
int vvadd( real alpha, vector x, vector y ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real* yv = y.data;
	real  x0,x1,x2,x3;

	if( size != y.size ) {
		return( -1 );
	}
	if( alpha == 0.0 ) {
		return( 0 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		yv[i]   += x0;
		yv[i+1] += x1;
		yv[i+2] += x2;
		yv[i+3] += x3;
	}
	for(;i<size;i++) {
		yv[i] += alpha*xv[i];
	}

	return( 0 );
}

/* y = alpha*x + beta*y */
int vvaddb( real alpha, vector x, real beta, vector y ) {
	int i,m;
	int size = x.size;
	real* xv = x.data;
	real* yv = y.data;
	real  x0,x1,x2,x3;

	if( size != y.size ) {
		return( -1 );
	}
	if( alpha == 0.0 ) {
		vscale( beta, y );
		return( 0 );
	}

	m = size - size%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*xv[i];
		x1 = alpha*xv[i+1];
		x2 = alpha*xv[i+2];
		x3 = alpha*xv[i+3];
		yv[i]   = x0 + beta*yv[i];
		yv[i+1] = x1 + beta*yv[i+1];
		yv[i+2] = x2 + beta*yv[i+2];
		yv[i+3] = x3 + beta*yv[i+3];
	}
	for(;i<size;i++) {
		yv[i] = alpha*xv[i] + beta*yv[i];
	}

	return( 0 );
}

/* vector data files have the following format:			*/
/* first line: "Vb#: size\n" indicates 'V' vector file, 'b' big	*/
/*   endian numbers, size of each number is '#' bytes, the 	*/
/*   vector is 'size' numbers long				*/
/* undefined data until location 128, where the data starts	*/
int vecread( char* filename, vector v1 ) {
	int r = 0;
	int sz;
	char buffer[128];
	FILE* fp;
	char* cp;
	int q;
	int gl_sz = DomainGlobalSize( v1.dtype );

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'V' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(real) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );
	if( sz != gl_sz ) {
		return( -5 );
	}

	/* skip over header info */
	fseek( fp, 128, SEEK_SET );

	r = ReadGlobalData( v1.data, v1.size, TYPE_REAL, 
		v1.dtype, fp );

	fclose( fp );

	q = 0x1234;
	cp = (char*)(&q);
	if( cp[0] == 0x34 ) {
		/* this machine is little endian */
		if( toupper(buffer[1]) != 'L' ) {
			swapendian( v1.data, sizeof(real), v1.size );
		}
	} else {
		/* this machine is big endian */
		if( toupper(buffer[1]) != 'B' ) {
			swapendian( v1.data, sizeof(real), v1.size );
		}
	}

	return( r );
}

int ivecread( char* filename, ivector v1 ) {
	int r = 0;
	int sz;
	char buffer[128];
	FILE* fp;
	int q;
	char* cp;
	int gl_sz = DomainGlobalSize( v1.dtype );

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'I' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(int) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );
	if( sz != gl_sz ) {
		return( -5 );
	}

	/* skip over header info */
	fseek( fp, 128, SEEK_SET );

	r = ReadGlobalData( v1.data, v1.size, TYPE_INT, 
		v1.dtype, fp );

	fclose( fp );

	q = 0x1234;
	cp = (char*)(&q);
	if( cp[0] == 0x34 ) {
		/* this machine is little endian */
		if( toupper(buffer[1]) != 'L' ) {
			swapendian( v1.data, sizeof(int), v1.size );
		}
	} else {
		/* this machine is big endian */
		if( toupper(buffer[1]) != 'B' ) {
			swapendian( v1.data, sizeof(int), v1.size );
		}
	}

	return( r );
}

int bvecread( char* filename, bvector v1 ) {
	int r = 0;
	int sz;
	char buffer[128];
	FILE* fp;
	int q;
	char* cp;
	int gl_sz = DomainGlobalSize( v1.dtype );

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'B' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(byte) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );
	if( sz != gl_sz ) {
		return( -5 );
	}

	/* skip over header info */
	fseek( fp, 128, SEEK_SET );

	r = ReadGlobalData( v1.data, v1.size, TYPE_BYTE, 
		v1.dtype, fp );

	fclose( fp );

	q = 0x1234;
	cp = (char*)(&q);
	if( cp[0] == 0x34 ) {
		/* this machine is little endian */
		if( toupper(buffer[1]) != 'L' ) {
			swapendian( v1.data, sizeof(byte), v1.size );
		}
	} else {
		/* this machine is big endian */
		if( toupper(buffer[1]) != 'B' ) {
			swapendian( v1.data, sizeof(byte), v1.size );
		}
	}

	return( r );
}

int vecreadinfo( char* filename, int* ssz ) {
	int   sz;
	FILE* fp;
	char  buffer[128];

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );
	fclose( fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'V' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(real) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );

	*ssz = sz;

	return( 0 );
}

int ivecreadinfo( char* filename, int* ssz ) {
	int   sz;
	FILE* fp;
	char  buffer[128];

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );
	fclose( fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'I' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(int) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );

	*ssz = sz;

	return( 0 );
}

int bvecreadinfo( char* filename, int* ssz ) {
	int   sz;
	FILE* fp;
	char  buffer[128];

	fp = fopen( filename, "r" );
	if( fp == NULL ) {
		return( -1 );
	}
	fgets( buffer, 128, fp );
	fclose( fp );

	/* check text IDs */
	if( toupper(buffer[0]) != 'B' ) {
		/* not a vector file */
		return( -2 );
	}
	if( (buffer[2]-'0') != sizeof(byte) ) {
		/* wrong data-type size - could convert! */
		return( -4 );
	}
	sz = atoi( buffer+4 );

	*ssz = sz;

	return( 0 );
}

int vecwrite( char* filename, vector v1 ) {
	FILE* fp;
	int q,sz;
	char* cp;

	if( SelfPE == 0 ) {
		fp = fopen( filename, "w" );
		if( fp == NULL ) {
			WriteGlobalData( NULL, 0, TYPE_REAL, -1, NULL );
			return( -1 );
		}

		sz = DomainGlobalSize( v1.dtype );

		q = 0x1234;
		cp = (char*)(&q);
		if( cp[0] == 0x34 ) {
			/* this machine is little endian */
			fprintf( fp, "VL%i: %i\n\014\n", 
				sizeof(real), sz );
		} else {
			/* this machine is big endian */
			fprintf( fp, "VB%i: %i\n\014\n", 
				sizeof(real), sz );
		}

		fseek( fp, 128, SEEK_SET );

		WriteGlobalData( v1.data, v1.size, TYPE_REAL, v1.dtype, fp );

		fclose( fp );
	} else {
		WriteGlobalData( v1.data, v1.size, TYPE_REAL, v1.dtype, fp );
	}

	MPI_Barrier( MPI_COMM_WORLD );

	return( 0 );
}

int ivecwrite( char* filename, ivector v1 ) {
	FILE* fp;
	int q,sz;
	char* cp;

	if( SelfPE == 0 ) {
		fp = fopen( filename, "w" );
		if( fp == NULL ) {
			WriteGlobalData( NULL, 0, TYPE_INT, -1, NULL );
			return( -1 );
		}

		sz = DomainGlobalSize( v1.dtype );

		q = 0x1234;
		cp = (char*)(&q);
		if( cp[0] == 0x34 ) {
			/* this machine is little endian */
			fprintf( fp, "IL%i: %i\n\014\n", 
				sizeof(int), sz );
		} else {
			/* this machine is big endian */
			fprintf( fp, "IB%i: %i\n\014\n", 
				sizeof(int), sz );
		}
		fseek( fp, 128, SEEK_SET );

		WriteGlobalData( v1.data, v1.size, TYPE_INT, v1.dtype, fp );

		fclose( fp );
	} else {
		WriteGlobalData( v1.data, v1.size, TYPE_INT, v1.dtype, fp );
	}

	MPI_Barrier( MPI_COMM_WORLD );

	return( 0 );
}

int bvecwrite( char* filename, bvector v1 ) {
	FILE* fp;
	int q,sz;
	char* cp;

	if( SelfPE == 0 ) {
		fp = fopen( filename, "w" );
		if( fp == NULL ) {
			WriteGlobalData( NULL, 0, TYPE_BYTE, -1, NULL );
			return( -1 );
		}

		sz = DomainGlobalSize( v1.dtype );

		q = 0x1234;
		cp = (char*)(&q);
		if( cp[0] == 0x34 ) {
			/* this machine is little endian */
			fprintf( fp, "BL%i: %i\n\014\n", 
				sizeof(byte), sz );
		} else {
			/* this machine is big endian */
			fprintf( fp, "BB%i: %i\n\014\n", 
				sizeof(byte), sz );
		}
		fseek( fp, 128, SEEK_SET );

		WriteGlobalData( v1.data, v1.size, TYPE_BYTE, v1.dtype, fp );

		fclose( fp );
	} else {
		WriteGlobalData( v1.data, v1.size, TYPE_BYTE, v1.dtype, fp );
	}

	MPI_Barrier( MPI_COMM_WORLD );

	return( 0 );
}

void vecprint( vector V ) {
	int i;
	real* vd = V.data;
	MPI_Status mpistat;
	MPI_Request mpireq;
	int NextPE,PrevPE;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* enter a sequential region */
	i = SelfPE;
	if( (NumPEs==1) or ((SequentialFirsttime==true) and (SelfPE==0)) ) {
		/* first time, PE0 passes through */
	} else {
		MPI_Recv( &i, 1, MPI_INT, PrevPE, TAG_SEQ,
			MPI_COMM_WORLD, &mpistat );
	}

	if( SelfPE == 0 ) {
		printf( "%s", prnt_hdr );
	}

	for(i=0;i<V.size;i++) {
		printf( prnt_fmt, vd[i] );
	}

	if( SelfPE == (NumPEs-1) ) {
		printf( "%s", prnt_ftr );
	}

	/* leave the sequential region */
	i = SelfPE;
	MPI_Isend( &i, 1, MPI_INT, NextPE, TAG_SEQ,
		MPI_COMM_WORLD, &mpireq );
	MPI_Request_free( &mpireq );
	SequentialFirsttime = false;

	return;
}

void ivecprint( ivector V ) {
	int i;
	int* vd = V.data;
	MPI_Status mpistat;
	MPI_Request mpireq;
	int NextPE,PrevPE;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* enter a sequential region */
	i = SelfPE;
	if( (NumPEs==1) or ((SequentialFirsttime==true) and (SelfPE==0)) ) {
		/* first time, PE0 passes through */
	} else {
		MPI_Recv( &i, 1, MPI_INT, PrevPE, TAG_SEQ,
			MPI_COMM_WORLD, &mpistat );
	}

	if( SelfPE == 0 ) {
		printf( "%s", prnt_hdr );
	}

	for(i=0;i<V.size;i++) {
		printf( prnt_ifmt, vd[i] );
	}

	if( SelfPE == (NumPEs-1) ) {
		printf( "%s", prnt_ftr );
	}

	/* leave the sequential region */
	i = SelfPE;
	MPI_Isend( &i, 1, MPI_INT, NextPE, TAG_SEQ,
		MPI_COMM_WORLD, &mpireq );
	MPI_Request_free( &mpireq );
	SequentialFirsttime = false;

	return;
}

void bvecprint( bvector V ) {
	int i;
	byte* vd = V.data;
	MPI_Status mpistat;
	MPI_Request mpireq;
	int NextPE,PrevPE;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* enter a sequential region */
	i = SelfPE;
	if( (NumPEs==1) or ((SequentialFirsttime==true) and (SelfPE==0)) ) {
		/* first time, PE0 passes through */
	} else {
		MPI_Recv( &i, 1, MPI_INT, PrevPE, TAG_SEQ,
			MPI_COMM_WORLD, &mpistat );
	}

	if( SelfPE == 0 ) {
		printf( "%s", prnt_hdr );
	}

	for(i=0;i<V.size;i++) {
		printf( prnt_bfmt, vd[i] );
	}

	if( SelfPE == (NumPEs-1) ) {
		printf( "%s", prnt_ftr );
	}

	/* leave the sequential region */
	i = SelfPE;
	MPI_Isend( &i, 1, MPI_INT, NextPE, TAG_SEQ,
		MPI_COMM_WORLD, &mpireq );
	MPI_Request_free( &mpireq );
	SequentialFirsttime = false;

	return;
}

void swapendian( void* ptr, int dsize, int num ) {
	char* cp = ptr;
	char b0,b1,b2,b3,b4,b5,b6,b7;
	int i,j;

	switch( dsize ) {
		case 8:
			for(i=0;i<num;i++) {
				b0 = cp[0];
				b1 = cp[1];
				b2 = cp[2];
				b3 = cp[3];
				b4 = cp[4];
				b5 = cp[5];
				b6 = cp[6];
				b7 = cp[7];
				cp[0] = b7;
				cp[1] = b6;
				cp[2] = b5;
				cp[3] = b4;
				cp[4] = b3;
				cp[5] = b2;
				cp[6] = b1;
				cp[7] = b0;
				cp += 8;
			}
			break;
		case 4:
			for(i=0;i<num;i++) {
				b0 = cp[0];
				b1 = cp[1];
				b2 = cp[2];
				b3 = cp[3];
				cp[0] = b3;
				cp[1] = b2;
				cp[2] = b1;
				cp[3] = b0;
				cp += 4;
			}
			break;
		case 2:
			for(i=0;i<num;i++) {
				b0 = cp[0];
				b1 = cp[1];
				cp[0] = b1;
				cp[1] = b0;
				cp += 2;
			}
			break;
		case 1:
			/* no work to do */
			break;
	}

	return;
}


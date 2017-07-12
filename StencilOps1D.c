/*	Basic Sparse matrix operations for C	*/

/*	John Pormann, Duke University */
/*		jbp1@duke.edu             */

/*	jcoef: Ndim [xdim ydim zdim] NumBCs [bc-list] [bc-flags] */
/*		0-D:  0 0			*/
/*		1-D:  1 xd 0		*/
/*		2-D:  2 xd yd 0		*/
/*		3-D:  3 xd yd zd 0	*/

/*	0-D is essentially just a diagonal matrix */

/* Copyright John B. Pormann, 2000-2007, all rights reserved */



#include "CardioWave.h"

int GetBandwidth1D( domain_t d );

char* SparseOps_RCSID = "$Id: StencilOps1D.c 14 2007-05-11 14:57:55Z jbp $";

int sprmpy0( real alpha, sparse A, vector x, real beta, vector y );
int sprmpy1( real alpha, sparse A, vector x, real beta, vector y );
int sprmpy2( real alpha, sparse A, vector x, real beta, vector y );
int sprmpy3( real alpha, sparse A, vector x, real beta, vector y );

int coefsize( int ndim ) {
	int r = 0;
	switch( ndim ) {
		case 3: r = 27;
			break;
		case 2: r = 9;
			break;
		case 1: r = 3;
			break;
		case 0: r = 1;
			break;
	}
	return( r );
}

/* allocate a sparse matrix of the given type */
int spralloc( sparse* M, domain_t dtp, int type, int rows,
		int cols, int maxnz ) {
	if( maxnz < 1024 ) {
		maxnz = 1024;
	}

	M->type   = type;
	M->rows   = rows;
	M->cols   = cols;
	M->maxnz  = maxnz;
	M->csep   = 0;
	M->dtype  = dtp;
	/* grab 2 message tags */
	M->msgtag = NextMsgTag();
	NextMsgTag();

	/* jcoef array is of size = maxnz */
	M->jcoef = (int*)malloc( maxnz*sizeof(int) );
	if( M->jcoef == NULL ) {
		return( -1 );
	}
	memset( (void*)(M->jcoef), 0, maxnz*sizeof(int) );

	/* coef array is of size = maxnz */
	/* : caller should set maxnz >= 3^Ndim */
	M->coef  = (real*)malloc( maxnz*sizeof(real) );
	if( M->coef == NULL ) {
		return( -2 );
	}
	memset( (void*)(M->coef), 0, maxnz*sizeof(real) );

	return( 0 );
}

int sprfree( sparse* M ) {
	if( M->jcoef != NULL ) {
		free( M->jcoef );
	}
	if( M->coef != NULL ) {
		free( M->coef );
	}

	M->type   = 0;
	M->rows   = 0;
	M->cols   = 0;
	M->maxnz  = 0;
	M->jcoef  = NULL;
	M->coef   = NULL;
	M->dtype  = -1;
	M->csep   = 0;
	M->msgtag = MPI_ANY_TAG;

	return( 0 );
}

int sswap( sparse* Ap, sparse* Bp ) {
	sparse C;

	C   = *Ap;
	*Ap = *Bp;
	*Bp = C;

	return( 0 );
}

int scopy( sparse A, sparse B ) {
	if( (A.type!=B.type) or (A.rows!=B.rows) or (A.cols!=B.cols)
	    or (A.maxnz!=B.maxnz) ) {
		return( -1 );
	}

	memcpy( B.jcoef, A.jcoef, A.maxnz*sizeof(int) );
	memcpy( B.coef, A.coef, A.maxnz*sizeof(real) );

	return( 0 );
}

int szero( sparse A ) {
	int sz;

	sz = coefsize( A.jcoef[0] );
	memset( (void*)(A.coef), 0, sz*sizeof(real) );

	return( 0 );
}

int sscale( real beta, sparse M1 ) {
	int i,m,sz;
	real* av = M1.coef;

	sz = coefsize( M1.jcoef[0] );
	m  = sz - sz%4;

	for(i=0;i<m;i+=4) {
		av[i]   *= beta;
		av[i+1] *= beta;
		av[i+2] *= beta;
		av[i+3] *= beta;
	}
	for(;i<sz;i++) {
		av[i] *= beta;
	}

	return( 0 );
}

int ssadd( real alpha, sparse A, sparse B ) {
	int i;
	int m,sz;
	int rrr  = A.rows;
	int ccc  = A.cols;
	real* av = A.coef;
	real* bv = B.coef;
	real  x0,x1,x2,x3;

	if( (rrr!=B.rows) or (ccc!=B.cols) or (A.type!=B.type) ) {
		return( -1 );
	}

	sz = coefsize( A.jcoef[0] );
	m  = sz - sz%4;

	for(i=0;i<m;i+=4) {
		x0 = alpha*av[i];
		x1 = alpha*av[i+1];
		x2 = alpha*av[i+2];
		x3 = alpha*av[i+3];
		bv[i]   += x0;
		bv[i+1] += x1;
		bv[i+2] += x2;
		bv[i+3] += x3;
	}
	for(;i<sz;i++) {
		bv[i] += alpha*av[i];
	}

	return( 0 );
}

int ssadddiag( real alpha, vector V1, sparse M2 ) {
	/* assume V1 is size=1 */
	switch( M2.jcoef[0] ) {
		case 3: M2.coef[13] += alpha*V1.data[0];
			break;
		case 2: M2.coef[4] += alpha*V1.data[0];
			break;
		case 1: M2.coef[1] += alpha*V1.data[0];
			break;
		case 0: M2.coef[0] += alpha*V1.data[0];
			break;
	}

	return( 0 );
}

void sprnonz( sparse A, int i, int* n, int* lst ) {
	int nd;
	int xd = 0, yd = 0, zd = 0;
	int x  = 0, y  = 0, z  = 0;
	int c;
	int tissuesize;

	/* unpack the dimensions from A.jcoef[] */
	nd = A.jcoef[0];

	tissuesize = DomainGlobalSize( A.dtype );

	/* : note there is fall-thru for case 3 and 2 */
	switch( nd ) {
		case 3: xd = A.jcoef[1];
			yd = A.jcoef[2];
			if( NumPEs > 1 ) {
				zd = tissuesize/(xd*yd);
			} else {
				zd = A.jcoef[3];
			}
			break;
		case 2: xd = A.jcoef[1];
			if( NumPEs > 1 ) {
				yd = tissuesize/xd;
			} else {
				yd = A.jcoef[2];
			}
			break;
		case 1: if( NumPEs > 1 ) {
				xd = tissuesize;
			} else {
				xd = A.jcoef[1];
			}
			break;
	}

	/* get coordinates */
	/* : note there is fall-thru for case 3 and 2 */
	/* ::old code:: c = i + TissueSplit[SelfPE]; */
	c = i + Local2Global( A.dtype, 0 );
	switch( nd ) {
		case 3: z = c/(xd*yd);
			c = c%(xd*yd);
		case 2: y = c/xd;
			c = c%xd;
		case 1: x = c;
	}

	/* store neighbors into lst */
	c = 0;
	lst[c] = i;
	c++;
	switch( nd ) {
		case 3: /* needs work! */
			if( z > 0 ) {
				lst[c] = i-(xd*yd);
				c++;
			}
			if( z < (zd-1) ) {
				lst[c] = i+(xd*yd);
				c++;
			}
		case 2: if( y > 0 ) {
				if( x > 0 ) {
					lst[c] = i-xd-1;
					c++;
				}
				lst[c] = i-xd;
				c++;
				if( x < (xd-1) ) {
					lst[c] = i-xd+1;
					c++;
				}
			}
			if( y < (yd-1) ) {
				if( x > 0 ) {
					lst[c] = i+xd-1;
					c++;
				}
				lst[c] = i+xd;
				c++;
				if( x < (xd-1) ) {
					lst[c] = i+xd+1;
					c++;
				}
			}
		case 1: if( x > 0 ) {
				lst[c] = i-1;
				c++;
			}
			if( x < (xd-1) ) {
				lst[c] = i+1;
				c++;
			}
	}

	*n = c;

	return;
}

int sprput( sparse m1, int i, int j, real v ) {
	int ndim = m1.jcoef[0];
	int xd,yd,zd,x1,x2,y1,y2,z1,z2;
	int xdiff,ydiff,zdiff;
	int n,flag;
	int r = -1;

	switch( ndim ) {
		case 3: xd = m1.jcoef[1];
			yd = m1.jcoef[2];
			zd = m1.jcoef[3];
			z1 = i/(xd*yd);
			y1 = (i/xd)%yd;
			x1 = i%xd;
			z2 = j/(xd*yd);
			y2 = (j/xd)%yd;
			x2 = j%xd;
			zdiff = z2 - z1;
			ydiff = y2 - y1;
			xdiff = x2 - x1;
			n = 13;
			flag = 0;
			if( (xdiff>=-1) and (xdiff<=1) ) {
				n += xdiff;
				flag++;
			}
			if( (ydiff>=-1) and (ydiff<=1) ) {
				n += ydiff*3;
				flag++;
			}
			if( (zdiff>=-1) and (zdiff<=1) ) {
				n += zdiff*9;
				flag++;
			}
			if( flag == 3 ) {
				m1.coef[n] = v;
				r = 0;
			}
			break;
		case 2: xd = m1.jcoef[1];
			yd = m1.jcoef[2];
			y1 = (i/xd)%yd;
			x1 = i%xd;
			y2 = (j/xd)%yd;
			x2 = j%xd;
			ydiff = y2 - y1;
			xdiff = x2 - x1;
			n = 4;
			flag = 0;
			if( (xdiff>=-1) and (xdiff<=1) ) {
				n += xdiff;
				flag++;
			}
			if( (ydiff>=-1) and (ydiff<=1) ) {
				n += ydiff*3;
				flag++;
			}
			if( flag == 2 ) {
				m1.coef[n] = v;
				r = 0;
			}
			break;
		case 1: xdiff = (i-j);
			if( (xdiff>=-1) and (xdiff<=1) ) {
				m1.coef[1+xdiff] = v;
				r = 0;
			}
			break;
		case 0: m1.coef[0] = v;
			r = 0;
			break;
	}

	return( r );
}

real sprget( sparse m1, int i, int j ) {
	int ndim = m1.jcoef[0];
	int diff = (i-j);
	int xd,yd,zd,x1,y1,z1,x2,y2,z2;
	int xdiff,ydiff,zdiff;
	int n,flag;

	switch( ndim ) {
		case 3: xd = m1.jcoef[1];
			yd = m1.jcoef[2];
			zd = m1.jcoef[3];
			z1 = i/(xd*yd);
			y1 = (i/xd)%yd;
			x1 = i%xd;
			z2 = j/(xd*yd);
			y2 = (j/xd)%yd;
			x2 = j%xd;
			zdiff = z2 - z1;
			ydiff = y2 - y1;
			xdiff = x2 - x1;
			n = 13;
			flag = 0;
			if( (xdiff>=-1) and (xdiff<=1) ) {
				n += xdiff;
				flag++;
			}
			if( (ydiff>=-1) and (ydiff<=1) ) {
				n += ydiff*3;
				flag++;
			}
			if( (zdiff>=-1) and (zdiff<=1) ) {
				n += zdiff*9;
				flag++;
			}
			if( flag == 3 ) {
				return( m1.coef[n] );
			}
			break;
		case 2: xd = m1.jcoef[1];
			yd = m1.jcoef[2];
			y1 = (i/xd)%yd;
			x1 = i%xd;
			y2 = (j/xd)%yd;
			x2 = j%xd;
			ydiff = y2 - y1;
			xdiff = x2 - x1;
			n = 4;
			flag = 0;
			if( (xdiff>=-1) and (xdiff<=1) ) {
				n += xdiff;
				flag++;
			}
			if( (ydiff>=-1) and (ydiff<=1) ) {
				n += ydiff*3;
				flag++;
			}
			if( flag == 2 ) {
				return( m1.coef[n] );
			}
			break;
		case 1: xdiff = (i-j);
			if( xdiff == 0 ) {
				return( m1.coef[1] );
			} else if( diff == 1 ) {
				return( m1.coef[2] );
			} else if( diff == -1 ) {
				return( m1.coef[0] );
			}
			break;
		case 0: return( m1.coef[0] );
			break;
	}

	return( 0.0 );
}

int sprvec( real alpha, sparse A, vector x, real beta, vector y ) {
	int ndim;
	int rtn = 0;

	ndim = A.jcoef[0];
	switch( ndim ) {
		case 0: rtn = sprmpy0( alpha, A, x, beta, y );
			break;
		case 1: rtn = sprmpy1( alpha, A, x, beta, y );
			break;
		case 2: rtn = sprmpy2( alpha, A, x, beta, y );
			break;
		case 3: rtn = sprmpy3( alpha, A, x, beta, y );
			break;
		default:
			rtn = -1;			
	}

	return( rtn );
}

int sprmpy0( real alpha, sparse A, vector x, real beta, vector y ) {
	int i,r,n;
	real* yv = y.data;
	real* xv = x.data;
	int ne = A.jcoef[1];

	if( ne > 0 ) {
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+2];
			/* save old values of the x/y-vec at the bc location */
			Workspace.data[2*i]   = xv[n];
			Workspace.data[2*i+1] = yv[n];
			if( A.jcoef[2+ne+i]&(_BC_COL|_BC_FULL) ) {
					xv[n] = 0.0;
					break;
			}
		}
	}

	r = vvaddb( alpha*A.coef[0], x, beta, y );

	if( ne > 0 ) {
		/* take care of DBCs */
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+2];
			if( A.jcoef[2+ne+i]&_BC_ROW ) {
				/* restore y-vec (undo calculations) */
				yv[n] = (1.0+beta)*Workspace.data[2*i+1];
			} else if( A.jcoef[2+ne+i]&_BC_COL ) {
				/* restore x-vec (had been set to 0) */
				xv[n] = Workspace.data[2*i];
			} else if( A.jcoef[2+ne+i]&_BC_FULL ) {
				/* account for the 1 on the diagonal */
				xv[n] = Workspace.data[2*i];
				yv[n] = alpha*xv[n] + beta*Workspace.data[2*i+1];
			}
		}
	}

	return( r );
}

int sprmpy1( real alpha, sparse A, vector x, real beta, vector y ) {
	int   i,n;
	int   sz = A.rows;
	int   xd;
	int   ne = A.jcoef[2];
	real* xv = x.data;
	real* yv = y.data;
	real  mp0,mp1,mp2;
	MPI_Request recvreqU,recvreqD;
	MPI_Request sendreqU,sendreqD;
	MPI_Status  recvstatU,recvstatD;
	int NextPE,PrevPE;
	int wkofs;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* note that global x-dim is stored in jcoef[1] */
	/* but we want local x-dim, so use sz instead   */
	xd = sz;
	
	/* where should comm-buffers start in workspace? */
	/* : we need space for original values of x-vec where DBCs are */
	wkofs = 2*ne;

	if( ne > 0 ) {
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+3];
			/* save old values of the x/y-vec at the bc location */
			Workspace.data[3*i]   = xv[n];
			Workspace.data[3*i+1] = yv[n];
			if( A.jcoef[3+ne+i]&(_BC_COL|_BC_FULL) ) {
					xv[n] = 0.0;
					break;
			}
		}
	}

	if( SelfPE != 0 ) {
		MPI_Irecv( Workspace.data+wkofs, 1, MPI_SSREAL, PrevPE, A.msgtag,
			MPI_COMM_WORLD, &recvreqD );
		MPI_Isend( xv, 1, MPI_SSREAL, PrevPE, A.msgtag+1,
			MPI_COMM_WORLD, &sendreqU );
	}
	if( SelfPE != (NumPEs-1) ) {
		MPI_Irecv( Workspace.data+wkofs+1, 1, MPI_SSREAL, NextPE, 
			A.msgtag+1, MPI_COMM_WORLD, &recvreqU );
		MPI_Isend( xv+sz-1, 1, MPI_SSREAL, NextPE, A.msgtag,
			MPI_COMM_WORLD, &sendreqD );
	}

	mp0 = alpha*A.coef[0];
	mp1 = alpha*A.coef[1];
	mp2 = alpha*A.coef[2];

	for(i=1;i<(xd-1);i++) {
		yv[i] = mp0*xv[i-1] + mp1*xv[i] + mp2*xv[i+1] + beta*yv[i];
	}

	/* do the boundaries - either by comm or with NBCs */
	if( SelfPE == 0 ) {
		yv[0] = (mp0+mp1)*xv[0] + mp2*xv[1] + beta*yv[0];
	} else {
		MPI_Wait( &recvreqD, &recvstatD );
		MPI_Request_free( &sendreqU );
		yv[0] = mp0*Workspace.data[wkofs+0] + mp1*xv[0] + mp2*xv[1]
			+ beta*yv[0];
	}
	if( SelfPE == (NumPEs-1) ) {
		i = xd-1;
		yv[i] = mp0*xv[i-1] + (mp1+mp2)*xv[i] + beta*yv[i];
	} else {
		MPI_Wait( &recvreqU, &recvstatU );
		MPI_Request_free( &sendreqD );
		i = xd-1;
		yv[i] = mp0*xv[i-1] + mp1*xv[i] + mp2*Workspace.data[wkofs+1]
			+ beta*yv[i];
	}

	if( ne > 0 ) {
		/* take care of DBCs */
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+3];
			if( A.jcoef[3+ne+i]&_BC_ROW ) {
				/* restore y-vec (undo calculations) */
				yv[n] = (1.0+beta)*Workspace.data[2*i+1];
			} else if( A.jcoef[3+ne+i]&_BC_COL ) {
				/* restore x-vec (had been set to 0) */
				xv[n] = Workspace.data[2*i];
			} else if( A.jcoef[3+ne+i]&_BC_FULL ) {
				/* account for the 1 on the diagonal */
				xv[n] = Workspace.data[2*i];
				yv[n] = alpha*xv[n] + beta*Workspace.data[2*i+1];
			}
		}
	}

	return( 0 );
}

int sprmpy2( real alpha, sparse A, vector x, real beta, vector y ) {
	int   i,j,k,o;
	int   sz = A.rows;
	int   xd = A.jcoef[1];
	int   yd;
	int   ne = A.jcoef[3];
	real* xv = x.data;
	real* yv = y.data;
	real* zv = NULL;
	real  mp0,mp1,mp2,mp3,mp4,mp5,mp6,mp7,mp8,mpx;
	MPI_Request recvreqU,recvreqD;
	MPI_Request sendreqU,sendreqD;
	MPI_Status  recvstatU,recvstatD;
	int NextPE,PrevPE;
	int wkofs,n;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* note that global y-dim is stored in jcoef[2] */
	/* but we want local y-dim = global-ydim/NumPEs */
	/* also local y-dim = sz/xd                     */
	yd = sz/xd;

	/* where should comm-buffers start in workspace? */
	/* : we need space for original values of x-vec where DBCs are */
	wkofs = 2*ne;

	if( ne > 0 ) {
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+4];
			/* save old values of the x/y-vec at the bc location */
			Workspace.data[4*i]   = xv[n];
			Workspace.data[4*i+1] = yv[n];
			if( A.jcoef[4+ne+i]&(_BC_COL|_BC_FULL) ) {
					xv[n] = 0.0;
					break;
			}
		}
	}

	if( SelfPE != 0 ) {
		MPI_Irecv( Workspace.data+wkofs, xd, MPI_SSREAL, PrevPE, A.msgtag,
			MPI_COMM_WORLD, &recvreqD );
		MPI_Isend( xv, xd, MPI_SSREAL, PrevPE, A.msgtag+1,
			MPI_COMM_WORLD, &sendreqU );
	}
	if( SelfPE != (NumPEs-1) ) {
		MPI_Irecv( Workspace.data+wkofs+xd, xd, MPI_SSREAL, NextPE, 
			A.msgtag+1, MPI_COMM_WORLD, &recvreqU );
		MPI_Isend( xv+sz-xd, xd, MPI_SSREAL, NextPE, A.msgtag,
			MPI_COMM_WORLD, &sendreqD );
	}

	mp0 = alpha*A.coef[0];
	mp1 = alpha*A.coef[1];
	mp2 = alpha*A.coef[2];
	mp3 = alpha*A.coef[3];
	mp4 = alpha*A.coef[4];
	mp5 = alpha*A.coef[5];
	mp6 = alpha*A.coef[6];
	mp7 = alpha*A.coef[7];
	mp8 = alpha*A.coef[8];

	/* do the bulk of the domain - interior points */
	k = xd;
	for(j=1;j<(yd-1);j++) {
		yv[k] =   beta*yv[k]
			+ mp1*xv[k-xd]            + mp2*xv[k-xd+1] 
			+ (mp0+mp3+mp4+mp6)*xv[k] + mp5*xv[k+1]
			+ mp7*xv[k+xd]            + mp8*xv[k+xd+1];
		k++;
		for(i=1;i<(xd-1);i++) {
			yv[k] =   beta*yv[k]
				+ mp0*xv[k-xd-1] + mp1*xv[k-xd] + mp2*xv[k-xd+1]
				+ mp3*xv[k-1]    + mp4*xv[k]    + mp5*xv[k+1]
				+ mp6*xv[k+xd-1] + mp7*xv[k+xd] + mp8*xv[k+xd+1];
			k++;
		}
		yv[k] =   beta*yv[k]
			+ mp0*xv[k-xd-1] + mp1*xv[k-xd]
			+ mp3*xv[k-1]    + (mp2+mp4+mp5+mp8)*xv[k] 
			+ mp6*xv[k+xd-1] + mp7*xv[k+xd];
		k++;
	}

	/* top edge of domain - either use comm values or NBC */
	if( SelfPE == (NumPEs-1) ) {
		k = sz - xd;
		yv[k] =   beta*yv[k]
			+ mp1*xv[k-xd] + mp2*xv[k-xd+1] 
			+ (mp3+mp4+mp0+mp6+mp7+mp8)*xv[k] + mp5*xv[k+1];
		k++;
		mpx = mp4 + mp6 + mp7 + mp8;
		for(i=1;i<(xd-1);i++) {
			yv[k] =   beta*yv[k]
				+ mp0*xv[k-xd-1] + mp1*xv[k-xd] + mp2*xv[k-xd+1]
				+ mp3*xv[k-1]    + mpx*xv[k]    + mp5*xv[k+1];
			k++;
		}
		yv[k] =   beta*yv[k]
			+ mp0*xv[k-xd-1] + mp1*xv[k-xd]
			+ mp3*xv[k-1]    + (mp4+mp5+mp2+mp6+mp7+mp8)*xv[k];
	} else {
		MPI_Wait( &recvreqU, &recvstatU );
		MPI_Request_free( &sendreqD );
		zv = Workspace.data + wkofs + xd;

		k = sz - xd;
		o = -k;
		yv[k] =   beta*yv[k]
			+ mp1*xv[k-xd]            + mp2*xv[k-xd+1] 
			+ (mp0+mp3+mp4+mp6)*xv[k] + mp5*xv[k+1]
			+ mp7*zv[k+o]             + mp8*zv[k+o+1];
		k++;
		for(i=1;i<(xd-1);i++) {
			yv[k] =   beta*yv[k]
				+ mp0*xv[k-xd-1] + mp1*xv[k-xd] + mp2*xv[k-xd+1]
				+ mp3*xv[k-1]    + mp4*xv[k]    + mp5*xv[k+1]
				+ mp6*zv[k+o-1]  + mp7*zv[k+o]  + mp8*zv[k+o+1];
			k++;
		}
		yv[k] =   beta*yv[k]
			+ mp0*xv[k-xd-1] + mp1*xv[k-xd]
			+ mp3*xv[k-1]    + (mp2+mp4+mp5+mp8)*xv[k] 
			+ mp6*zv[k+o-1]  + mp7*zv[k+o];
	}
	/* bottom edge of domain - either use comm values or NBC */
	if( SelfPE == 0 ) {
		k = 0;
		yv[k] =   beta*yv[k]
			+ (mp0+mp1+mp2+mp3+mp4+mp6)*xv[k] + mp5*xv[k+1]
			+ mp7*xv[k+xd]                    + mp8*xv[k+xd+1];
		k++;
		mpx = mp4 + mp0 + mp1 + mp2;
		for(i=1;i<(xd-1);i++) {
			yv[k] =   beta*yv[k]
				+ mp3*xv[k-1]    + mpx*xv[k]    + mp5*xv[k+1]
				+ mp6*xv[k+xd-1] + mp7*xv[k+xd] + mp8*xv[k+xd+1];
			k++;
		}
		yv[k] =   beta*yv[k]
			+ mp3*xv[k-1]    + (mp0+mp1+mp2+mp4+mp5+mp8)*xv[k] 
			+ mp6*xv[k+xd-1] + mp7*xv[k+xd];
	} else {
		MPI_Wait( &recvreqD, &recvstatD );
		MPI_Request_free( &sendreqU );
		zv = Workspace.data + wkofs;

		k = 0;
		yv[k] =   beta*yv[k]
			+ mp1*zv[k]               + mp2*zv[k+1] 
			+ (mp0+mp3+mp4+mp6)*xv[k] + mp5*xv[k+1]
			+ mp7*xv[k+xd]            + mp8*xv[k+xd+1];
		k++;
		for(i=1;i<(xd-1);i++) {
			yv[k] =   beta*yv[k]
				+ mp0*zv[k-1]    + mp1*zv[k]    + mp2*zv[k+1]
				+ mp3*xv[k-1]    + mp4*xv[k]    + mp5*xv[k+1]
				+ mp6*xv[k+xd-1] + mp7*xv[k+xd] + mp8*xv[k+xd+1];
			k++;
		}
		yv[k] =   beta*yv[k]
			+ mp0*zv[k-1]    + mp1*zv[k]
			+ mp3*xv[k-1]    + (mp2+mp4+mp5+mp8)*xv[k] 
			+ mp6*xv[k+xd-1] + mp7*xv[k+xd];
	}

	if( ne > 0 ) {
		/* take care of DBCs */
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+4];
			if( A.jcoef[4+ne+i]&_BC_ROW ) {
				/* restore y-vec (undo calculations) */
				yv[n] = (1.0+beta)*Workspace.data[2*i+1];
			} else if( A.jcoef[4+ne+i]&_BC_COL ) {
				/* restore x-vec (had been set to 0) */
				xv[n] = Workspace.data[2*i];
			} else if( A.jcoef[4+ne+i]&_BC_FULL ) {
				/* account for the 1 on the diagonal */
				xv[n] = Workspace.data[2*i];
				yv[n] = alpha*xv[n] + beta*Workspace.data[2*i+1];
			}
		}
	}

	return( 0 );
}

int sprmpy3( real alpha, sparse A, vector x, real beta, vector y ) {
	int   i,j,k,m,o;
	real* xv = x.data;
	real* yv = y.data;
	real* zv;
	int   xd = A.jcoef[1];
	int   yd = A.jcoef[2];
	int   zd;
	int   ne = A.jcoef[4];
	int   xy;
	int   sz = A.rows;
	real  mpx,mpy,mpz,mpd,mpdd;
	MPI_Request recvreqU,recvreqD;
	MPI_Request sendreqU,sendreqD;
	MPI_Status  recvstatU,recvstatD;
	int NextPE,PrevPE;
	int pp,wkofs,n;
	static int counter = 0;

	counter++;

	/* compute z as sz/(xd*yd) */
	xy = xd*yd;
	zd = sz/xy;

	NextPE = (SelfPE+1)%NumPEs;
	PrevPE = (SelfPE+NumPEs-1)%NumPEs;

	/* where should comm-buffers start in workspace? */
	/* : we need space for original values of x-vec where DBCs are */
	wkofs = 2*ne;

	if( ne > 0 ) {
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+5];
			/* save old values of the x/y-vec at the bc location */
			Workspace.data[5*i]   = xv[n];
			Workspace.data[5*i+1] = yv[n];
			if( A.jcoef[5+ne+i]&(_BC_COL|_BC_FULL) ) {
					xv[n] = 0.0;
					break;
			}
		}
	}

	if( SelfPE != 0 ) {
		MPI_Irecv( Workspace.data+wkofs, xy, MPI_SSREAL, PrevPE, A.msgtag,
			MPI_COMM_WORLD, &recvreqD );
		MPI_Isend( xv, xy, MPI_SSREAL, PrevPE, A.msgtag+1,
			MPI_COMM_WORLD, &sendreqU );
/* printf("PE%i: %i : sent %lf %lf %lf\n",SelfPE,counter,xv[0],xv[1],xv[2]);
fflush(stdout); */
	}
	if( SelfPE != (NumPEs-1) ) {
		MPI_Irecv( Workspace.data+wkofs+xy, xy, MPI_SSREAL, NextPE, 
			A.msgtag+1, MPI_COMM_WORLD, &recvreqU );
		MPI_Isend( xv+sz-xy, xy, MPI_SSREAL, NextPE, A.msgtag,
			MPI_COMM_WORLD, &sendreqD );
/* printf("PE%i: %i : sent %lf %lf %lf\n",SelfPE,counter,xv[sz-xy+0],xv[sz-xy+1],xv[sz-xy+2]);
fflush(stdout); */
	}

	mpx = alpha*A.coef[12];
	mpy = alpha*A.coef[10];
	mpz = alpha*A.coef[4];
	mpd = alpha*A.coef[13];

	/* : interior planes */
	m = xy;
	for(k=1;k<(zd-1);k++) {
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m-xy] + xv[m+xy] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] =  beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m+xd]
				+ mpz*( xv[m+xy] + xv[m-xy] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m-xy] + xv[m+xy] );
		m++;
		/* interior of plane */
		mpdd = mpd+mpx;
		for(j=1;j<(yd-1);j++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*xv[m+1]
				+ mpy*( xv[m+xd] + xv[m-xd] )
				+ mpz*( xv[m+xy] + xv[m-xy] );
			m++;
			for(i=1;i<(xd-1);i++) {
				yv[m] = beta*yv[m] + mpd*xv[m]
					+ mpx*( xv[m-1] + xv[m+1] )
					+ mpy*( xv[m-xd] + xv[m+xd] )
					+ mpz*( xv[m+xy] + xv[m-xy] );
				m++;
			}
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*xv[m-1]
				+ mpy*( xv[m-xd] + xv[m+xd] )
				+ mpz*( xv[m+xy] + xv[m-xy] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m+xy] + xv[m-xy] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m-xd] 
				+ mpz*( xv[m+xy] + xv[m-xy] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m+xy] + xv[m-xy] );
		m++;
	}

	/* bottom plane of domain - either use comm values or NBC */
	if( SelfPE == 0 ) {
		m = 0;
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m+xd]
			+ mpz*xv[m+xy];
		m++;
		mpdd = mpd+mpy+mpz;
		for(i=1;i<(xd-1);i++) {
			yv[m] =  beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m+xd]
				+ mpz*xv[m+xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m+xd]
			+ mpz*xv[m+xy];
		m++;
		mpdd = mpd+mpz;
		for(j=1;j<(yd-1);j++) {
			yv[m] = beta*yv[m] + (mpd+mpx+mpz)*xv[m]
				+ mpx*xv[m+1]
				+ mpy*( xv[m+xd] + xv[m-xd] )
				+ mpz*xv[m+xy];
			m++;
			for(i=1;i<(xd-1);i++) {
				yv[m] = beta*yv[m] + mpdd*xv[m]
					+ mpx*( xv[m-1] + xv[m+1] )
					+ mpy*( xv[m-xd] + xv[m+xd] )
					+ mpz*xv[m+xy];
				m++;
			}
			yv[m] = beta*yv[m] + (mpd+mpx+mpz)*xv[m]
				+ mpx*xv[m-1]
				+ mpy*( xv[m-xd] + xv[m+xd] )
				+ mpz*xv[m+xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m-xd]
			+ mpz*xv[m+xy];
		m++;
		mpdd = mpd+mpy+mpz;
		for(i=1;i<(xd-1);i++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m-xd] 
				+ mpz*xv[m+xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m-xd]
			+ mpz*xv[m+xy];
		m++;
	} else {
		MPI_Wait( &recvreqD, &recvstatD );
		MPI_Request_free( &sendreqU );
		zv = Workspace.data + wkofs;
/* printf("PE%i: %i : recv %lf %lf %lf\n",SelfPE,counter,zv[0],zv[1],zv[2]);
fflush(stdout); */

		m = 0;
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m+xy] + zv[m] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] =  beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m+xd]
				+ mpz*( xv[m+xy] + zv[m] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m+xy] + zv[m] );
		m++;
		for(j=1;j<(yd-1);j++) {
			yv[m] = beta*yv[m] + (mpd+mpx)*xv[m]
				+ mpx*xv[m+1]
				+ mpy*( xv[m+xd] + xv[m-xd] )
				+ mpz*( xv[m+xy] + zv[m] );
			m++;
			for(i=1;i<(xd-1);i++) {
				yv[m] = beta*yv[m] + mpd*xv[m]
					+ mpx*( xv[m-1] + xv[m+1] )
					+ mpy*( xv[m-xd] + xv[m+xd] )
					+ mpz*( xv[m+xy] + zv[m] );
				m++;
			}
			yv[m] = beta*yv[m] + (mpd+mpx)*xv[m]
				+ mpx*xv[m-1]
				+ mpy*( xv[m-xd] + xv[m+xd] )
				+ mpz*( xv[m+xy] + zv[m] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m+xy] + zv[m] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m-xd] 
				+ mpz*( xv[m+xy] + zv[m] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m+xy] + zv[m] );
		m++;
	}
	if( SelfPE == (NumPEs-1) ) {
		m = sz - xy;
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m+xd]
			+ mpz*xv[m-xy];
		m++;
		mpdd = mpd+mpy+mpz;
		for(i=1;i<(xd-1);i++) {
			yv[m] =  beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m+xd]
				+ mpz*xv[m-xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m+xd]
			+ mpz*xv[m-xy];
		m++;
		mpdd = mpd+mpz;
		for(j=1;j<(yd-1);j++) {
			yv[m] = beta*yv[m] + (mpd+mpx+mpz)*xv[m]
				+ mpx*xv[m+1]
				+ mpy*( xv[m+xd] + xv[m-xd] )
				+ mpz*xv[m-xy];
			m++;
			for(i=1;i<(xd-1);i++) {
				yv[m] = beta*yv[m] + mpdd*xv[m]
					+ mpx*( xv[m-1] + xv[m+1] )
					+ mpy*( xv[m-xd] + xv[m+xd] )
					+ mpz*xv[m-xy];
				m++;
			}
			yv[m] = beta*yv[m] + (mpd+mpx+mpz)*xv[m]
				+ mpx*xv[m-1]
				+ mpy*( xv[m-xd] + xv[m+xd] )
				+ mpz*xv[m-xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m-xd]
			+ mpz*xv[m-xy];
		m++;
		mpdd = mpd+mpy+mpz;
		for(i=1;i<(xd-1);i++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m-xd] 
				+ mpz*xv[m-xy];
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy+mpz)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m-xd]
			+ mpz*xv[m-xy];
		m++;
	} else {
		MPI_Wait( &recvreqU, &recvstatU );
		MPI_Request_free( &sendreqD );
		zv = Workspace.data + wkofs + xy;
/* printf("PE%i: %i : recv %lf %lf %lf\n",SelfPE,counter,zv[0],zv[1],zv[2]);
fflush(stdout); */

		m = sz - xy;
		o = -m;
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m-xy] + zv[m+o] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] =  beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m+xd]
				+ mpz*( xv[m-xy] + zv[m+o] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m+xd]
			+ mpz*( xv[m-xy] + zv[m+o] );
		m++;
		for(j=1;j<(yd-1);j++) {
			yv[m] = beta*yv[m] + (mpd+mpx)*xv[m]
				+ mpx*xv[m+1]
				+ mpy*( xv[m+xd] + xv[m-xd] )
				+ mpz*( xv[m-xy] + zv[m+o] );
			m++;
			for(i=1;i<(xd-1);i++) {
				yv[m] = beta*yv[m] + mpd*xv[m]
					+ mpx*( xv[m-1] + xv[m+1] )
					+ mpy*( xv[m-xd] + xv[m+xd] )
					+ mpz*( xv[m-xy] + zv[m+o] );
				m++;
			}
			yv[m] = beta*yv[m] + (mpd+mpx)*xv[m]
				+ mpx*xv[m-1]
				+ mpy*( xv[m-xd] + xv[m+xd] )
				+ mpz*( xv[m-xy] + zv[m+o] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m+1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m-xy] + zv[m+o] );
		m++;
		mpdd = mpd+mpy;
		for(i=1;i<(xd-1);i++) {
			yv[m] = beta*yv[m] + mpdd*xv[m]
				+ mpx*( xv[m-1] + xv[m+1] )
				+ mpy*xv[m-xd] 
				+ mpz*( xv[m-xy] + zv[m+o] );
			m++;
		}
		yv[m] = beta*yv[m] + (mpd+mpx+mpy)*xv[m]
			+ mpx*xv[m-1]
			+ mpy*xv[m-xd]
			+ mpz*( xv[m-xy] + zv[m+o] );
		m++;
	}

	if( ne > 0 ) {
		/* take care of DBCs */
		for(i=0;i<ne;i++) {
			n = A.jcoef[i+5];
			if( A.jcoef[5+ne+i]&_BC_ROW ) {
				/* restore y-vec (undo calculations) */
				yv[n] = (1.0+beta)*Workspace.data[2*i+1];
			} else if( A.jcoef[5+ne+i]&_BC_COL ) {
				/* restore x-vec (had been set to 0) */
				xv[n] = Workspace.data[2*i];
			} else if( A.jcoef[5+ne+i]&_BC_FULL ) {
				/* account for the 1 on the diagonal */
				xv[n] = Workspace.data[2*i];
				yv[n] = alpha*xv[n] + beta*Workspace.data[2*i+1];
			}
		}
	}

	return( 0 );
}


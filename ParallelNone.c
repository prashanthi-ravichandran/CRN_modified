/* dummy file to handle non-parallel cases - duplicate the functionality */
/* of MPI (Message Passing Interface) even if no MPI library is found    */

/* see the README files for information on how to obtain an MPI */
/* library for free, for a wide range of computers              */

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */



#include "CardioWave.h"

static char* RCSID = "$Id: ParallelNone.c 14 2007-05-11 14:57:55Z jbp $";

static int bandw = -1, bath_bandw = -1;

int InitParallel( char** res ) {
	if( (DebugLevel>0) and (SelfPE==0) ) {
		printf("ParallelNone: single CPU\n");
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("ParallelNone: RCSID: %s\n",RCSID);
	}

	return( 0 );
}

int CreateRegularSplitting( domain_t dtp, int xd, int yd, int zd ) {
	int nd;

	/* dimensionality of the problem */
	nd = (xd>1)+(yd>1)+(zd>1);

	switch( dtp ) {
		case Tissue:
		case Intra:
		case Extra:
			TissueSplit = (int*)malloc( 2*sizeof(int) );
			if( TissueSplit == NULL ) {
				return( -1 );
			}
			TissueSplit[0]   = 0;
			TissueSplit[1]   = xd*yd*zd;
			TissueGlobalSize = xd*yd*zd;
			TissueLocalSize  = xd*yd*zd;
			switch( nd ) {
				case 3: bandw = xd*yd;
					break;
				case 2: bandw = xd;
					break;
				case 1: bandw = 1;
					break;
				default: bandw = 0;
			}
			break;
		case Bath:
			BathSplit = (int*)malloc( 2*sizeof(int) );
			if( BathSplit == NULL ) {
				return( -1 );
			}
			BathSplit[0]   = 0;
			BathSplit[1]   = xd*yd*zd;
			BathGlobalSize = xd*yd*zd;
			BathLocalSize  = xd*yd*zd;
			switch( nd ) {
				case 3: bath_bandw = xd*yd;
					break;
				case 2: bath_bandw = xd;
					break;
				case 1: bath_bandw = 1;
					break;
				default: bath_bandw = 0;
			}
			break;
	}

	return( 0 );
}

int CreateIrregularSplitting( domain_t dtp, int sz, int bw ) {
	switch( dtp ) {
		case Tissue:
		case Intra:
		case Extra:
			TissueSplit = (int*)malloc( 2*sizeof(int) );
			if( TissueSplit == NULL ) {
				return( -1 );
			}
			TissueSplit[0]   = 0;
			TissueSplit[1]   = sz;
			TissueGlobalSize = sz;
			TissueLocalSize  = sz;
			bandw = bw;
			break;
		case Bath:
			BathSplit = (int*)malloc( 2*sizeof(int) );
			if( BathSplit == NULL ) {
				return( -1 );
			}
			BathSplit[0]   = 0;
			BathSplit[1]   = sz;
			BathGlobalSize = sz;
			BathLocalSize  = sz;
			bath_bandw = bw;
			break;
	}

	return( 0 );
}

int GetBandwidth1D( domain_t d ) {
	switch( d ) {
		case Tissue:
		case Intra:
		case Extra:
			return( bandw );
			break;
		case Bath:
			return( bath_bandw );
			break;
	}
	return( -1 );
}

/* this is a utility routine to write data from all PEs to 1 file */
/* : all data is collected by PE0, then written by PE0 only       */
int WriteGlobalData( void* data, int num, datatype_t dtp, 
		domain_t dt, FILE* fp ) {
	int r = 0;
	int dsz;

	if( (data==NULL) or (fp==NULL) or (num==0) ) {
		return( -1 );
	}

	switch( dtp ) {
		case TYPE_BYTE:
			dsz = 1;
			break;
		case TYPE_INT:
			dsz = sizeof(int);
			break;
		case TYPE_REAL:
			dsz = sizeof(real);
			break;
	}

	fwrite( data, dsz, num, fp );
	r = 0;

	return( r );
}

/* this assumes the file has already been 'fseek'ed to the */
/* start of the first piece of global data                 */
/* - eg. for vector files, already fseek'd past the 128    */
/*   byte header area                                      */
int ReadGlobalData( void* data, int num, datatype_t dtp, 
		domain_t dt, FILE* fp ) {
	int r = 0;
	int dsz;

	if( (data==NULL) or (fp==NULL) or (num==0) ) {
		return( -1 );
	}

	switch( dtp ) {
		case TYPE_BYTE:
			dsz = 1;
			break;
		case TYPE_INT:
			dsz = sizeof(int);
			break;
		case TYPE_REAL:
			dsz = sizeof(real);
			break;
	}

	fread( data, dsz, num, fp );
	r = 0;

	return( r );
}

int Global2Local( domain_t dt, int g ) {
	return( g );
}

int Local2Global( domain_t dt, int l ) {
	return( l );
}

int OwnerOf( domain_t dt, int g ) {
	return( 0 );
}

int OwnerOfLocal( domain_t dt, int g ) {
	return( 0 );
}

int OwnerOfGlobal( domain_t dt, int g ) {
	return( 0 );
}

/* note that for Ve-only, nblks=1 */
/*   for Vi-Ve, nblks=2; for Ve-Vb, nblks=2 */
/*   for Vi-Ve-Vb, nblks=3, but we don't use the (1,3) & (3,1) corners */
int MatAlloc( bigmatrix* M, int nblks, int* blkspl, domain_t* dtypes,
		int type, int mnz, int mnz2 ) {
	int r,i,j,c;
	int localsz;

	localsz = blkspl[nblks];
        M->Rows = localsz;
        M->Cols = localsz;
	M->NumBlocks = nblks;
	M->BlockSplit = blkspl;
	M->MatrixInfo = (sparse*)malloc( nblks*nblks*sizeof(sparse) );
	if( M->MatrixInfo == NULL ) {
		return( -1 );
	}

	c = 0;
	for(i=0;i<nblks;i++) {
	  for(j=0;j<nblks;j++) {
		if( i == j ) {
			r = spralloc( &(M->MatrixInfo[c]), dtypes[i], type, 
				blkspl[i+1]-blkspl[i], blkspl[j+1]-blkspl[j], 
				mnz );
		} else if( abs(i-j) == 1 ) {
			r = spralloc( &(M->MatrixInfo[c]), dtypes[i], -type,
				blkspl[i+1]-blkspl[i], blkspl[j+1]-blkspl[j], 
				mnz2 );
		} else {
			M->MatrixInfo[c].type = _ZERO;
		}
		if( r < 0 ) {
			break;
		}
		c++;
	  }
	  if( r < 0 ) {
		break;
	  }
	}

	return( r );
}

int MatCopy( bigmatrix A, bigmatrix B ) {
	int r = 0,i;

	if( (A.Rows!=B.Rows) 
	    or (A.Cols!=B.Cols) 
	    or (A.NumBlocks!=B.NumBlocks) ) {
		return( -1 );
	}

	for(i=0;i<(A.NumBlocks*A.NumBlocks);i++) {
		r |= scopy( A.MatrixInfo[i], B.MatrixInfo[i] );
	}

	return( r );
}

int MatScale( real alpha, bigmatrix A ) {
	int r = 0,i;

	for(i=0;i<(A.NumBlocks*A.NumBlocks);i++) {
		r |= sscale( alpha, A.MatrixInfo[i] );
	}

	return( r );
}

int MatAdd( real alpha, bigmatrix A, bigmatrix B ) {
	int r = 0,i;

	for(i=0;i<(A.NumBlocks*A.NumBlocks);i++) {
		r |= ssadd( alpha, A.MatrixInfo[i], B.MatrixInfo[i] );
	}

	return( r );
}

int MatVec( real alpha, bigmatrix A, vector x, real beta, vector y ) {
	int r = 0,i,j,c;
	int nblks = A.NumBlocks;
	int* split = A.BlockSplit;
	vector w,z;

	/* do diagonal blocks first & incorporate beta*y term */
        for(i=0;i<nblks;i++) {
                z.size = split[i+1] - split[i];
                z.data = y.data + split[i];
                z.dtype = A.MatrixInfo[i+i*nblks].dtype;
                w.size = split[i+1] - split[i];
                w.data = x.data + split[i];
                w.dtype = A.MatrixInfo[i+i*nblks].dtype;
                r |= sprvec( alpha, A.MatrixInfo[i+i*nblks], w, beta, z );
        }

	/* now do off-diagonal blocks */
	c = 0;
	for(i=0;i<nblks;i++) {
	  z.size = split[i+1] - split[i];
	  z.data = y.data + split[i];
	  z.dtype = A.MatrixInfo[i].dtype;
	  for(j=0;j<nblks;j++) {
		if( (i!=j) and (A.MatrixInfo[c].type!=_ZERO) ) {
			w.size = split[j+1] - split[j];
			w.data = x.data + split[j];
			w.dtype = A.MatrixInfo[c].dtype;
			r |= sprvec( alpha, A.MatrixInfo[c], w, 1.0, z );
		}
		c++;
	  }
	}

        return( r );
}

real MatGet( bigmatrix A, int i, int j ) {
	int c,k,ii,jj;
	real r = 0.0;
	int* split = A.BlockSplit;

	if( (i<0) or (i>=A.Rows) ) {
		return( r );
	}

	ii = -1;
	jj = -1;
	for(k=0;k<A.NumBlocks;k++) {
		if( (i>=split[k]) and (i<split[k+1]) ) {
			ii = k;
		}
		if( (j>=split[k]) and (j<split[k+1]) ) {
			jj = k;
		}
	}

	c = ii + jj*A.NumBlocks;
	r = sprget( A.MatrixInfo[c], i-split[ii], j-split[jj] );

        return( r );
}

int MatPut( bigmatrix A, int i, int j, real v ) {
	int c,k,ii,jj;
	int r = -1;
	int* split = A.BlockSplit;

	if( (i<0) or (i>=A.Rows) ) {
		return( -1 );
	}

	ii = -1;
	jj = -1;
	for(k=0;k<A.NumBlocks;k++) {
		if( (i>=split[k]) and (i<split[k+1]) ) {
			ii = k;
		}
		if( (j>=split[k]) and (j<split[k+1]) ) {
			jj = k;
		}
	}

	c = ii + jj*A.NumBlocks;
	r = sprput( A.MatrixInfo[c], i-split[ii], j-split[jj], v );

	return( r );
}

/* returns the list of non-zeros for row-i */
/* : uses the Workspace vector's data area */
void MatGetNonzeros( bigmatrix A, int i, int* n, int* lst ) {
	int  n2,ntot,ii,jj,c;

	if( (i<0) or (i>=A.Rows) ) {
		*n   = -1;
		return;
	}

	ntot = 0;
	c = 0;
	for(ii=0;ii<A.NumBlocks;ii++) {
	  for(jj=0;jj<A.NumBlocks;jj++) {
		sprnonz( A.MatrixInfo[c], i, &n2, lst+ntot );
		if( n2 > 0 ) {
			ntot += n2;
		}
		c++;
	  }
	}

	*n = ntot;

	return;
}

/**********************************************************************/
/*                     fake MPI commands follow                       */
/**********************************************************************/

int MPI_Init( int* ac, char*** av ) {
	return( MPI_SUCCESS );
}

int MPI_Finalize( void ) {
	return( MPI_SUCCESS );
}

double MPI_Wtime( void ) {
	return( 0.0 );
}

int MPI_Comm_size( MPI_Comm c, int* np ) {
	*np = 1;
	return( MPI_SUCCESS );
}

int MPI_Comm_rank( MPI_Comm c, int* sp ) {
	*sp = 0;
	return( MPI_SUCCESS );
}

int MPI_Send(void* p, int i, MPI_Datatype d, int j, int k, MPI_Comm c) {
	return( MPI_SUCCESS );
}

int MPI_Recv( void* p, int i, MPI_Datatype d, int j, int k, MPI_Comm c, 
		MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Isend( void* p, int i, MPI_Datatype d, int j, int k, MPI_Comm c, 
		MPI_Request* r ) {
	return( MPI_SUCCESS );
}
int MPI_Issend( void* p, int i, MPI_Datatype d, int j, int k, MPI_Comm c, 
		MPI_Request* r ) {
	return( MPI_SUCCESS );
}
int MPI_Irecv(void* p, int i, MPI_Datatype d, int j, int k, MPI_Comm c, 
		MPI_Request* r ) {
	return( MPI_SUCCESS );
}
int MPI_Wait( MPI_Request* r, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Test( MPI_Request* r, int* i, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Request_free( MPI_Request* r ) {
	return( MPI_SUCCESS );
}
int MPI_Waitany( int i, MPI_Request* r, int* j, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Testany(int i, MPI_Request* r, int* j, int* k, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Waitall(int i, MPI_Request* r, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Testall(int i, MPI_Request* r, int* j, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Waitsome(int i, MPI_Request* r, int* j, int* k, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Testsome(int i, MPI_Request* r, int* j, int* k, MPI_Status* s ) {
	return( MPI_SUCCESS );
}
int MPI_Iprobe(int i, int j, MPI_Comm c, int* flag, MPI_Status* s ) {
	*flag = 0;
	return( MPI_SUCCESS );
}
int MPI_Probe(int i, int j, MPI_Comm c, MPI_Status* s) {
	return( MPI_SUCCESS );
}
int MPI_Cancel(MPI_Request* r) {
	return( MPI_SUCCESS );
}
int MPI_Barrier( MPI_Comm c ) {
	return( MPI_SUCCESS );
}
int MPI_Alltoall(void* p, int num, MPI_Datatype d, void* q, int j, 
		MPI_Datatype e, MPI_Comm c ) {
	int i;
	switch( d ) {
		case MPI_DOUBLE:
			{ double* x = p;
			double* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_FLOAT:
			{ float* x = p;
			float* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_INT:
			{ int* x = p;
			int* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_LONG:
			{ long* x = p;
			long* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
	}

	return( MPI_SUCCESS );
}
int MPI_Allgather( void* p, int num, MPI_Datatype d, void* q, int num2,
		MPI_Datatype e, MPI_Comm c ) {
	int i;
	switch( d ) {
		case MPI_DOUBLE:
			{ double* x = p;
			double* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_FLOAT:
			{ float* x = p;
			float* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_INT:
			{ int* x = p;
			int* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_LONG:
			{ long* x = p;
			long* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
	}

	return( MPI_SUCCESS );
}
int MPI_Allreduce( void* p, void* q, int num, MPI_Datatype d, MPI_Op op, 
		MPI_Comm c ) {
	int i;
	switch( d ) {
		case MPI_DOUBLE:
			{ double* x = p;
			double* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_FLOAT:
			{ float* x = p;
			float* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_INT:
			{ int* x = p;
			int* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_LONG:
			{ long* x = p;
			long* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
	}

	return( MPI_SUCCESS );
}
int MPI_Gather( void* p, int num, MPI_Datatype d, void* q, int num2,
		MPI_Datatype d2, int root, MPI_Comm c ) {
	int i;
	switch( d ) {
		case MPI_DOUBLE:
			{ double* x = p;
			double* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_FLOAT:
			{ float* x = p;
			float* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_INT:
			{ int* x = p;
			int* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_LONG:
			{ long* x = p;
			long* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
	}

	return( MPI_SUCCESS );
}
int MPI_Gatherv( void* p, int num, MPI_Datatype d, void* q, int* num2,
		int* displ, MPI_Datatype d2, int root, MPI_Comm c ) {
	int i;
	switch( d ) {
		case MPI_DOUBLE:
			{ double* x = p;
			double* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_FLOAT:
			{ float* x = p;
			float* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_INT:
			{ int* x = p;
			int* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
		case MPI_LONG:
			{ long* x = p;
			long* y = q;
				for(i=0;i<num;i++) {
					y[i] = x[i]; 
				}
			}
			break;
	}

	return( MPI_SUCCESS );
}
int MPI_Comm_dup( MPI_Comm c1, MPI_Comm* c2 ) {
	*c2 = c1;
	return( MPI_SUCCESS );
}
int MPI_Type_struct( int n, int* dl, MPI_Aint* df, MPI_Datatype *dt, 
	MPI_Datatype* new ) {
	return( MPI_SUCCESS );
}
int MPI_Type_commit( MPI_Datatype* new ) {
	return( MPI_SUCCESS );
}
int MPI_Get_count( MPI_Status* stat, MPI_Datatype dt, int* cnt ) {
	*cnt = 0;
	return( MPI_SUCCESS );
}



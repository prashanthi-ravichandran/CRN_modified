/* a fudged MPI header file for CardioWave */

/* Note that you can get MPICH and LAM/MPI free of charge from: */
/*    MPICH:   http://www-unix.mcs.anl.gov/mpi/mpich/           */
/*    LAM/MPI: http://www.mpi.nd.edu/lam/                       */

/* RCSID: $Id: mpi.h 14 2007-05-11 14:57:55Z jbp $ */


#ifndef _MPI_INCLUDE
#define _MPI_INCLUDE

#define MPI_SUCCESS 0

typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Request;
typedef int MPI_Aint;
typedef int MPI_Comm;

#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_BYTE           ((MPI_Datatype)2)
#define MPI_INT            ((MPI_Datatype)3)
#define MPI_LONG           ((MPI_Datatype)4)
#define MPI_FLOAT          ((MPI_Datatype)5)
#define MPI_DOUBLE         ((MPI_Datatype)6)

#define MPI_SUM            ((MPI_Op)100)
#define MPI_MAX            ((MPI_Op)101)
#define MPI_MIN            ((MPI_Op)102)

#define MPI_COMM_WORLD     ((MPI_Comm)777)

#define MPI_PROC_NULL      (-1)
#define MPI_ANY_SOURCE 	   (-2)
#define MPI_ANY_TAG        (-3)

#define MPI_REQUEST_NULL   (0)

/* 
   Status object.  It is the only user-visible MPI data-structure 
   The "count" field is PRIVATE; use MPI_Get_count to access it. 
 */
typedef struct { 
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Wait(MPI_Request *, MPI_Status *);
int MPI_Test(MPI_Request *, int *, MPI_Status *);
int MPI_Request_free(MPI_Request *);
int MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Waitall(int, MPI_Request *, MPI_Status *);
int MPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Iprobe(int, int, MPI_Comm, int *flag, MPI_Status *);
int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
int MPI_Cancel(MPI_Request *);
int MPI_Test_cancelled(MPI_Status *, int *);
int MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
int MPI_Type_commit(MPI_Datatype *);
int MPI_Type_free(MPI_Datatype *);
int MPI_Barrier(MPI_Comm );
int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
int MPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm); 
int MPI_Gatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, int, MPI_Comm); 
int MPI_Scatter(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int MPI_Scatterv(void* , int *, int *,  MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int MPI_Allgather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int MPI_Allgatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int MPI_Alltoall(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int MPI_Alltoallv(void* , int *, int *, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
int MPI_Comm_size(MPI_Comm, int *);
int MPI_Comm_rank(MPI_Comm, int *);
int MPI_Get_processor_name(char *, int *);
int MPI_Get_version(int *, int *);
double MPI_Wtime(void);
double MPI_Wtick(void);
int MPI_Init(int *, char ***);
int MPI_Init_thread( int *, char ***, int, int * );
int MPI_Finalize(void);
int MPI_Initialized(int *);
int MPI_Abort(MPI_Comm, int);
int MPI_Type_struct( int n, int* dl, MPI_Aint* df, MPI_Datatype *dt, MPI_Datatype* new );
int MPI_Type_commit( MPI_Datatype* new );

#endif


/* CardioWave - main include file */

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */

/* RCSID: $Id: CardioWave.h 14 2007-05-11 14:57:55Z jbp $ */

#ifndef _SIMSYSTEM_H
#define _SIMSYSTEM_H

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <malloc/malloc.h>
#include <memory.h>
#include <time.h>
#include "mpi.h"

#if defined(_SGL)||defined(_SINGLE)
	typedef float real;
	#define MPI_SSREAL MPI_FLOAT
#else
	typedef double real;
	#define MPI_SSREAL MPI_DOUBLE
#endif

/* compute offset of a member in a structure */
#define OffsetOf(type,member) (offsetof(type,member)/sizeof(real))

/* some other handy data types */
typedef enum { false = 0, true = -1, no = 0, yes = -1 } logical;
typedef unsigned char byte;
#define and &&
#define or  ||
#define not !
#define mod %

/* simulation type */
typedef enum {
	Error = -1,
	Monodomain = 1,
	ReducedBidomain = 2,
	FullBidomain = 3,
	ReducedBidomainBath = 4,
	FullBidomainBath = 5
} sim_t;

/* domain type */
typedef enum {
	Undefd = -1,
	Tissue = 0,
	Intracellular = 1,
	Intra = 1,
	Extracellular = 2,
	Extra = 2,
	Bath = 3,
	Solver = 4,
	State = 5,
	StateVar = 5,
	Aux = 6,
	AuxVar = 6,
	Memparam = 7,
	Param = 7
} domain_t;

typedef enum {
	TYPE_BYTE = 1,
	TYPE_CHAR = 1,
	TYPE_INT  = 2,
	TYPE_REAL = 3
} datatype_t;

/* this header is for the full-matrix, vector, and sparse-matrix operations */
#include "VectorOps.h"

/* the big matrix structure - for the solver */
/* : other arrays are size NumComm           */
typedef struct {
        int     Rows,Cols;
	int     NumBlocks;
	int*    BlockSplit;
	sparse* MatrixInfo;
} bigmatrix;
#define _INIT_BIGMATRIX {0,0,0,NULL,NULL}

/* some defined constants */
#define MAX_PROCS     256
#define MAX_RESOURCES 256
#define MAX_FILENAME  128
#define MAX_PATCHSIZE  64
#define MAX_STIM       32
#define MAX_MEMBRANE    8
#define MAX_OFFSETS   256
#define G2L_REMOTE     -1
#define GETJ_ALL        0
#define GETJ_DIDV       1

/* for the byte type, -1 won't work */
#define ByteError     255

/* for the AssertAny/All functions */
typedef enum {
	LOGOP_LT = 1,
	LOGOP_GT = 2,
	LOGOP_EQ = 4,
	LOGOP_LE = 5,
	LOGOP_GE = 6,
	LOGOP_NE = 3
} logop_t;

/* for parsing the resource strings */
typedef struct {
	char* txt;
	int   cmd;
} rword;

/* function prototypes */
int  InitSimulation( char** resources );
void ExitSimulation( vector Vm, vector Vx, vector Q, vector Av );
int  GetNodeInfo( void );
int  ComputeSizes( char** resources );
int  DomainLocalSize( domain_t dt );
int  DomainGlobalSize( domain_t dt );
int  GetInitialValues( vector* Vm, vector* Vx, vector* Q, vector* Av );
int  Checkpoint( real t, vector Vm, vector Vx, vector Q, vector Av );
int  NextMsgTag( void );
real GetMemUsed( void );
/* : defined in SimSystem.c - calls ??_GetF in Membrane models */
int  GetF( real t,real dt, vector Vm,vector Q,vector Fv,vector Fq,vector Av );
/* : for the membrane patches */
int  InitMembrane( char** resources );
void SetPatches( vector Vm, vector Q );
void GetJ( int flag, real dt, real theta, vector Vm, vector Q, vector Gv, 
	vector Gq, vector JJ, vector ZZ );
void GetDQ( real dt, real theta, vector Vm, vector Q, vector Gv, vector Gq, 
	vector Dv, vector Dq );
void ExitMembrane( void );
/* : for the time-steppers */
int  InitTimeStepper( char** res, real* t );
int  Update( real* t, vector Vm, vector Vx, vector Q, vector Av );
void ExitTimeStepper( void );
/* :: defined in time-stepper but called as an interrupt handler */
void CatchSignal( int i );
/* : for the domain model */
int  InitDomain( char** resources );
void ExitDomain( void );
/* : for the stimulus model */
int  InitStimulus( char** resources );
void ApplyStimulus( int IntraOrExtra, real t, real alpha,vector Fv,vector Fq );
void ExitStimulus( void );
/* : for the linear solvers */
int  InitLinearSolver( char** resources );
int  Solve( bigmatrix A, vector x, vector b );
void ExitLinearSolver( void );
/* : for the preconditioners */
int  InitPreconditioner( char** resources );
void ComputeFactors( bigmatrix A );
int  PreconditionSystem( bigmatrix A, vector x, vector b );
void ExitPreconditioner( void );
/* : for the boundary conditions */
int  InitBCs( char** resources );
void ApplyBCmatrix( bigmatrix A );
void ApplyBCvector( real t, bigmatrix A, vector x, vector y );
void ExitBCs( void );
/* : Parallel kernel */
int  InitParallel( char** resources );
int  CreateRegularSplitting( domain_t dt, int xd, int yd, int zd );
int  CreateIrregularSplitting( domain_t dt, int sz, int bw );
int  WriteGlobalData( void* ptr, int num, datatype_t dtp, 
	domain_t d, FILE* fp );
int  ReadGlobalData( void* ptr, int num, datatype_t dtp, 
	domain_t d, FILE* fp );
int  Global2Local( domain_t dt, int gl );
int  Local2Global( domain_t dt, int lc );
int  OwnerOfLocal( domain_t dt, int lc );
int  OwnerOfGlobal( domain_t dt, int gl );
/* defined in Parallel kernels */
int  MatAlloc( bigmatrix* M, int nblks, int* blksz, domain_t* dtps,
	int type, int mnz, int mnz2 );
int  MatCommSetup( bigmatrix M );
int  MatVec( real alpha, bigmatrix A, vector x, real beta, vector y );
int  MatAdd( real alpha, bigmatrix A, bigmatrix B );
int  MatScale( real alpha, bigmatrix A );
int  MatCopy( bigmatrix A, bigmatrix B );
real MatGet( bigmatrix A, int i, int j );
int  MatPut( bigmatrix A, int i, int j, real v );
void MatGetNonzeros( bigmatrix A, int i, int* n, int* lst );
/* :: used in catching signals */
real ProcessSignal( real x );
/* : alias into the Vx vector */
int  AliasVxVi( vector Vx, vector* Vi );
int  AliasVxVe( vector Vx, vector* Ve );
int  AliasVxVb( vector Vx, vector* Vb );

/* helper functions */
int   AppendResources( char** res, int ac, char** av );
int   NumResources( char** res );
int   NumWords( rword* list );
real  GetRealValue( char* txt );
int   GetIntValue( char* txt );
logical GetTFValue( char* txt );
byte  GetByteValue( char* txt );
int   GetNumValues( char* txt );
int*  GetIntArray( char* txt );
real* GetRealArray( char* txt );
char* GetStringValue( char* txt );
int   FindCommand( rword* list, char* res );
int   FindNum( char* txt );
int   RegisterOffset( int ofs, int nt, domain_t dt, char* text );
int   FindOffset( char* text, int* ofs, int* nt, domain_t* dt );

/* trace/debug info */
int  InitDebug( char** res );
void ExitDebug( void );
void CWINT_DebugEnter( char* txt, char* file, const char* func, int line );
void CWINT_DebugLeave( char* txt, char* file, const char* func, int line );
void CWINT_DebugMark( char* txt, char* file, const char* func, int line );
#if defined(NO_DEBUG)
#define DebugEnter(X)
#define DebugLeave(X)
#define DebugMark(X)
#else
#define DebugEnter(X) CWINT_DebugEnter(X,__FILE__,__FUNCTION__,__LINE__)
#define DebugLeave(X) CWINT_DebugLeave(X,__FILE__,__FUNCTION__,__LINE__)
#define DebugMark(X) CWINT_DebugMark(X,__FILE__,__FUNCTION__,__LINE__)
#endif
void CWINT_Assert( int x, char* file, const char* func, int line );
void CWINT_AssertAny( vector x, logop_t lop, real value, char* file, const char* func, int line );
void CWINT_AssertAll( vector x, logop_t lop, real value, char* file, const char* func, int line );
#if defined(TEST_ASSERT)
#define Assert(X) CWINT_Assert(X,__FILE__,__FUNCTION__,__LINE__)
#define AssertAny(X,Y,Z) CWINT_AssertAny(X,Y,Z,__FILE__,__FUNCTION__,__LINE__)
#define AssertAll(X,Y,Z) CWINT_AssertAll(X,Y,Z,__FILE__,__FUNCTION__,__LINE__)
#else
#define Assert(X)
#define AssertAny(X,Y,Z)
#define AssertAll(X,Y,Z)
#endif

/* global variables */
/* : defined in SimSystem.c */
extern int     NumPEs,SelfPE,NextPE,PrevPE;
extern int     NumMem;
/* : MaxComm is defined in Parallel module */
extern int     MaxComm;
/* : filled in by Mem modules */
extern int     PatchSize[MAX_MEMBRANE];
extern logical UseAuxvars;
extern int     AuxiliarySize[MAX_MEMBRANE];
extern int     MemParamSize[MAX_MEMBRANE];
extern int     (*FunctionTable[MAX_MEMBRANE])();
extern void    (*SetpatchTable[MAX_MEMBRANE])();
/* : in SimSystem kernel */
extern int     DebugLevel;
extern logical ShowVersion;
extern logical LoadState,SaveState;
extern FILE*   FpResources;
extern char*   LoadStateFilename;
extern char*   SaveStateFilename;
extern real    Beta,Cm;
extern vector  Workspace;
/* : defined in the Time-stepper */
extern sim_t   SimType;
extern logical NeedVolume, NeedInvVolume;
/* : defined in the kernel, but initialized by Grid routines */
extern vector  Volume;
extern vector  InvVolume;
extern vector  MemParams;
extern bvector NodeType;
extern int     TissueGlobalSize,TissueLocalSize;
extern int     BathGlobalSize,BathLocalSize;
extern int     StateGlobalSize,StateLocalSize;
extern int     AuxiliaryGlobalSize,AuxiliaryLocalSize;
extern int     SolverGlobalSize,SolverLocalSize;
extern int     MemparamGlobalSize,MemparamLocalSize;
extern int*    TissueSplit;
extern int*    BathSplit;
extern int*    SolverSplit;
extern int*    StateSplit;
extern int*    AuxiliarySplit;
extern int*    MemparamSplit;
extern sparse  MatrixINT,MatrixEXT,MatrixBATH;
extern sparse  MatrixE2B,MatrixB2E;

#endif

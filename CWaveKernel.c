/* CardioWave - main kernel routines */

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */



#include "CardioWave.h"
#include <signal.h>

#if defined(CWAVEARCH_Linux)
#  include <unistd.h>
#  include <sys/types.h>
#  define _use_procstat
#elif defined(CWAVEARCH_AIX) || defined(CWAVEARCH_SunOS)
#  include <sys/resource.h>
#  define _use_getrusage
#elif defined(CWAVEARCH_Darwin)
#  include <malloc/malloc.h>
#  define _use_mallocstat
#else
#  include <unistd.h>
#  define _use_sbrk
#endif

typedef struct {
	int step,sig;
} sigdata_t;

typedef struct {
	int      ofs;
	int      ntype;   /* could use type 'byte' here */
	domain_t dtype;
	char     name[MAX_FILENAME];
} offset_t;

void CatchSignal( int i );

int  NumPEs = 1, SelfPE = 0;
int  DebugLevel = 0;
logical ShowVersion = false;
logical LoadState  = false;
logical SaveState  = false;
logical CheckpointState = false;
logical UseAuxvars = false;
FILE*   FpResources  = NULL;
real Beta = 2000.0;
real Cm   = 1.0;
char* LoadStateFilename = NULL;
char* SaveStateFilename = NULL;
char* CheckpointFilename = NULL;

int   NumMem = 1;
int   (*FunctionTable[MAX_MEMBRANE])();
void  (*SetpatchTable[MAX_MEMBRANE])();
int   PatchSize[MAX_MEMBRANE];
int   AuxiliarySize[MAX_MEMBRANE];
int   MemParamSize[MAX_MEMBRANE];

logical NeedVolume    = false;
logical NeedInvVolume = false;

bvector NodeType   = _INIT_VECTOR;
vector  Volume     = _INIT_VECTOR;
vector  InvVolume  = _INIT_VECTOR;
vector  Workspace  = _INIT_VECTOR;
vector  MemParams  = _INIT_VECTOR;
sparse  MatrixINT  = _INIT_SPARSE;
sparse  MatrixEXT  = _INIT_SPARSE;
sparse  MatrixBATH = _INIT_SPARSE;
sparse  MatrixE2B  = _INIT_SPARSE;
sparse  MatrixB2E  = _INIT_SPARSE;

int    TissueGlobalSize = 0, TissueLocalSize = 0;
int    BathGlobalSize = 0, BathLocalSize = 0;
int    SolverGlobalSize = 0, SolverLocalSize = 0;
int    StateGlobalSize = 0, StateLocalSize = 0;
int    AuxiliaryGlobalSize = 0, AuxiliaryLocalSize = 0;
int    MemparamGlobalSize = 0, MemparamLocalSize = 0;
int*   TissueSplit;
int*   BathSplit;
int*   SolverSplit;
int*   StateSplit;
int*   AuxiliarySplit;
int*   MemparamSplit;

extern char* VectorOps_RCSID;
extern char* SparseOps_RCSID;

static char* RCSID = "$Id: CWaveKernel.c 14 2007-05-11 14:57:55Z jbp $";

static char  CkptVmFile[MAX_FILENAME],CkptVmTemp[MAX_FILENAME];
static char  CkptViFile[MAX_FILENAME],CkptViTemp[MAX_FILENAME];
static char  CkptVeFile[MAX_FILENAME],CkptVeTemp[MAX_FILENAME];
static char  CkptVbFile[MAX_FILENAME],CkptVbTemp[MAX_FILENAME];
static char  CkptQFile[MAX_FILENAME],CkptQTemp[MAX_FILENAME];
static int   CheckpointInterval = -1, Ckpt_count = 0;
static int   CheckpointIntTime  = -1;
static int   CkptMsgTag;

static char* InitCondFileVm = NULL;
static char* InitCondFile[MAX_PATCHSIZE];

static int   CW_SignalRecvd = 0;
static int   CW_NumSteps = 0;
static int   CW_AlmostDone = 0;
static MPI_Comm     SignalComm;
static sigdata_t*   SignalData = NULL;
static MPI_Request* SignalReqs = NULL;
static MPI_Status*  SignalStats = NULL;

static char* MemparamFile = NULL;

static int NextTagNum = 1;

static int NumOffsets = 0;
static offset_t OffsetList[MAX_OFFSETS];

static rword resources[] = {
	{ "auxvar",	9044 },
	{ "auxvars",	9044 },
	{ "beta",	9010 },
	{ "catchsig",	9099 },
	{ "catchsignal",9099 },
	{ "checkpoint",	9100 },
	{ "checkpointfile", 9101 },
	{ "checkpointinterval", 9102 },
	{ "ckpt",	9100 },
	{ "ckptfile",	9101 },
	{ "ckptint",	9102 },
	{ "ckptinterval", 9102 },
	{ "checkpointtime", 9103 },
	{ "ckpttime",	9103 },
	{ "ckptinttime", 9103 },
	{ "ckptintervaltime", 9103 },
	{ "cm",		9011 },
	{ "debug",	9002 },
	{ "debuglevel",	9002 },
	{ "debugflush",	9003 },
	{ "debug_flush",9003 },
	{ "loadfile",   9040 },
	{ "loadstate",	9030 },
	{ "memfile",	9033 },
	{ "memparam",	9033 },
	{ "memparamfile", 9033 },
	{ "savefile",	9041 },
	{ "savestate",  9031 },
	{ "statefile",	9021 },
	{ "statfile",	9021 },
	{ "statusfile",	9021 },
	{ "forceauxvars",9045 },
	{ "forceauxvar",9045 },
	{ "forceaux",	9045 },
	{ "useauxvar",	9044 },
	{ "useauxvars",	9044 },
	{ "useaux",	9044 },
	{ "initvm",	9201 },
	{ "initialvm",	9201 },
	{ "initcondvm",	9201 },
	{ "initcondfilevm", 9201 },
	{ "initcondvmfile", 9201 },
	{ "initialconditionvm", 9201 },
	{ "initcond",	9200 },
	{ "initcondfile", 9200 },
	{ "initialcondition", 9200 },
	{ "showversion", 9300 },
	{ "show_version", 9300 },
	{ "showver", 9300 },
	{ "show_ver", 9300 },
	{ "showrev", 9300 },
	{ "show_rev", 9300 },
	{ NULL, 0 }
};

int InitSimulation( char** res ) {
	int i;
	int max_npes = -1;
	int cmd,num;
	char fname[MAX_FILENAME];
	char* av[2] = { NULL, NULL };
	char* statefilename = NULL;
	logical catchsig = false;
	logical DebugFlush = false;
	logical forceauxvars = false;
	MPI_Request sendreq;
	sigdata_t sigd;

	InitCondFileVm = NULL;
	for(i=0;i<MAX_PATCHSIZE;i++) {
		InitCondFile[i] = NULL;
	}

	/* go through resource list and look for word/val pairs */
	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		switch( cmd ) {
			case 9002:
				DebugLevel = GetIntValue( res[i] );
				break;
			case 9003:
				DebugFlush = GetTFValue( res[i] );
				break;
			case 9010:
				Beta = GetRealValue( res[i] );
				break;
			case 9011:
				Cm = GetRealValue( res[i] );
				break;
			case 9021:
				statefilename = GetStringValue( res[i] );
				break;
			case 9030:
				LoadState = GetTFValue( res[i] );
				break;
			case 9031:
				SaveState = GetTFValue( res[i] );
				break;
			case 9033:
				MemparamFile = GetStringValue( res[i] );
				break;
			case 9040:
				LoadStateFilename = GetStringValue( res[i] );
				break;
			case 9041:
				SaveStateFilename = GetStringValue( res[i] );
				break;
			case 9044:
				UseAuxvars = GetTFValue( res[i] );
				break;
			case 9045:
				/* you can only force auxvars to true */
				/* you cannot un-force them to false */
				if( GetTFValue(res[i]) == true ) {
					forceauxvars = true;
				}
				break;
			case 9099:
				catchsig = GetTFValue( res[i] );
				break;
			case 9100:
				CheckpointState = GetTFValue( res[i] );
				break;
			case 9101:
				CheckpointFilename = GetStringValue( res[i] );
				break;
			case 9102:
				CheckpointInterval = GetIntValue( res[i] );
				break;
			case 9103:
				CheckpointIntTime = GetTimeValue( res[i] );
				break;
			case 9200:
				num = FindNum( res[i] );
				InitCondFile[num] = GetStringValue( res[i] );
				break;
			case 9201:
				InitCondFileVm = GetStringValue( res[i] );
				break;
			case 9300:
				ShowVersion = GetTFValue( res[i] );
				break;
		}
		i++;
	}

	/* see who we are */
	MPI_Comm_size( MPI_COMM_WORLD, &NumPEs );
	MPI_Comm_rank( MPI_COMM_WORLD, &SelfPE );

	/* create a duplicate communicator so that signal processing */
	/* does not interfere with normal processing */
	MPI_Comm_dup( MPI_COMM_WORLD, &SignalComm );
	SignalReqs = (MPI_Request*)malloc( NumPEs*sizeof(MPI_Request) );
	if( SignalReqs == NULL ) {
		return( -1 );
	}
	SignalData = (sigdata_t*)malloc( NumPEs*sizeof(sigdata_t) );
	if( SignalData == NULL ) {
		return( -2 );
	}
	for(i=0;i<NumPEs;i++) {
		SignalData[i].step = 0;
		SignalData[i].sig  = 0;
	}
	SignalStats = (MPI_Status*)malloc( NumPEs*sizeof(MPI_Status) );
	if( SignalStats == NULL ) {
		return( -3 );
	}
	/* post open Irecv's from each remote pe */
	for(i=0;i<NumPEs;i++) {
		if( i == SelfPE ) {
			SignalReqs[i] = MPI_REQUEST_NULL;
		} else {
			MPI_Irecv( &SignalData[i], 2, MPI_INT, i, MPI_ANY_TAG,
				SignalComm, &SignalReqs[i] );
		}
	}
	for(i=0;i<NumPEs;i++) {
		sigd.step = 0;
		sigd.sig  = 0;
		if( i != SelfPE ) {
			MPI_Isend( &sigd, 2, MPI_INT, i, 123,
				SignalComm, &sendreq );
			MPI_Request_free( &sendreq );
		}
	}

	/* do we need to use auxvars, even if user says not to? */
	if( forceauxvars == true ) {
		UseAuxvars = true;
	}

	if( (LoadStateFilename==NULL) and (statefilename!=NULL) ) {
		LoadStateFilename = statefilename;
	}
	if( (SaveStateFilename==NULL) and (statefilename!=NULL) ) {
		SaveStateFilename = statefilename;
	}
	if( (CheckpointFilename==NULL) and (statefilename!=NULL) ) {
		CheckpointFilename = statefilename;
	}

	if( SaveState or CheckpointState or catchsig ) {
		/* attach to the CPU Time Limit signal(s) */
		#if defined(SIGCPULIM)
			signal( SIGCPULIM, CatchSignal );
		#endif
		#if defined(SIGXCPU)
			signal( SIGXCPU, CatchSignal );
		#endif
		/* also catch the terminate/quit/interrupt signals */
		signal( SIGHUP, CatchSignal );
		signal( SIGINT, CatchSignal );
		signal( SIGQUIT, CatchSignal );
		signal( SIGTERM, CatchSignal );
		/* use ALARM signal for watchdog timer? */
		/* : or just as a signal to force a checkpoint */
		signal( SIGALRM, CatchSignal );
		/* NOTE: MPICH uses SIGUSR1, so we can't touch it */
		/* NOTE: LAM uses SIGUSR2, so we can't touch it */
	}

	if( SaveState or CheckpointState ) {
		/* reopen file for possible restart of this run */
		if( SelfPE == 0 ) {
			FpResources = fopen( "restart.in", "a" );
		} else {
			FpResources = NULL;
		}
	}

	CkptMsgTag = NextMsgTag();
	if( CheckpointState ) {
		/* : Vm first */
		strncpy( CkptVmFile, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVmFile, ".vm.vec" );
		strncpy( CkptVmTemp, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVmTemp, ".vm_tmp.vec" );

		/* : Ve/Vi/Vb next */
		strncpy( CkptViFile, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptViFile, ".vi.vec" );
		strncpy( CkptViTemp, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptViTemp, ".vi_tmp.vec" );
		strncpy( CkptVeFile, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVeFile, ".ve.vec" );
		strncpy( CkptVeTemp, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVeTemp, ".ve_tmp.vec" );
		strncpy( CkptVbFile, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVbFile, ".vb.vec" );
		strncpy( CkptVbTemp, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptVbTemp, ".vb_tmp.vec" );

		/* : Q last */
		strncpy( CkptQFile, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptQFile, ".q.vec" );
		strncpy( CkptQTemp, CheckpointFilename, MAX_FILENAME );
		strcat(  CkptQTemp, ".q_tmp.vec" );
	} else {
		CheckpointInterval = -1;
		CheckpointIntTime  = -1;
	}

	if( (DebugLevel>0) and (SelfPE==0) ) {
#if defined(TEST_ASSERT)
		i = 1;
#else
		i = 0;
#endif
		printf("CardioWave: debug=%i/%c, load=%c, save=%c, ckpt=%c/%i/%i, aux=%c\n",
			DebugLevel, (i?'Y':'N'),
			(LoadState?'Y':'N'),(SaveState?'Y':'N'),
			(CheckpointState?'Y':'N'),CheckpointInterval,
			CheckpointIntTime,(UseAuxvars?'Y':'N') );
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("CWaveKernel: RCSID: %s\n",RCSID);
		printf("VectorOps: RCSID: %s\n",VectorOps_RCSID);
		printf("SparseOps: RCSID: %s\n",SparseOps_RCSID);
	}

	return( 0 );
}

void ExitSimulation( vector Vm, vector Vx, vector Q, vector Av ) {
	char  fname[MAX_FILENAME];
	vector vtmp;

	if( SaveState ) {
		/* dump state information */

		/* : Vm first */
		strncpy( fname, SaveStateFilename, MAX_FILENAME );
		strcat(  fname, ".vm.vec" );
		vecwrite( fname, Vm );

		/* : Ve/Vi/Vb next */
		if( (SimType==FullBidomain) or (SimType==FullBidomainBath) ) {
			strncpy( fname, SaveStateFilename, MAX_FILENAME );
			strcat(  fname, ".vi.vec" );
			AliasVxVi( Vx, &vtmp );
			vecwrite( fname, vtmp );
		}
		if( SimType != Monodomain ) {
			strncpy( fname, SaveStateFilename, MAX_FILENAME );
			strcat(  fname, ".ve.vec" );
			AliasVxVe( Vx, &vtmp );
			vecwrite( fname, vtmp );
		}
		if( (SimType==ReducedBidomainBath) 
		    or (SimType==FullBidomainBath) ) {
			strncpy( fname, SaveStateFilename, MAX_FILENAME );
			strcat(  fname, ".vb.vec" );
			AliasVxVb( Vx, &vtmp );
			vecwrite( fname, vtmp );
		}

		/* : Q last */
		if( StateGlobalSize != 0 ) {
			strncpy( fname, SaveStateFilename, MAX_FILENAME );
			strcat(  fname, ".q.vec" );
			vecwrite( fname, Q );
		}

		/* : MemParams too ? */
		if( MemparamGlobalSize != 0 ) {
			strncpy( fname, SaveStateFilename, MAX_FILENAME );
			strcat(  fname, ".mp.vec" );
			vecwrite( fname, MemParams );
		}

		/* note that this will be the LAST function called */
		if( SelfPE == 0 ) {
			fclose( FpResources );
		}
	}

	return;
}

int Checkpoint( real t, vector Vm, vector Vx, vector Q, vector Av ) {
	vector Vi,Ve,Vb;
	int p,ex_flag,ck_flag,e,sig;
	MPI_Request sendreq;
	sigdata_t senddata;
	logical do_ckpt;

	/* determine if we need to ckpt */
	/* : if so, the ckpt will occur at next t-step */
	sig = CW_SignalRecvd;
	Ckpt_count++;
	if( sig == SIGALRM ) {
		/* reset the watchdog timer */
		CW_SignalRecvd = 0;
		if( CheckpointIntTime > 0 ) {
			if( SelfPE == 0 ) {
				alarm( CheckpointIntTime );
			}
		}
	} else if( sig > 0 ) {
		/* reset the signal flag */
		CW_SignalRecvd = 0;
	} else if( (CheckpointInterval>0)
			and (Ckpt_count>CheckpointInterval) ) {
		/* signal a ckpt */
printf("PE%i: ckpt\n",SelfPE); fflush(stdout);
		sig = SIGALRM;
		Ckpt_count = 0;
	}

	/* post sends for this global-step */
	senddata.step = CW_NumSteps;
	if( sig ) {
		senddata.sig = sig;
	} else {
		senddata.sig = 0;
	}
	for(p=0;p<NumPEs;p++) {
		if( p != SelfPE ) {
			MPI_Isend( &senddata, 2, MPI_INT, p, 123,
				SignalComm, &sendreq );
			MPI_Request_free( &sendreq );
		}
	}

	/* wait for recv from previous global-step */
	MPI_Waitall( NumPEs, SignalReqs, SignalStats );
	/* : all PEs must now be at least at NumSteps */

	/* now check for error/ckpt conditions */
	do_ckpt = false;
	ck_flag = 0;
	ex_flag = 0;
	for(p=0;p<NumPEs;p++) {
		if( SignalData[p].sig == SIGALRM ) {
			ck_flag = 1;
		} else if( SignalData[p].sig ) {
			ex_flag = 1;
		}
	}

	if( ex_flag ) {
		/* force an early exit */
		return( -1 );
	} else if( ck_flag ) {
		/* force an early checkpoint */
		do_ckpt = true;
	}

	/* post new recvs for this global-step */
	for(p=0;p<NumPEs;p++) {
		if( p == SelfPE ) {
			SignalReqs[p] = MPI_REQUEST_NULL;
		} else {
			MPI_Irecv( &SignalData[p], 2, MPI_INT, p, MPI_ANY_TAG,
				SignalComm, &SignalReqs[p] );
		}
	}

	CW_NumSteps++;

	if( CheckpointState and do_ckpt ) {
printf("PE%i: ckpt system\n",SelfPE);
fflush(stdout);
		/* first, we'll hit a barrier to make */
		/* sure everyone is still alive       */
		e = MPI_Barrier( SignalComm );
		if( e < 0 ) {
			return( -100 );
		}

		/* initially, we write to temp files */
		vecwrite( CkptVmTemp, Vm );

		if( (SimType==FullBidomain) 
		    or (SimType==FullBidomainBath) ) {
			AliasVxVi( Vx, &Vi );
			vecwrite( CkptViTemp, Vi );
		}
		if( SimType != Monodomain ) {
			AliasVxVe( Vx, &Ve );
			vecwrite( CkptVeTemp, Ve );
		}
		if( (SimType==ReducedBidomainBath)
		    or (SimType==FullBidomainBath) ) {
			AliasVxVb( Vx, &Vb );
			vecwrite( CkptVbTemp, Vb );
		}

		if( StateGlobalSize > 0 ) {
			vecwrite( CkptQTemp, Q );
		}

		/* now, we can rename the files */
		if( SelfPE == 0 ) {
			rename( CkptVmTemp, CkptVmFile );
			if( (SimType==FullBidomain)
			    or (SimType==FullBidomainBath) ) {
				rename( CkptViTemp, CkptViFile );
			}
			if( SimType != Monodomain ) {
				rename( CkptVeTemp, CkptVeFile );
			}
			if( (SimType==ReducedBidomainBath)
			    or (SimType==FullBidomainBath) ) {
				rename( CkptVbTemp, CkptVbFile );
			}
			if( StateGlobalSize > 0 ) {
				rename( CkptQTemp, CkptQFile );
			}
		}

		/* finally, we write out the time-step that */
		/* we just completed */
		if( SelfPE == 0 ) {
			rewind( FpResources );
			fprintf( FpResources, "t0=%le\n", t );
			fflush( FpResources );
		}

	}

	/* mimic a send of the the data to self */
	/* : thus, NEXT t-step we will detect error along with other PEs */
	SignalData[SelfPE] = senddata;

	return( 0 );
}

int ComputeSizes( char** res ) {
	int i;
	byte* nt = NodeType.data;
	int* itmp;

	/* allocate space for a work or temporary vector */
	if( BathLocalSize > TissueLocalSize ) {
		if( vecalloc(&Workspace,Bath) < 0 ) {
			return( -8 );
		}
	} else {
		if( vecalloc(&Workspace,Tissue) < 0 ) {
			return( -9 );
		}
	}

	itmp = (int*)Workspace.data;

	/* add up the state sizes */
	StateSplit = (int*)malloc( (NumPEs+1)*sizeof(int) );
	if( StateSplit == NULL ) {
		return( -1 );
	}
	StateLocalSize = 0;
	for(i=0;i<TissueLocalSize;i++) {
		StateLocalSize += PatchSize[ nt[i] ];
	}
	MPI_Allgather( &StateLocalSize, 1, MPI_INT, itmp, 
		1, MPI_INT, MPI_COMM_WORLD );
	StateGlobalSize = 0;
	for(i=0;i<NumPEs;i++) {
		StateSplit[i] = StateGlobalSize;
		StateGlobalSize += itmp[i];
	}
	StateSplit[i] = StateGlobalSize;

	/* add up the auxiliary sizes */
	AuxiliarySplit = (int*)malloc( (NumPEs+1)*sizeof(int) );
	if( AuxiliarySplit == NULL ) {
		return( -2 );
	}
	AuxiliaryLocalSize = 0;
	for(i=0;i<TissueLocalSize;i++) {
		AuxiliaryLocalSize += AuxiliarySize[ nt[i] ];
	}
	MPI_Allgather( &AuxiliaryLocalSize, 1, MPI_INT, itmp, 
		1, MPI_INT, MPI_COMM_WORLD );
	AuxiliaryGlobalSize = 0;
	for(i=0;i<NumPEs;i++) {
		AuxiliarySplit[i] = AuxiliaryGlobalSize;
		AuxiliaryGlobalSize += itmp[i];
	}
	AuxiliarySplit[i] = AuxiliaryGlobalSize;

	/* add up the mem param sizes */
	MemparamSplit = (int*)malloc( (NumPEs+1)*sizeof(int) );
	if( MemparamSplit == NULL ) {
		return( -3 );
	}
	MemparamLocalSize = 0;
	for(i=0;i<TissueLocalSize;i++) {
		MemparamLocalSize += MemParamSize[ nt[i] ];
	}
	MPI_Allgather( &MemparamLocalSize, 1, MPI_INT, itmp, 
		1, MPI_INT, MPI_COMM_WORLD );
	MemparamGlobalSize = 0;
	for(i=0;i<NumPEs;i++) {
		MemparamSplit[i] = MemparamGlobalSize;
		MemparamGlobalSize += itmp[i];
	}
	MemparamSplit[i] = MemparamGlobalSize;

	SolverSplit = (int*)malloc( (NumPEs+1)*sizeof(int) );
	if( SolverSplit == NULL ) {
		return( -4 );
	}
	switch( SimType ) {
		case Monodomain:
		case ReducedBidomain:
			SolverGlobalSize = TissueGlobalSize;
			SolverLocalSize = TissueLocalSize;
			for(i=0;i<(NumPEs+1);i++) {
				SolverSplit[i] = TissueSplit[i];
			}
			break;
		case FullBidomain:
			SolverGlobalSize = 2*TissueGlobalSize;
			SolverLocalSize = 2*TissueLocalSize;
			for(i=0;i<(NumPEs+1);i++) {
				SolverSplit[i] = 2*TissueSplit[i];
			}
			break;
		case ReducedBidomainBath:
			SolverGlobalSize = TissueGlobalSize + BathGlobalSize;
			SolverLocalSize = TissueLocalSize + BathLocalSize;
			for(i=0;i<(NumPEs+1);i++) {
				SolverSplit[i] = TissueSplit[i] + BathSplit[i];
			}
			break;
		case FullBidomainBath:
			SolverGlobalSize = 2*TissueGlobalSize + BathGlobalSize;
			SolverLocalSize = 2*TissueLocalSize + BathLocalSize;
			for(i=0;i<(NumPEs+1);i++) {
				SolverSplit[i] = 2*TissueSplit[i]
						+ BathSplit[i];
			}
			break;
	}

	return( 0 );
}

int DomainLocalSize( domain_t dt ) {
        switch( dt ) {
                case Tissue:
                case Intra:
                case Extra:
                        return( TissueLocalSize );
                        break;
                case Bath:
                        return( BathLocalSize );
                        break;
                case Solver:
                        return( SolverLocalSize );
                        break;
                case State:
                        return( StateLocalSize );
                        break;
                case Aux:
                        return( AuxiliaryLocalSize );
                        break;
		case Param:
			return( MemparamLocalSize );
			break;
        }
        return( -1 );
}

int DomainGlobalSize( domain_t dt ) {
        switch( dt ) {
                case Tissue:
                case Intra:
                case Extra:
                        return( TissueGlobalSize );
                        break;
                case Bath:
                        return( BathGlobalSize );
                        break;
                case Solver:
                        return( SolverGlobalSize );
                        break;
                case State:
                        return( StateGlobalSize );
                        break;
                case Aux:
                        return( AuxiliaryGlobalSize );
                        break;
		case Param:
			return( MemparamGlobalSize );
			break;
        }
        return( -1 );
}

int AliasVxVi( vector Vx, vector* Vi ) {
	switch( SimType ) {
		case ReducedBidomain:
		case ReducedBidomainBath:
			Vi->size = 0;
			Vi->data = NULL;
			Vi->dtype = Tissue;
			break;
		case Monodomain:
			/* fudge it so that Vi looks like Vm */
		case FullBidomain:
		case FullBidomainBath:
			Vi->size = TissueLocalSize;
			Vi->data = Vx.data;
			Vi->dtype = Tissue;
			break;
	}

	return( 0 );
}

int AliasVxVe( vector Vx, vector* Ve ) {
	switch( SimType ) {
		case Monodomain:
			Ve->size = 0;
			Ve->data = NULL;
			Ve->dtype = -1;
			break;
		case ReducedBidomain:
		case ReducedBidomainBath:
			Ve->size = TissueLocalSize;
			Ve->data = Vx.data;
			Ve->dtype = Tissue;
			break;
		case FullBidomain:
		case FullBidomainBath:
			Ve->size = TissueLocalSize;
			Ve->data = Vx.data + TissueLocalSize;
			Ve->dtype = Tissue;
			break;
	}

	return( 0 );
}

int AliasVxVb( vector Vx, vector* Vb ) {
	switch( SimType ) {
		case Monodomain:
		case FullBidomain:
		case ReducedBidomain:
			Vb->size = 0;
			Vb->data = NULL;
			Vb->dtype = -1;
			break;
		case FullBidomainBath:
			Vb->size = TissueLocalSize;
			Vb->data = Vx.data + 2*TissueLocalSize;
			Vb->dtype = Bath;
			break;
		case ReducedBidomainBath:
			Vb->size = TissueLocalSize;
			Vb->data = Vx.data + TissueLocalSize;
			Vb->dtype = Bath;
			break;
	}

	return( 0 );
}

/* allocate the vectors and load the initial values */
int GetInitialValues( vector* Vm, vector* Vx, vector* Q, vector* Av ) {
	int   i,j,k,r,n;
	int   status = 0;
	int   initflag = 0;
	real* rp;
	char  fnamevm[MAX_FILENAME];
	char  fnamevi[MAX_FILENAME];
	char  fnameve[MAX_FILENAME];
	char  fnamevb[MAX_FILENAME];
	char  fnameq[MAX_FILENAME];
	vector Vi,Ve,Vb;
	real* wksp;
	real* qp;
	byte* nt = NodeType.data;

	/* allocate or load in the state vectors */
	if( LoadState ) {
		/* attempt to load in the state file vectors */
		/* : load vectors will return a negative indicator on error */
		/* first, check files to make sure they are correct size */
		strncpy( fnamevm, LoadStateFilename, MAX_FILENAME );
		strcat(  fnamevm, ".vm.vec" );
		r = vecreadinfo( fnamevm, &i );
		if( (r<0) or (i!=TissueGlobalSize) ) {
			if( SelfPE == 0 ) {
				printf("Problems loading vm restart file [%s]\n",
					fnamevm );
			}
			status = -1;
		}

		if( SimType != Monodomain ) {
			strncpy( fnameve, LoadStateFilename, MAX_FILENAME );
			strcat(  fnameve, ".ve.vec" );
			r = vecreadinfo( fnameve, &i );
			if( (r<0) or (i!=TissueGlobalSize) ) {
			  if( SelfPE == 0 ) {
				printf("Problems loading ve restart file [%s]\n",
					fnameve );
			  }
			  status = -1;
			}
		}
		if( (SimType==FullBidomain) or (SimType==FullBidomainBath) ) {
			strncpy( fnamevi, LoadStateFilename, MAX_FILENAME );
			strcat(  fnamevi, ".vi.vec" );
			r = vecreadinfo( fnamevi, &i );
			if( (r<0) or (i!=TissueGlobalSize) ) {
			  if( SelfPE == 0 ) {
				printf("Problems loading vi restart file [%s]\n",
					fnamevi );
			  }
			  status = -1;
			}
		}
		if( (SimType==ReducedBidomainBath) 
		    or (SimType==FullBidomainBath) ) {
			strncpy( fnamevb, LoadStateFilename, MAX_FILENAME );
			strcat(  fnamevb, ".vb.vec" );
			r = vecreadinfo( fnamevb, &i );
			if( (r<0) or (i!=BathGlobalSize) ) {
			  if( SelfPE == 0 ) {
				printf("Problems loading vb restart file [%s]\n",
					fnamevb );
			  }
			  status = -1;
			}
		}
		if( StateGlobalSize > 0 ) {
			strncpy( fnameq, LoadStateFilename, MAX_FILENAME );
			strcat(  fnameq, ".q.vec" );
			r = vecreadinfo( fnameq, &i );
			if( (r<0) or (i!=StateGlobalSize) ) {
			  if( SelfPE == 0 ) {
				printf("Problems loading q restart file [%s]\n",
					fnameq );
			  }
			  status = -1;
			}
		}
	}

	/* allocate the vectors */
	if( vecalloc(Vm,Tissue) < 0 ) {
		return( -1 );
	}
	if( vecalloc(Vx,Solver) < 0 ) {
		return( -2 );
	}
	if( StateGlobalSize > 0 ) {
		if( vecalloc(Q,State) < 0 ) {
			return( -3 );
		}
	}
	/* allocate space for the AuxVars */
	if( (AuxiliaryGlobalSize>0) and UseAuxvars ) {
		if( vecalloc(Av,Aux) < 0 ) {
			return( -4 );
		}
	}

	if( MemparamGlobalSize > 0 ) {
		if( vecalloc(&MemParams,Memparam) < 0 ) {
			return( -5 );
		}
		r = vecreadinfo( MemparamFile, &i );
		if( (r<0) or (i!=MemparamGlobalSize) ) {
		  for(i=0;i<MAX_MEMBRANE;i++) {
			if( SetpatchTable[i] != NULL ) {
				SetpatchTable[i]( Memparam, *Vm, *Q, *Av );
			}
		  }
		  if( SelfPE == 0 ) {
		    if( MemparamFile != NULL ) {
		      printf("Problems loading membrane parameter file [%s]\n",
				MemparamFile );
		    } else {
		      printf("No membrane parameter file given\n");
		    }
		    printf("Using modules' Setpatch functions\n");
		  }
		} else {
		  vecread( MemparamFile, MemParams );
		  if( SelfPE == 0 ) {
			printf("Loaded membrane parameters from file\n");
		  }
		}
	}

	if( LoadState and (status>=0) ) {
		/* get Vm */
		vecread( fnamevm, *Vm );

		/* get Ve & Vi */
		if( (SimType==FullBidomain) 
		    or (SimType==FullBidomainBath) ) {
			AliasVxVi( *Vx, &Vi );
			vecread( fnamevi, Vi );
		}
		if( SimType != Monodomain ) {
			AliasVxVe( *Vx, &Ve );
			vecread( fnameve, Ve );
		}
		if( (SimType==ReducedBidomainBath)
		    or (SimType==FullBidomainBath) ) {
			AliasVxVb( *Vx, &Vb );
			vecread( fnamevb, Vb );
		}

		/* get Q */
		if( StateGlobalSize > 0 ) {
			vecread( fnameq, *Q );
		}

		if( (SelfPE==0) and (DebugLevel>0) ) {
			printf("Loaded initial values from file\n");
		}
	} else {
		/* now call InitPatches() */
		for(i=0;i<MAX_MEMBRANE;i++) {
			if( SetpatchTable[i] != NULL ) {
				SetpatchTable[i]( Tissue, *Vm, *Q, *Av );
			}
		}
		for(i=0;i<MAX_MEMBRANE;i++) {
			if( SetpatchTable[i] != NULL ) {
				SetpatchTable[i]( AuxVar, *Vm, *Q, *Av );
			}
		}

		if( (SelfPE==0) and (DebugLevel>0) ) {
			printf("Computed initial values as resting conditions\n");
		}
	}

	/* now overwrite with any 'initcond[]=file' settings */
	n = -1;
	for(i=0;i<MAX_MEMBRANE;i++) {
		if( PatchSize[i] > n ) {
			n = PatchSize[i];
		}
	}
	if( InitCondFileVm ) {
		if( (SelfPE==0) and (DebugLevel>0) ) {
		   printf("Overwriting initial voltage with file [%s]\n",
				InitCondFileVm );
		}
		r = vecreadinfo( InitCondFileVm, &j );
		if( (r==0) and (j==TissueGlobalSize) ) {
			vecread( InitCondFileVm, Workspace );
			wksp = Workspace.data;
			qp = Vm->data;
			for(k=0;k<TissueLocalSize;k++) {
				qp[i] = wksp[i];
				qp++;
				wksp++;
			}
		}
	}
	initflag = 0;
	for(i=0;i<n;i++) {
		if( InitCondFile[i] ) {
			initflag = 1;
			if( (SelfPE==0) and (DebugLevel>0) ) {
			   printf("Overwriting state var %d with file [%s]\n",
					i, InitCondFile[i] );
			}
			r = vecreadinfo( InitCondFile[i], &j );
			if( (r==0) and (j==TissueGlobalSize) ) {
				vecread( InitCondFile[i], Workspace );
				wksp = Workspace.data;
				qp = Q->data;
				for(k=0;k<TissueLocalSize;k++) {
					qp[i] = wksp[i];
					qp += PatchSize[ nt[k] ];
					wksp += PatchSize[ nt[k] ];
				}
			}
		}
	} 
	if( initflag and (SelfPE==0) and (DebugLevel>0) ) {
		printf("Found %i state vars max\n",n);
		if( NumMem > 1 ) {
			printf("Warning: multiple membrane models found [%d]\n",
				NumMem );
		}
	}

	/* this is the last fcn called before main sim is started */
	/* so we'll start the watchdog timer here */
	if( CheckpointIntTime > 0 ) {
		if( SelfPE == 0 ) {
			alarm( CheckpointIntTime );
		}
	}

	return( 0 );
}

/* GetF is just a loop over the MEM*_GetF functions */
#ifndef _NEURO
int GetF( real t, real dt, vector Vm, vector Q, 
		vector Fv, vector Fq, vector Av ) {
	int i,r;

	DebugEnter( "GetF" );

	r = 0;

	for(i=0;i<MAX_MEMBRANE;i++) {
		if( FunctionTable[i] != NULL ) {
			r |= FunctionTable[i]( t, dt, Vm, Q, Fv, Fq, Av );
		}
	}

	DebugLeave( "GetF" );

	return( r );
}
#endif

void CatchSignal( int i ) {
	CW_SignalRecvd = i;
	return;
}

int NextMsgTag( void ) {
	int r = NextTagNum;
	NextTagNum++;
	return( r );
}

int AppendResources( char** res, int ac, char** av ) {
	FILE* fp;
	int   i,k,c;
	long  l;
	char* buffer;
	char* ptr;

	c = NumResources( res );
	for(i=1;i<ac;i++) {
		if( (av[i][0]=='-') or (av[i][0]=='+') ) {
			/* input file */
			fp = fopen( av[i]+1, "r" );
			if( fp != NULL ) {
				/* get file size */
				l = 0;
				fseek( fp, l, SEEK_END );
				l = ftell( fp );
				/* allocate memory to hold the file */
				buffer = (char*)malloc( l+1 );
				/* load file into malloc'd memory */
				rewind( fp );
				fread( buffer, 1, l, fp );
				fclose( fp );
				buffer[l] = 0;
				/* parse each line */
				k = 0;
				ptr = buffer;
				while( k < l ) {
					if( buffer[k] == '\n' ) {
						buffer[k] = 0;
						if( (ptr[0]!='#') 
						  and (ptr[0]!='%') 
						  and (ptr[0]!='!')
						  and (ptr[0]!='/')
						  and (ptr[0]!=0) ) {
							res[c] = ptr;
							c++;
						}
						/* set up for next one */
						ptr = buffer+k+1;
					} else if( (buffer[k]=='\\') 
						   and (buffer[k+1]=='\n') ) {
						/* catenate with */
						/* previous line */
						buffer[k] = ' ';
						buffer[k+1] = ' ';
					}
					k++;
				}
				/* check to see if we forgot to put a return */
				/* after the last line of the input file     */
				if( ptr[0] != 0 ) {
					res[c] = ptr;
					c++;
					if( c >= MAX_RESOURCES ) {
						return( -1 );
					}
				}
			} else {
				printf("Error opening file %s\n",av[i]+1);
			}
		} else {
			/* resource spec:  stuff=value */
			res[c] = av[i];
			c++;
			if( c >= MAX_RESOURCES ) {
				return( -1 );
			}
		}
	}
	res[c] = NULL;

	return( 0 );
}

int NumResources( char** res ) {
	int i;

	i = 0;
	while( res[i] != NULL ) {
		i++;
	}

	return( i );
}

int NumWords( rword* list ) {
	int i;

	i = 0;
	while( list[i].txt != NULL ) {
		i++;
	}

	return( i );
}

int GetTimeValue( char* txt ) {
	int i,j,len,sec;
	logical colon;
	char* cp;

	/* move to the '=' */
	len = strlen( txt );
	i   = strcspn( txt, "=" );
	if( i == 0 ) {
		return( -1 );
	}
	i++;

	/* look for ':hms' */
	colon = false;
	sec = 0;
	while( i < len ) {
		j = strcspn( txt+i, ":hms" );
		if( j == 0 ) {
			break;
		}
		switch( txt[i+j] ) {
			case ':':
				colon = true;
				sec = 60*sec + atoi( txt+i );
				break;
			case 'H':
			case 'h':
				sec += 3600*atoi( txt+i );
				break;
			case 'M':
			case 'm':
				sec += 60*atoi( txt+i );
				break;
			case 'S':
			case 's':
				sec += atoi( txt+i );
				break;
			case 0:
				if( colon ) {
					sec = 60*sec + atoi( txt+i );
				} else {
					sec += atoi( txt+i );
				}
		}
		i += j + 1;
	}

	return( sec );
}

real GetRealValue( char* txt ) {
	int i;
	real r;

	/* move to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}

	if( txt[i] == 0 ) {
		return( -9.9e9 );
	}

	return( atof(txt+i+1) );
}

int GetIntValue( char* txt ) {
	int i;

	/* move to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}

	if( txt[i] == 0 ) {
		return( -1 );
	}

	return( atoi(txt+i+1) );
}

byte GetByteValue( char* txt ) {
	int i;
	byte b;

	/* move to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}

	if( txt[i] == 0 ) {
		return( -1 );
	}

	b = (byte)(atoi(txt+i+1));
	return( b );
}

logical GetTFValue( char* txt ) {
	int i;

	/* move to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}
	i++;

	if( txt[i] == 0 ) {
		/* this is really an error condition, but no real way */
		/* to indicate that with only 'true' or 'false'       */
		return( false );
	}

	/* look for the 't' in 'true' or the 'y' in 'yes' */
	if( (txt[i]=='t') or (txt[i]=='T') or (txt[i]=='1') 
		or (txt[i]=='y') or (txt[i]=='Y') ) {
		return( true );
	}	
	/* look for the word 'on' */
	if( (txt[i]=='o') or (txt[i]=='O') ) {
		if( (txt[i+1]=='n') or (txt[i+1]=='N') ) {
			return( true );
		}
	}

	return( false );
}

int GetNumValues( char* txt ) {
	int i,j,k,s,nc;
	char* p;
	char* cpy;
	int count;
	static char* lasttxt   = NULL;
	static int   lastcount = -1;

	if( txt == lasttxt ) {
		return( lastcount );
	}

	/* go to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}
	i++;

	/* go thru each token in the text */
	count = 0;
	cpy = strdup( txt+i );
	p   = strtok( cpy, "," );
	while( p != NULL ) {
		/* look for commas and colons */
		i = 0;
		nc = 0;
		while( p[i] != 0 ) {
			if( p[i] == ':' ) {
				if( nc == 0 ) {
					j = i;
				} else if( nc == 1 ) {
					k = i;
				}
				nc++;
			}
			i++;
		}
		switch( nc ) {
			case 0: /* no colons, just a single number */
				count++;
				break;
			case 1: /* 1 colon, assume unit stride */
				k = atoi( p+j+1 );
				j = atoi( p );
				count += (k-j)+1;
				break;
			case 2: /* 2 colon, non-unit stride */
				k = atoi( p+k+1 );
				s = atoi( p+j+1 );
				j = atoi( p );
				count += ((k-j)/s) + 1;
				break;
		}
		p = strtok( NULL, "," );
	}

	lasttxt = txt;
	lastcount = count;

	return( count );
}

int* GetIntArray( char* txt ) {
	int  i,j,k,s,nc,l;
	int* list = NULL;
	char* p   = NULL;
	char* cpy = NULL;
	int count;

	count = GetNumValues( txt );
	list = (int*)malloc( count*sizeof(int) );

	/* go to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}
	i++;

	/* go thru each token in the text */
	count = 0;
	cpy = strdup( txt+i );
	p   = strtok( cpy, "," );
	while( p != NULL ) {
		/* look for commas and colons */
		i = 0;
		nc = 0;
		while( p[i] != 0 ) {
			if( p[i] == ':' ) {
				if( nc == 0 ) {
					j = i;
				} else if( nc == 1 ) {
					k = i;
				}
				nc++;
			}
			i++;
		}
		switch( nc ) {
			case 0: /* no colons, just a single number */
				list[count] = atoi( p );
				count++;
				break;
			case 1: /* 1 colon, assume unit stride */
				k = atoi( p+j+1 );
				j = atoi( p );
				for(l=j;l<=k;l++) {
					list[count] = l;
					count++;
				}
				break;
			case 2: /* 2 colons, non-unit stride */
				k = atoi( p+k+1 );
				s = atoi( p+j+1 );
				j = atoi( p );
				for(l=j;l<=k;l+=s) {
					list[count] = l;
					count++;
				}
				break;
		}
		p = strtok( NULL, "," );
	}

	return( list );
}

real* GetRealArray( char* txt ) {
	int  i,j,k,s,nc;
	char* p   = NULL;
	char* cpy = NULL;
	int count;
	real  l;
	real* list = NULL;

	count = GetNumValues( txt );
	list = (real*)malloc( count*sizeof(real) );

	/* go to the equals sign */
	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}
	i++;

	/* go thru each token in the text */
	count = 0;
	cpy = strdup( txt+i );
	p   = strtok( cpy, "," );
	while( p != NULL ) {
		/* look for commas and colons */
		i = 0;
		nc = 0;
		while( p[i] != 0 ) {
			if( p[i] == ':' ) {
				if( nc == 0 ) {
					j = i;
				} else if( nc == 1 ) {
					k = i;
				}
				nc++;
			}
			i++;
		}
		switch( nc ) {
			case 0: /* no colons, just a single number */
				list[count] = atof( p );
				count++;
				break;
			case 1: /* 1 colon, assume unit stride */
				k = atof( p+j+1 );
				j = atof( p );
				for(l=j;l<=k;l+=1.0) {
					list[count] = l;
					count++;
				}
				break;
			case 2: /* 2 colons, non-unit stride */
				k = atof( p+k+1 );
				s = atof( p+j+1 );
				j = atof( p );
				for(l=j;l<=k;l+=s) {
					list[count] = l;
					count++;
				}
				break;
		}
		p = strtok( NULL, "," );
	}

	return( list );
}

char* GetStringValue( char* txt ) {
	int i;

	i = 0;
	while( txt[i] != '=' ) {
		i++;
	}

	if( txt[i+1] == 0 ) {
		return( 0 );
	}

	return( txt+i+1 );
}

int wordcompare( const void* v1, const void* v2 ) {
	const rword* w1 = (rword*)v1;
	const rword* w2 = (rword*)v2;

	return( strcasecmp(w1->txt,w2->txt) );
}

int FindCommand( rword* list, char* txt ) {
	int    i;
	int    cmd = -1;
	rword* wf;
	static char keytxt[128];
	rword v1,v2;

	i = 0;
	while( (txt[i]!='(') and (txt[i]!='[') 
		and (txt[i]!=':') and (txt[i]!='=') 
		and (txt[i]!=' ') and (txt[i]!='\t') ) {
		keytxt[i] = txt[i];
		i++;
	}
	keytxt[i] = 0;

	wf = list;
	while( (wf->txt) != NULL ) {
		i = strcasecmp( wf->txt, keytxt );
		if( i == 0 ) {
			/* match found */
			cmd = wf->cmd;
			break;
		} else {
			wf++;
		}
	}

	return( cmd );
}

int FindNum( char* txt ) {
	int i;

	/* move forward to a bracket/paren/colon */
	i = 0;
	while( (txt[i]!='(') and (txt[i]!='[') 
		and (txt[i]!=':') and (txt[i]!=0) ) {
		i++;
	}
	if( txt[i] == 0 ) {
		return( -1 );
	}
	i++;

	if( (txt[i]<'0') or (txt[i]>'9') ) {
		if( txt[i] < 'a' ) {
			return( (int)txt[i]+32 );
		} else {
			return( (int)txt[i] );
		}
	}
	return( atoi(txt+i) );
}

int FindNum2( char* txt ) {
	int i;

	/* move forward to a bracket/paren/colon */
	i = 0;
	while( (txt[i]!='(') and (txt[i]!='[') 
		and (txt[i]!=':') and (txt[i]!=0) ) {
		i++;
	}
	if( txt[i] == 0 ) {
		return( -1 );
	}
	i++;
	/* move forward to the second bracket/paren */
	i = 0;
	while( (txt[i]!='(') and (txt[i]!='[') 
		and (txt[i]!=':') and (txt[i]!=0) ) {
		i++;
	}
	if( txt[i] == 0 ) {
		return( -1 );
	}
	i++;

	if( (txt[i]<'0') or (txt[i]>'9') ) {
		if( txt[i] < 'a' ) {
			return( (int)txt[i]+32 );
		} else {
			return( (int)txt[i] );
		}
	}
	return( atoi(txt+i) );
}

int RegisterOffset( int oo, int nt, domain_t tt, char* text ) {
	if( NumOffsets >= MAX_OFFSETS ) {
		return( -1 );
	}

	OffsetList[NumOffsets].ofs   = oo;
	OffsetList[NumOffsets].ntype = nt;
	OffsetList[NumOffsets].dtype = tt;
	strncpy( OffsetList[NumOffsets].name, text, MAX_FILENAME );

	NumOffsets++;

	return( 0 );
}

int FindOffset( char* text, int* oo, int* nt, domain_t* dt  ) {
	int i;

	for(i=0;i<NumOffsets;i++) {
		if( strncasecmp(OffsetList[i].name,text,MAX_FILENAME) == 0 ) {
			*oo = OffsetList[i].ofs;
			*nt = OffsetList[i].ntype;
			*dt = OffsetList[i].dtype;
			return( 0 );
		}
	}

	return( -1 );
}


/* module test functionality */
void CWINT_Assert( int x, char* file, const char* func, int line ) {
	if( x == 0 ) {
		printf( "    PE%i: assertion failed, file %s, func %s, line %i\n",
			SelfPE, file, func, line );
	}
	return;
}

/* : not quite ready for parallelism yet */
void CWINT_AssertAny( vector x, logop_t op, real val, 
			char* file, const char* func, int line ) {
	int i,ltmp,f;
	real xx;
	f = 0;
	for(i=0;i<x.size;i++) {
		xx = x.data[i];
		ltmp = (xx<val) + 2*(xx>val) + 4*(xx==val);
		if( (ltmp&op) > 0 ) {
			f = 1;
			break;
		}
	}
	/* **** Need to check across all PEs **** */
	if( f ) {
		printf( "    PE%i: assert-any failed, file %s, func %s, line %i\n",
				SelfPE, file, func, line );
	}
	return;
}

/* : not quite ready for parallelism yet */
void CWINT_AssertAll( vector x, logop_t op, real val, 
			char* file, const char* func, int line ) {
	int i,ltmp,s;
	real xx;
	s = 0;
	for(i=0;i<x.size;i++) {
		xx = x.data[i];
		ltmp = (xx<val) + 2*(xx>val) + 4*(xx==val);
		if( (ltmp&op) > 0 ) {
			s++;
		}
	}
	/* **** Need to sum across all PEs **** */
	if( s == x.size ) {
		printf( "    PE%i: assert-all failed, file %s, func %s, line %i\n",
			SelfPE, file, func, line );
	}
	return;
}


/* amount of memory used by the program in KB */
#if defined(_use_sbrk )
real GetMemUsed( void ) {
	char* ptr;
	int   imem;
	real  mem;

	ptr = (char*)sbrk(0);
	imem = (int)(ptr);
	mem = (real)(imem)/1024.0;

	return( mem );
}
#elif defined(_use_getrusage)
real GetMemUsed( void ) {
	struct rusage ru;
	getrusage( RUSAGE_SELF, &ru );
	return( (real)(ru.ru_maxrss) );
}
#elif defined(_use_procstat)
real GetMemUsed( void ) {
	FILE* fp;
	char buffer[512];
	int i,c,l;
	long long npgs;
	size_t pgsz;

	pgsz = getpagesize();

	sprintf( buffer, "/proc/%i/stat", getpid() );
	fp = fopen( buffer, "r" );
	fgets( buffer, 512, fp );
	fclose( fp );

	npgs = -1;
	c = 0;
	for(i=0;i<512;i++) {
		if( buffer[i] == ' ' ) {
			l = i;
			c++;
			if( c == 23 ) {
				npgs = atoll( buffer+l+1 );
			}
		}
	}

	return( (real)(npgs*pgsz)/1024.0 );
}
#elif defined(_use_mallocstat)
real GetMemUsed( void ) {
	malloc_statistics_t mstat = {0,0,0,0};

	malloc_zone_statistics( NULL, &mstat );
	/* convert to KB */
	mstat.size_allocated >>= 10;

	return( (real)(mstat.size_allocated) );
}
#else
real GetMemUsed( void ) {
	return( 0.0 );
}
#endif

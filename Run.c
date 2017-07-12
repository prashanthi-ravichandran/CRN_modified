/* CardioWave - main program, generate by cw_run.pl */

/* MemCRN.o Domain123.o TstepEMrk23_tol.o StimTrain.o OutputDump.o ParallelNone.o DebugNone.o CWaveKernel.o VectorOps.o StencilOps1D.o */

#include "CardioWave.h"
#include <sys/times.h>
#include <limits.h>

static char* resources[MAX_RESOURCES] = {
  "xdim=1",  "ydim=1",  "zdim=1",  "outputfile=default.out",
  "beta=2000.0",  "debug=0",  "loadstate=no",  "savestate=no",
  "checkpoint=no",  "defaultnodetype=0",
  NULL
};


/* main program */
int main( int argc, char** argv ) {
  int i,err,ckpt;
  real t;
  vector Vm = _INIT_VECTOR;
  vector Vx = _INIT_VECTOR;
  vector Q  = _INIT_VECTOR;
  vector Av = _INIT_VECTOR;
  clock_t t0,t1,t2;
  real   m0,m1;
  struct tms tmsbuf;

  t0 = times( &tmsbuf );
  m0 = GetMemUsed();

  MPI_Init( &argc, &argv );

  /* set up for multiple membrane models */
  for(i=0;i<MAX_MEMBRANE;i++) {
    PatchSize[i]     = 0;
    AuxiliarySize[i] = 0;
    MemParamSize[i]  = 0;
    FunctionTable[i] = NULL;
    SetpatchTable[i] = NULL;
  }
  NumMem = 1;

  /* append any command line settings/files to the list */
  err = AppendResources( resources, argc, argv );
  if( err < 0 ) {
    printf("Error in AppendResources() = %i\n",err);
    return( -1 );
  }
  err = InitSimulation(resources);
  if( err < 0 ) {
    printf("Error in InitSimulation() = %i\n",err);
    return( -2 );
  }
  err = InitParallel(resources);
  if( err < 0 ) {
    printf("Error in InitParallel() = %i\n",err);
    return( -3 );
  }
  err = InitDebug(resources);
  if( err < 0 ) {
    printf("Error in InitDebug() = %i\n",err);
    return( -4 );
  }
  err = InitDomain(resources);
  if( err < 0 ) {
    printf("Error in InitDomain() = %i\n",err);
    return( -5 );
  }
  err = InitMembrane_CRN(resources);
  if( err < 0 ) {
    printf("Error in InitMembrane_CRN() = %i\n",err);
    return( -6 );
  }
  err = ComputeSizes(resources);
  if( err < 0 ) {
    printf("Error in ComputeSizes() = %i\n",err);
    return( -7 );
  }
  err = InitStimulus(resources);
  if( err < 0 ) {
    printf("Error in InitStimulus() = %i\n",err);
    return( -8 );
  }
  err = InitOutput_Dump(resources);
  if( err < 0 ) {
    printf("Error in InitOutput_Dump() = %i\n",err);
    return( -9 );
  }
  err = InitTimeStepper(resources,&t);
  if( err < 0 ) {
    printf("Error in InitTimeStepper() = %i\n",err);
    return( -10 );
  }

  err = GetInitialValues(&Vm,&Vx,&Q,&Av);
  if( err < 0 ) {
    printf("Error in GetInitialValues() = %i\n",err);
    return( -11 );
  }

  m1 = GetMemUsed();
  if( DebugLevel > 1 ) {
    printf("PE%i: memory detail: %10.3lf - %10.3lf = %10.3lf KB\n",
           SelfPE,m1,m0,m1-m0);
  }
  if( DebugLevel > 0 ) {
    m1 = m1 - m0;
    MPI_Allreduce( &m1, &m0, 1, MPI_SSREAL, MPI_SUM, MPI_COMM_WORLD );
    if( SelfPE == 0 ) {
      printf("Memory Used = %10.3lf MB\n",m0/1024.0);
    }
  }

  if( DebugLevel > 0 ) {
    fflush( stdout );
  }

  MPI_Barrier( MPI_COMM_WORLD );
  t1 = times( &tmsbuf );

  /* output first time-step */
  if( t == 0.0 ) {
    Output_Dump( t, Vm, Vx, Q, Av );
  }

  /* main loop */
  err = 0;
  ckpt = 0;
  while( (err>=0) and (ckpt>=0) ) {
    /* update the variables & check for errors */
    err = Update( &t, Vm, Vx, Q, Av );

    Output_Dump( t, Vm, Vx, Q, Av );

    /* AdaptGrid( t, Vm, Vx, Q, Av ); */

    ckpt = Checkpoint( t, Vm, Vx, Q, Av );
  }

  /* all done, call the exit routines */
  ExitDebug();
  ExitDomain();
  ExitMembrane_CRN();
  ExitStimulus();
  ExitOutput_Dump();
  ExitTimeStepper();
  ExitSimulation( Vm, Vx, Q, Av );

  MPI_Barrier( MPI_COMM_WORLD );
  t2 = times( &tmsbuf );

  MPI_Finalize();

  if( SelfPE == 0 ) {
    printf("Init Time = %10.3lf\n",(real)(t1-t0)/(real)(CLK_TCK) );
    printf("Run Time  = %10.3lf\n",(real)(t2-t1)/(real)(CLK_TCK) );
  }

  return( 0 );
}



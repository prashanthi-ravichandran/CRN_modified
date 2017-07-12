/* output module which dumps the entire domain, in binary, to a file */

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */



#include "CardioWave.h"

#if defined(_WT)
	extern char* Active;
#endif

/* local function prototypes */
void ExitOutput_Dump( void );

static char* RCSID = "$Id: OutputDump.c 14 2007-05-11 14:57:55Z jbp $";

static char* filename = NULL;
static int   num_outputs = 0;
static logical print_vm = true;
static FILE* fp_vm = NULL;
static char  fname_vm[MAX_FILENAME];
static logical print_ve = true;
static FILE* fp_ve = NULL;
static char  fname_ve[MAX_FILENAME];
static logical print_vi = false;
static FILE* fp_vi = NULL;
static char  fname_vi[MAX_FILENAME];
static logical print_vb = false;
static FILE* fp_vb  = NULL;
static char  fname_vb[MAX_FILENAME];
static logical print_q = false;
static FILE* fp_q  = NULL;
static char  fname_q[MAX_FILENAME];
static logical print_aux = false;
static FILE* fp_a  = NULL;
static char  fname_aux[MAX_FILENAME];
static logical print_wt = false;
static FILE* fp_wt = NULL;
static logical print_time = false;
static FILE* fp_time = NULL;
static logical print_hdr = true;
static char  fname_hdr[MAX_FILENAME];
static FILE* fp_hdr = NULL;
static real  Tspacing = 0.1;
static real  Tnext    = 0.0;
static logical AppendOutput = false;
static logical SavePEfiles = false;
static real* statedataarea = NULL;
static real* auxvardataarea   = NULL;

static rword resources[] = {
	{ "appendout",	6040 },
	{ "appendoutput",6040 },
	{ "dump_a",	6011 },
	{ "dump_aux",	6011 },
	{ "dump_b",	6024 },
	{ "dump_bath",	6024 },
	{ "dump_q",	6010 },
	{ "dump_spacing",6001 },
	{ "dump_save",	6030 },
	{ "dump_savepe",6030 },
	{ "dump_time",	6025 },
	{ "dump_tnext",	6099 },
	{ "dump_vb",	6024 },
	{ "dump_ve",	6023 },
	{ "dump_vi",	6022 },
	{ "dump_vm",	6021 },
	{ "dump_wt",	6026 },
	{ "dumpa",	6011 },
	{ "dumpaux",	6011 },
	{ "dumpb",	6024 },
	{ "dumpbath",	6024 },
	{ "dumpq",	6010 },
	{ "dumpspacing",6001 },
	{ "dumptime",	6025 },
	{ "dumptnext",	6099 },
	{ "dumpvb",	6024 },
	{ "dumpve",	6023 },
	{ "dumpvi",	6022 },
	{ "dumpvm",	6021 },
	{ "dumpwt",	6026 },
	{ "dump_hdr",	6050 },
	{ "dumphdr",	6050 },
	{ "dump_header",6050 },
	{ "dumpheader",	6050 },
	{ "outfile",	6003 },
	{ "outfilename",6003 },
	{ "outputfile",	6003 },
	{ "outputfilename",6003 },
	{ "outputspacing",6001 },
	{ "outspacing",	6001 },
	{ NULL, 0 }
};

int InitOutput_Dump( char** res ) {
	int i;
	int cmd;
	char fn[MAX_FILENAME];

	DebugEnter( "InitOutput_Dump" );

	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		switch( cmd ) {
			case 6001:
				Tspacing = GetRealValue( res[i] );
				break;
			case 6003:
				filename = GetStringValue( res[i] );
				break;
			case 6010:
				print_q = GetTFValue( res[i] );
				break;
			case 6011:
				print_aux = GetTFValue( res[i] );
				break;
			case 6021:
				print_vm = GetTFValue( res[i] );
				break;
			case 6022:
				print_vi = GetTFValue( res[i] );
				break;
			case 6023:
				print_ve = GetTFValue( res[i] );
				break;
			case 6024:
				print_vb = GetTFValue( res[i] );
				break;
			case 6025:
				print_time = GetTFValue( res[i] );
				break;
			case 6026:
				print_wt = GetTFValue( res[i] );
				break;
			case 6030:
				SavePEfiles = GetTFValue( res[i] );
				break;
			case 6040:
				AppendOutput = GetTFValue( res[i] );
				break;
			case 6050:
				print_hdr = GetTFValue( res[i] );
				break;
			case 6099:
				Tnext = GetRealValue( res[i] );
				break;
		}
		i++;
	}

	if( filename == NULL ) {
		return( -1 );
	}

	/* make sure we don't try to print out non-existant values */
	#if !defined(_WT)
		print_wt = false;
	#endif
	if( SimType == Monodomain ) {
		print_vi = false;
		print_ve = false;
	}
	if( (SimType==ReducedBidomain) or (SimType==ReducedBidomainBath) ) {
		print_vi = false;
	}
	if( (SimType!=ReducedBidomainBath) and (SimType!=FullBidomainBath) ) {
		print_vb = false;
	}

	if( print_time and (SelfPE!=0) ) {
		print_time = false;
	}

	if( print_hdr and (SelfPE!=0) ) {
		print_hdr = false;
	}

	/* we always overwrite the old temp files */
	/* then look at AppendOutput in ExitOutput() */
	if( print_hdr ) {
		sprintf( fname_hdr, "%s.dhdr", filename );
		fp_hdr = fopen( fname_hdr, "w" );
		if( fp_hdr == NULL ) {
			printf("Error opening file [%s]\n",fname_hdr);
			return( -2 );
		}
	}
	if( print_vm ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_vm, "%s.vm", filename );
			fp_vm = fopen( fname_vm, "w" );
			fclose( fp_vm );
		}
		sprintf( fname_vm, "%s.%i.vm", filename, SelfPE );
		fp_vm = fopen( fname_vm, "w+" );
		if( fp_vm == NULL ) {
			printf("Error opening file [%s]\n",fname_vm);
			return( -2 );
		}
	}
	if( print_vi ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_vi, "%s.vi", filename );
			fp_vi = fopen( fname_vi, "w" );
			fclose( fp_vi );
		}
		sprintf( fname_vi, "%s.%i.vi", filename, SelfPE );
		fp_vi = fopen( fname_vi, "w+" );
		if( fp_vi == NULL ) {
			printf("Error opening file [%s]\n",fname_vi);
			return( -3 );
		}
	}
	if( print_ve ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_ve, "%s.ve", filename );
			fp_ve = fopen( fname_ve, "w" );
			fclose( fp_ve );
		}
		sprintf( fname_ve, "%s.%i.ve", filename, SelfPE );
		fp_ve = fopen( fname_ve, "w+" );
		if( fp_ve == NULL ) {
			printf("Error opening file [%s]\n",fname_ve);
			return( -3 );
		}
	}
	if( print_vb ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_vb, "%s.vb", filename );
			fp_vb = fopen( fname_vb, "w" );
			fclose( fp_vb );
		}
		sprintf( fname_vb, "%s.%i.vb", filename, SelfPE );
		fp_vb = fopen( fname_vb, "w+" );
		if( fp_vb == NULL ) {
			printf("Error opening file [%s]\n",fname_vb);
			return( -3 );
		}
	}
	if( print_q ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_q, "%s.q", filename );
			fp_q = fopen( fname_q, "w" );
			fclose( fp_q );
		}
		sprintf( fname_q, "%s.%i.q", filename, SelfPE );
		fp_q = fopen( fname_q, "w+" );
		if( fp_q == NULL ) {
			printf("Error opening file [%s]\n",fname_q);
			return( -3 );
		}
	}
	if( print_aux ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fname_aux, "%s.aux", filename );
			fp_a = fopen( fname_aux, "w" );
			fclose( fp_a );
		}
		sprintf( fname_aux, "%s.%i.aux", filename, SelfPE );
		fp_a = fopen( fname_aux, "w+" );
		if( fp_a == NULL ) {
			printf("Error opening file [%s]\n",fname_aux);
			return( -3 );
		}
	}
	if( print_wt ) {
		if( (SelfPE==0) and (AppendOutput==false) ) {
			sprintf( fn, "%s.wt", filename );
			fp_wt = fopen( fn, "w" );
			fclose( fp_wt );
		}
		sprintf( fn, "%s.%i.wt", filename, SelfPE );
		fp_wt = fopen( fn, "w+" );
		if( fp_wt == NULL ) {
			printf("Error opening file [%s]\n",fn);
			return( -3 );
		}
	}
	if( print_time ) {
		sprintf( fn, "%s.time", filename );
		if( AppendOutput ) {
			fp_time = fopen( fn, "a" );
		} else {
			fp_time = fopen( fn, "w" );
		}
		if( fp_time == NULL ) {
			printf("Error opening file [%s]\n",fn);
			return( -4 );
		}
	}

	num_outputs = 0;

	MPI_Barrier( MPI_COMM_WORLD );

	if( (DebugLevel>0) and (SelfPE==0) ) {
		printf("OutputDump: spacing = %lf, dump = [%c,%c,%c,%c,%c,%c,%c]\n",
			Tspacing, (print_vm?'Y':'N'), (print_vi?'Y':'N'),
			(print_ve?'Y':'N'), (print_vb?'Y':'N'),
			(print_q?'Y':'N'), (print_aux?'Y':'N'), 
			(print_time?'Y':'N') );
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("OutputDump: RCSID: %s\n",RCSID);
	}

	DebugLeave( "InitOutput_Dump" );

	return( 0 );
}

void ExitOutput_Dump( void ) {
	int i,n;
	char fn[MAX_FILENAME];
	FILE* wr_vm;
	FILE* wr_vi;
	FILE* wr_ve;
	FILE* wr_vb;
	FILE* wr_q;
	FILE* wr_a;
	char* cp;

	DebugEnter( "ExitOutput_Dump" );

	if( SaveState and (SelfPE==0) ) {
		fprintf( FpResources, "dump_tnext=%le\n", Tnext );
	}


	/* stitch data together then close files */
	if( print_vm ) {
		rewind( fp_vm );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.vm", filename );
			if( AppendOutput ) {
				wr_vm = fopen( fn, "a" );
			} else {
				wr_vm = fopen( fn, "w" );
			}
		}
	}
	if( print_vi ) {
		rewind( fp_vi );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.vi", filename );
			if( AppendOutput ) {
				wr_vi = fopen( fn, "a" );
			} else {
				wr_vi = fopen( fn, "w" );
			}
		}
	}
	if( print_ve ) {
		rewind( fp_ve );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.ve", filename );
			if( AppendOutput ) {
				wr_ve = fopen( fn, "a" );
			} else {
				wr_ve = fopen( fn, "w" );
			}
		}
	}
	if( print_vb ) {
		rewind( fp_vb );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.vb", filename );
			if( AppendOutput ) {
				wr_vb = fopen( fn, "a" );
			} else {
				wr_vb = fopen( fn, "w" );
			}
		}
	}
	if( print_q ) {
		rewind( fp_q );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.q", filename );
			if( AppendOutput ) {
				wr_q = fopen( fn, "a" );
			} else {
				wr_q = fopen( fn, "w" );
			}
		}
	}
	if( print_aux ) {
		rewind( fp_a );
		if( SelfPE == 0 ) {
			sprintf( fn, "%s.aux", filename );
			if( AppendOutput ) {
				wr_a = fopen( fn, "a" );
			} else {
				wr_a = fopen( fn, "w" );
			}
		}
	}

	/* stitch data together then close files */
	for(n=0;n<num_outputs;n++) {
		if( print_vm ) {
			fread( Workspace.data, sizeof(real), TissueLocalSize,
				fp_vm );
			WriteGlobalData( Workspace.data, TissueLocalSize,
				TYPE_REAL, Tissue, wr_vm );
		}
		if( print_vi ) {
			fread( Workspace.data, sizeof(real), TissueLocalSize,
				fp_vi );
			WriteGlobalData( Workspace.data, TissueLocalSize,
				TYPE_REAL, Tissue, wr_vi );
		}
		if( print_ve ) {
			fread( Workspace.data, sizeof(real), TissueLocalSize,
				fp_ve );
			WriteGlobalData( Workspace.data, TissueLocalSize,
				TYPE_REAL, Tissue, wr_ve );
		}
		if( print_vb ) {
			fread( Workspace.data, sizeof(real), BathLocalSize,
				fp_vb );
			WriteGlobalData( Workspace.data, BathLocalSize,
				TYPE_REAL, Bath, wr_vb );
		}
		if( print_q ) {
			fread( statedataarea, sizeof(real), StateLocalSize,
				fp_q );
			WriteGlobalData( statedataarea, StateLocalSize,
				TYPE_REAL, State, wr_q );
		}
		if( print_aux ) {
			fread( auxvardataarea, sizeof(real), AuxiliaryLocalSize,
				fp_a );
			WriteGlobalData( auxvardataarea, AuxiliaryLocalSize,
				TYPE_REAL, AuxVar, wr_a );
		}
	}

	if( print_vm ) {
		fclose( fp_vm );
		if( SavePEfiles == false ) {
			remove( fname_vm );
		}
		if( SelfPE == 0 ) {
			fclose( wr_vm );
		}
	}
	if( print_vi ) {
		fclose( fp_vi );
		if( SavePEfiles == false ) {
			remove( fname_vi );
		}
		if( SelfPE == 0 ) {
			fclose( wr_vi );
		}
	}
	if( print_ve ) {
		fclose( fp_ve );
		if( SavePEfiles == false ) {
			remove( fname_ve );
		}
		if( SelfPE == 0 ) {
			fclose( wr_ve );
		}
	}
	if( print_vb ) {
		fclose( fp_vb );
		if( SavePEfiles == false ) {
			remove( fname_vb );
		}
		if( SelfPE == 0 ) {
			fclose( wr_vb );
		}
	}
	if( print_q ) {
		fclose( fp_q );
		if( SavePEfiles == false ) {
			remove( fname_q );
		}
		if( SelfPE == 0 ) {
			fclose( wr_q );
		}
	}
	if( print_aux ) {
		fclose( fp_a );
		if( SavePEfiles == false ) {
			remove( fname_aux );
		}
		if( SelfPE == 0 ) {
			fclose( wr_a );
		}
	}

	if( (SelfPE==0) and print_hdr ) {
		fprintf( fp_hdr, "data_size: %i\n", sizeof(double) );
		n = 0x1234;
		cp = (char*)(&n);
		if( cp[0] == 0x34 ) {
			fprintf( fp_hdr, "endianness: little\n" );
		} else {
			fprintf( fp_hdr, "endianness: big\n" );
		}
		fprintf( fp_hdr, "domain_size: %i\n", TissueGlobalSize );
		fprintf( fp_hdr, "state_size: %i\n", StateGlobalSize );
		fprintf( fp_hdr, "aux_size: %i\n", AuxiliaryGlobalSize );
		fprintf( fp_hdr, "num_tsteps: %i\n", num_outputs );
		n = TissueGlobalSize*num_outputs*sizeof(double);
		fprintf( fp_hdr, "vm_file_size: %i\n", n );
		n = 0;
		for(i=0;i<MAX_MEMBRANE;i++) {
			if( PatchSize[i] > n ) {
				n = PatchSize[i];
			}
		}
		fprintf( fp_hdr, "max_mem_size: %i\n", n );
		n = StateGlobalSize*num_outputs*sizeof(double);
		fprintf( fp_hdr, "state_file_size: %i\n", n );
		n = 0;
		for(i=0;i<MAX_MEMBRANE;i++) {
			if( AuxiliarySize[i] > n ) {
				n = AuxiliarySize[i];
			}
		}
		fprintf( fp_hdr, "max_aux_size: %i\n", n );
		n = AuxiliaryGlobalSize*num_outputs*sizeof(double);
		fprintf( fp_hdr, "aux_file_size: %i\n", n );
	}

	DebugLeave( "ExitOutput_Dump" );

	return;
}

void Output_Dump( real t, vector Vm, vector Vx, vector Q, vector Av ) {
	char* cp;
	vector Vi,Ve,Vb;

	DebugEnter( "Output_Dump" );

	if( t >= Tnext ) {
		DebugMark( "Output_mark" );

		num_outputs++;

		if( print_vm ) {
			fwrite( Vm.data, sizeof(real), Vm.size, fp_vm );
			fflush( fp_vm );
		}
		if( print_vi ) {
			AliasVxVi( Vx, &Vi );
			fwrite( Vi.data, sizeof(real), Vi.size, fp_vi );
			fflush( fp_vi );
		}
		if( print_ve ) {
			AliasVxVe( Vx, &Ve );
			fwrite( Ve.data, sizeof(real), Ve.size, fp_ve );
			fflush( fp_ve );
		}
		if( print_vb ) {
			AliasVxVb( Vx, &Vb );
			fwrite( Vb.data, sizeof(real), Vb.size, fp_vb );
			fflush( fp_vb );
		}
		if( print_q ) {
			fwrite( Q.data, sizeof(real), Q.size, fp_q );
			fflush( fp_q );

			/* save the state data area for use in Exit routine */
			statedataarea = Q.data;
		}
		if( print_aux ) {
			fwrite( Av.data, sizeof(real), Av.size, fp_a );
			fflush( fp_a );

			/* save the auxvar data area for use in Exit routine */
			auxvardataarea = Av.data;
		}

		#if defined(_WT)
			if( print_wt ) {
				cp = Active;
				fwrite( cp, sizeof(char), LocalSize, fp_wt );
				fflush( fp_wt );
			}
		#endif

		/* recall print_time can only be true for PE0 */
		if( print_time ) {
			fprintf( fp_time, "%lf\n", t );
			fflush( fp_time );
		}

		Tnext += Tspacing;
	}

	DebugLeave( "Output_Dump" );

	return;
}


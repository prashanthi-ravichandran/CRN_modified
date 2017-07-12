/* stimulus functions for the patch files */

/* : allows a train of stimuli to be specified easily */

/* : we will overwrite the Start arrays with the start of the NEXT train */
/*   (same for Finish array) */

/* Copyright John B. Pormann, 5 March 2001, all rights reserved */



#include "CardioWave.h"

#if defined(_WT)
	extern char* Active;
#endif

static char* RCSID = "$Id: StimTrain.c 14 2007-05-11 14:57:55Z jbp $";

static int  NumStim = 0;
static real Start[MAX_STIM];
static real Duration[MAX_STIM];
static real Period[MAX_STIM];
static real Finish[MAX_STIM];
static real Strength[MAX_STIM];
static int  NumNodes[MAX_STIM];
static int* NodeList[MAX_STIM] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
static int  IntOrExt[MAX_STIM] = {0,0,0,0,0,0,0,0};

static rword resources[] = {
	{ "dur",	3005 },
	{ "duration",	3005 },
	{ "is",		3003 },
	{ "istim",	3003 },
	{ "list",	3002 },
	{ "nlist",	3002 },
	{ "node",	3002 },
	{ "nodelist",	3002 },
	{ "nodes",	3002 },
	{ "nstim",	3001 },
	{ "numstim",	3001 },
	{ "onset",	3004 },
	{ "start",	3004 },
	{ "strength",	3003 },
	{ "tdur",	3005 },
	{ "ts",		3004 },
	{ "tstart",	3004 },
	{ "type",	3006 },
	{ "period",	3007 },
	{ "frequency",	3008 },
	{ "finish",	3009 },
	{ NULL, 0 }
};

int InitStimulus( char** res ) {
	int i,j,k,c,n;
	int cmd;
	char* bf;
	int*  nbrs;

	DebugEnter( "InitStimulus_Train" );

	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		switch( cmd ) {
			case 3001:
				NumStim = GetIntValue( res[i] );
				break;
			case 3002:
				k = FindNum( res[i] );
				c = GetNumValues( res[i] );
				if( NumNodes[k] != c ) {
					NumNodes[k] = c;
				}
				if( NodeList[k] != NULL ) {
					free( NodeList[k] );
				}
				NodeList[k] = GetIntArray( res[i] );
				break;
			case 3003:
				k = FindNum( res[i] );
				Strength[k] = GetRealValue( res[i] );
				break;
			case 3004:
				k = FindNum( res[i] );
				Start[k] = GetRealValue( res[i] );
				break;
			case 3005:
				/* just store durations for now */
				k = FindNum( res[i] );
				Duration[k] = GetRealValue( res[i] );
				break;
			case 3006:
				/* type = intra- or extra-cellular */
				k = FindNum( res[i] );
				bf = GetStringValue( res[i] );
				if( (bf[0]=='i') or (bf[0]=='I') ) {
					/* intracellular */
					IntOrExt[k] = Intracellular;
				} else {
					/* extracellular */
					IntOrExt[k] = Extracellular;
				}
				break;
			case 3007:
				k = FindNum( res[i] );
				Period[k] = GetRealValue( res[i] );
				break;
			case 3008:
				k = FindNum( res[i] );
				Period[k] = 1.0/GetRealValue( res[i] );
				break;
			case 3009:
				k = FindNum( res[i] );
				Finish[k] = GetRealValue( res[i] );
				break;
		}
		i++;
	}

	/* remove global node numbers that fall outside this PE's local area */
	/* : just converts from node numbers to local node numbers	*/
	/*   (not to matrix row numbers!)				*/
	for(i=0;i<NumStim;i++) {
		/* see how many stim nodes this PE has in its local domain */
		c = 0;
		for(j=0;j<NumNodes[i];j++) {
			n = Global2Local( Tissue, NodeList[i][j] );
			if( (n>=0) and (n<TissueLocalSize) ) {
				c++;
			}
		}
		/* now collect the local nodes and compact the list */
		k = 0;
		for(j=0;j<c;j++) {
		  while( true ) {
			/* check the k-th nodes in original list */
			n = Global2Local( Tissue, NodeList[i][k] );
			if( (n>=0) and (n<TissueLocalSize) ) {
				/* store it in j-th position */
				NodeList[i][j] = n;
				break;
			}
			k++;
		  }
		  k++;
		}
		NumNodes[i] = c;
	}

	/* NodeList[][] now holds the matrix rows to apply the given */
	/* stimuli to.  Thus, ApplyStimulus() can just run through   */
	/* this list without doing any conversions on-the-fly.       */

	#if defined(_WT)
		/* need to set up the Active list to reflect the initial stim */
		for(i=0;i<NumStim;i++) {
			if( Start[i] == 0.0 ) {
				for(j=0;j<NumNodes[i];j++) {
					/* get neighbors to node c */
					n = NodeList[i][j];
					MatGetNonzeros( MatrixINT, n, &c, &nbrs );
					for(k=0;k<c;k++) {
						Active[nbrs[k]] |= 1;
					}
				}
			}
		}
	#endif

	if( (DebugLevel>0) and (SelfPE==0) ) {
		if( NumStim < MAX_STIM ) {
			if( DebugLevel > 2 ) {
				printf("StimTrain: nstim = %i: ",NumStim);
				for(i=0;i<NumStim;i++) {
					printf("%i ",NumNodes[i]);
				}
				printf("\n");
			} else {
				printf("StimTrain: nstim = %i\n",NumStim);
			}
		} else {
			printf("StimTrain: error: too many stimuli (%i,%i)\n",
				NumStim,MAX_STIM);
		}
	}
	if( DebugLevel > 3 ) {
		for(i=0;i<NumStim;i++) {
			printf("StimTrain: PE%i: stim %i: numnodes=%i, strt=%lf, fnsh=%lf, extra=%i\n",
				SelfPE,i,NumNodes[i],Start[i],
				Finish[i],IntOrExt[i] );
		}
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("StimTrain: RCSID: %s\n",RCSID);
	}

	DebugLeave( "InitStimulus_Train" );

	return( 0 );
}

void ExitStimulus( void ) {
	return;
}

void ApplyStimulus( int IorE, real t, real alpha, vector fv, vector fq ) {
	int  i,j,k,c,n;
	int* list;
	int* nbrs;
	real t2;

	DebugEnter( "ApplyStimulus_Train" );

	for(i=0;i<NumStim;i++) {
		if( (IntOrExt[i]==IorE) and (t>=Start[i]) and (t<=Finish[i]) ) {
		   /* check if we are in the "on" phase of the train */
		   t2 = fmod( (t-Start[i]), Period[i] );
		   if( t2 <= Duration[i] ) {
			list = NodeList[i];
			for(j=0;j<NumNodes[i];j++) {
				n = list[j];
				fv.data[n] += alpha*Strength[i];

				#if defined(_WT)
					MatGetNonzeros( MatrixINT, n, &c, &nbrs );
					for(k=0;k<c;k++) {
						Active[nbrs[k]] |= 1;
					}
				#endif
			}
		   }
		}
	}

	DebugLeave( "ApplyStimulus_Train" );

	return;
}

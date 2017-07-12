/* Tstep Module for the Explicit Monodomain case	*/
/*							*/
/* uses a 2nd-3rd order Runge-Kutta-Fehlberg method	*/
/*	converted from the Matlab code 'ode23.m'	*/
/*	(which is copyright The Mathworks, Inc., 1984)	*/

/* RCSID: $Id: TstepEMrk23.c 14 2007-05-11 14:57:55Z jbp $ */


#include "CardioWave.h"

/* global variables */
sim_t SimType = Monodomain;

/* local variables */
static real   TT = 0.0;
static real   Tf = 5.000;
static real   dT = 0.001;
static vector Fv1 = _INIT_VECTOR;
static vector Fq1 = _INIT_VECTOR;
static vector Fv2 = _INIT_VECTOR;
static vector Fq2 = _INIT_VECTOR;
static vector Fv3 = _INIT_VECTOR;
static vector Fq3 = _INIT_VECTOR;
static vector Vtmp = _INIT_VECTOR;
static vector Qtmp = _INIT_VECTOR;
static real   DT_min = -1.0;
static real   DT_max = -1.0;
static FILE*  fp_time = NULL;

static rword resources[] = {
	{ "deltat",	2001 },
	{ "dt",		2001 },
	{ "dt0",	2001 },
	{ "dt_max",	2011 },
	{ "dt_min",	2010 },
	{ "dtmax",	2011 },
	{ "dtmin",	2010 },
	{ "t0",		2003 },
	{ "tf",		2002 },
	{ "tfile",	2020 },
	{ "tfinal",	2002 },
	{ "timefile",	2020 },
	{ "tinit",	2003 },
	{ "tinitial",	2003 },
	{ NULL, 0 }
};

int InitTimeStepper( char** res, real* t ) {
	int   i,cmd;
	real* x;
	real* y;
	char* filename = NULL;

	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		switch( cmd ) {
			case 2001:
				dT = GetRealValue( res[i] );
				break;
			case 2002:
				Tf = GetRealValue( res[i] );
				break;
			case 2003:
				TT = GetRealValue( res[i] );
				break;
			case 2010:
				DT_min = GetRealValue( res[i] );
				break;
			case 2011:
				DT_max = GetRealValue( res[i] );
				break;
			case 2020:
				filename = GetStringValue( res[i] );
				break;
		}
		i++;
	}

	if( DT_max < 0.0 ) {
		DT_max = (1.0e2)*dT;
	}
	if( DT_min < 0.0 ) {
		DT_min = (1.0e-2)*dT;
	}

	CreateMatrix();

	if( vecalloc(&Fv1,Tissue) < 0 ) {
		return( -1 );
	}
	if( vecalloc(&Fq1,State) < 0 ) {
		return( -2 );
	}
	if( vecalloc(&Fv2,Tissue) < 0 ) {
		return( -3 );
	}
	if( vecalloc(&Fq2,State) < 0 ) {
		return( -4 );
	}
	if( vecalloc(&Fv3,Tissue) < 0 ) {
		return( -5 );
	}
	if( vecalloc(&Fq3,State) < 0 ) {
		return( -6 );
	}
	if( vecalloc(&Vtmp,Tissue) < 0 ) {
		return( -7 );
	}
	if( vecalloc(&Qtmp,State) < 0 ) {
		return( -8 );
	}

	if( (filename!=NULL) and (SelfPE==0) ) {
		if( LoadState ) {
			fp_time = fopen( filename, "a" );
		} else {
			fp_time = fopen( filename, "w" );
		}
	}
	
	Assert(0);

	*t = TT;

	if( (DebugLevel>0) and (SelfPE==0) ) {
		printf("TstepEMrk23: dt = [%le,%le,%le], tf = %le\n",
			dT,DT_min,DT_max,Tf);
	}

	return( 0 );
}

/* what to do at exit of program */
void ExitTimeStepper( void ) {
	if( fp_time != NULL ) {
		fclose( fp_time );
	}

	if( SaveState and (SelfPE==0) ) {
		fprintf( FpResources, "t0=%le\ndt0=%le\n", TT, dT );
	}

	return;
}

int CreateMatrix( void ) {
	/* : for Explicit-Monodomain, no work needs be done   */
	/*   on any matrices, but we do need to compute 1/Vol */
	if( InvVolume.size == 0 ) {
		/* : since we don't really need Volume directly we'll */
		/*   just swap instead of allocating and copying      */
		vswap( &Volume, &InvVolume );
		vrecip( 1.0, InvVolume );
	}

	return( 0 );
}

int Update( real* t, vector Vm, vector Vx, vector Q, vector Av ) {
	int  flag;
	real delta_v,delta_q,tau_v,tau_q;
	real delta,tau;

	flag = 1;
	while( flag == 1 ) {
		/* make sure we don't overrun the endpoint */
		if( (TT+dT) > Tf ) {
			dT = Tf - TT;
		}

		/* get F1=F(X) */
		/* : compute the diffusion currents */
		sprvec( 1.0, MatrixINT, Vm, 0.0, Fv1 );
		/* : adjust for voronoi volumes/areas */
		vvscale( 1.0/Beta, InvVolume, Fv1 );
		/* : now, add the ionic currents into Fv & Fq */
		GetF( TT, dT, Vm, Q, Fv1, Fq1, Av );
		/* : this applies 1.0*Istim */
		ApplyStimulus( Intracellular, TT, 1.0, Fv1, Fq1 );

		/* compute F2=F(X+dt*F1) */
		vcopy( Vm, Vtmp );
		vcopy( Q,  Qtmp );
		vvadd( dT/Cm, Fv1, Vtmp );
		vvadd( dT, Fq1, Qtmp );
		/* : compute the diffusion currents */
		sprvec( 1.0, MatrixINT, Vtmp, 0.0, Fv2 );
		/* : adjust for voronoi volumes/areas */
		vvscale( 1.0/Beta, InvVolume, Fv2 );
		/* : now, add the ionic currents into Fv & Fq */
		GetF( TT+dT, dT, Vtmp, Qtmp, Fv2, Fq2, Av );
		/* : this applies 1.0*Istim */
		ApplyStimulus( Intracellular, TT+dT, 1.0, Fv2, Fq2 );

		/* compute F3=F(X+dt*(F1+F2)/4) */
		vcopy( Vm, Vtmp );
		vcopy( Q,  Qtmp );
		vvadd( dT/(Cm*4.0), Fv1, Vtmp );
		vvadd( dT/4.0, Fq1, Qtmp );
		vvadd( dT/(Cm*4.0), Fv2, Vtmp );
		vvadd( dT/4.0, Fq2, Qtmp );
		/* : compute the diffusion currents */
		sprvec( 1.0, MatrixINT, Vtmp, 0.0, Fv3 );
		/* : adjust for voronoi volumes/areas */
		vvscale( 1.0/Beta, InvVolume, Fv3 );
		/* : now, add the ionic currents into Fv & Fq */
		GetF( TT+0.250*dT, dT, Vtmp, Qtmp, Fv3, Fq3, Av );
		/* : this applies 1.0*Istim */
		ApplyStimulus( Intracellular, TT+0.250*dT, 1.0, Fv3, Fq3 );

		/* estimate error */
		/* : set Vtmp = dt*(F1-2*F3+F2)/2 */
		vcopy( Fv1, Vtmp );
		vcopy( Fq1, Qtmp );
		vvaddb( dT/(Cm*2.0), Fv2, dT/(Cm*2.0), Vtmp );
		vvaddb( dT/2.0, Fq2, dT/2.0, Qtmp );
		vvadd( -2.0*dT/(Cm*2.0), Fv3, Vtmp );
		vvadd( -2.0*dT/2.0, Fq3, Qtmp );
		delta_v = normI( Vtmp );
		delta_q = normI( Qtmp );
		if( delta_v < delta_q ) {
			delta = delta_q;
		} else {
			delta = delta_v;
		}
		tau_v = normI( Vm );
		tau_q = normI( Q );
		if( tau_v < tau_q ) {
			tau = tau_q;
		} else {
			tau = tau_v;
		}
		if( tau < 1.0 ) {
			tau = 1.0;
		}
		tau *= 1.0e-4;

		if( delta <= tau ) {
			/* we're ok - error is small enough */
			/* : X(t+1) = X + dT*(F1+F2+4*F3)/6 */
			vvadd( dT/(Cm*6.0), Fv1, Vm );
			vvadd( dT/6.0, Fq1, Q );
			vvadd( dT/(Cm*6.0), Fv2, Vm );
			vvadd( dT/6.0, Fq2, Q );
			vvadd( 4.0*dT/(Cm*6.0), Fv3, Vm );
			vvadd( 4.0*dT/6.0, Fq3, Q );

			/* update time count */
			TT += dT;
			*t = TT;

			flag = 0;
		}

		if( delta != 0.0 ) {
			/* taken straight from Matlab's 'ode23.m' code */
			dT = 0.9*dT*pow(tau/delta,1.0/3.0);
			if( dT > DT_max ) {
				dT = DT_max;
			}
			if( dT < DT_min ) {
				dT = DT_min;
			}
		}
	}

	if( (fp_time!=NULL) and (SelfPE==0) ) {
		fprintf( fp_time, "%lf\t%lf\n", TT, dT );
	}

	if( TT >= Tf ) {
		return( -1 );
	}

	return( 0 );
}


/********************************************************************************
Action Potential Model of a Human Atrial Cell

Program Written By Joseph Tranquillo
First Draft on 4/29/99

Program corrected and benchmarked against below reference
by Robert Oliver on 9/14/00


From The Journal Article : "Ionic Mechanisms Underlying Human Atrial Action
	Properties: Insights from a Mathematical Model"

Authors: Courtemanche,M. , Ramirez,R. and Nattel,S.

Am. J. Physiol. 275 (Heart Circ. Physiol. 44); pages H301-H321. 1998.

Note: has been converted to standard units

Program revised by Chad Johnson on 10/11/02:
	Changes were minor: Scale Iion by membrane surface area, but not
		      	    individual conductances or other constants.
			    This was done to avoid slight truncation errors.
	This membrane module was verified by comparing patch output to the
	output from code known to be accurate. (crj)

Written by Joe Tranquillo

********************************************************************************/



#include "CardioWave.h"

void SetPatches_CRN( domain_t tp, vector Vm, vector Q );
int GetF_CRN( real t, real dt, vector Vm, vector Qv, vector Fv, vector Fq, vector Av );

/* LMCG patch */
typedef struct {
	real C1_RyR;
	real O1_RyR;
	real C2_RyR;
	real O2_RyR;
	real C0;
	real C1;
	real C2;
	real C3;
	real C4;
	real Open;
	real CCa0;
	real CCa1;
	real CCa2;
	real CCa3;
	real CCa4;
	real yCa;
	real HTRPNCa;
	real LTRPNCa;
	real m,h,j,sa,si,ua,ui,xr,xs,d,f,fca,u,v,w;
	real Nai,Ki;
	real Cai;
	real CaNSR;
	real CaSS;
	real CaJSR;
} CRN_patch;

typedef struct{
	real Inai,Ik1i,Itoi,Ikuri,Ikri,Iksi;
	real Icali,Inaki,Inacai,Ibcai,Ibnai;
	real Ipcai,Jreli,Jtri,Jupi,Jupleaki;
	real Iioni;
	real Caii,CaNSRi,CaSSi,CaJSRi;
} CRN_auxvar;

/* globals */
int  CRN_PatchSize     = sizeof(CRN_patch)/sizeof(real);
int  CRN_AuxiliarySize = sizeof(CRN_auxvar)/sizeof(real);
int  CRN_MemParamSize  = 0;
byte CRN_NodeType = ByteError;

static char* RCSID = "$Id: MemCRN.c 14 2007-05-11 14:57:55Z jbp $";

/* constants */
static real cellLength = 100.0;		/* um */
static real cellDiameter = 16.0;	/* um */
static real R = 8.3143;			/* J/(Kmol) */
static real T = 310.0;			/* K */
static real F = 96.4867;		/* C/mmol */
static real Vcell = 20100.0;		/* um^3 */
static real Vi = 13668.0;		/* um^3 */
static real Vup = 1109.52;		/* um^3 */
static real Vrel = 96.48;		/* um^3 */
static real Ko = 5.4;			/* mM */
static real Nao = 140.0;		/* mM */
static real Cao = 1.8;			/* mM */
static real gna = 7.8e-4;		/* mS */
static real gk1 = 0.09e-4;		/* mS */
static real gto = 0.1652e-4;		/* mS */
static real gkr = 0.0294e-4;		/* mS */
static real gks = 0.129e-4;		/* mS */
static real gcal = 0.1238e-4;		/* mS */
static real gbca = 0.00113e-4;		/* mS */
static real gbna = 0.000674e-4;		/* mS */
static real Inakmax = 0.60e-4;		/* uA */
static real Inacamax = 1600.0e-4;	/* uA */
static real Ipcamax = 0.275e-4;		/* uA */
static real Jupmax = 0.005;		/* mM/msec */
static real Kq10 = 3.0;			/* unitless */
static real lambda = 0.35;		/* unitless */
static real Kmnai = 10.0;		/* mM */
static real Kmko = 1.5;			/* mM */
static real Kmna = 87.5;		/* mM */
static real Kmca = 1.38;		/* mM */
static real ksat = 0.1;			/* unitless */
static real krel = 30.0;		/* msec^(-1) */
static real kup = 0.00092;		/* mM */
static real CaNSRmax = 15.0;		/* mM */
static real Cmdnmax = 0.05;		/* mM */
static real Trpnmax = 0.07;		/* mM */
static real Csqnmax = 10.0;		/* mM */
static real KmCmdn = 0.00238;		/* mM */
static real KmTrpn = 0.0005;		/* mM */
static real KmCsqn = 0.8;		/* mM */
static real taufca = 2.0;		/* msec */
static real tautr = 180.0;		/* msec */
static real tauu = 8.0;			/* msec */

//IMW PARAMETERS

//General Parameters
static real IMW_R    = 8314.0;
static real IMW_frdy = 96485.0;
static real IMW_temp = 310.0;
static real RT_over_F = 26.7156;

//Volumes
static real Acap = 1.534e-4;             /*  cm2  */
static real Vmyo = 25.84e-6;               /*  cm  */
static real VJSR = 0.16e-6;
static real VNSR = 2.1e-6;
static real VSS  = 1.2e-9;

//Ical
static real	fL				= 0.3;
static real gL				= 4.0;
static real bL				= 2.0;
static real aL				= 2.0;
static real omega			= 0.0025;
static real Pscale			= 7.0;
static real PCa				= 2.469e-4; //Pscale*

//Buffering
static real LTRPNtot		= 70.0e-3;
static real HTRPNtot		= 140.e-3;
static real khtrpn_plus		= 20.0;
static real khtrpn_minus	= 0.066e-3;
static real kltrpn_plus		= 40.0;
static real kltrpn_minus	= 40.0e-3;
static real CMDNtot			= 50.0e-3;
static real CSQNtot			= 15.0;
static real EGTAtot			= 0.0;
static real KmCMDN			= 2.38e-3;
static real KmCSQN			= 0.8;
static real KmEGTA			= 1.5e-4;

//Ca Handling parameters
static real v1				= 1.8;
static real Kfb				= 0.000168;
static real Krb				= 3.29;
static real KSR				= 1.2;
static real Nfb				= 1.2;
static real Nrb				= 1.0;
static real vmaxf			= 0.0748e-3;
static real vmaxr			= 0.318e-3;
//static real tautr			= 0.5747; //IMW tautr
static real tauxfer			= 26.7;
static real kaplus			= 0.01215;
static real kaminus			= 0.576;
static real kbplus			= 0.00405;
static real kbminus			= 1.93;
static real kcplus			= 0.1;
static real kcminus			= 0.0008;
static real ncoop			= 4.0;
static real mcoop			= 3.0;

static real CRN_RestVoltage = -81.2;	/* mV */

static CRN_patch CRN_RestPatch = {
		0.4448038946E+00,
		0.6601783287E-03,
		0.5545359240E+00,
		0.3973920060E-08,
		0.8621934054E+00,
		0.1141398211E-01,
		0.5666325225E-04,
		0.1250213570E-06,
		0.1034426748E-09,
		0.7757981563E-11,
		0.1198547081E+00,
		0.6346794302E-02,
		0.1260326488E-03,
		0.1112315238E-05,
		0.3681315892E-08,
		0.9997157074E+00,
		0.9772152481E+00,
		0.8046584973E-01,
		2.91e-3,0.965,0.978,3.04e-2,0.999,
	  4.96e-3,0.999,3.29e-5,1.87e-2,1.37e-4,
	  0.999,0.775,0.0,1.0,0.999,11.2,139.0,
	  1.02e-4, //Cai CRN
		1.49, //CaNSR CRN
		1.02e-4, //CaSS CRN
		1.49 //CaJSR CRN
		/*0.8601192016E-04, //Cai IMW
		0.2855294915E+00, //CaNSR IMW
		0.1420215245E-03, //CaSS
		0.2852239446E+00 //CaJSR*/
	};


static rword resources[] = {
	{ "CRN_IV",	  1001 },
	{ "CRN_Patch",	  1001 },
	{ "CRN_Vr",	  1002 },
	{ "CRN_Vrest",	  1002 },
	{ "CRN_gbca",	  1003 },
	{ "CRN_gbna",	  1004 },
	{ "CRN_gcal",	  1005 },
	{ "CRN_gk1",	  1006 },
	{ "CRN_gkr",	  1007 },
	{ "CRN_gks",	  1008 },
	{ "CRN_gna",	  1009 },
	{ "CRN_gto",	  1010 },
	{ "CRN_nodetype", 1100 },
	{ "CRN_node",	  1100 },
	{ NULL, 0 }
};



int InitMembrane_CRN( char** res ) {
	int i,j,c;
	int cmd;
	real* iv;
	real* p;

	DebugEnter( "InitMembrane_CRN" );

	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		switch( cmd ) {
			case 1001:
				iv = GetRealArray( res[i] );
				p = (real*)(&CRN_RestPatch);
				c  = GetNumValues( res[i] );
				if( c > CRN_PatchSize ) {
					c = CRN_PatchSize;
				}
				for(j=0;j<c;j++) {
					p[j] = iv[j];
				}
				break;
			case 1002:
				CRN_RestVoltage = GetRealValue( res[i] );
				break;
			case 1003:
				gbca = GetRealValue( res[i] );
				break;
			case 1004:
				gbna = GetRealValue( res[i] );
				break;
			case 1005:
				gcal = GetRealValue( res[i] );
				break;
			case 1006:
				gk1 = GetRealValue( res[i] );
				break;
			case 1007:
				gkr = GetRealValue( res[i] );
				break;
			case 1008:
				gks = GetRealValue( res[i] );
				break;
			case 1009:
				gna = GetRealValue( res[i] );
				break;
			case 1010:
				gto = GetRealValue( res[i] );
				break;
			case 1100:
				CRN_NodeType = GetByteValue( res[i] );
				break;
		}
		i++;
	}

        /* required!! */
	if( CRN_NodeType == ByteError ) {
		return( -1 );
	}
	PatchSize[CRN_NodeType]     = CRN_PatchSize;
	if( UseAuxvars ) {
        	AuxiliarySize[CRN_NodeType] = CRN_AuxiliarySize;
	} else {
        	AuxiliarySize[CRN_NodeType] = 0;
	}
	MemParamSize[CRN_NodeType]  = CRN_MemParamSize;
	FunctionTable[CRN_NodeType] = GetF_CRN;
	SetpatchTable[CRN_NodeType] = SetPatches_CRN;

	RegisterOffset( OffsetOf(CRN_patch,m), CRN_NodeType, StateVar, "CRN_M" );
	RegisterOffset( OffsetOf(CRN_patch,h), CRN_NodeType, StateVar, "CRN_H" );
	RegisterOffset( OffsetOf(CRN_patch,j), CRN_NodeType, StateVar, "CRN_J" );
	RegisterOffset( OffsetOf(CRN_patch,sa), CRN_NodeType, StateVar, "CRN_SA" );
	RegisterOffset( OffsetOf(CRN_patch,si), CRN_NodeType, StateVar, "CRN_SI" );
	RegisterOffset( OffsetOf(CRN_patch,ua), CRN_NodeType, StateVar, "CRN_ua" );
	RegisterOffset( OffsetOf(CRN_patch,ui), CRN_NodeType, StateVar, "CRN_ui" );
	RegisterOffset( OffsetOf(CRN_patch,xr), CRN_NodeType, StateVar, "CRN_xr" );
	RegisterOffset( OffsetOf(CRN_patch,xs), CRN_NodeType, StateVar, "CRN_xs" );
	RegisterOffset( OffsetOf(CRN_patch,d), CRN_NodeType, StateVar, "CRN_d" );
	RegisterOffset( OffsetOf(CRN_patch,f), CRN_NodeType, StateVar, "CRN_f" );
	RegisterOffset( OffsetOf(CRN_patch,fca), CRN_NodeType, StateVar, "CRN_fca" );
	RegisterOffset( OffsetOf(CRN_patch,u), CRN_NodeType, StateVar, "CRN_u" );
	RegisterOffset( OffsetOf(CRN_patch,v), CRN_NodeType, StateVar, "CRN_v" );
	RegisterOffset( OffsetOf(CRN_patch,w), CRN_NodeType, StateVar, "CRN_w" );
	RegisterOffset( OffsetOf(CRN_patch,Nai), CRN_NodeType, StateVar, "CRN_Nai" );
	RegisterOffset( OffsetOf(CRN_patch,Ki), CRN_NodeType, StateVar, "CRN_Ki" );
	RegisterOffset( OffsetOf(CRN_patch,Cai), CRN_NodeType, StateVar, "CRN_Cai" );
	RegisterOffset( OffsetOf(CRN_patch,CaNSR), CRN_NodeType, StateVar, "CRN_CaNSR" );
	RegisterOffset( OffsetOf(CRN_patch,CaSS), CRN_NodeType, StateVar, "CRN_CaSS" );
	RegisterOffset( OffsetOf(CRN_patch,CaJSR), CRN_NodeType, StateVar, "CRN_CaJSR" );


	RegisterOffset( OffsetOf(CRN_auxvar,Inai), CRN_NodeType, AuxVar, "CRN_Inai" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ik1i), CRN_NodeType, AuxVar, "CRN_Ik1i" );
	RegisterOffset( OffsetOf(CRN_auxvar,Itoi), CRN_NodeType, AuxVar, "CRN_Itoi" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ikuri), CRN_NodeType, AuxVar, "CRN_Ikuri" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ikri), CRN_NodeType, AuxVar, "CRN_Ikri" );
	RegisterOffset( OffsetOf(CRN_auxvar,Iksi), CRN_NodeType, AuxVar, "CRN_Iksi" );
	RegisterOffset( OffsetOf(CRN_auxvar,Icali), CRN_NodeType, AuxVar, "CRN_Icali" );
	RegisterOffset( OffsetOf(CRN_auxvar,Inaki), CRN_NodeType, AuxVar, "CRN_Inaki" );
	RegisterOffset( OffsetOf(CRN_auxvar,Inacai), CRN_NodeType, AuxVar, "CRN_Inacai" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ibcai), CRN_NodeType, AuxVar, "CRN_Ibcai" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ibnai), CRN_NodeType, AuxVar, "CRN_Ibnai" );
	RegisterOffset( OffsetOf(CRN_auxvar,Ipcai), CRN_NodeType, AuxVar, "CRN_Ipcai" );
	RegisterOffset( OffsetOf(CRN_auxvar,Jreli), CRN_NodeType, AuxVar, "CRN_Jreli" );
	RegisterOffset( OffsetOf(CRN_auxvar,Jtri), CRN_NodeType, AuxVar, "CRN_Jtri" );
	RegisterOffset( OffsetOf(CRN_auxvar,Jupi), CRN_NodeType, AuxVar, "CRN_Jupi" );
	RegisterOffset( OffsetOf(CRN_auxvar,Jupleaki), CRN_NodeType, AuxVar, "CRN_Jupleaki" );
	RegisterOffset( OffsetOf(CRN_auxvar,Iioni), CRN_NodeType, AuxVar, "CRN_Iioni" );
	RegisterOffset( OffsetOf(CRN_auxvar,Caii), CRN_NodeType, AuxVar, "CRN_Caii" );
	RegisterOffset( OffsetOf(CRN_auxvar,CaNSRi), CRN_NodeType, AuxVar, "CRN_CaNSRi" );
	RegisterOffset( OffsetOf(CRN_auxvar,CaSSi), CRN_NodeType, AuxVar, "CRN_CaSSi" );
	RegisterOffset( OffsetOf(CRN_auxvar,CaJSRi), CRN_NodeType, AuxVar, "CRN_CaJSRi" );

	if( (DebugLevel>0) and (SelfPE==0) ) {
		printf("MemCRN: CRN_patchsize = %i\n",CRN_PatchSize);
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("MemCRN: RCSID: %s\n",RCSID);
	}

	DebugLeave( "InitMembrane_CRN" );

	return( 0 );
}


void ExitMembrane_CRN( void ) {
	return;
}


void SetPatches_CRN( domain_t tp, vector Vm, vector Q ) {
        int i,j;
        real* vm = Vm.data;
        real* qp = Q.data;
        real* rp = (real*)(&CRN_RestPatch);
        byte* nt = NodeType.data;

	DebugEnter( "SetPatches_CRN" );

	if( tp == Tissue ) {
          for(i=0;i<TissueLocalSize;i++) {
                if( nt[i] == CRN_NodeType ) {
                        vm[i] = CRN_RestVoltage;
                        for(j=0;j<CRN_PatchSize;j++) {
                                qp[j] = rp[j];
                        }
                }
                qp += PatchSize[ nt[i] ];
	  }
        }

	DebugLeave( "SetPatches_CRN" );

	return;
}


/* GetF() should fill in vectors Fv and Fq with the new currents/updates */
/* NOTE: this should be from the form:  dV/dt + Fv = 0, dq/dt + Fq = 0 */
int GetF_CRN( real t, real dt, vector Vm, vector Qv,
		vector Fv, vector Fq, vector Av ) {
  int    i;
  real*  vmptr = Vm.data;
  real*  fv = Fv.data;
  real*  qm = Qv.data;
  real*  fq = Fq.data;
  real*  av = Av.data;
  CRN_patch* qp;
  CRN_patch* fp;
  CRN_auxvar* ap;
  byte*  nt = NodeType.data;
  real vm;
  real m,h,j,sa,si,ua,ui,xr,xs,d,f,fca,u,v,w;
  real Nai,Ki,Cai,CaNSR,CaSS,CaJSR;
  real Ina,Ik1,Ito,Ikur,Ikr,Iks,Ical,Iion;
  real Inak,Inaca,Ibca,Ibna,Ipca,Jupleak;
  real Ena,Ek,Eca,gkur,sigma,fnak;
  real am,bm,aj,bj,ah,bh,taum,tauh,tauj,infm,infh,infj;
  real asa,bsa,asi,bsi,tausa,infsa,tausi,infsi;
  real aui,bui,aua,bua,tauui,infui,tauua,infua;
  real axr,bxr,infxr,tauxr,axs,bxs,tauxs,infxs;
  real taud,infd,tauf,inff,inffca;
  real infu,tauv,infv,tauw,infw,Fn,B1,B2;

	//IMW PARAMETERS
	real a1,a2,a3,a4;
	real VF_over_RT,VFsq_over_RT,ICamax;
	real fb,rb,Jup,Jrel,Jxfer,Jtr,dLTRPNCa,dHTRPNCa,Jtrpn,beta_SS,beta_JSR,beta_i;
	real dNai,dKi,dCai,dCaSS,dCaNSR,dCaJSR,dC1_RyR,dC2_RyR,dO1_RyR,dO2_RyR;
	real vi,alpha,beta,gamma,omega,alpha_prime,beta_prime;
	real C1_RyR,O1_RyR,C2_RyR,O2_RyR;
	real C0,C1,C2,C3,C4,Open,CCa0,CCa1,CCa2,CCa3,CCa4,yCa;
	real HTRPNCa,LTRPNCa;
	real C0_to_C1,C1_to_C2,C2_to_C3,C3_to_C4,CCa0_to_CCa1,CCa1_to_CCa2,CCa2_to_CCa3,CCa3_to_CCa4;
	real C1_to_C0,C2_to_C1,C3_to_C2,C4_to_C3,CCa1_to_CCa0,CCa2_to_CCa1,CCa3_to_CCa2,CCa4_to_CCa3,C0_to_CCa0;
	real C1_to_CCa1,C2_to_CCa2,C3_to_CCa3,C4_to_CCa4,CCa0_to_C0,CCa1_to_C1,CCa2_to_C2,CCa3_to_C3,CCa4_to_C4;
	real dC0,dC1,dC2,dC3,dC4,dOpen,dCCa0,dCCa1,dCCa2,dCCa3,dCCa4,yCa_inf,tau_yCa,dyCa;

  DebugEnter( "GetF_CRN" );

  for(i=0;i<TissueLocalSize;i++) {
    if( nt[i] == CRN_NodeType ) {
      /* map pointers to patch structures */
      qp = (CRN_patch*)qm;
      fp = (CRN_patch*)fq;
      ap = (CRN_auxvar*)av;

      /* get local variables */
      vm = vmptr[i];
			C1_RyR			= qp->C1_RyR;
			O1_RyR			= qp->O1_RyR;
			C2_RyR			= qp->C2_RyR;
			O2_RyR			= qp->O2_RyR;
			C0				= qp->C0;
			C1				= qp->C1;
			C2				= qp->C2;
			C3				= qp->C3;
			C4				= qp->C4;
			Open			= qp->Open;
			CCa0			= qp->CCa0;
			CCa1			= qp->CCa1;
			CCa2			= qp->CCa2;
			CCa3			= qp->CCa3;
			CCa4			= qp->CCa4;
			yCa				= qp->yCa;
			HTRPNCa			= qp->HTRPNCa;
			LTRPNCa			= qp->LTRPNCa;
      m  = qp->m;
      h  = qp->h;
      j  = qp->j;
      sa = qp->sa;
      si = qp->si;
      ua = qp->ua;
      ui = qp->ui;
      xr = qp->xr;
      xs = qp->xs;
      d  = qp->d;
      f  = qp->f;
      fca = qp->fca;
      u  = qp->u;
      v  = qp->v;
      w  = qp->w;
      Nai = qp->Nai;
      Ki  = qp->Ki;
			Cai				= qp->Cai;
			CaNSR			= qp->CaNSR;
			CaSS			= qp->CaSS;
			CaJSR			= qp->CaJSR;

//Driving Force et al. Calculations
Ena = (R*T/F)*log(Nao/Nai);
Ek = (R*T/F)*log(Ko/Ki);
Eca = (R*T/(2.0*F))*log(Cao/Cai);
gkur = 1.0e-4*(0.005+0.05/(1.0+exp(-1.0*(vm-15.0)/13.0)));
sigma = (1.0/7.0)*(exp(Nao/67.3)-1.0);
fnak = 1.0/(1.0+0.1245*exp(-0.1*F*vm/(R*T))+0.0365*sigma*exp(-1.0*F*vm/(R*T)));

//CRN Currents

Ina = gna*m*m*m*h*j*(vm-Ena);
Ik1 = gk1*(vm-Ek)/(1.0+exp(0.07*(vm+80.0)));
Ito = gto*sa*sa*sa*si*(vm-Ek);
Ikur = gkur*ua*ua*ua*ui*(vm-Ek);
Ikr = gkr*xr*(vm-Ek)/(1.0 + exp((vm+15.0)/22.4));
Iks = gks*xs*xs*(vm-Ek);
Ical = gcal*d*f*fca*(vm-65.0);
Inak = Inakmax*fnak*(1.0/(1.0 + pow(Kmnai/Nai,1.5))) * Ko/(Ko+Kmko);
Inaca = Inacamax*(exp(lambda*F*vm/(R*T))*Nai*Nai*Nai*Cao-
			exp((lambda-1)*F*vm/(R*T))*Nao*Nao*Nao*Cai)/
	((Kmna*Kmna*Kmna + Nao*Nao*Nao)*(Kmca + Cao)*(1.0 + ksat*exp((lambda-1)*F*vm/(R*T))));
Ibca = gbca*(vm-Eca);
Ibna = gbna*(vm-Ena);
Ipca = Ipcamax*Cai/(0.0005 + Cai);

Jrel = krel*u*u*v*w*(CaJSR-Cai);
Jtr = (CaNSR-CaJSR)/tautr;
Jup = Jupmax/(1.0+(kup/Cai));
Jupleak = CaNSR*Jupmax/CaNSRmax;

//CRN Ion Fluxes

B1 = 1.0e6*((2.0*Inaca-Ipca-Ical-Ibca)/(2.0*F*Vi))+
	(Vup*(Jupleak-Jup)+Jrel*Vrel)/Vi;
B2 = 1.0+Trpnmax*KmTrpn/((Cai+KmTrpn)*(Cai+KmTrpn))+Cmdnmax*KmCmdn/
	((Cai+KmCmdn)*(Cai+KmCmdn));

dNai = 1.0e6*((-3.0*Inak-3.0*Inaca-Ibna-Ina)/(F*Vi));
dKi = 1.0e6*((2.0*Inak-Ik1-Ito-Ikur-Ikr-Iks)/(F*Vi));
dCai = B1/B2;
dCaNSR = Jup-Jupleak-Jtr*Vrel/Vup;
dCaSS = dCai;
dCaJSR = (Jtr-Jrel)*(1.0/(1.0+(Csqnmax*KmCsqn)/
				   ((CaJSR+KmCsqn)*(CaJSR+KmCsqn))));

//IMW Currents & Calcium Handling

//Ical Current Phase

	VF_over_RT=vm/RT_over_F;
	VFsq_over_RT=(IMW_frdy)*VF_over_RT;
	a1 =  1.0e-3*exp(2.0*VF_over_RT)-Cao*0.341;
	a2 =  exp(2.0*VF_over_RT)-1.0;
	ICamax = Pscale*PCa*4.0*VFsq_over_RT*(a1/a2);
	//Ical = ICamax*yCa*Open;

//Calcium Flux
	fb			= pow((Cai/Kfb),Nfb);
	rb			= pow((CaNSR/Krb),Nrb);
	//Jup			= KSR*(vmaxf*fb - vmaxr*rb)/(1.0 + fb + rb);
	//Jrel			= v1*(O1_RyR+O2_RyR)*(CaJSR-CaSS);
	//Jtr			= (CaNSR - CaJSR)/tautr;
	//Jxfer		= (CaSS-Cai)/tauxfer;

//Buffers
		a1			= kltrpn_minus * LTRPNCa;
		dLTRPNCa		= kltrpn_plus*Cai*(1.0 - LTRPNCa) - a1;
		a1			= khtrpn_minus * HTRPNCa;
		dHTRPNCa		= khtrpn_plus*Cai*(1.0 - HTRPNCa) - a1;
		Jtrpn		= LTRPNtot*dLTRPNCa+HTRPNtot*dHTRPNCa;
		a1			= CMDNtot*KmCMDN/(pow((CaSS+KmCMDN),2.0));
		a2			= EGTAtot*KmEGTA/(pow((CaSS+KmEGTA),2.0));
		beta_SS		= 1.0/(1.0+a1+a2);
		a1			= CSQNtot*KmCSQN/(pow((CaJSR+KmCSQN),2.0));
		beta_JSR		= 1.0/(1.0+a1);

		a1			= CMDNtot*KmCMDN/(pow((Cai+KmCMDN),2.0));
		a2			= EGTAtot*KmEGTA/(pow((Cai+KmEGTA),2.0));
		beta_i		= 1.0/(1.0+a1+a2);

//RyR
		a1			= pow((CaSS*1000.0),mcoop);
		a2			= pow((CaSS*1000.0),ncoop);
		dC1_RyR		= -kaplus*a2*C1_RyR+kaminus*O1_RyR;
		dO2_RyR		=  kbplus*a1*O1_RyR - kbminus*O2_RyR;
		dC2_RyR		=  kcplus*O1_RyR - kcminus*C2_RyR;
		dO1_RyR		= -(dC1_RyR + dO2_RyR + dC2_RyR);

//Ical Compute Phase

		alpha		= 4.0*1.20*0.4160* exp(0.0120*(vm-35.0));
		beta			= 4.0*0.450*0.0490* exp(-0.0650*(vm-22.0));
		alpha_prime	= aL*alpha;
		beta_prime	= beta/bL;
		C0_to_C1		= 4.0*alpha;
		C1_to_C2		= 3.0*alpha;
		C2_to_C3		= 2.0*alpha;
		C3_to_C4		= alpha;
		CCa0_to_CCa1 = 4.0*alpha_prime;
		CCa1_to_CCa2 = 3.0*alpha_prime;
		CCa2_to_CCa3 = 2.0*alpha_prime;
		CCa3_to_CCa4 = alpha_prime;
		C1_to_C0		= beta;
		C2_to_C1		= 2.0*beta;
		C3_to_C2		= 3.0*beta;
		C4_to_C3		= 4.0*beta;
		CCa1_to_CCa0 = beta_prime;
		CCa2_to_CCa1 = 2.0*beta_prime;
		CCa3_to_CCa2 = 3.0*beta_prime;
		CCa4_to_CCa3 = 4.0*beta_prime;

		gamma		= 0.6*0.092330*CaSS;
		C0_to_CCa0	= gamma;				// = gamma
		C1_to_CCa1	= aL*C0_to_CCa0;		// = gamma*aL
		C2_to_CCa2	= aL*C1_to_CCa1;		// = gamma*aL^2
		C3_to_CCa3	= aL*C2_to_CCa2;		// = gamma*aL^3
		C4_to_CCa4	= aL*C3_to_CCa3;		// = gamma*aL^4
		CCa0_to_C0	= omega;				// = omega
		CCa1_to_C1	= CCa0_to_C0/bL;		// = omega/bL
		CCa2_to_C2	= CCa1_to_C1/bL;		// = omega/bL^2
		CCa3_to_C3	= CCa2_to_C2/bL;		// = omega/bL^3
		CCa4_to_C4	= CCa3_to_C3/bL;		// = omega/bL^4

		a1			= (C0_to_C1+C0_to_CCa0)*C0;
		a2			= C1_to_C0*C1 + CCa0_to_C0*CCa0;
		dC0			= a2 - a1;
		a1			= (C1_to_C0+C1_to_C2+C1_to_CCa1)*C1;
		a2			= C0_to_C1*C0 + C2_to_C1*C2 + CCa1_to_C1*CCa1;
		dC1			= a2 - a1;
		a1			= (C2_to_C1+C2_to_C3+C2_to_CCa2)*C2;
		a2			= C1_to_C2*C1 + C3_to_C2*C3 + CCa2_to_C2*CCa2;
		dC2			= a2 - a1;
		a1			= (C3_to_C2+C3_to_C4+C3_to_CCa3)*C3;
		a2			= C2_to_C3*C2 + C4_to_C3*C4 + CCa3_to_C3*CCa3;
		dC3			= a2 - a1;
		a1			= (C4_to_C3+fL+C4_to_CCa4)*C4;
		a2			= C3_to_C4*C3 + gL*Open + CCa4_to_C4*CCa4;
		dC4			= a2 - a1;
		dOpen		= fL*C4 - gL*Open;
		a1			= (CCa0_to_CCa1+CCa0_to_C0)*CCa0;
		a2			= CCa1_to_CCa0*CCa1 + C0_to_CCa0*C0;
		dCCa0		= a2 - a1;
		a1			= (CCa1_to_CCa0+CCa1_to_CCa2+CCa1_to_C1)*CCa1;
		a2			= CCa0_to_CCa1*CCa0 + CCa2_to_CCa1*CCa2 + C1_to_CCa1*C1;
		dCCa1		= a2 - a1;
		a1			= (CCa2_to_CCa1+CCa2_to_CCa3+CCa2_to_C2)*CCa2;
		a2			= CCa1_to_CCa2*CCa1 + CCa3_to_CCa2*CCa3 + C2_to_CCa2*C2;
		dCCa2		= a2 - a1;
		a1			= (CCa3_to_CCa2+CCa3_to_CCa4+CCa3_to_C3)*CCa3;
		a2			= CCa2_to_CCa3*CCa2 + CCa4_to_CCa3*CCa4 + C3_to_CCa3*C3;
		dCCa3		= a2 - a1;
		a1			= (CCa4_to_CCa3+CCa4_to_C4)*CCa4;
		a2			= CCa3_to_CCa4*CCa3 + C4_to_CCa4*C4;
		dCCa4		= a2 - a1;

		a1			= 0.82;
		yCa_inf		= a1/(1.0+exp((vm + 28.50)/(7.80))) + (1.0-a1);
		tau_yCa		= 1.0/(0.00336336209452/(0.50+exp(vm/(-5.53899874036055))) + 0.00779046570737*exp(vm/ (-49.51039631160386))  );
		//		   tau_yCa		= 1.0/(0.00653/(0.50+exp(vm/(7.1))) + 0.00512*exp(vm/ (-39.8))  );
		dyCa			= (yCa_inf-yCa)/tau_yCa;

// COMPUTE Iion

/* Iion: scale by membrane surface area in cm^2 */
Iion = Ina + Ik1 + Ito + Ikur + Ikr + Iks + Ical + Ipca + Inak + Inaca + Ibna + Ibca;
Iion /= (M_PI*cellDiameter*cellLength*1.0e-8);

//Gating Variable Calculations
if(vm==-47.13){
am = 3.2;
}
 else{
am = 0.32*(vm+47.13)/(1.0-exp(-0.1*(vm+47.13)));
 }
bm = 0.08*exp(-vm/11.0);
taum = 1.0/(am+bm);
infm = am*taum;


if(vm>=-40.0){
ah = 0.0;
bh = 1.0/(0.13*(1.0+exp(-1.0*(vm+10.66)/11.1)));
aj = 0.0;
bj = 0.3*exp((-2.535e-7)*vm)/(1.0+exp(-0.1*(vm+32.0)));
}
else{
ah = 0.135*exp(-1.0*(vm+80.0)/6.8);
bh = 3.56*exp(0.079*vm)+(3.1e5)*exp(0.35*vm);
aj = (-127140.0*exp(0.2444*vm)-(3.474e-5)*exp(-0.04391*vm))*(vm+37.78)/
	  (1.0+exp(0.311*(vm+79.23)));
bj = 0.1212*exp(-0.01052*vm)/(1.0+exp(-0.1378*(vm+40.14)));
}

tauh = 1.0/(ah+bh);
infh = ah*tauh;

tauj = 1.0/(aj+bj);
infj = aj*tauj;

asa = 0.65/((exp(-1.0*(vm+10.0)/8.5))+(exp(-1.0*(vm-30.0)/59.0)));
bsa = 0.65/(2.5 + exp((vm+82.0)/17.0) );
infsa = 1.0/(1.0+exp(-1.0*(vm+20.47)/17.54));
tausa = (1.0/(asa+bsa))/Kq10;

asi = 1.0/(18.53+exp((vm+113.7)/10.95));
bsi = 1.0/(35.56+exp(-1.0*(vm+1.26)/7.44));
infsi = 1.0/(1.0+exp((vm+43.1)/5.3));
tausi = (1.0/(asi+bsi))/Kq10;

aua = 0.65*1.0/(exp(-1.0*(vm+10.0)/8.5)+exp(-1.0*(vm-30.0)/59.0));
bua = 0.65*1.0/(2.5+exp((vm+82.0)/17.0));
infua = 1.0/(1.0+exp(-1.0*(vm+30.3)/9.6));
tauua = (1.0/(aua+bua))/Kq10;

aui = 1.0/(21.0+exp(-1.0*(vm-185.0)/28.0));
bui = exp((vm-158.0)/16.0);
infui = 1.0/(1.0+exp((vm-99.45)/27.48));
tauui = (1.0/(aui+bui))/Kq10;

axr = 0.0003*(vm+14.1)/(1.0-exp(-1.0*(vm+14.1)/5.0));
bxr = (7.3898e-5)*(vm-3.3328)/(exp((vm-3.3328)/5.1237)-1.0);
infxr = 1.0/(1.0+exp(-1.0*(vm+14.1)/6.5));
tauxr = (1.0/(axr+bxr));

axs = (4.0e-5)*(vm-19.9)/(1.0-exp(-1.0*(vm-19.9)/17.0));
bxs = (3.5e-5)*(vm-19.9)/(exp((vm-19.9)/9.0) - 1.0);
tauxs = 0.5/(axs+bxs);
infxs = 1.0/sqrt(1.0+exp(-1.0*(vm-19.9)/12.7));

taud = (1.0-exp(-1.0*(vm+10.0)/6.24))/(0.035*(vm+10.0)*(1.0+exp(-1.0*(vm+10.0)/6.24)));
infd = 1.0/(1.0+exp(-1.0*(vm+10.0)/8.0));

tauf = 9.0/(0.0197*exp(-1.0*0.0337*0.0337*(vm+10.0)*(vm+10.0)) + 0.02);
inff = 1.0/(1.0+exp((vm+28.0)/6.9));

inffca = 1.0/(1.0+Cai/0.00035);

Fn = 1.0e-12*Vrel*Jrel-(5.0e-7/F)*(0.5*Ical-0.2*Inaca);

infu = 1.0/(1.0+exp(-1.0*(Fn-3.4175e-13)/13.67e-16));

tauv = 1.91+2.09/(1.0+exp(-1.0*(Fn-3.4175e-13)/13.67e-16));
infv = 1.0-1.0/(1.0+exp(-1.0*(Fn-6.835e-14)/13.67e-16));

tauw = 6.0*(1.0-exp(-1.0*(vm-7.9)/5.0))/((1.0+0.3*exp(-1.0*(vm-7.9)/5.0))*(vm-7.9));
infw = 1.0-1.0/(1.0+exp(-1.0*(vm-40.0)/17.0));

/* update differential equations */
fv[i] += -1.0*Iion;

fp->Nai = dNai;
fp->Ki = dKi;
fp->Cai = dCai;
fp->CaNSR = dCaNSR;
fp->CaJSR = dCaJSR;

fp->m = (infm-m)/taum;
fp->h = (infh-h)/tauh;
fp->j = (infj-j)/tauj;
fp->sa = (infsa-sa)/tausa;
fp->si = (infsi-si)/tausi;
fp->ua = (infua-ua)/tauua;
fp->ui = (infui-ui)/tauui;
fp->xr = (infxr-xr)/tauxr;
fp->xs = (infxs-xs)/tauxs;
fp->d = (infd-d)/taud;
fp->f = (inff-f)/tauf;
fp->fca = (inffca-fca)/taufca;
fp->u = (infu-u)/tauu;
fp->v = (infv-v)/tauv;
fp->w = (infw-w)/tauw;

 if( UseAuxvars ) {
	ap->Inai = Ina;
	ap->Ik1i = Ik1;
	ap->Itoi = Ito;
	ap->Ikuri = Ikur;
	ap->Ikri = Ikr;
	ap->Iksi = Iks;
	ap->Icali = Ical;
	ap->Inaki = Inak;
	ap->Inacai = Inaca;
	ap->Ibcai = Ibca;
	ap->Ibnai = Ibna;
	ap->Ipcai = Ipca;
	ap->Jreli = Jrel;
	ap->Jtri = Jtr;
	ap->Jupi = Jup;
	ap->Jupleaki = Jupleak;
	ap->Iioni = Iion ;
	ap->Caii = Cai;
	ap->CaNSRi = CaNSR;
	ap->CaSSi = CaSS;
	ap->CaJSRi = CaJSR;
 }

    } /* end-if */

    /* go to next patch */
    qm += PatchSize[ nt[i] ];
    fq += PatchSize[ nt[i] ];
    av += AuxiliarySize[ nt[i] ];

  } /* end-for */

  DebugLeave( "GetF_CRN" );

  return( 0 );
}

/*	Grid123 Module					*/
/*							*/
/*	used for 1,2,3-D domains, especially for	*/
/*	implicit codes which need storage of the main	*/
/*	diagonal					*/
/*							*/
/*		sigma is in mS/cm			*/
/*		dx & dy are in cm			*/
/*		dz or thickness is in cm		*/
/*		volume is in cm^3 (if used)		*/
/*							*/
/*	uses the sparse banded matrix format		*/

/* Copyright John B. Pormann, 21 July 2000, all rights reserved */



#include "CardioWave.h"

void MakeMatrix( sparse M, domain_t dt, int lsz, int Nd,
        int xd, int yd, int zd,
        real sx, real sy, real sz, real th, real ph,
        real dx, real dy, real dz );

static real __volume = 1.0;
static real __invvolume = 1.0;
static char* RCSID = "$Id: Domain123.c 14 2007-05-11 14:57:55Z jbp $";

static rword resources[] = {
	{ "bathdx",	4111 },
	{ "bathdy",	4112 },
	{ "bathdz",	4113 },
	{ "bathxd",	4101 },
	{ "bathxdim",	4101 },
	{ "bathyd",	4102 },
	{ "bathydim",	4102 },
	{ "bathzd",	4103 },
	{ "bathzdim",	4103 },
	{ "bdx",	4111 },
	{ "bdy",	4112 },
	{ "bdz",	4113 },
	{ "defaultnodetype",	4201 },
	{ "default_nodetype",	4201 },
	{ "nodefile",	4200 },
	{ "dx",		4001 },
	{ "dy",		4002 },
	{ "dz",		4003 },
	{ "scale",	4040 },
	{ "scalefactor",4040 },
	{ "scaling",	4040 },
	{ "sig",	4999 },
	{ "sigma",	4999 },
	{ "sigma_0",	4030 },
	{ "sigma_1",	4031 },
	{ "sigma_2",	4032 },
	{ "sigma_x",	4030 },
	{ "sigma_y",	4031 },
	{ "sigma_z",	4032 },
	{ "sigma_th",	4036 },
	{ "sigma_theta",4036 },
	{ "sigma_angle",4036 },
	{ "sigma_ph",	4037 },
	{ "sigma_phi",	4037 },
	{ "sigma_int_0",4030 },
	{ "sigma_int_1",4031 },
	{ "sigma_int_2",4032 },
	{ "sigma_int_x",4030 },
	{ "sigma_int_y",4031 },
	{ "sigma_int_z",4032 },
	{ "sigma_int_th",4036 },
	{ "sigma_int_theta",4036 },
	{ "sigma_int_angle",4036 },
	{ "sigma_int_ph",4037 },
	{ "sigma_int_phi",4037 },
	{ "sigma_ext_0",4033 },
	{ "sigma_ext_1",4034 },
	{ "sigma_ext_2",4035 },
	{ "sigma_ext_x",4033 },
	{ "sigma_ext_y",4034 },
	{ "sigma_ext_z",4035 },
	{ "sigma_ext_th",4036 },
	{ "sigma_ext_theta",4038 },
	{ "sigma_ext_angle",4038 },
	{ "sigma_ext_ph",4039 },
	{ "sigma_ext_phi",4039 },
	{ "sigma_b",	4131 },
	{ "sigma_bath",	4131 },
	{ "vol",	4051 },
	{ "volume",	4051 },
	{ "xd",		4004 },
	{ "xdim",	4004 },
	{ "yd",		4005 },
	{ "ydim",	4005 },
	{ "zd",		4006 },
	{ "zdim",	4006 },
	{ NULL, 0 }
};

int InitDomain( char** res ) {
	int  i,j;
	int  cmd;
	real sex = 1.0, sey = 1.0, sez = 1.0;
	real six = 1.0, siy = 1.0, siz = 1.0, sss;
	real theta_int = 0.0, phi_int = 0.0;
	real theta_ext = 0.0, phi_ext = 0.0;
	real scale = 1.0;
	real dx = 1.0, dy = 1.0, dz = 1.0;
	real volume = 1.0;
	int  Ndim = 0, Xdim = 1, Ydim = 1, Zdim = 1;
	int  Flags;
	int  BathXdim = 1, BathYdim = 1, BathZdim = 1;
	int  BathBW;
	real sb = 1.0, sbx, sby, sbz, sb2;
	real bdx = 1.0, bdy = 1.0, bdz = 1.0;
	byte defaultnodetype = ByteError;
	char* nodefile = NULL;

	DebugEnter( "InitDomain_123" );

	i = 0;
	while( res[i] != NULL ) {
		cmd = FindCommand( resources, res[i] );
		if( cmd == 4999 ) {
			/* assume a resource like:	*/
			/*   sigma[e][y]=1.0		*/
			/*   sigma(i)=1.0		*/
			/*   sigma:bath=1.0		*/
			j = FindNum( res[i] );
			switch( j ) {
				case 'E':
				case 'e':
					cmd = 4033;
					break;
				case 'I':
				case 'i':
					cmd = 4030;
					break;
				case 'B':
				case 'b':
					cmd = 4111;
					break;
			}
			/* FindNum2() will return -1 if no second paren is found */
			j = FindNum2( res[i] );
			if( j == 'Y' ) {
				/* they specified the y-dir conductivity */
				cmd++;
			}
			if( j == 'A' ) {
				/* they specified the angle for the fiber */
				cmd += 2;
			}
		}
		switch( cmd ) {
			case 4001:
				dx = GetRealValue( res[i] );
				break;
			case 4002:
				dy = GetRealValue( res[i] );
				break;
			case 4003:
				dz = GetRealValue( res[i] );
				break;
			case 4004:
				Xdim = GetIntValue( res[i] );
				break;
			case 4005:
				Ydim = GetIntValue( res[i] );
				break;
			case 4006:
				Zdim = GetIntValue( res[i] );
				break;
			case 4030:
				six = GetRealValue( res[i] );
				break;
			case 4031:
				siy = GetRealValue( res[i] );
				break;
			case 4032:
				siz = GetRealValue( res[i] );
				break;
			case 4033:
				sex = GetRealValue( res[i] );
				break;
			case 4034:
				sey = GetRealValue( res[i] );
				break;
			case 4035:
				sez = GetRealValue( res[i] );
				break;
			case 4036:
				theta_int = GetRealValue( res[i] );
				break;
			case 4037:
				phi_int = GetRealValue( res[i] );
				break;
			case 4038:
				theta_ext = GetRealValue( res[i] );
				break;
			case 4039:
				phi_ext = GetRealValue( res[i] );
				break;
			case 4040:
				scale = GetRealValue( res[i] );
				break;
			case 4051:
				volume = GetRealValue( res[i] );
				break;
			/* bath properties */
			case 4101:
				BathXdim = GetIntValue( res[i] );
				break;
			case 4102:
				BathYdim = GetIntValue( res[i] );
				break;
			case 4103:
				BathZdim = GetIntValue( res[i] );
				break;
			case 4111:
				bdx = GetRealValue( res[i] );
				break;
			case 4112:
				bdy = GetRealValue( res[i] );
				break;
			case 4113:
				bdz = GetRealValue( res[i] );
				break;
			case 4131:
				sb = GetRealValue( res[i] );
				break;
			case 4200:
				nodefile = GetStringValue( res[i] );
				break;
			case 4201:
				defaultnodetype = GetByteValue( res[i] );
				break;
		}
		i++;
	}

	/* set up the system size */
	Ndim  = (Xdim!=1)+(Ydim!=1)+(Zdim!=1);
	CreateRegularSplitting( Tissue, Xdim, Ydim, Zdim );

	/* now create the matrices */
	i = spralloc( &MatrixINT, Intra, _STENCIL, TissueLocalSize, 
		TissueLocalSize, 1024 );
	if( i < 0 ) {
		return( -1 );
	}
	/* create the matrix */
	MatrixINT.jcoef[0] = Ndim;
	switch( Ndim ) {
		case 3: MatrixINT.jcoef[3] = Zdim;
		case 2: MatrixINT.jcoef[2] = Ydim;
		case 1: MatrixINT.jcoef[1] = Xdim;
	}
	MakeMatrix( MatrixINT, Intra, TissueLocalSize, Ndim, Xdim, Ydim, Zdim,
		scale*six, scale*siy, scale*siz, theta_int, phi_int,
		dx, dy, dz );
	DebugMark( "MakeMatrix-X" );

	if( SimType != Monodomain ) {
		i = spralloc( &MatrixEXT, Extra, _STENCIL, 
			TissueLocalSize, TissueLocalSize, 1024 );
		if( i < 0 ) {
			return( -3 );
		}
		MatrixEXT.jcoef[0] = Ndim;
		switch( Ndim ) {
			case 3: MatrixEXT.jcoef[3] = Zdim;
			case 2: MatrixEXT.jcoef[2] = Ydim;
			case 1: MatrixEXT.jcoef[1] = Xdim;
		}
		MakeMatrix( MatrixEXT, Extra, TissueLocalSize, Ndim, Xdim,
			Ydim, Zdim, scale*sex, scale*sey, scale*sez,
			theta_ext, phi_ext, dx, dy, dz );
	}

	if( (SimType==ReducedBidomainBath) or (SimType==FullBidomainBath) ) {
		if( SelfPE == 0 ) {
			printf("Stencil matrix cannot be used with Bath!\n");
		}
		return( -100 );
	}

	/* fudge the volume vector */
	/* : we usually set volume=1 */
	if( volume < 0.0 ) {
		volume = 1.0;
		switch( Ndim ) {
			case 3:
				volume *= dz;
			case 2:
				volume *= dy;
			case 1:
				volume *= dx;
		}
	}
	/* fudge the volume vector */
	__volume = volume;
	__invvolume = 1.0/volume;
	Volume.size  = 1;
	Volume.dtype = Tissue;
	Volume.data  = &__volume;
	InvVolume.size  = 1;
	InvVolume.dtype = Tissue;
	InvVolume.data  = &__invvolume;

	/* now load in the node info */
	if( nodefile != NULL ) {
		i = bvecreadinfo( nodefile, &j );
		if( i < 0 ) {
			return( -15 );
		}
		if( j != TissueGlobalSize ) {
			return( -16 );
		}
	}
	bvecalloc( &NodeType, Tissue );
	if( nodefile != NULL ) {
		bvecread( nodefile, NodeType );
	} else if( defaultnodetype == ByteError ) {
		bvfill( 0, NodeType );
	} else {
		bvfill( defaultnodetype, NodeType );
	}

	if( (DebugLevel>0) and (SelfPE==0) ) {
		printf("Domain123: %iD: %i x %i x %i = %i\n",
			Ndim, Xdim, Ydim, Zdim, TissueGlobalSize );
	}
	if( ShowVersion and (SelfPE==0) ) {
		printf("Domain123: RCSID: %s\n",RCSID);
	}

	DebugLeave( "InitDomain_123" );

	return( 0 );
}

void ExitDomain( void ) {
	return;
}

void MakeMatrix( sparse MM, domain_t tp, int localsz, int Ndim, int Xdim, 
		int Ydim, int Zdim, real sx, real sy, real sz,
		real th, real ph, real dx, real dy, real dz ) {
	int  m,g,f,x,y,z;
	real ssdd,ssd;
	real mx,my,mz;
	real RX,RY,RZ,RXY,RXZ,RYZ,sn,cs,sn2,cs2;
	int  XYdim = Xdim*Ydim;
	int  Xd1   = Xdim-1;
	int  Yd1   = Ydim-1;
	int  Zd1   = Zdim-1;

	/* compute the main diagonal term */
	ssd = 0.0;
	switch( Ndim ) {
		case 3: ssd = 0.0;
			cs  = cos( th );
			sn  = sin( th );
			cs2 = cos( ph );
			sn2 = sin( ph );
			RX  = ( sx*cs*cs + sy*sn*sn )/dx/dx;
			RY  = ( sx*cs2*cs2*sn*sn + sy*cs2*cs2*cs*cs
				+ sz*sn2*sn2 )/dy/dy;
			RZ  = ( sx*sn2*sn2*sn*sn + sy*sn2*sn2*cs*cs
				+ sz*cs2*cs2 )/dz/dz;
			RXY = ( (sy - sx)*cs2*cs*sn )/2.0/dx/dy;
			RXZ = ( (sx - sy)*sn2*cs*sn )/2.0/dx/dz;
			RYZ = ( -1.0*sn2*cs2*( sx*sn*sn - sy*cs*cs )
				+ sz*sn2*cs2 )/2.0/dy/dz;
			m = XYdim + Xdim + 1;
			break;
		case 2: ssd = 0.0;
			cs = cos( th );
			sn = sin( th );
			RX  = ( sx*cs*cs + sy*sn*sn )/dx/dx;
			RXY = cs*sn*( sy - sx )/2.0/dx/dy;
			RY  = ( sx*sn*sn + sy*cs*cs )/dy/dy;
			m = Xdim + 1;
			break;
		case 1: ssd = 0.0;
			RX = sx/dx/dx;
			m = 1;
			break;
		case 0: ssd = 0.0;
			m = 0;
			break;
	}

	/* m is the local node number */
	/* : for Stencil matrices, we fudged it above */
	/* decompose node number to (x,y,z) */
	z = 1;
	y = 1;
	x = 1;

		/* in case of 0D domain, we'll set ssdd to 0.0 */
		ssdd = ssd;

		/* first dimension */
		if( Ndim > 0 ) {
			if( x < Xd1 ) {
				sprput(MM,m,m+1,RX);
				ssdd += RX;
			}
			if( x > 0 ) {
				sprput(MM,m,m-1,RX);
				ssdd += RX;
			}
		}

		/* second dimension */
		if( Ndim > 1 ) {
			if( y < Yd1 ) {
				sprput(MM,m,m+Xdim,RY);
				ssdd += RY;
				if( x < Xd1 ) {
					sprput(MM,m,m+Xdim+1,RXY);
					ssdd += RXY;
				}
				if( x > 0 ) {
					sprput(MM,m,m+Xdim-1,-RXY);
					ssdd -= RXY;
				}
			}
			if( y > 0 ) {
				sprput(MM,m,m-Xdim,RY);
				ssdd += RY;
				if( x < Xd1 ) {
					sprput(MM,m,m-Xdim+1,-RXY);
					ssdd -= RXY;
				}
				if( x > 0 ) {
					sprput(MM,m,m-Xdim-1,RXY);
					ssdd += RXY;
				}
			}
		}

		/* third dimension */
		if( Ndim > 2 ) {
			if( z < Zd1 ) {
				sprput(MM,m,m+XYdim,RZ);
				ssdd += RZ;
				if( y < Yd1 ) {
					sprput(MM,m,m+XYdim+Xdim,RYZ);
					ssdd += RYZ;
				}
				if( y > 0 ) {
					sprput(MM,m,m+XYdim-Xdim,-RYZ);
					ssdd -= RYZ;
				}
				if( x < Xd1 ) {
					sprput(MM,m,m+XYdim+1,RXZ);
					ssdd += RXZ;
				}
				if( x > 0 ) {
					sprput(MM,m,m+XYdim-1,-RXZ);
					ssdd -= RXZ;
				}
			}
			if( z > 0 ) {
				sprput(MM,m,m-XYdim,RZ);
				ssdd += RZ;
				if( y < Yd1 ) {
					sprput(MM,m,m-XYdim+Xdim,-RYZ);
					ssdd -= RYZ;
				}
				if( y > 0 ) {
					sprput(MM,m,m-XYdim-Xdim,RYZ);
					ssdd += RYZ;
				}
				if( x < Xd1 ) {
					sprput(MM,m,m-XYdim+1,-RXZ);
					ssdd -= RXZ;
				}
				if( x > 0 ) {
					sprput(MM,m,m-XYdim-1,RXZ);
					ssdd += RXZ;
				}
			}
		}

		/* now correct the main diagonal term */
		sprput(MM,m,m,-ssdd);

	return;
}



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fdm2d.h"

void Ndiv2Nx(int *Ndiv,int ityp, int *Nx){
	switch(ityp)
	{
		case 0:	// v1-grid
			Nx[0]=Ndiv[0]+1;
			Nx[1]=Ndiv[1];
			break;
		case 1: // v2-grid
			Nx[0]=Ndiv[0];
			Nx[1]=Ndiv[1]+1;
			break;
		case 2: // s11-grid
			Nx[0]=Ndiv[0];
			Nx[1]=Ndiv[1];
			break;
		case 3: // s12-grid
			Nx[0]=Ndiv[0]+1;
			Nx[1]=Ndiv[1]+1;
			break;
		case 4: // s22-grid
			Nx[0]=Ndiv[0];
			Nx[1]=Ndiv[1];
			break;
	}	
}

//   Index (i,j) <--> Index l 
int ij2l(int i, int j, int *Nx){
	return i*Nx[1]+j;		
}
int ij2l(int i, int j, int ityp, int *Ndiv){
	int Nx[2];
	Ndiv2Nx(Ndiv,ityp,Nx);
	return ij2l(i,j,Nx);
}

void l2ij(int l, int *Nx, int *i, int *j){
	*i=l/Nx[1];
	*j=l%Nx[1];
}

void l2ij(int l, int ityp, int *Ndiv, int *i, int *j){
	int Nx[2];
	Ndiv2Nx(Ndiv,ityp,Nx);
	l2ij(l,Nx,i,j);
}



//   Index (i,j) --> Coordinate :(x,y) 
void indx2cod(int i, int j, int ityp, double *Xa, double *dx, double *xcod)
{
	int i0[2];

	indx_ofst(ityp,i0);
	xcod[0]=Xa[0]+dx[0]*(i+i0[0]*0.5);
	xcod[1]=Xa[1]+dx[1]*(j+i0[1]*0.5);
}

void indx_ofst(int ityp,int *i0){

	i0[0]=0;
	i0[1]=0;
	switch(ityp){
	case 0:	// v1-grid
		i0[1]=1;
		break;
	case 1: // v2-grid
		i0[0]=1;
		break;
	case 2: // s11-grid
		i0[0]=1;
		i0[1]=1;
		break;
	case 3: // s12-grid
		break;
	case 4: // s22-grid
		i0[0]=1;
		i0[1]=1;
		break;
	default:
		puts("invalid grid type!");
		puts("--> process terminated");
		exit(-1);
	};
		
}
//   Coordinate :(x,y) --> Index (i,j)
void cod2indx(double *xcod, int ityp, double *Xa, double *dx, int *indx){
	int i0[2];
	indx_ofst(ityp,i0);

//	indx[0]=(int)( (xcod[0]-Xa[0])/dx[0]-i0[0]*0.5 );
//	indx[1]=(int)( (xcod[1]-Xa[1])/dx[1]-i0[1]*0.5 );
	indx[0]=(floor)( (xcod[0]-Xa[0])/dx[0]+(1-i0[0])*0.5 );
	indx[1]=(floor)( (xcod[1]-Xa[1])/dx[1]+(1-i0[1])*0.5 );
}

// Distance
double dist2d(double *x1, double *x2){
	double xx=x1[0]-x2[0];
	double yy=x1[1]-x2[1];

	return sqrt(xx*xx+yy*yy);
}


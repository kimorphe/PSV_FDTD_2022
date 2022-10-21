#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fdm2d.h"

//using namespace std;

double plane_wave(
	int ig, int jg,	// 2D grid index 	
	int igrd,	// grid type (0: v1,1:v2, 2:s11,3:s12,4:s22)
	int it,		// time step
	int ityp,	// wave type (L-wave=0, T-wave=1)
	double ain,	// incident angle [deg]
	double *xp0,	// zero phase point
	InWv inwv,	// wave data(class)
	Dom2D dom,	// domain data
	double *ginx, 	// incident field || (para)
	double *giny	// incident field _|_ (vert)
){

	double pin[2],din[2],xcod[2];
	double PI=4.0*atan(1.0);
	double fin,fin1,fin2,cin;
	double t0,dt=inwv.dt;
	double arg,al,am,rm,rc,pd;
	double vin[2],sij[2][2];
	int Nt=inwv.Nt; 
	int i,j;
	int ismp1,ismp2;

	t0=it*dt;
	if( igrd >1) t0=(it-0.5)*dt;

	ain*=(PI/180.0);
	pin[0]=cos(ain);
	pin[1]=sin(ain);

	switch(ityp){
	case 0:	// L-wave 
		din[0]=pin[0];
		din[1]=pin[1];	
		cin=dom.cL;
		break;
	case 1:	// T-wave 
		din[0]=-pin[1];
		din[1]=pin[0];	
		cin=dom.cT;
		break;
	default: 
		puts("Illegal wave type is specified in plane_wave . ");
		puts("--> process terminated");
		exit(-1);
		break;
	};
	pd=pin[0]*din[0]+pin[1]*din[1];

	am=dom.rho*dom.cT*dom.cT;
	al=dom.rho*(dom.cL*dom.cL-2.0*dom.cT*dom.cT);
	rm=al+2.*am;
	rc=dom.rho*cin;
	

	indx2cod(ig,jg,igrd,dom.Xa,dom.dx,xcod);
	xcod[0]-=xp0[0];
	xcod[1]-=xp0[1];
	arg=(t0-(xcod[0]*pin[0]+xcod[1]*pin[1])/cin);
	ismp1=floor(arg/dt);
	ismp2=ismp1++;

	fin1=0.0; fin2=0.0;

	if(ismp1>=0 && ismp1 < Nt) fin1=inwv.amp[ismp1];
	if(ismp2>=0 && ismp2 < Nt) fin2=inwv.amp[ismp2];

	fin=fin1*(ismp2-arg/dt)+fin2*(arg/dt-ismp1);
	vin[0]=din[0]*fin;
	vin[1]=din[1]*fin;	


	for(i=0;i<2;i++){	
	for(j=0;j<2;j++){	
		sij[i][j]=am*(pin[i]*din[j]+pin[j]*din[i]);	
		if(i==j) sij[i][j]+=(al*pd);
		sij[i][j]*=(-fin/cin);
	}
	}

	switch(igrd){
	case 0:	// v1-grid
		(*ginx)=-sij[0][0]*pin[0]/rc;
		(*giny)=-sij[0][1]*pin[1]/rc;
		break;
	case 1:	// v2-grid
		(*ginx)=-sij[0][1]*pin[0]/rc;
		(*giny)=-sij[1][1]*pin[1]/rc;
		break;
	case 2:	// s11-grid
		(*ginx)=-vin[0]*rm*pin[0]/cin;
		(*giny)=-vin[1]*al*pin[0]/cin;
		break;
	case 3:	// s12-grid
		(*ginx)=-vin[1]*am*pin[0]/cin;
		(*giny)=-vin[0]*am*pin[1]/cin;
		break;
	case 4:	// s22-grid
		(*ginx)=-vin[0]*al*pin[0]/cin;
		(*giny)=-vin[1]*rm*pin[1]/cin;
		break;
	};
	return((*ginx)+(*giny));
}

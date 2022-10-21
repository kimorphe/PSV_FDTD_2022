#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdm2d.h"


//-----------------------------------------------------
//
//	INITIALIZE PML
//
void PML :: setup(Dom2D dom, char *ficon){

	int i,j,ndim=2;
//,i0,j0;
	double W2,A0; 
	double ss,xx,yy;
//	FILE *fp;
//	char fname[]="icon.inp";

	for(i=0;i<ndim;i++){
		ia[i]=(floor)((Xa[i]-dom.Xa[i])/dom.dx[i]+0.5);
		ib[i]=(floor)((Xb[i]-dom.Xa[i])/dom.dx[i]+0.5);

		if(ia[i] > dom.Ndiv[i]) ia[i]=dom.Ndiv[i];
		if(ib[i] > dom.Ndiv[i]) ib[i]=dom.Ndiv[i];

		if(ia[i]<0) ia[i]=0;
		if(ib[i]<0) ib[i]=0;

		nx[i]=abs(ib[i]-ia[i])+1;	// number of grids

		Xa[i]=dom.Xa[i]+ia[i]*dom.dx[i];
		Xb[i]=dom.Xa[i]+ib[i]*dom.dx[i];
		dx[i]=0.0;
		if(nx[i]>1) dx[i]=(Xb[i]-Xa[i])/(nx[i]-1); //dx can be negative


	}
		printf("ia=%d %d\n",ia[0],ia[1]);
		printf("ib=%d %d\n",ib[0],ib[1]);


	mem_alloc();
	set_IC(dom,ficon);

	mexp=2;
	gmm=1.e-05;

	for(i=0;i<nx[0]-1;i++){
		xx=Xa[0]+dx[0]*(i+0.5);
		W2=0.0; A0=0.0;

		ss=dom.xa[0]-xx;
		if(ss > 0.0){
			W2=2.0*fabs(dom.Wa[0]);
			A0=-(mexp+1)*dom.cL*logf(gmm)/pow(W2,mexp+1);	
			dmpx[i]=A0*pow(ss,mexp);
		}
		ss=xx-dom.xb[0];
		if(ss > 0.0){
			W2=2.0*fabs(dom.Wb[0]);
			A0=-(mexp+1)*dom.cL*logf(gmm)/pow(W2,mexp+1);	
			dmpx[i]=A0*pow(ss,mexp);
		}
	}

	
	for(j=0;j<nx[1]-1;j++){
		yy=Xa[1]+dx[1]*(j+0.5);
		W2=0.0; A0=0.0;

		ss=dom.xa[1]-yy;
		if(ss > 0.0){
			W2=2.0*fabs(dom.Wa[1]);
			A0=-(mexp+1)*dom.cL*logf(gmm)/pow(W2,mexp+1);	
			dmpy[j]=A0*pow(ss,mexp);
		}
		ss=yy-dom.xb[1];
		if(ss > 0.0){
			W2=2.0*fabs(dom.Wb[1]);
			A0=-(mexp+1)*dom.cL*logf(gmm)/pow(W2,mexp+1);	
			dmpy[j]=A0*pow(ss,mexp);
		}
	}		

}

//-----------------------------------------------------
//		DYNAMIC MEMORY ALLOCATION
//
void PML :: mem_alloc(){
	
	int j,n1,n2;
//	double *pt;
	
	n1=nx[0]-1;
	n2=nx[1]-1;
	mem_alloc2D(n1,n2,&s11p);
	mem_alloc2D(n1,n2,&s11v);
	mem_alloc2D(n1,n2,&s22p);
	mem_alloc2D(n1,n2,&s22v);
	mem_alloc2D(n1,n2,&s11);
	mem_alloc2D(n1,n2,&s22);

	dmpx=(double *)malloc(sizeof(double)*(n1));
	dmpy=(double *)malloc(sizeof(double)*(n2));

	for(j=0;j<n1;j++) dmpx[j]=0.0;
	for(j=0;j<n2;j++) dmpy[j]=0.0;

	n1=nx[0];
	n2=nx[1]-1;
	mem_alloc2D(n1,n2,&v1p);
	mem_alloc2D(n1,n2,&v1v);
	mem_alloc2D(n1,n2,&v1);

	n1=nx[0]-1;
	n2=nx[1];
	mem_alloc2D(n1,n2,&v2p);
	mem_alloc2D(n1,n2,&v2v);
	mem_alloc2D(n1,n2,&v2);

	n1=nx[0];
	n2=nx[1];
	mem_alloc2D(n1,n2,&s12p);
	mem_alloc2D(n1,n2,&s12v);
	mem_alloc2D(n1,n2,&s12);

};
int PML::set_IC(Dom2D dom, char *ficon){

	FILE *fp=fopen(ficon,"r");
	char cbff[128],fname[128];
	int i,j,l,i0,j0,it0=0;
	int init,inwv_typ,inwv_num;
	double ain,xp0[2];
//	double gin,ginx,giny;
	double ginx,giny;
	
	if(fp == NULL){
		s110=0.0;
		s220=0.0;
		puts("Can't open icon.inp.");
		puts(" --> initial field in PML set to zero.");
		return(-1);
	}	

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&init);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&ain);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&inwv_typ);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf \n",xp0,xp0+1);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&inwv_num);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",&s110,&s220);

	fclose(fp);
	
	sprintf(fname,"inwv%d.dat",inwv_num);
	InWv inwv=InWv(fname);

	for(l=0;l<5;l++){
		gridNum(l);
//		printf("l=%d\n",l);
//		printf("Nx=%d %d\n",Nx[0],Nx[1]);
	for(i=0;i<Nx[0];i++){
		i0=ia[0];
	for(j=0;j<Nx[1];j++){
		j0=ia[1];
		ginx=0.0; giny=0.0;
		if(init ==1) plane_wave(i+i0,j+j0,l,it0,inwv_typ,ain,xp0,inwv,dom,&ginx,&giny);
		switch(l){
		case 0:
			v1p[i][j]=ginx;
			v1v[i][j]=giny;
			v1[i][j]=v1p[i][j]+v1v[i][j];
			break;
		case 1: 
			v2p[i][j]=ginx;
			v2v[i][j]=giny;
			v2[i][j]=v2p[i][j]+v2v[i][j];
			break;
		case 2: 
			s11p[i][j]=ginx+s110;
			s11v[i][j]=giny;
			s11[i][j]=s11p[i][j]+s11v[i][j];
			break;
		case 3:
			s12p[i][j]=ginx;
			s12v[i][j]=giny;
			s12[i][j]=s12p[i][j]+s12v[i][j];
			break;
		case 4:
			s22p[i][j]=ginx+s220;
			s22v[i][j]=giny;
			s22[i][j]=s22p[i][j]+s22v[i][j];
			break;
		}
	}
	}
	}
	return(0);
}
void PML ::gridNum(int ityp){
	switch(ityp){
	case 0:	// v1-grid
		Nx[0]=nx[0];
		Nx[1]=nx[1]-1;
		break;
	case 1:	// v2-grid
		Nx[0]=nx[0]-1;
		Nx[1]=nx[1];
		break;
	case 2:	// s11-grid
		Nx[0]=nx[0]-1;
		Nx[1]=nx[1]-1;
		break;
	case 3:	// s12-grid
		Nx[0]=nx[0];
		Nx[1]=nx[1];
		break;
	case 4:	// s22-grid
		Nx[0]=nx[0]-1;
		Nx[1]=nx[1]-1;
		break;
	};
}
//-----------------------------------------------------
//
//		TIME-STEPPING OPERATION 1
//		  ( Velocity --> Stress)
//
void PML :: v2s(double amu, double almb, double dt){
	int i,j,ib,jb;
//	double dmpX,dmpY,alph,beta,dv;
	double dmpX,dmpY,alph,beta;
	double dv1,dv2;
	double dtr=1.0/dt;
	double a2m=almb+2.0*amu;

//		Normal Stress: s11, s22
	for(i=0; i<nx[0]-1 ;i++){
		dmpX=dmpx[i];
		alph=(dtr-dmpX)/(dtr+dmpX);
	for(j=0; j<nx[1]-1 ;j++){
		dmpY=dmpy[j];
		beta=(dtr-dmpY)/(dtr+dmpY);

		dv1=(v1[i+1][j]-v1[i][j])/(dtr+dmpX);
		dv2=(v2[i][j+1]-v2[i][j])/(dtr+dmpY);

		s11p[i][j]*=alph;
		s11v[i][j]*=beta;
		s11p[i][j]+=(dv1*a2m/dx[0]);
		s11p[i][j]+=(s110*2.*dmpX/(dtr+dmpX));	// <--------------
		s11v[i][j]+=(dv2*almb/dx[1]);
		s11[i][j]=s11p[i][j]+s11v[i][j];

		s22p[i][j]*=alph;
		s22v[i][j]*=beta;
		s22p[i][j]+=(dv1*almb/dx[0]);
		s22p[i][j]+=(s220*2.*dmpX/(dtr+dmpX)); // <--------------
		s22v[i][j]+=(dv2*a2m/dx[1]);
		s22[i][j]=s22p[i][j]+s22v[i][j];
	}
	}

//		Shear Stress: s12
	for(i=1; i<nx[0]-1 ;i++){
		dmpX=.5*(dmpx[i]+dmpx[i-1]);
		alph=(dtr-dmpX)/(dtr+dmpX);
	for(j=1; j<nx[1]-1 ;j++){
		dmpY=.5*(dmpy[j]+dmpy[j-1]);
		beta=(dtr-dmpY)/(dtr+dmpY);
		
		dv2=(v2[i][j]-v2[i-1][j])*amu/(dtr+dmpX);
		dv1=(v1[i][j]-v1[i][j-1])*amu/(dtr+dmpY);

		s12p[i][j]*=alph;
		s12v[i][j]*=beta;
		s12p[i][j]+=(dv2/dx[0]);
		s12v[i][j]+=(dv1/dx[1]);
		s12[i][j]=s12p[i][j]+s12v[i][j];
	}
	}

//		Stress Free B.C.	
	for(ib=0;ib<2;ib++){
		i=0;
		if(ib==1) i=nx[0]-1;
	for(j=0;j<nx[1];j++){
		s12p[i][j]=0.0;
		s12v[i][j]=0.0;
		s12[i][j]=0.0;
	}
	}

	for(jb=0;jb<2;jb++){
		j=0;
		if(jb==1) j=nx[1]-1;
	for(i=1;i<nx[0]-1;i++){
		s12p[i][j]=0.0;
		s12v[i][j]=0.0;
		s12[i][j]=0.0;
	}
	}
};

//-----------------------------------------------------
//
//		TIME-STEPPING OPERATION 2
//		  ( Stress --> Velocity )
//
void PML :: s2v(double rho, double dt){
	int i,j,ib,jb;
	double dmpX,dmpY,alph,beta,dS11,dS22,dS12;
	double dtr=1.0/dt,sgn;
	double rdx=dx[0]*rho, rdy=dx[1]*rho;

//	################ Velocity v1 ##################
	for(i=1; i<nx[0]-1 ;i++){
		dmpX=.5*(dmpx[i]+dmpx[i-1]);
		alph=(dtr-dmpX)/(dtr+dmpX);
	for(j=0; j<nx[1]-1; j++){
		dmpY=dmpy[j];
		beta=(dtr-dmpY)/(dtr+dmpY);
		dS11=(s11[i][j]-s11[i-1][j])/((dtr+dmpX)*rdx);
		dS12=(s12[i][j+1]-s12[i][j])/((dtr+dmpY)*rdy);

		v1p[i][j]*=alph;
		v1v[i][j]*=beta;
		v1p[i][j]+=dS11;
		v1v[i][j]+=dS12;
		v1[i][j]=v1p[i][j]+v1v[i][j];

	}
	}
//		Boundary Nodes (stress-free B.C.)
	i=0;
	sgn=1.0;
	dmpX=dmpx[i];
	for(ib=0;ib<2;ib++){
		if(ib==1){
			i=nx[0]-1;
			dmpX=dmpx[i-1];
			sgn=-1.0;
		}
		alph=(dtr-dmpX)/(dtr+dmpX);
	for(j=0;j<nx[1]-1;j++){
		dmpY=dmpy[j];
		beta=(dtr-dmpY)/(dtr+dmpY);
		dS11=sgn*2.*(s11[i-ib][j]-s110)/((dtr+dmpX)*rdx);
		dS12=(s12[i][j+1]-s12[i][j])/((dtr+dmpY)*rdy);
		v1p[i][j]*=alph;
		v1v[i][j]*=beta;
		v1p[i][j]+=dS11;
		v1v[i][j]+=dS12;
		v1[i][j]=v1p[i][j]+v1v[i][j];
	}	
	}

//	################ Velocity v2 ##################
	for(i=0; i<nx[0]-1; i++){
		dmpX=dmpx[i];
		alph=(dtr-dmpX)/(dtr+dmpX);
	for(j=1; j<nx[1]-1; j++){
		dmpY=.5*(dmpy[j]+dmpy[j-1]);
		beta=(dtr-dmpY)/(dtr+dmpY);

		dS12=(s12[i+1][j]-s12[i][j])/((dtr+dmpX)*rdx);
		dS22=(s22[i][j]-s22[i][j-1])/((dtr+dmpY)*rdy);

		v2p[i][j]*=alph;
		v2v[i][j]*=beta;
		v2p[i][j]+=dS12;
		v2v[i][j]+=dS22;
		v2[i][j]=v2p[i][j]+v2v[i][j];
	}
	}

//		Boundary Nodes (stress-free B.C.)
	j=0;
	sgn=1.0;
	dmpY=dmpy[j];
	for(jb=0; jb<2; jb++){
		if(jb==1){
			j=nx[1]-1;
			sgn=-1.0;
			dmpY=dmpy[j-1];
		}	
		beta=(dtr-dmpY)/(dtr+dmpY);
	for(i=0; i<nx[0]-1; i++){
		dmpX=dmpx[i];
		alph=(dtr-dmpX)/(dtr+dmpX);

		dS12=(s12[i+1][j]-s12[i][j])/((dtr+dmpX)*rdx);
		dS22=2.0*sgn*(s22[i][j-jb]-s220)/((dtr+dmpY)*rdy);

		v2p[i][j]*=alph;
		v2v[i][j]*=beta;
		v2p[i][j]+=dS12;
		v2v[i][j]+=dS22;
		v2[i][j]=v2p[i][j]+v2v[i][j];
	}
	}
	
};

//----------------------------------------------------

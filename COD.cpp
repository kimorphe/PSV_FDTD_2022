#include <stdio.h>
#include <math.h>
#include "fdm2d.h"
#include "crack.h"
//using namespace std;

void Crack :: COD_hrz(int l, double dt, double rho, double *dx, double s22b_now){
	
	double dt2=dt*dt,dtx,dte;
	double a2p,a2m,a2;	// accelerations
	double u2p_tmp,u2m_tmp; // displacements
	double s22b_pre=s22b[l];// crack-opening stress
	double cod,cov,coa;	// crack opening disp., vel. and acc. 
	double xi,xi1,xi2;	// non-dimensional time
	double vxip,vxim,vxi; 	// velocities
	double EPS=1.e-10;	// small number used to avoide division by 0
	double dxr[2], ds12;

	dxr[0]=dx[0]*rho;
	dxr[1]=dx[1]*rho;
	ds12=(s12[l+1]-s12[l])/dxr[0];

	switch(iop[l]){
	case 0:		// OPEN PART
		a2p=( 2.0*s22p[l]/dxr[1]+ds12);
		a2m=(-2.0*s22m[l]/dxr[1]+ds12);
		u2p_tmp=u2p[l]+v2p[l]*dt+0.5*a2p*dt2;
		u2m_tmp=u2m[l]+v2m[l]*dt+0.5*a2m*dt2;

		cod=u2p_tmp-u2m_tmp;
		if(cod+cod0 >= 0.0){		// OPEN --> OPEN
			u2p[l]=u2p_tmp;	
			u2m[l]=u2m_tmp;
			v2p[l]+=(a2p*dt);
			v2m[l]+=(a2m*dt);
		}else{			// OPEN --> CLOSE
			iop[l]=1;

			cod=u2p[l]-u2m[l];
			cov=v2p[l]-v2m[l];
			coa=a2p-a2m;
			
			xi=sqrt(cov*cov-2.*cod*coa);
			xi1=(-cov+xi)/(coa*dt+EPS);
			xi2=(-cov-xi)/(coa*dt+EPS);

			xi=xi1;
			if(xi <= -EPS) xi=xi2;
			if(xi-1.0 >= EPS) xi=xi2;

			if(xi <0.0) puts("xi_c <0 !");
			if(xi >1.0) puts("xi_c >1 !");
			dtx=dt*xi;
			u2p[l]+=(v2p[l]*dtx+0.5*a2p*dtx*dtx);
			u2m[l]=u2p[l];

			vxip=v2p[l]+a2p*dtx;
			vxim=v2m[l]+a2m*dtx;
			vxi=0.5*(vxip+vxim);

			a2=(s22p[l]-s22m[l])/dxr[1]+ds12;
			
			dte=dt*(1.0-xi);
			u2p[l]+=(vxi*dte+0.5*a2*dte*dte);
			u2m[l]=u2p[l];
			v2p[l]=vxi+a2*dte;
			v2m[l]=v2p[l];
			
		}
		break;
	case 1:		// CLOSED PART
		a2=(s22p[l]-s22m[l])/dxr[1]+ds12;
		if(s22b_now <= 0.0){	// CLOSE --> CLOSE
			u2p[l]+=(v2p[l]*dt+0.5*a2*dt2);
			u2m[l]=u2p[l];
			v2p[l]+=(a2*dt);
			v2m[l]=v2p[l];
		}else{			// CLOSE --> OPEN
			iop[l]=0;
			xi=fabs(s22b_pre)/(fabs(s22b_now)+fabs(s22b_pre)+EPS)-0.5;
			dtx=dt*xi;
			dte=dt*(1.0-xi);

			if(xi <0.0) xi=0.0;
			if(xi >1.0) puts("xi_o >1 !");

			u2p[l]+=(v2p[l]*dtx+0.5*a2*dtx*dtx);
			u2m[l]=u2p[l];

			v2p[l]+=(a2*dtx);
			v2m[l]=v2p[l];
			
			a2p=( 2.0*s22p[l]/dxr[1]+ds12);
			a2m=(-2.0*s22m[l]/dxr[1]+ds12);
			if(a2p < a2m ) dte=0.0;
			u2p[l]+=(v2p[l]*dte+0.5*a2p*dte*dte);
			u2m[l]+=(v2m[l]*dte+0.5*a2m*dte*dte);
			if( u2p[l] < u2m[l]) printf("u2p=%lf,  u2m=%lf\n",u2p[l],u2m[l]);
			v2p[l]+=(a2p*dte);
			v2m[l]+=(a2m*dte);
		}
		break;
	};
}
void Crack :: COD_vrt(int l, double dt, double rho, double *dx, double s11b_now){
	
	double dt2=dt*dt,dtx,dte;
	double a1p,a1m,a1;	// accelerations
	double u1p_tmp,u1m_tmp; // displacements 
	double s11b_pre=s11b[l];// crack opening stress (at previous step) 
	double cod,cov,coa;	// crack opening disp., vel. and acc. 
	double xi,xi1,xi2;	// non-dimensional times
	double vxip,vxim,vxi;	// interim velocities
	double EPS=1.e-10;	// small number used to avoid division by 0 
	double dxr[2],ds12;

	dxr[0]=dx[0]*rho;
	dxr[1]=dx[1]*rho;
	ds12=(s12[l+1]-s12[l])/dxr[1];

	switch(iop[l]){
	case 0:		// OPEN PART
		a1p=( 2.0*s11p[l]/dxr[0]+ds12);
		a1m=(-2.0*s11m[l]/dxr[0]+ds12);

		u1p_tmp=u1p[l]+v1p[l]*dt+0.5*a1p*dt2;
		u1m_tmp=u1m[l]+v1m[l]*dt+0.5*a1m*dt2;

		cod=u1p_tmp-u1m_tmp;
		if(cod+cod0 >= 0.0){		// OPEN --> OPEN
			u1p[l]=u1p_tmp;	
			u1m[l]=u1m_tmp;
			v1p[l]+=(a1p*dt);
			v1m[l]+=(a1m*dt);
		}else{			// OPEN --> CLOSE
			iop[l]=1;

			cod=u1p[l]-u1m[l];
			cov=v1p[l]-v1m[l];
			coa=a1p-a1m;
			
			xi=sqrt(cov*cov-2.*cod*coa);
			xi1=(-cov+xi)/(coa*dt+EPS);
			xi2=(-cov-xi)/(coa*dt+EPS);

			xi=xi1;
			if(xi <= -EPS) xi=xi2;
			if(xi-1.0 >= EPS) xi=xi2;

			if(xi <0.0){
				puts("xi_c <0 !");
				printf("xi1,xi2=%lf %lf\n",xi1,xi2);
				printf("cod=%lf, cov=%lf, coa=%lf\n",cod,cov,coa);
				printf("l=%d\n",l);
				printf("u1p,u1m=%lf %lf\n",u1p[l],u1m[l]);
			}
			if(xi >1.0) puts("xi_c >1 !");
			dtx=dt*xi;
			u1p[l]+=(v1p[l]*dtx+0.5*a1p*dtx*dtx);
			u1m[l]=u1p[l];

			vxip=v1p[l]+a1p*dtx;
			vxim=v1m[l]+a1m*dtx;
			vxi=0.5*(vxip+vxim);

			a1=(s11p[l]-s11m[l])/dxr[0]+ds12;
			
			dte=dt*(1.0-xi);
			u1p[l]+=(vxi*dte+0.5*a1*dte*dte);
			u1m[l]=u1p[l];
			v1p[l]=vxi+a1*dte;
			v1m[l]=v1p[l];
		}
		break;
	case 1:		// CLOSED PART
		a1=(s11p[l]-s11m[l])/dxr[0]+ds12;
		if(s11b_now <= 0.0){	// CLOSE --> CLOSE
			u1p[l]+=(v1p[l]*dt+0.5*a1*dt2);
			u1m[l]=u1p[l];
			v1p[l]+=(a1*dt);
			v1m[l]=v1p[l];
		}else{			// CLOSE --> OPEN
			iop[l]=0;
			xi=fabs(s11b_pre)/(fabs(s11b_now)+fabs(s11b_pre)+EPS)-0.5;
			dtx=dt*xi;
			dte=dt*(1.0-xi);

			if(xi <0.0) xi=0.0;
			if(xi >1.0) puts("xi_o >1 !");

			u1p[l]+=(v1p[l]*dtx+0.5*a1*dtx*dtx);
			u1m[l]=u1p[l];

			v1p[l]+=(a1*dtx);
			v1m[l]=v1p[l];
			
			a1p=( 2.0*s11p[l]/dxr[0]+ds12);
			a1m=(-2.0*s11m[l]/dxr[0]+ds12);

			if(a1p < a1m ) dte=0.0;
			u1p[l]+=(v1p[l]*dte+0.5*a1p*dte*dte);
			u1m[l]+=(v1m[l]*dte+0.5*a1m*dte*dte);
			if( u1p[l] < u1m[l]) printf("u1p=%lf,  u1m=%lf\n",u1p[l],u1m[l]);
			v1p[l]+=(a1p*dte);
			v1m[l]+=(a1m*dte);
		}
		break;
	};
}

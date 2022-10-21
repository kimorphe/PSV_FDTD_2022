#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdm2d.h"
#include "crack.h"
//using namespace std;

Crack::Crack(){ 	// Class constructor
}

void Crack::gen_Crack(double *Y1, double *Y2, Dom2D dom){
	int i,j,jtyp;
	int indx1[2],indx2[2];

	int ndat, *ktmp;

	for(i=0;i<2;i++) {
		X1[i]=Y1[i];
		X2[i]=Y2[i];
	}

//	Count number of grids;
	Ng=0;
	switch(idir){
		case 0:	//horizontal crack
			jtyp=1;
			cod2indx(X1,jtyp,dom.Xa,dom.dx,indx1);
			cod2indx(X2,jtyp,dom.Xa,dom.dx,indx2);
			indx2cod(indx1[0],indx1[1],jtyp,dom.Xa,dom.dx,X1);
			indx2cod(indx2[0],indx2[1],jtyp,dom.Xa,dom.dx,X2);
			Len=dist2d(X1,X2)+dom.dx[0];
			j=indx1[1];
			if(j<0) break;
			if(j>=dom.Ndiv[1]) break;
			ndat=ceil(Len/dom.dx[0])+1;
			ktmp=(int *)malloc(sizeof(int)*ndat);
			for(i=indx1[0];i<=indx2[0];i++){
				if(i < 0) continue;
				if(i>=dom.Ndiv[0]) continue;
				if(dom.kcell[i][j] <0) continue;
				ktmp[Ng]=ij2l(i,j,jtyp,dom.Ndiv);	
				Ng++;
			}
			Len=dom.dx[0]*Ng;
			break;
		case 1:	// vertical crack
			jtyp=0;
			cod2indx(X1,jtyp,dom.Xa,dom.dx,indx1);
			cod2indx(X2,jtyp,dom.Xa,dom.dx,indx2);
			indx2cod(indx1[0],indx1[1],jtyp,dom.Xa,dom.dx,X1);
			indx2cod(indx2[0],indx2[1],jtyp,dom.Xa,dom.dx,X2);

			Len=dist2d(X1,X2)+dom.dx[1];

			i=indx1[0];
			if(i<0) break;
			if(i>=dom.Ndiv[0]) break;
			ndat=ceil(Len/dom.dx[1])+1;
			ktmp=(int *)malloc(sizeof(int)*ndat);
			for(j=indx1[1];j<=indx2[1];j++){
				if(j>=dom.Ndiv[1]) continue;
				if(j<0) continue;
				if(dom.kcell[i][j] <0) continue;
				ktmp[Ng]=ij2l(i,j,jtyp,dom.Ndiv);
				Ng++;
			}
			Len=dom.dx[1]*Ng;
			break;	
	};

	stat=(int *)malloc(sizeof(int)*Ng);
	iop=(int *)malloc(sizeof(int)*Ng);
	c0=(double *)malloc(sizeof(double)*Ng);
	c1=(double *)malloc(sizeof(double)*Ng);
	
	for(i=0;i<Ng;i++){
		stat[i]=0;	// status set to "stick"
		iop[i]=0;	// status set to "open"
		c0[i]=0.0;
		c1[i]=0.0;
	}
	iop[0]=1;
	iop[Ng-1]=1;

	s11p=(double *)malloc(sizeof(double)*Ng);
	s11m=(double *)malloc(sizeof(double)*Ng);
	s22p=(double *)malloc(sizeof(double)*Ng);
	s22m=(double *)malloc(sizeof(double)*Ng);
	s12=(double *)malloc(sizeof(double)*(Ng+1));
	for(i=0;i<Ng;i++){
		s11p[i]=0.0; s11m[i]=0.0;
		s22p[i]=0.0; s22m[i]=0.0;
		s12[i]=0.0;
	}
	s12[Ng]=0.0;

	switch(idir){
	case 0:		//Horizontal Crack
		v2p=(double *)malloc(sizeof(double)*Ng);
		v2m=(double *)malloc(sizeof(double)*Ng);
		u2p=(double *)malloc(sizeof(double)*Ng);
		u2m=(double *)malloc(sizeof(double)*Ng);
		s22b=(double *)malloc(sizeof(double)*Ng);
		kv2=(int *)malloc(sizeof(int)*Ng);
		for(i=0;i<Ng;i++){
			v2p[i]=0.0; v2m[i]=0.0;
			u2p[i]=0.0; u2m[i]=0.0;
			s22b[i]=0.0;
			kv2[i]=ktmp[i];
		}
		break;
	case 1:		//Vertical Crack
		v1p=(double *)malloc(sizeof(double)*Ng);
		v1m=(double *)malloc(sizeof(double)*Ng);
		u1p=(double *)malloc(sizeof(double)*Ng);
		u1m=(double *)malloc(sizeof(double)*Ng);
		s11b=(double *)malloc(sizeof(double)*Ng);
		kv1=(int *)malloc(sizeof(int)*Ng);
		for(i=0;i<Ng;i++){
			v1p[i]=0.0; v1m[i]=0.0;
			u1p[i]=0.0; u1m[i]=0.0;
			s11b[i]=0.0;
			kv1[i]=ktmp[i];
		}
		break;
	};
	
};

void Crack :: print_grid(int num,Dom2D dom){
	char fname[128];
	sprintf(fname,"crack%d.grd",num);
	FILE *fp=fopen(fname,"w");
	int i,j,l,ityp;
	double xcod[2];

	switch(idir){
		case 0: //horizontal crack
			ityp=1;
			for(l=0;l<Ng;l++){
				l2ij(kv2[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0],xcod[1]);
			}
			break;
		case 1: //vertical crack
			ityp=0;
			for(l=0;l<Ng;l++){
				l2ij(kv1[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0],xcod[1]);
			}
			break;
	};

	fclose(fp);
}
void Crack :: out_cod(FILE *fp,int it, int Nt,double dt, Dom2D dom){
	int i,j,l,ityp;
	double xcod[2];
	if(it==0){
		fprintf(fp,"## Nt (time steps), Ng (number of nodes) \n");
		fprintf(fp,"%d %d \n",Nt,Ng);
		fprintf(fp,"## dt (time increment) \n");
		fprintf(fp,"%lf \n",dt);
		fprintf(fp,"## (x1 x2): nodal coordinates \n");

		switch(idir){
		case 0	:	// horizontal crack
			ityp=1;
			for(l=0;l<Ng;l++){
				l2ij(kv2[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0],xcod[1]);
			}	
			break;
		case 1	:	// vertical crack
			ityp=0;
			for(l=0;l<Ng;l++){
				l2ij(kv1[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0],xcod[1]);
			}
			break;
		};
		fprintf(fp,"## vp, vm, up, um \n");

	}

	switch(idir){
	case 0:		//horizontal crack
		for(l=0;l<Ng;l++){
			fprintf(fp,"%lf %lf %lf %lf\n",v2p[l],v2m[l],u2p[l],u2m[l]);
		}
		break;
	case 1:		//vertical crack
		for(l=0;l<Ng;l++){
			fprintf(fp,"%lf %lf %lf %lf\n",v1p[l],v1m[l],u1p[l],u1m[l]);
		}
		break;
	};
}
void Crack :: out_csd(FILE *fp,int it, int Nt,double dt, Dom2D dom, Fld2D fld){
	int i,j,l,ityp;
	double xcod[2],vp,vm;
	if(it==0){
		fprintf(fp,"## Nt (time steps), Ng (number of nodes) \n");
		fprintf(fp,"%d %d \n",Nt,Ng-1);
		fprintf(fp,"## dt (time increment) \n");
		fprintf(fp,"%lf \n",dt);
		fprintf(fp,"## (x1 x2): nodal coordinates \n");

		switch(idir){
		case 0	:	// horizontal crack
			ityp=1;
			for(l=0;l<Ng-1;l++){
				l2ij(kv2[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0]+0.5*dom.dx[0],xcod[1]);
			}	
			break;
		case 1	:	// vertical crack
			ityp=0;
			for(l=0;l<Ng-1;l++){
				l2ij(kv1[l],ityp,dom.Ndiv,&i,&j);
				indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
				fprintf(fp,"%lf %lf \n",xcod[0],xcod[1]+0.5*dom.dx[1]);
			}
			break;
		};
		fprintf(fp,"## vp, vm \n");

	}

	switch(idir){
	case 0:		//horizontal crack
		ityp=1;
		for(l=0;l<Ng-1;l++){
			l2ij(kv2[l],ityp,dom.Ndiv,&i,&j);
			vp=fld.v1[i+1][j];
			vm=fld.v1[i+1][j-1];
			fprintf(fp,"%lf %lf\n",vp,vm);
		}
		break;
	case 1:		//vertical crack
		ityp=0;
		for(l=0;l<Ng-1;l++){
			l2ij(kv1[l],ityp,dom.Ndiv,&i,&j);
			vp=fld.v2[i][j+1];
			vm=fld.v2[i-1][j+1];
			fprintf(fp,"%lf %lf\n",vp,vm);
		}
		break;
	};
}


void Crack :: s2v(Fld2D fld, Dom2D dom, double dt){

	int i,j,l,jtyp;
	double dt2=dt*0.5;
	double s22b_now,s11b_now;
//	double hx,hy,dtr;

//	hx=dom.dx[0];
//	hy=dom.dx[1];
//	dtr=dt/dom.rho;

	switch(idir){
	case 0:		// Horizontal crack
		jtyp=1;	// node type set to v2

		u2p[0]+=(v2p[0]*dt2);
		u2p[Ng-1]+=(v2p[Ng-1]*dt2);
		
		for(l=1;l<Ng-1;l++){
			l2ij(kv2[l],jtyp,dom.Ndiv,&i,&j);
			s22b_now=(fld.s22[i][j]+fld.s22[i][j-1])*0.5;
			COD_hrz(l,dt,dom.rho,dom.dx,s22b_now);
			fld.v2[i][j]=v2p[l];
		}

		l2ij(kv2[0],jtyp,dom.Ndiv,&i,&j);
		v2p[0]=fld.v2[i][j];

		l2ij(kv2[Ng-1],jtyp,dom.Ndiv,&i,&j);
		v2p[Ng-1]=fld.v2[i][j];

		u2p[0]+=(v2p[0]*dt2);
		u2p[Ng-1]+=(v2p[Ng-1]*dt2);

		v2m[0]=v2p[0];
		v2m[Ng-1]=v2p[Ng-1];
		u2m[0]=u2p[0];
		u2m[Ng-1]=u2p[Ng-1];

		break;
	case 1:		// Vertical crack
		jtyp=0;	// node type set to v1 

		u1p[0]+=(v1p[0]*dt2);
		u1p[Ng-1]+=(v1p[Ng-1]*dt2);

		for(l=1;l<Ng-1;l++){
			l2ij(kv1[l],jtyp,dom.Ndiv,&i,&j);
			s11b_now=(fld.s11[i-1][j]+fld.s11[i][j])*0.5;
			COD_vrt(l,dt,dom.rho,dom.dx,s11b_now);
			fld.v1[i][j]=v1p[l];
		}

		l2ij(kv1[0],jtyp,dom.Ndiv,&i,&j);
		v1p[0]=fld.v1[i][j];

		l2ij(kv1[Ng-1],jtyp,dom.Ndiv,&i,&j);
		v1p[Ng-1]=fld.v1[i][j];

		u1p[0]+=(v1p[0]*dt2);
		u1p[Ng-1]+=(v1p[Ng-1]*dt2);

		v1m[0]=v1p[0];
		v1m[Ng-1]=v1p[Ng-1];
		u1m[0]=u1p[0];
		u1m[Ng-1]=u1p[Ng-1];

		break;
	};
}
void Crack :: v2s(Fld2D fld, Dom2D dom, double dt){

	int i,j,l,jtyp;
	double hx,hy,amu,almb;
	double Ax0,Ax1;
	double Ay0,Ay1;
	double dv1p,dv1m,dv2p,dv2m;

	hx=dom.dx[0];
	hy=dom.dx[1];
	amu=dom.amu;
	almb=dom.almb;

	Ax0=(almb+2.*amu)*dt/hx;
	Ay0=almb*dt/hy;

	Ax1=almb*dt/hx;
	Ay1=(almb+2.*amu)*dt/hy;

	switch(idir){
	case 0:		// horizontal crack
		jtyp=1;
		for(l=0;l<Ng;l++){
			l2ij(kv2[l],jtyp,dom.Ndiv,&i,&j);

			dv1p=fld.v1[i+1][j]-fld.v1[i][j];
			dv1m=fld.v1[i+1][j-1]-fld.v1[i][j-1];
			dv2p=fld.v2[i][j+1]-v2p[l];
			dv2m=v2m[l]-fld.v2[i][j-1];

			s22b[l]=0.5*(s22p[l]+s22m[l]);

			s11p[l]+=Ax0*dv1p+Ay0*dv2p;
			s11m[l]+=Ax0*dv1m+Ay0*dv2m;
			s22p[l]+=Ax1*dv1p+Ay1*dv2p;
			s22m[l]+=Ax1*dv1m+Ay1*dv2m;

			fld.s11[i][j]=s11p[l];
			fld.s11[i][j-1]=s11m[l];
			fld.s22[i][j]=s22p[l];
			fld.s22[i][j-1]=s22m[l];
			
		}
		for(l=1;l<Ng;l++){
			l2ij(kv2[l],jtyp,dom.Ndiv,&i,&j);
			switch(iop[l-1]+iop[l]){
			case 2:	// close
				s12[l]=fld.s12[i][j];
//				s12[l]=0.0;
//				fld.s12[i][j]=0.0;
				break;	
			default :
				s12[l]=0.0;
				fld.s12[i][j]=s12[l];
				break;
			};
		}
		l2ij(kv2[0],jtyp,dom.Ndiv,&i,&j);
		s12[0]=fld.s12[i][j];
		l2ij(kv2[Ng-1],jtyp,dom.Ndiv,&i,&j);
		s12[Ng]=fld.s12[i+1][j];
		break;
	case 1:		// vertical crack
		jtyp=0;
		for(l=0;l<Ng;l++){
			l2ij(kv1[l],jtyp,dom.Ndiv,&i,&j);

			dv1p=fld.v1[i+1][j]-v1p[l];
			dv1m=v1m[l]-fld.v1[i-1][j];

			dv2p=fld.v2[i][j+1]-fld.v2[i][j];
			dv2m=fld.v2[i-1][j+1]-fld.v2[i-1][j];

			s11b[l]=0.5*(s11p[l]+s11m[l]);

			s11p[l]+=Ax0*dv1p+Ay0*dv2p;
			s11m[l]+=Ax0*dv1m+Ay0*dv2m;
			s22p[l]+=Ax1*dv1p+Ay1*dv2p;
			s22m[l]+=Ax1*dv1m+Ay1*dv2m;
		
			fld.s11[i-1][j]=s11m[l];
			fld.s11[i][j]=s11p[l];
			fld.s22[i-1][j]=s22m[l];
			fld.s22[i][j]=s22p[l];

//			if(l < Ng-1) fld.s12[i][j+1]=0.0;
		}

		for(l=1;l<Ng;l++){
			l2ij(kv1[l],jtyp,dom.Ndiv,&i,&j);
			switch(iop[l-1]+iop[l]){
			case 2:	// close
				s12[l]=fld.s12[i][j];
//				s12[l]=0.0;
//				fld.s12[i][j]=0.0;
				break;	
			default :
				s12[l]=0.0;
				fld.s12[i][j]=s12[l];
				break;
			};
		}
		l2ij(kv1[0],jtyp,dom.Ndiv,&i,&j);
		s12[0]=fld.s12[i][j];
		l2ij(kv1[Ng-1],jtyp,dom.Ndiv,&i,&j);
		s12[Ng]=fld.s12[i][j+1];

		break;
	};
}

double Crack :: IsSlip(int l, double str_tmp){
	
	double  str,sgn;

	sgn=1.0;
	if(str_tmp<0.0) sgn=-1.0;

	switch(stat[l]){
	case 0 : // stick 
		if(fabs(str_tmp) > c0[l]){
			str=sgn*c0[l];	
			stat[l]=1;
		}else{
			str=str_tmp;
		}
		break;
	case 1 : // slip
		if(fabs(str_tmp) >= c1[l]){
			str=sgn*c1[l];	
		}else{
			str=str_tmp;
			stat[l]=0;
		}
		break;
	};

	return str;
}
void Crack :: set_c01(int nseg, double *xi, double *c_0, double *c_1){

	int i,iseg;
	double ss,ds,s2;

	ds=(X2[0]-X1[0])/Ng;
	if(idir==1) ds=(X2[1]-X1[1])/Ng;

	ss=0.5*ds;
	s2=xi[0]*Len;
	iseg=0;
	for(i=0;i<Ng;i++){
		if(ss>s2){
			if(iseg == nseg-1) continue;
			iseg++;
			s2+=xi[iseg]*Len;
		}
		c0[i]=c_0[iseg];
		c1[i]=c_1[iseg];
//		printf("ss=%lf, c0=%lf, c1=%lf\n",ss,c0[i],c1[i]);
		ss+=ds;
	}

}

void Crack :: setIC(Dom2D dom, Fld2D fld){
	
	int i,j,l,jtyp;

	switch(idir){
	case 0:		//Horizontal Crack
		jtyp=1; 
		for(l=0;l<Ng;l++){
			l2ij(kv2[l],jtyp,dom.Ndiv,&i,&j);
			v2p[l]=fld.v2[i][j];
			v2m[l]=fld.v2[i][j];
			u2p[l]=0.0;
			u2m[l]=0.0;

			s11p[l]=fld.s11[i][j];
			s11m[l]=fld.s11[i][j-1];
			s22p[l]=fld.s22[i][j];
			s22m[l]=fld.s22[i][j-1];

			s12[l]=fld.s12[i][j];

			s22b[l]=0.5*(s22p[l]+s22m[l]);
		}
		s12[Ng]=fld.s12[i+1][j];
		break;
	case 1:		//Vertical Crack
		jtyp=0; 
		for(l=0;l<Ng;l++){
			l2ij(kv1[l],jtyp,dom.Ndiv,&i,&j);
			v1p[l]=fld.v1[i][j];
			v1m[l]=fld.v1[i][j];
			u1p[l]=0.0;
			u1m[l]=0.0;

			s11p[l]=fld.s11[i][j];
			s11m[l]=fld.s11[i-1][j];
			s22p[l]=fld.s22[i][j];
			s22m[l]=fld.s22[i-1][j];

			s12[l]=fld.s12[i][j];

			s11b[l]=0.5*(s11p[l]+s11m[l]);
		}
		s12[Ng]=fld.s12[i][j+1];
		break;
	}
	
}


int make_cracks(char *fname, Dom2D dom, Crack **crk){
	FILE *fp=fopen(fname,"r");
	int i,j,ncrk,idir,nseg;
	char cbff[128];
	double Y1[2],Y2[2];
	double *cdat0,*cdat1,*xi;
	double cod0;

	ncrk=0;
	while( fgets(cbff,8,fp)!=NULL){
		if(strcmp(cbff,"##Crack")==0){
			fgets(cbff,2,fp);
			fgets(cbff,128,fp);
			fscanf(fp,"%d\n",&ncrk);
			
			(*crk)=(Crack *)malloc(sizeof(Crack)*ncrk);
			for(i=0;i<ncrk;i++){
				fgets(cbff,128,fp);
				fscanf(fp,"%d %lf %lf %lf %lf\n",&idir,Y1,Y1+1,Y2,Y2+1);
				(*crk)[i].idir=idir;
				(*crk)[i].gen_Crack(Y1,Y2, dom);

				fgets(cbff,128,fp);
				fscanf(fp,"%d\n",&nseg);
//					printf("nseg=%d\n",nseg);
				cdat0=(double *)malloc(sizeof(double)*nseg);
				cdat1=(double *)malloc(sizeof(double)*nseg);
				xi=(double *)malloc(sizeof(double)*nseg);

				fgets(cbff,128,fp);
				for(j=0;j<nseg;j++){
					fscanf(fp,"%lf %lf %lf\n",xi+j,cdat0+j,cdat1+j);
//					printf("%lf %lf %lf\n",xi[j],cdat0[j],cdat1[j]);
				}
				(*crk)[i].set_c01(nseg,xi,cdat0,cdat1);
				free(xi);
				free(cdat0);
				free(cdat1);
				fgets(cbff,128,fp);
				fscanf(fp,"%lf\n",&cod0);
				(*crk)[i].cod0=cod0;
				printf("cod0=%lf\n",cod0);
			}			
		}
	};
	fclose(fp);
	return ncrk;
}


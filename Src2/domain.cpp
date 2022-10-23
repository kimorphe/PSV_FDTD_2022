#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"fdm2d.h"

double* DOMAIN::ij2xf(int i, int j, int type){

	if(type==0){	// v3 (midpoint) grid 
		xf[0]=Xa[0]+(i+0.5)*dx[0];
		xf[1]=Xa[1]+(j+0.5)*dx[1];
		return(xf);
	};

	if(type==1){ // q1-grid  (leftword shifted grid)
		xf[0]=Xa[0]+i*dx[0];
		xf[1]=Xa[1]+(j+0.5)*dx[1];
		return(xf);
	};

	if(type==2){ // q2-grid (upward shifted grid)
		xf[0]=Xa[0]+(i+0.5)*dx[0];
		xf[1]=Xa[1]+j*dx[1];
		return(xf);
	};

	if(type==3){ // vertices 
		xf[0]=Xa[0]+i*dx[0];
		xf[1]=Xa[1]+j*dx[1];
		return(xf);
	};
	return(NULL);
};	
void DOMAIN::setup(double *xll, double *wdt, double *dh){
	Xa=xll;
	Wd=wdt;
	dx=dh;
};
void DOMAIN::print_prms(){
	printf("---------- DOMAIN SETTING-------------\n");
	printf("Xa=%lf, %lf\n",Xa[0],Xa[1]);
	printf("Wd=%lf, %lf\n",Wd[0],Wd[1]);
	printf("Ndiv=%d, %d\n",Ndiv[0],Ndiv[1]);
	printf("dx=%lf, %lf\n",dx[0],dx[1]);
	printf("--------------------------------------\n");
};
void DOMAIN::fwrite(){
	FILE *fp=fopen("domain.out","w");

	fprintf(fp,"  Xa=(%lf, %lf)\n",Xa[0],Xa[1]);
	fprintf(fp,"  Wd=(%lf, %lf)\n",Wd[0],Wd[1]);
	fprintf(fp,"Ndiv=(%d, %d)\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"  dx=(%lf, %lf)\n",dx[0],dx[1]);
	fprintf(fp,"---------  PML ----------\n");
	fprintf(fp,"  Ha=(%lf, %lf)\n",Ha[0],Ha[1]);
	fprintf(fp,"  Hb=(%lf, %lf)\n",Hb[0],Hb[1]);
	fprintf(fp," NHa=(%d, %d)\n",NHa[0],NHa[1]);
	fprintf(fp," NHb=(%d, %d)\n",NHa[0],NHa[1]);
	fprintf(fp,"---Material Consntants---\n");
	fprintf(fp,"  cT=%lf [km/s]\n",cT);
	fprintf(fp," rho=%lf [g/cm3]\n",rho);
	fprintf(fp," amu=%lf [GPa]\n",amu);
	fclose(fp);
	printf(" --> domain.out\n");
};
void DOMAIN::init(int *nxy){
	Ndiv=nxy;
	ndat=Ndiv[0]*Ndiv[1];

	int *p=(int *)malloc(sizeof(int)*ndat);
	kcell=(int **)malloc(sizeof(int*)*Ndiv[0]);

	int k=0;
	for( int i=0; i<Ndiv[0];i++){
		kcell[i]=p+i*Ndiv[1];
	for( int j=0; j<Ndiv[1];j++){
		p[k++]=0;
	}
	}
};
void DOMAIN::print_kcell(){
	for( int i=0; i<Ndiv[0];i++){
	for( int j=0; j<Ndiv[1];j++){
		printf("%d ", kcell[i][j]);
	}
		printf("\n");
	}
}
//  ----------- DOMAIN PERFORATION  --------------
void DOMAIN::perfo_ellip(char *fname){
	int i,j;
	double xcod[2],xc[2],a,b;
	double xx,yy;

	FILE *fp;
	char cbff[8];
	fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname); 
	i=0;
	while(fgets(cbff,8,fp) !=NULL){
		if(strcmp(cbff,"##Ellip")==0){
			fscanf(fp,"%lf, %lf, %lf, %lf\n",xc,xc+1,&a,&b);
			printf("%lf %lf %lf %lf\n",xc[0],xc[1],a,b);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				xx=(xcod[0]-xc[0])/a;
				yy=(xcod[1]-xc[1])/b;
				if(xx*xx+yy*yy <1.0) kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
};
void DOMAIN::slit(char *fname){
	int i,j;
	double xs[2],xe[2],xc[2],wx,wy;
	int ixs[2],ixe[2];

	FILE *fp;
	char cbff[7];
	fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname); 
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Rect")==0){
			printf("Flag ##Rect found !!\n");
			fscanf(fp,"%lf, %lf\n",xc,xc+1); // center
			fscanf(fp,"%lf, %lf\n",&wx,&wy); // width and height
			//fscanf(fp,"%lf, %lf\n",xs,xs+1);
			//fscanf(fp,"%lf, %lf\n",xe,xe+1);
			xs[0]=xc[0]-0.5*wx;
			xe[0]=xc[0]+0.5*wx;
			xs[1]=xc[1]-0.5*wy;
			xe[1]=xc[1]+0.5*wy;

			for(i=0;i<2;i++){
				ixs[i]=floor((xs[i]-Xa[i])/dx[i]);
				ixe[i]=floor((xe[i]-Xa[i])/dx[i]);
				if( ixs[i] <0) ixs[i]=0;
				if( ixe[i] >= Ndiv[i]) ixe[i]=Ndiv[i]-1;
			}

			for(i=ixs[0]; i<=ixe[0]; i++){
			for(j=ixs[1]; j<=ixe[1]; j++){
				 kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
}
//-------------EXPORT GEOMETRY DATA  ------------------
void DOMAIN:: out_kcell(){

	FILE *fp=fopen("kcell.dat","w");
	int i,j;
	double Xb[2];
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];


	fprintf(fp,"# Xa[2], Xb[2] (computational domian)\n");
	fprintf(fp," %lf, %lf\n %lf, %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d, %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# kcell[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp, "%d\n",kcell[i][j]);
	}
	}

	fflush(fp);
};
void DOMAIN::topography(char *fname){
	FILE *fp=fopen(fname,"r");
	int typ=1;	//piecewise linear function
	double xx,yy;
	double *xknot,*yknot;
	int nseg;
	char cbff[9];

	if(fp==NULL) show_msg(fname);

	int i,j,jy;
	int isgn=1;
	while(fgets(cbff,9,fp) !=NULL){
	if(strcmp(cbff,"##Tpgrph")==0){
		fscanf(fp,"%d\n",&typ);
		printf("type=%d\n",typ);
		if(typ<0) isgn=-1;
		typ=abs(typ);

		if(typ<3){
			fscanf(fp,"%d\n",&nseg);
			printf("nseg=%d\n",nseg);
			xknot=(double *)malloc(sizeof(double)*nseg);
			yknot=(double *)malloc(sizeof(double)*nseg);
			for(i=0;i<nseg;i++) fscanf(fp,"%lf, %lf\n",xknot+i,yknot+i);

			for(i=0;i<Ndiv[0];i++){
				xx=Xa[0]+(i+0.5)*dx[0];
				yy=pwfun(typ,xx,xknot,yknot,nseg,1.0);
				jy=int((yy-Xa[1])/dx[1]);
				if(isgn==1){
					for(j=jy;j<Ndiv[1];j++) kcell[i][j]=1;
				}
				if(isgn==-1){
					for(j=0;j<jy;j++) kcell[i][j]=1;
				}
			};
		free(xknot);
		free(yknot);
		}	
	}
	};	// end of while loop
	fclose(fp);

};
double pwfun(
	int typ, 	// piecewise constant(0), linear(1)
	double x, 	// x-coordinate
	double *xknot, double *yknot,  // discontinuous/broken  points
	int ndiv, // number of intervals
	double y0 // magnicifation
){

	int idiv=0;

	for(idiv=0; idiv<ndiv; idiv++){
		if(xknot[idiv]>x) break;
	};
	if(idiv==0) return(yknot[0]);
	if(idiv==ndiv) return(yknot[ndiv-1]);

	double y;
	switch(typ){
	case 1: // piecsewise linear function
		double xl,xr,xi;
		xl=xknot[idiv-1];
		xr=xknot[idiv];
		xi=(x-xl)/(xr-xl);
		y=(1.0-xi)*yknot[idiv-1]+xi*yknot[idiv];
		break;
	case 2: //piecewise constant function
		return(yknot[idiv]*y0);
	};

	return(y*y0);
};
int DOMAIN::find_q1bnd(int nx, double y){

	int jy=int((y-Xa[1])/dx[1]);
	int i,ibnd,il,ir;
	switch(nx){ 
	case 1:
		ir=1;
		ibnd=-1;
		for(i=Ndiv[0]-1;i>=0;i--){
			il=kcell[i][jy];
			if((ir-il)==1){
				ibnd=i+1;
				break;
			};
			ir=il;
		};
		break;
	case -1:
		il=1;
		ibnd=-1;
		for(i=0;i<Ndiv[0];i++){
			ir=kcell[i][jy];
			if((il-ir)==1){
				ibnd=i;
				break;
			};
			il=ir;
		};
		break;
	};
	return(ibnd);
};
int DOMAIN::find_q2bnd(int ny, double x){
	int ix=int((x-Xa[0])/dx[0]);
	int j,jbnd,iu,id;

	switch(ny){
	case 1:
		iu=1;
		jbnd=-1;
		for(j=Ndiv[1]-1;j>=0;j--){
			id=kcell[ix][j];
			if((iu-id)==1){
				jbnd=j+1;
				break;
			};
			iu=id;
		};
		break;
	case -1:
		id=1;
		jbnd=-1;
		for(j=0; j<Ndiv[1]; j++){
			iu=kcell[ix][j];
			if((id-iu)==1){
				jbnd=j;
				break;
			};
			id=iu;
		};
		break;
	};
	return(jbnd);
};
void DOMAIN::PML_setup(double gmm){
	double L;
	int m=2;
	double eps=1.e-10;
	for(int i=0; i<2; i++){
		//L=2.*Ha[i];
		L=Ha[i];
		A0[i]=0.0; 
		if(L>eps) A0[i]=-0.5*(m+1)*cT/pow(L,m+1)*logf(gmm);

		//L=2.*Hb[i];
		L=Hb[i];
		B0[i]=0.0; 
		if(L>eps) B0[i]=-0.5*(m+1)*cT/pow(L,m+1)*logf(gmm);
	};
};
double DOMAIN::PML_dcy(int idir, double xy){ // direction: idir=0(x), 1(y)
	double r;
	double s=0.0,C0=0.0;
	r=Xa[idir]+Ha[idir]-xy;
	if(r > 0.0){
		s=r;
		C0=A0[idir];
	};
	r=xy-(Xa[idir]+Wd[idir]-Hb[idir]);
	if(r > 0.0){
		s=r;
		C0=B0[idir];
	}
	return(C0*s*s);
};


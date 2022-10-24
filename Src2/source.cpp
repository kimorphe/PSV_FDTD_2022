#include<stdio.h>
#include<stdlib.h>
#include"fdm2d.h"

void SOURCE::print(){
	printf("----------- Source Parameters -----------\n");
	printf("ng=%d\n",ng);
	double x,y;
	printf("wv ID=%d\n",iwv);

	double xofst=0.5, yofst=0.5;
	if(type==1) xofst=0.0;
	if(type==2) yofst=0.0;

	for(int i=0;i<ng;i++){
		printf("(i,j)=(%d, %d) ",isrc[i],jsrc[i]);
		x=Xa[0]+(isrc[i]+xofst)*dx[0];
		y=Xa[1]+(jsrc[i]+yofst)*dx[1];
		printf("(x,y)=(%lf, %lf)\n",x,y);
	};
};
void SOURCE::set_center(){
	x0=dx[0]*(isrc[0]+isrc[ng-1]+1.)*0.5+Xa[0];
	y0=dx[1]*(jsrc[0]+jsrc[ng-1]+1.)*0.5+Xa[1];
};
void SOURCE::set_inc_ang(double th, double ct){
	double PI=4.0*atan(1.0);
	th_in=th;
	th=th/180.*PI;	
	Ctd=sin(th)/ct*dx[type-1];
};
void SOURCE::mem_alloc(int n){
	ng=n;
	isrc=(int *)malloc(sizeof(int)*ng);
	jsrc=(int *)malloc(sizeof(int)*ng);
	ksrc=(int *)malloc(sizeof(int)*ng);
	for(int i=0;i<ng;i++){
		isrc[i]=-1;
		jsrc[i]=-1;
		ksrc[i]=-1;
	};
};
void TRNSDCR::init_bwv(int nn, double dtau){
	Nt=nn;
	dt=dtau;
	double *p=(double *)malloc(sizeof(double)*Nt*ng);
	bwv=(double **)malloc(sizeof(double *)*ng);
	int i;
	for(int i=0; i<Nt*ng; i++) p[i]=0.0;
	for(i=0;i<ng;i++) bwv[i]=p+i*Nt;
};
double TRNSDCR::val(int ig, double tt){
	if(tt<0.0) return(0.0);
	int it=int(tt/dt);
	if(it>Nt-1) return(0.0);
	if(it==Nt-1) return(bwv[ig][it]);
	double xi=tt/dt-it;
	return((1.-xi)*bwv[ig][it]+xi*bwv[ig][it+1]);
};
double TRNSDCR::amp_synth(int it){
	int ig;
	double amp=0.0,tt;
	for(ig=0;ig<ng;ig++){
		tt=it*dt-ig*Ctd;
		if(Ctd<0.0) tt=it*dt+(ng-1-ig)*Ctd;
		amp+=val(ig,tt);
	};
	return(amp/ng);
};
void TRNSDCR::record(int jt, double **fld){
	int ig,i,j;
	for(ig=0; ig<ng; ig++){
		i=isrc[ig];
		j=jsrc[ig];
		if(nml==1){ 
			if(type==1) i-=nml;
			if(type==2) j-=nml;
		};
		bwv[ig][jt]=fld[i][j]; 
	};
};
void TRNSDCR::fwrite_setting(char *fn, char *mode, int num){
	FILE *fp=fopen(fn,mode);
	if(fp==NULL) show_msg(fn);

	fprintf(fp,"#---------- T/R=%d -----------\n",ID);
	fprintf(fp,"# type=%d, nrml=%d\n",type,nml);
	fprintf(fp,"# th_in=%lf[deg], Ctd=%lf\n",th_in,Ctd);
	double x,y;
	double xofst=0.5, yofst=0.5;
	if(type==1) xofst=0.0;
	if(type==2) yofst=0.0;
	fprintf(fp,"# waveform ID=%d\n",iwv);
	fprintf(fp,"# ng=%d\n",ng);
	fprintf(fp,"# index (i,j),  coordinate (x,y)\n");
	for(int i=0;i<ng;i++){
		//printf("(i,j)=(%d, %d) ",isrc[i],jsrc[i]);
		x=Xa[0]+(isrc[i]+xofst)*dx[0];
		y=Xa[1]+(jsrc[i]+yofst)*dx[1];
		fprintf(fp,"%d, %d, %lf, %lf\n",isrc[i],jsrc[i],x,y);
		//printf("(x,y)=(%lf, %lf)\n",x,y);
	};

	fclose(fp);
};
void TRNSDCR::fwrite(){
	int i,j;
	char fname[128];
	sprintf(fname,"tr%d.out",ID);
	FILE *fp=fopen(fname,"w");
	fprintf(fp,"# Nt, dt\n");
	fprintf(fp,"%d, %lf\n",Nt,dt);
	double x,y;
	for(i=0;i<ng;i++){
		x=Xa[0]+dx[0]*(isrc[i]+0.5);
		y=Xa[1]+dx[1]*(jsrc[i]+0.5);
		fprintf(fp,"# %lf, %lf\n",x,y);
	for(j=0;j<Nt;j++){
		fprintf(fp,"%lf\n",bwv[i][j]);
	};
	};
	fclose(fp);
};
double TRNSDCR::mean_amp(int it){
	if(it<0) return(0.0);
	if(it>=Nt) return(0.0);
	double amp=0.0;
	for(int j=0;j<ng;j++) amp+=bwv[j][it];
	return(amp/ng);
};
void TRNSDCR::clear(){
	int i,j;
	for(i=0;i<ng;i++){
	for(j=0;j<Nt;j++){
		bwv[i][j]=0.0;
	}
	}	
};
void ARRAY::init(int nn, int nm){
	nele=nn;
	nmeas=nm;
	actv=(int *)malloc(sizeof(int)*nele*nmeas);
	a0=(double *)malloc(sizeof(double)*nele*nmeas);
	tdly=(double *)malloc(sizeof(double)*nele*nmeas);

	for(int i=0;i<nele*nmeas;i++){
		actv[i]=0;
		a0[i]=1.0;
		tdly[i]=0.0;
	};
	i0=0;
};
void ARRAY::print(){
	int k=0;
	printf("nmeas=%d\n",nmeas);
	for(int j=0;j<nmeas;j++){
	for(int i=0;i<nele;i++){
		printf("i=%d, actv=%d, a0=%lf, tdly=%lf\n",i,actv[k],a0[k], tdly[k]);
		k++;
	}
	printf("\n");
	}
};
/*
void ARRAY::fwrite(){
	char fname[128];
	int i,j,k;
	FILE *fp;
	int Nt=trs[0].Nt;
	double dt=trs[0].dt;
	sprintf(fname,"ary%d.out\n",round);
	fp=fopen(fname,"w");
	for(j=0;j<nele;j++){
		fprintf(fp,"# e=%d\n",j);
		for(k=0;k<Nt;k++){
			fprintf(fp,"%lf, %lf\n",dt*k,trs[j].mean_amp(k));
		}
	};
	fclose(fp);
};
*/

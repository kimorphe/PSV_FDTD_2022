#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fdm2d.h"

using namespace std;

int set_srcs(char *fname, Dom2D dom, Src **src){
//
//		Define Source Class Array
//		--> source parameters are read from src.inp
//
	double xs[2],xe[2];
	FILE *fp=fopen(fname,"r");
	char cbff[128];
	int nsrc,i,iseq,is_sync,dly_seq;
	int ifld,ityp;
	double pi=4.0*atan(1.0);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&is_sync);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&dly_seq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsrc);

	(*src)=(Src *)malloc(sizeof(Src)*nsrc);
	for(i=0;i<nsrc;i++){
		fgets(cbff,128,fp);
		fscanf(fp,"%d %d\n",&ityp,&ifld);
		printf("%d %d\n",ityp,ifld);

		fgets(cbff,128,fp);
		fscanf(fp,"%lf %lf\n",xs,xs+1);
		fscanf(fp,"%lf %lf\n",xe,xe+1);
		fgets(cbff,128,fp);
		fscanf(fp,"%d \n",&iseq);
	
//		(*src)[i].ityp=1;	// source type set to s22
		(*src)[i].ityp=ityp;	// source grid type 
		(*src)[i].ifld=ifld;	// B.C. type (Dirichlet/Neumann)

		(*src)[i].Xsa[0]=xs[0];
		(*src)[i].Xsa[1]=xs[1];
		(*src)[i].Xsb[0]=xe[0];
		(*src)[i].Xsb[1]=xe[1];
		(*src)[i].inmb=iseq;
		(*src)[i].dly_seq=dly_seq;
		(*src)[i].is_sync=is_sync;
		
		(*src)[i].gen_Grid(dom);
		(*src)[i].set_Wvfm();
	}


	fclose(fp);
	return nsrc;
}

int sync_src(char *fname, Dom2D dom,Src **src){
	int i,wvtyp,ndat=1;
	int is_sync; // 0:Yes, 1: No
	int dly_seq; // 0:steering, 1:focusing
	int nsrc;    // number of array elements
	double xc[2],x0[2],pin[2],xf[2];
	char cbff[128];
	double *dlt,dlt_min,dlt_max,ain;
	double pi=4.0*atan(1.0);

	FILE *fp=fopen(fname,"r");

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&is_sync);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&dly_seq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsrc);

	dlt=(double *)malloc(sizeof(double)*nsrc);	// element-wise delay time
	fclose(fp);

	ndat=nsrc;	// number of beam-forming parameters  
	if(is_sync==0) ndat=1;
	switch(dly_seq){
	case 0:	// beam steering
		fp=fopen(fname,"r");
		while(fgets(cbff,8,fp) != NULL){
			if(strcmp(cbff,"##InDir")==0){	
			fgets(cbff,128,fp);
			for(i=0;i<ndat;i++){
				fscanf(fp,"%lf %d\n",&ain,&wvtyp);
				ain=ain/180.0*pi;
				(*src)[i].ain=ain;
				(*src)[i].wvtyp=wvtyp;

			switch(wvtyp){
				case 0:
					(*src)[i].cin=dom.cL;
					break;
				case 1:
					(*src)[i].cin=dom.cT;
					break;
				case 2:
					(*src)[i].cin=dom.cR;
					break;
				default:
					puts("Illeagal inc-wave type specified.!");
					printf("wvtyp=%d\n",wvtyp);
					puts("-->process terminated.");
					exit(-1);
					break;
				};

//				(*src)[i].cin=dom.cL;
//				if(wvtyp==1) (*src)[i].cin=dom.cT;
//				printf("ain=%lf\n",ain);
			}
			}
		};
		fclose(fp);
		break;
	case 1:  // beam focusing
		fp=fopen(fname,"r");	
		while(fgets(cbff,8,fp)!=NULL){
			if(strcmp(cbff,"##Focus")==0){
			fgets(cbff,128,fp);
			for(i=0;i<ndat;i++){
				fscanf(fp,"%lf %lf %d\n",xf,xf+1,&wvtyp);
				(*src)[i].xf[0]=xf[0];
				(*src)[i].xf[1]=xf[1];
				(*src)[i].cin=dom.cL;
				(*src)[i].wvtyp=wvtyp;
				if(wvtyp==1) (*src)[i].cin=dom.cT;
				printf("xf=%lf %lf\n",xf[0],xf[1]);
				printf("i,ndat=%d %d\n",i,ndat);
			}
			}
		}
		fclose(fp);
		break;
	};

	switch(is_sync){
		case 0:	// synchronize array elements
		if(dly_seq==0){	// beam steering
			pin[0]=cos(ain);
			pin[1]=sin(ain);
			dlt_min=0.0;
			for(i=0;i<nsrc;i++){
				xc[0]=0.5*((*src)[i].Xsa[0]+(*src)[i].Xsb[0]);
				xc[1]=0.5*((*src)[i].Xsa[1]+(*src)[i].Xsb[1]);
				if(i==0){
					x0[0]=xc[0]; 
					x0[1]=xc[1];
				}
				
				dlt[i]=(xc[0]-x0[0])*pin[0]+(xc[1]-x0[1])*pin[1];
				dlt[i]/=((*src)[i].cin);
				if(dlt[i]<dlt_min) dlt_min=dlt[i];
			}
			for(i=0;i<nsrc;i++) (*src)[i].delay=dlt[i]-dlt_min;
//			for(i=0;i<nsrc;i++) printf("delay=%lf\n",(*src)[i].delay);
		}
		if(dly_seq==1){	//beam focusing
			for(i=0;i<nsrc;i++){
				xc[0]=0.5*((*src)[i].Xsa[0]+(*src)[i].Xsb[0]);
				xc[1]=0.5*((*src)[i].Xsa[1]+(*src)[i].Xsb[1]);
				dlt[i]=dist2d(xc,xf)/((*src)[i].cin);
				if(i==0) dlt_max=dlt[i];
				if(dlt[i]> dlt_max) dlt_max=dlt[i];
			}
			for(i=0;i<nsrc;i++) (*src)[i].delay=dlt_max-dlt[i];
			
		}
			break;
		case 1: //
			break;
	};
	
	return(0);
}

//	Source Class Constructor 
Src::Src(int ii){
	ityp=ii;
};
Src::Src(){
};

void Src::isType(){
	printf("source type is %d\n",ityp);
};
void Src::gen_Src(char *fname){
	FILE *fp;
	char ch[128];

	fp=fopen(fname,"r");
	fgets(ch,128,fp);
	fscanf(fp,"%lf %lf",Xsa,Xsa+1);	
	fscanf(fp,"%lf %lf\n",Xsb,Xsb+1);	
	fgets(ch,128,fp);
	fscanf(fp,"%d",&inmb);	

	fclose(fp);
};
double Src::set_Wvfm(){
	int i,ifile=inmb;
	char fname[128];
	inwvs=(InWv *)malloc(sizeof(InWv)*nwv);

	for(i=0;i<nwv;i++){
		if(inmb==-1) ifile=i;
		sprintf(fname,"inwv%d.dat",ifile);
		puts(fname);
		inwvs[i]=InWv(fname);
//		printf("dt=%lf\n",inwvs[i].dt);
//		inwvs[i].disp();
//		inwvs[i].gen_wvfm();
	}
	Nt=inwvs[0].Nt;
	dt=inwvs[0].dt;
	
	return inwvs[0].dt;
}
void Src::gen_Grid(Dom2D dom){
	int ind1[2],ind2[2];
	int iwv=0,ntmp;
	int i,j,il,id;

	cod2indx(Xsa,ityp,dom.Xa,dom.dx,ind1);
	cod2indx(Xsb,ityp,dom.Xa,dom.dx,ind2);
	
	switch(ityp){
	case 0:	// s11-source
		i=ind1[0];
		ngrd=ind2[1]-ind1[1]+1;
		ksrc=(int *)malloc(sizeof(int)*ngrd);
		wvID=(int *)malloc(sizeof(int)*ngrd);
		ngrd=0;
		for(j=ind1[1]; j<=ind2[1]; j++){
			il=1; 
			if(i>0) il=dom.kcell[i-1][j];
			ksrc[ngrd]=ij2l(i,j,ityp,dom.Ndiv);
			if(il==1) ksrc[ngrd]*=-1;
			wvID[ngrd]=iwv;
			if(inmb==-1) wvID[ngrd]=iwv++;
			ngrd++;
		}
		break;
	case 1: // s22-source
		j=ind1[1];
		ngrd=ind2[0]-ind1[0]+1;
		ksrc=(int *)malloc(sizeof(int)*ngrd);
		wvID=(int *)malloc(sizeof(int)*ngrd);
		ngrd=0;
		for(i=ind1[0]; i<=ind2[0];i++){
			id=1; 
			if(j>0) id=dom.kcell[i][j-1];
			ksrc[ngrd]=ij2l(i,j,ityp,dom.Ndiv);
			if(id==1) ksrc[ngrd]*=-1;
			wvID[ngrd]=iwv;
			if(inmb==-1) wvID[ngrd]=iwv++;
			ngrd++;
		}
		break;
	case 3:	// s12 source
		ngrd=ind2[0]-ind1[0]+1;
		ntmp=ind2[1]-ind1[1]+1;
		if(ntmp > ngrd){ // vertically aligned source grids
			ngrd=ntmp;
			i=ind1[0];
			ksrc=(int *)malloc(sizeof(int)*ngrd);
			wvID=(int *)malloc(sizeof(int)*ngrd);
			ngrd=0;
			for(j=ind1[1]; j<=ind2[1]; j++){
				ksrc[ngrd]=ij2l(i,j,ityp,dom.Ndiv);
				wvID[ngrd]=iwv;
				if(inmb==-1) wvID[ngrd]=iwv++;
				ngrd++;
			}
		}else{		 // horizontally aligned source grids
			j=ind1[1];
			ksrc=(int *)malloc(sizeof(int)*ngrd);
			wvID=(int *)malloc(sizeof(int)*ngrd);
			ngrd=0;
			for(i=ind1[0]; i<=ind2[0];i++){
				ksrc[ngrd]=ij2l(i,j,ityp,dom.Ndiv);
				wvID[ngrd]=iwv;
				if(inmb==-1) wvID[ngrd]=iwv++;
				ngrd++;
			}
		}
		break;
	};
//	nwv=wvID[ngrd-1]+1;
	nwv=iwv+1;
}
void Src::set_Delay(int *Nx, Dom2D dom){
	int j,indx[2];
	double xcod[2],x0[2],pin[2];
	double *dlt,dlt_min=0.0,dlt_max;

	dlt=(double *)malloc(sizeof(double)*ngrd);
	idly=(int *)malloc(sizeof(int)*ngrd);

	switch (is_sync){ 
	case 0:	// synchronize array 
		for(j=0;j<ngrd;j++){
			idly[j]=delay/dt;
		}
		break;
	case 1: // use array elements independently 
		if(dly_seq==0){ // beam steering
			pin[0]=cos(ain);
			pin[1]=sin(ain);
			for(j=0;j<ngrd;j++){
				l2ij(ksrc[j],Nx,indx,indx+1);
				indx2cod(indx[0],indx[1],ityp,dom.Xa,dom.dx,xcod);
				if(j==0){
					x0[0]=xcod[0];
					x0[1]=xcod[1];
				}
				dlt[j]=(xcod[0]-x0[0])*pin[0]+(xcod[1]-x0[1])*pin[1];
//				dlt[j]/=dom.cT;
				dlt[j]/=cin;
				if(dlt[j]<dlt_min) dlt_min=dlt[j];
			}

			for(j=0;j<ngrd;j++) idly[j]=(dlt[j]-dlt_min)/dt;
//			for(j=0;j<ngrd;j++) printf("tdly=%lf  idly=%d\n",dlt[j],idly[j]);
		}
		if(dly_seq==1){	// beam focusing
			for(j=0;j<ngrd;j++){
				l2ij(ksrc[j],Nx,indx,indx+1);
				indx2cod(indx[0],indx[1],ityp,dom.Xa,dom.dx,xcod);
//				dlt[j]=dist2d(xf,xcod)/dom.cT;
				dlt[j]=dist2d(xf,xcod)/cin;
				if(j==0) dlt_max=dlt[j];
				if(dlt[j] > dlt_max) dlt_max=dlt[j];
			}	
			for(j=0;j<ngrd;j++) idly[j]=(dlt_max-dlt[j])/dt;
		}

		break;
	};
}
void Src::out(int isrc, int *Nx,Dom2D dom){

//		isrc: source ID
//		Nx[2]: number of grids in x,y-directions
//		dom : domain class data 

	int j;
	int indx[2];
	double xcod[2];
	double t1,tt;
	char fname[128];
	sprintf(fname,"source%d.out",isrc);
	FILE *fout=fopen(fname,"w");

	fprintf(fout,"## Grid type(0:v1, 1:v2, 3:s12)\n");
	fprintf(fout,"%d\n",ityp);
	fprintf(fout,"## Incident wave type(P=0, S=1, R=2), delay pattern\n");
	fprintf(fout,"%d %d\n",wvtyp, dly_seq);
	if(dly_seq==0){		// steering
		fprintf(fout,"## icident angle\n");
		fprintf(fout,"%lf\n",ain);
	}else if(dly_seq==1){	// focusing
		fprintf(fout,"## focal point\n");
		fprintf(fout,"%lf %lf\n",xf[0],xf[1]);
	}
	fprintf(fout,"## Number of Grids\n");
	fprintf(fout,"%d\n",ngrd);
	fprintf(fout,"## Source Grids & waveform ID :\n");
	for(j=0;j<ngrd;j++){
		l2ij(abs(ksrc[j]),Nx,indx,indx+1);
		indx2cod(indx[0],indx[1],ityp,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf %d\n",xcod[0],xcod[1],wvID[j]);
	}
	fprintf(fout,"## dt, Nt\n");
	fprintf(fout,"%lf %d\n",dt,Nt);
	fprintf(fout,"## number of waveforms\n");
	fprintf(fout,"%d\n",nwv);
	for(j=0;j<nwv;j++){
		fprintf(fout,"## amp (wave ID=%d)\n",wvID[j]);
		t1=inwvs[j].t1;
		dt=inwvs[j].dt;
		Nt=inwvs[j].Nt;	
		for(int i=0;i<Nt;i++){
			tt=t1+dt*i;
//			fprintf(fout,"%lf %lf\n",tt,inwvs[j].amp[i]);
			fprintf(fout,"%lf\n",inwvs[j].amp[i]);
		}
//		fprintf(fout,"\n");
	}
//	fprintf(fout,"----------------------------------\n");
	fclose(fout);
}
//-----------------------------------------------------------

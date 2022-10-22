#include<iostream> 

#include<stdio.h>
#include<math.h> 
#include<string.h>
#include<stdlib.h>
#include<sys/stat.h>
#include "fdm2d.h"
#include "crack.h"
using namespace std;

int main(int argv, char *argc[]){

	int i;
	int ncrk;	// number of cracks
	int nsrc;	// number of source classes
	int nrec;	// number of receiver classes
	int Nt;		// number of time steps
	int ninc,nout,isum;
	double dt;
	FILE *fout,*fdat;
	FILE **fCOD,**fCSD;
	double xx;

//		Default Input Data Files
	char fgeom[128]="geom.inp";
	char fsrc[128]="src.inp";
	char ficon[128]="icon.inp";
	char fprnt[128]="mdata.out";

	puts("------ INPUT DATA FILES -------");
	if(argv > 1) strcpy(fgeom,argc[1]);
		printf(" 1. %s\n",fgeom);
	if(argv > 2) strcpy(fsrc,argc[2]);
		printf(" 2. %s\n",fsrc);
	if(argv > 3) strcpy(ficon,argc[3]);
		printf(" 2. %s\n",ficon);
	puts("-------------------------------");

	char fname[256];
	Src *src;	// Source class array
	Recs *rec;	// Receiver class array 
	Crack *crk;	// Crack class array


	fdat=fopen(fprnt,"w");
	
	
	
//		COMPUTATIONAL DOMAIN
	Dom2D dom1=Dom2D(fgeom);
	dom1.perfo(fgeom);// perforate computational domain
	dom1.perfo_ellip(fgeom);// perforate computational domain
	dom1.slit(fgeom);// perforate computational domain
	dom1.angled_slit(fgeom);// perforate computational domain
	dom1.perfo_para(fgeom);
	dom1.polygon(fgeom);
	dom1.Cut(fgeom);
	dom1.butt_weld(fgeom);
	//dom1.WireCut();
	dom1.out_kcell(); // export kcell data
	dom1.out_kcell_tight(); // export kcell data without PML

	ncrk=make_cracks(fgeom,dom1,&crk);
	fCOD=(FILE **)malloc(sizeof(FILE*)*ncrk);
	fCSD=(FILE **)malloc(sizeof(FILE*)*ncrk);
	for(i=0;i<ncrk;i++){
		crk[i].print_grid(i,dom1);
		sprintf(fname,"cod%d.out",i);
		fCOD[i]=fopen(fname,"w");
		sprintf(fname,"csd%d.out",i);
		fCSD[i]=fopen(fname,"w");
	}

//		WAVE FIELD DATA
	Fld2D fld1=Fld2D(dom1,ficon);
	fld1.gen_indx1(dom1);
	fld1.gen_indx2(dom1);
	fld1.gen_indx12(dom1);

	for(i=0;i<ncrk;i++) crk[i].setIC(dom1,fld1);

//		SETUP SOURCE
	nsrc=set_srcs(fsrc,dom1, &src);
	sync_src(fsrc,dom1,&src);
	dt=src[0].dt;
	Nt=src[0].Nt;

	dom1.CFL(dt);	// show CFL number 

	fld1.gridNum(src[0].ityp);
	for(i=0;i<nsrc;i++) src[i].out(i,fld1.Nx,dom1);
	for(i=0;i<nsrc;i++) src[i].set_Delay(fld1.Nx,dom1);
	for(i=0;i<nsrc;i++) fld1.del_src_grid(src[i]);
	
	nout=floor(dom1.tout_s/dt);
	ninc=floor( (dom1.tout_e-dom1.tout_s)/((dom1.Nout-1)*dt) );
	printf("ninc=%d, Nout=%d, %lf dT_out=%lf dt=%lf\n",ninc,dom1.Nout,ninc*(dom1.Nout-1)*dt,ninc*dt,dt);
	if(ninc < 1) ninc=1;
//		SETUP RECEIVER
	nrec=set_recs(fgeom,dom1,&rec,Nt,dt);
	printf("nrec=%d\n",nrec);

		fprintf(fdat,"# nsrc, nrec: number of sources & recs)\n");
		fprintf(fdat,"%d %d\n",nsrc,nrec);
		fprintf(fdat,"# ncrk: number of cracks\n");
		fprintf(fdat,"%d\n",ncrk);
		fprintf(fdat,"# rho: mass density\n");
		fprintf(fdat,"%lf\n",dom1.rho);
		fprintf(fdat,"# cL, cT: phase velocities\n");
		fprintf(fdat,"%lf %lf\n",dom1.cL,dom1.cT);
		fprintf(fdat,"# CFL: Counrant Number (x,y)\n");
		fprintf(fdat,"%lf %lf\n",dom1.cfl0,dom1.cfl1);
		fprintf(fdat,"# output times: ts, te, Nout\n");
		fprintf(fdat,"%lf %lf %d\n",dom1.tout_s,dom1.tout_e,dom1.Nout);
		fflush(fdat);

	int ndamp,ndinc,idsum=1;
	int ndout=10;	// Number of Anchoring(entry/exit) points M'
	ndinc=(Nt-1)/ndout;
	if((Nt-1)%ndout!=0){
		printf("Warning !! ndout(=%d) is not a factor of Nt(=%d)\n",ndout,Nt);
	};

	ndamp=ndinc;
	char dir_rst[128]="Rst";
	char dir_out[128]=".";
	int irst0,irst1,it0,it1;
	//int Kd=100;
	int Kd=250; // 2022 06 20
	int incK=(Nt-1)/ndout/Kd;

	if(dom1.irst==0) mkdir(dir_rst,0777);	
	it0=0;
	it1=Nt-1;
	isum=0;
	int j=0;	// loop counter
	if(dom1.irst==1){
		j++;
		isum=1;
		irst0=int(dom1.tr1/(ndinc*dt));
		irst1=int(dom1.tr2/(ndinc*dt));
		it0=irst0*ndinc+1;
		it1=irst1*ndinc;
		printf("Trst1=%lf, irst=%d, it0=%d \n",it0*dt,irst0,it0);
		printf("Trst2=%lf, irst=%d, it1=%d \n",it1*dt,irst1,it1);
		if(it1>=Nt) it1=Nt-1;
		fld1.import(irst0);	// load field data to restart

		ninc=incK;
		nout=incK;
		sprintf(dir_out,"%s",dir_rst);

		sprintf(fname,"%s/v%d.out",dir_out,0);
		printf("%s\n",fname);
		fld1.snap_out(0,fname,(it0-1)*dt,dom1);

		sprintf(fname,"%s/s%d.out",dir_out,0);
		printf("%s\n",fname);
		fld1.snap_out(2,fname,it0*dt,dom1);
	};

	printf("ninc(inck)=%d, nout(init)=%d\n",ninc,nout);
	printf("it0=%d, it1=%d\n",it0,it1);


//		TIME-STEPPING
	//for(int it=0;it<Nt;it++){
	for(int it=it0; it<=it1; it++){

		fld1.v2s(dom1,dt);
		
		for(i=0;i<ncrk;i++) crk[i].v2s(fld1,dom1,dt);
		for(i=0;i<nsrc;i++){
			fld1.apply_src(it,src[i],dom1);
		}
		fld1.s2v(dom1,dt);
		for(i=0;i<ncrk;i++) crk[i].s2v(fld1,dom1,dt);
		for(i=0;i<ncrk;i++) crk[i].out_cod(fCOD[i],it,Nt,dt,dom1);
		for(i=0;i<ncrk;i++) crk[i].out_csd(fCSD[i],it,Nt,dt,dom1,fld1);
		//if((it-it0)==nout){
		if(j==nout){
			nout+=ninc;
			if((it-0.5)*dt <= dom1.tout_e){
				sprintf(fname,"%s/v%d.out",dir_out,isum);
				printf("%s\n",fname);
				fld1.snap_out(0,fname,it*dt,dom1);

				//sprintf(fname,"w%d.out",isum);
				//printf("%s\n",fname);
				//fld1.snap_out_del(fname,it*dt,dom1);

				sprintf(fname,"%s/s%d.out",dir_out,isum);
				printf("%s\n",fname);
				fld1.snap_out(2,fname,it*dt,dom1);
				isum++;
			}
		}
		for(i=0;i<nrec;i++) fld1.out(it,rec[i],dom1);

		if(dom1.irst==0){ // if fresh start 
		if(it==ndamp){
			printf("it=%d, %d, t=%lf\n",it,ndinc,it*dt);
			fld1.damp(dir_rst,idsum);
			ndamp+=ndinc;
			idsum++;
		}
		}
		j++;
	}
	
	sprintf(fname,"%s/nfile.out",dir_out);
	//fout=fopen("nfile.out","w");	
	fout=fopen(fname,"w");	
	fprintf(fout,"# number of output files (field snapshot)\n");
	fprintf(fout,"%d\n",isum);
	fprintf(fout,"%s\n",fgeom);
	fprintf(fout,"%s\n",fsrc);
	fprintf(fout,"%s\n",ficon);
	fprintf(fout,"# Data in this file is read when executing gsaft.exe");
	fclose(fout);

	return 0;
}

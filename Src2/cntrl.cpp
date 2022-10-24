#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdm2d.h"

void show_msg(char *fn){
	printf("Cannot find file %s \n",fn);
	printf(" --> abort process ...\n");
	exit(-1);

};


void CNTRL::setup_domain(char *fname){
	int ndim=2;
	int i,j,k;

	double Ya[2],Yb[2];
	FILE *fp;
	char cbff[128],fkcell[128];

//		Read Input Domain Data 
	fp=fopen(fname,"r");	
	if(fp==NULL) show_msg(fname);

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fkcell);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf",Xa,Xa+1);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];

	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	int nchar_min=3;
	if(strlen(fkcell) > nchar_min){
		printf("read geometry from %s\n",fkcell);
		// read xa,xb,Ndiv & kcell from fkcell
//		Dom2D::load_kcell(fkcell);
	};	

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Ha,Ha+1);	// PML thickness
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Hb,Hb+1); // PML thickness

	dx[0]=Wd[0]/Ndiv[0];
	dx[1]=Wd[1]/Ndiv[1];

	// domain extention 
	for(i=0; i<ndim; i++){
		NHa[i]=ceil(Ha[i]/dx[i]);
		NHb[i]=ceil(Hb[i]/dx[i]);
		Ha[i]=dx[i]*NHa[i];
		Hb[i]=dx[i]*NHb[i];

		Xa[i]-=Ha[i];
		Wd[i]=Wd[i]+Ha[i]+Hb[i];
		Ndiv[i]=Ndiv[i]+NHa[i]+NHb[i];
	}

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf\n",&rho,&cL,&cT);

	amu=rho*cT*cT;
	almb=rho*(cL*cL-2.*cT*cT);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Ya,Ya+1);
	fscanf(fp,"%lf %lf\n",Yb,Yb+1);

	fclose(fp);

//		Setup Domain Object
	dm.setup(Xa,Wd,dx);
	dm.cT=cT; dm.cL=cL;
	dm.rho=rho; 
	dm.amu=amu; dm.almb=almb;

	dm.init(Ndiv);
	dm.NHa=NHa;
	dm.NHb=NHb;

	double gmm=1.e-05;	// expected decay 
	dm.Ha=Ha;
	dm.Hb=Hb;
	dm.PML_setup(gmm);

	FILE *ftmp=fopen("pml_val.out","w");
	double tmp;
	for(i=0; i<Ndiv[0]; i++){
		tmp=Xa[0]+dx[0]*(i+0.5);
		fprintf(ftmp,"%lf, %lf\n",tmp,dm.PML_dcy(0,tmp));
	}
	fprintf(ftmp,"\n");
	for(i=0; i<Ndiv[1]; i++){
		tmp=Xa[1]+dx[1]*(i+0.5);
		fprintf(ftmp,"%lf, %lf\n",tmp,dm.PML_dcy(1,tmp));
	}
	fclose(ftmp);
	printf(" -->pml_val.out\n");

	dm.perfo_ellip(fname);
	dm.slit(fname);
	dm.topography(fname);
	dm.out_kcell();		// write kcell data 
	printf(" -->kcell.dat\n");
	dm.fwrite();	// write domain setting

//		Setup Field Objects
	s11.init(Ndiv,0);	
	s22.init(Ndiv,0);	
	v1.init(Ndiv,1);	
	v2.init(Ndiv,2);	
	s12.init(Ndiv,3);

	s11.setup(Xa,Wd,dx);
	s22.setup(Xa,Wd,dx);
	v1.setup(Xa,Wd,dx);
	v2.setup(Xa,Wd,dx);
	s12.setup(Xa,Wd,dx);

	s11.gen_indx0(dm.kcell);
	s22.gen_indx0(dm.kcell);
	v1.gen_indx1(dm.kcell);
	v2.gen_indx2(dm.kcell);
	s12.gen_indx3(dm.kcell);

	char fnf[128]="field_setting.out";
	char md[6]="w",name[6];
	sprintf(name,"s11"); s11.fwrite_prms(fnf,md,name);
	sprintf(md,"a");
	sprintf(name,"v1"); v1.fwrite_prms(fnf,md,name);
	sprintf(name,"v2"); v2.fwrite_prms(fnf,md,name);
	sprintf(name,"s12"); s12.fwrite_prms(fnf,md,name);
	printf(" -->field_setting.out\n");
	sprintf(fnf,"v1bnd.dat"); 
	v1.fwrite_bnd(fnf);
	sprintf(fnf,"v2bnd.dat"); 
	v2.fwrite_bnd(fnf);
	sprintf(fnf,"s12bnd.dat"); 
	s12.fwrite_bnd(fnf);

	printf(" | -->v1bnd.dat\n");
	printf(" | -->v2bnd.dat\n");
	printf(" | -->s12bnd.dat\n");

/*		following needed updating if kcell read from a file 
	ptmp=(long int *)malloc(sizeof(long int)*Ndiv[0]*Ndiv[1]);
	kcell=(long int **)malloc(sizeof(long int*)*Ndiv[0]);
	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}

	if(strlen(fkcell) > nchar_min){
		fp=fopen(fkcell,"r");
		for(i=0;i<5;i++){
		       	fgets(cbff,128,fp);
		}
		for(i=0;i<ndiv[0];i++){
		for(j=0;j<ndiv[1];j++){
			fscanf(fp,"%d\n",&k);
			kcell[i+nwa[0]][j+nwa[1]]=k;
		}
		}
		fclose(fp);
	}
*/
};
bool CNTRL::out_time(int it){
	if(it==iout){
		printf("it/Nt=%d/%d (Nout=%d)\n",it,Nt,iout);
		iout+=Ninc;
		return(true);
	}
	return(false);
};

void CNTRL::time_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	if(fp==NULL) show_msg(fname);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %d\n",&Tf, &Nt);
	dt=Tf/(Nt-1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf,%d\n",&tout_s,&tout_e,&Nout);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&irst);
	if(irst==1){
		fgets(cbff,128,fp);
		fscanf(fp,"%lf, %lf\n",&tr1,&tr2);
		printf("(tr1,tr2)=(%lf, %lf)\n",tr1,tr2);
	}

	printf(" CFL=%lf\n",CNTRL::CFL());

	Ninc=(tout_e-tout_s)/((Nout-1)*dt);
	if(Ninc<1) Ninc=1;
	iout=floor(tout_s/dt);
	if(iout==0) iout+=Ninc;
	iout0=iout;

	/*
		If Necessary READ I.C. Here
	*/

	fclose(fp);

	fp=fopen("time_setting.out","w");
	fprintf(fp,"irst(restart 0:No, 1:Yes)=%d\n",irst);
	fprintf(fp,"tlim=[ %lf, %lf]\n",0.0,Tf);
	fprintf(fp,"  Nt=%d\n",Nt);
	fprintf(fp,"  dt=%lf\n",dt);
	fprintf(fp," CFL=%lf\n",CNTRL::CFL());
	fprintf(fp,"Nout=%d\n",Nout);
	fprintf(fp,"Ninc=%d,(tinc=%lf)\n",Ninc,Ninc*dt);
	fprintf(fp,"iout_start=%d,(tout_start=%lf)\n",iout0, iout0*dt);
	fclose(fp);
	printf(" -->time_setting.out\n");
	
};
double CNTRL::CFL(){
	double dh=dx[0];
	double Crt;
	if(dh >dx[1]) dh=dx[1];
	Crt=cL*dt/dh*sqrt(2.0);
	if(Crt>1.0){
		printf(" stability condition is not satisfied (CFL=%lf)!!\n --> abort proces\n",Crt);
		exit(-1);
	};
	return(Crt);
};

void CNTRL::wvfm_setting(){
	char fname[128];

	int i,j;
	wvs=(Wv1D *)malloc(sizeof(Wv1D)*nwv);
	for(i=0;i<nwv;i++){
		sprintf(fname,"inwv%d.dat",i);
		wvs[i].gen_wv(fname);
	}

	FILE *fp=fopen("inwvs.out","w");
	fprintf(fp,"# nwv, dt, Nt\n");
	fprintf(fp,"%d, %lf, %d\n",nwv,wvs[0].dt, wvs[0].Nt);
	for(i=0;i<nwv;i++){
		fprintf(fp,"#iwv=%d\n",i);
	for(j=0;j<wvs[i].Nt;j++){
		fprintf(fp,"%lf\n",wvs[i].amp[j]);
	}
	}
	fclose(fp);
	printf(" -->inwvs.out\n");
	
	//char fnout[128];
	//wv.FFT(1);
	//sprintf(fnout,"awvt.out");
	//wv.out_amp(fnout,',');
	//awv.Butterworth(16.0,10.0);
	//sprintf(fnout,"awvw.out");
	//wv.out_Amp(fnout,0);
};

int CNTRL::src_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	double xy,wdt,w2;
	int ityp;	// 1: v1-source, 2: v2-source
	int nml;	// 1,-1 (normal vector, nx/ny)
	int iwv;	// waveform No.
	double xsrc,ysrc,ain;
	int i0,i,j,i1,i2,j1,j2,ng;
	int isrc,jsrc,isum;

	char fout[128]="tr_elems.out";
	char mode[6]="w";

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsrc);
	srcs=(TRNSDCR *)malloc(sizeof(TRNSDCR)*nsrc);
	fgets(cbff,128,fp);
	nwv=0;
	for(int k=0; k<nsrc; k++){	// k-th source element
		fscanf(fp,"%d, %lf, %lf, %d, %d, %lf\n",&ityp, &xy, &wdt, &nml, &iwv, &ain);
		if(iwv>nwv) nwv=iwv;
		srcs[k].ID=k;
		srcs[k].nml=nml;
		srcs[k].type=ityp;
		srcs[k].iwv=iwv;
		srcs[k].Xa=Xa;
		srcs[k].dx=dx;
		srcs[k].set_inc_ang(ain,cT);
		w2=wdt*0.5;
		isum=0;
		switch(ityp){
		case 1:	// v1-source
			ng=ceil(wdt/dx[1]);
			if(ng==0) ng=1;
			if(ng%2==0){
				i0=int((xy-Xa[1])/dx[1]+0.5);
				xy=Xa[1]+i0*dx[1];
			}else{
				i0=int((xy-Xa[1])/dx[1]);
				xy=Xa[1]+(i0+0.5)*dx[1];
			};
			j1=int((xy-w2-Xa[1])/dx[1]);
			srcs[k].mem_alloc(ng);
			for(j=0; j<ng; j++){
				jsrc=j+j1;
				if(jsrc<0) continue;
				if(jsrc>=Ndiv[1]) break;
				ysrc=Xa[1]+(jsrc+0.5)*dx[1];
				isrc=dm.find_q1bnd(nml,ysrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		case 2: // v2-source
			ng=ceil(wdt/dx[0]);
			if(ng%2==0){
				i0=int((xy-Xa[0])/dx[0]+0.5);
				xy=Xa[0]+i0*dx[0];
			}else{
				i0=int((xy-Xa[0])/dx[0]);
				xy=Xa[0]+(i0+0.5)*dx[0];
			};
			i1=int((xy-w2-Xa[0])/dx[0]);
			if(ng==0) ng=1;
			srcs[k].mem_alloc(ng);
			for(i=0; i<ng; i++){
				isrc=i+i1;
				if(isrc<0) continue;
				if(isrc>=Ndiv[0]) break;
				xsrc=Xa[0]+(0.5+isrc)*dx[0];
				jsrc=dm.find_q2bnd(nml,xsrc);
				srcs[k].isrc[isum]=isrc;
				srcs[k].jsrc[isum]=jsrc;
				srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				isum++;
			}
			break;
		};
		ng=isum;
		//srcs[k].print();
		srcs[k].fwrite_setting(fout,mode,k);
		sprintf(mode,"a");
		srcs[k].set_center();
		srcs[k].init_bwv(Nt,dt);
	}
	nwv++;
	fclose(fp);
	printf(" -->tr_elems.out\n");
	return(nwv);
};
void CNTRL::array_setting(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];	
	if(fp==NULL) show_msg(fname);
	int nele,nmeas;
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nele);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nmeas);
	if(nele != nsrc){
		printf(" Error: nele(=%d) must be equal to nsrc(=%d) !\n",nele,nsrc);
		exit(-1);
	};
	ary.init(nsrc,nmeas);
	int i,j,k=0;
	for(j=0;j<nmeas;j++){ 
		fgets(cbff,128,fp);
		for(i=0;i<nsrc;i++){
			fscanf(fp,"%d, %lf, %lf\n",ary.actv+k, ary.a0+k, ary.tdly+k);
			k++;
		};
	};
	round=0;
	//ary.print();
	fclose(fp);
};
/*
void CNTRL::capture(int jt){
	int i,j,type;
	for(j=0;j<nsrc;j++) srcs[j].record(jt,v3.F);
};
void CNTRL::q2v(int itime){
	int i,j,k,l;
	int nx,ny;
	double **F1,**F2;
	double dFx, dFy;
	double beta,gmm, *xcod;

	nx=v3.Ng[0];
	ny=v3.Ng[1];
	F1=q1.F;
	F2=q2.F;

	for(k=0; k<v3.Nin; k++){
		l=v3.kint[k];
		v3.l2ij(l,&i,&j);
		xcod=dm.ij2xf(i,j,0);

		beta=0.5*dt*dm.PML_dcy(0,xcod[0]);
		gmm=dt/(rho*dx[0]);
		dFx=(F1[i+1][j]-F1[i][j])*gmm;
		v3x.F[i][j]=((1.-beta)*v3x.F[i][j]+dFx)/(1.+beta);
		//printf("dcy=%lf, xcod=%lf\n",(1.-beta)/(1.+beta),xcod[0]);

		beta=0.5*dt*dm.PML_dcy(1,xcod[1]);
		gmm=dt/(rho*dx[1]);
		dFy=(F2[i][j+1]-F2[i][j])*gmm;
		v3y.F[i][j]=((1.-beta)*v3y.F[i][j]+dFy)/(1.+beta);

		//v3.F[i][j]+=(dFx+dFy);
		v3.F[i][j]=v3x.F[i][j]+v3y.F[i][j];
	}
};
void CNTRL::v2q(int itime){
	int i,j,k,l;
	int nx,ny;
	double **F;
	double dFx,dFy,alph;

	F=v3.F;
	nx=q1.Ng[0]; ny=q1.Ng[1];

	for(k=0; k<q1.Nin; k++){
		l=q1.kint[k];
		q1.l2ij(l,&i,&j);
		alph=amu*dt/dx[0];
		//dFx=(F[i][j]-F[i-1][j])/dx[0]*dt/rho;
		dFx=alph*(F[i][j]-F[i-1][j]);
		q1.F[i][j]+=dFx;
	};

	for(k=0; k<q2.Nin; k++){
		l=q2.kint[k];
		q2.l2ij(l,&i,&j);
		alph=amu*dt/dx[1];
		//dFy=(F[i][j]-F[i][j-1])/dx[1]*dt/rho;
		dFy=alph*(F[i][j]-F[i][j-1]);
		q2.F[i][j]+=dFy;
	};

	//SOURCE src;
	TRNSDCR src;

	double bvl,tdly,tdly0,a0;
	int sgn,sft,iwv,i0;
	i0=nsrc*round;
	for(int isrc=0;isrc<nsrc;isrc++){
		if(ary.actv[isrc+i0]==0) continue;

		src=srcs[isrc]; // SOURCE
		iwv=src.iwv;

		a0=ary.a0[isrc+i0];
		tdly0=ary.tdly[isrc+i0];
		//bvl=wvs[iwv].amp[itime];
		sft=1;
		sgn=src.nml;
		if(sgn ==-1) sft=0;

		switch(src.type){
		case 1:	// v1 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				if(src.Ctd<0.0) tdly=(k-src.ng+1)*src.Ctd+tdly0;	// negative incident angle
				bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				//dFx=sgn*(bvl-F[i-sft][j])/dx[0]*dt/rho;
				//q1.F[i][j]+=dFx;
				q1.F[i][j]=bvl*sgn;
			}
			break;
		case 2: // v2 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				if(src.Ctd<0.0) tdly=(k-src.ng+1)*src.Ctd+tdly0;	// negative incident angle
				bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				//dFy=sgn*(bvl-F[i][j-sft])/dx[1]*dt/rho;
				//q2.F[i][j]+=dFy;
				q2.F[i][j]=bvl*sgn;
			}
			break;
		};
	}

};

void CNTRL::clear(){
	v3.clear();
	q1.clear();
	q2.clear();
	v3x.clear();
	v3y.clear();
	for(int i=0;i<nsrc;i++) srcs[i].clear();
	iout=iout0;
};
void CNTRL::fwrite_ary(){
	int j,k;
	char fname[128];
	FILE *fp;

	//	ARRAY WAVEFORMS 
	sprintf(fname,"T%d/ary.out",round);
	fp=fopen(fname,"w");
	fprintf(fp,"# nele\n");
	fprintf(fp,"%d\n",nsrc);
	fprintf(fp,"# Nt, dt\n");
	fprintf(fp,"%d, %lf\n",Nt,dt);
	for(j=0;j<nsrc;j++){
		fprintf(fp,"# e=%d, %lf, %lf\n",j,srcs[j].x0,srcs[j].y0);
		for(k=0;k<Nt;k++){
			fprintf(fp,"%lf\n",srcs[j].amp_synth(k));
		}
	};
	fclose(fp);
};
void CNTRL::snapshot(int n_meas, int isum, int it){
	v3.fwrite_trim(n_meas,isum,NHa,NHb,it*dt);
};




*/

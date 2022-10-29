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
				//srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				srcs[k].ksrc[isum]=isrc*Ndiv[1]+jsrc;
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
				//srcs[k].ksrc[isum]=isrc+(Ndiv[1]+1)*jsrc;
				srcs[k].ksrc[isum]=isrc*(Ndiv[1]+1)+jsrc;
				isum++;
			}
			break;
		};
		ng=isum;
		srcs[k].ng=isum;
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
void CNTRL::mark_src_grid(){
	int i,j,loc;
	int type;
	FIELD V;

	for(i=0; i<nsrc; i++){
		type=srcs[i].type;
		if(type==1){
			V=v1;
		}else if(type ==2){
			V=v2;
		}else{
			printf("invalid source type provided to mark_src_grid\n");
			exit(-1);
		}
		for(j=0;j<srcs[i].ng;j++){
			loc=CNTRL::find_src_index(V,srcs[i].ksrc[j]);
			V.ksrc[loc]=true;
			//printf("j=%d, loc=%d\n",j,loc);
		};
	};

	/*
	int isum=0;
	for(i=0; i<v2.Nbnd; i++){
		if(v2.ksrc[i]){
			printf("i=%d, %d\n",i,v2.ksrc[i]);
			isum++;
		}
	};
	printf("isum=%d\n",isum);
	*/
};
int CNTRL::find_src_index(FIELD V, int k){
	int il,im,ih;
	
	il=0;
	ih=V.Nbnd-1;
	if(abs(V.kbnd[il])==k) return(il);
	if(abs(V.kbnd[ih])==k) return(ih);
	im=int((il+ih)*0.5);

	while(abs(V.kbnd[im]) !=k){
		if(ih-il<=1){
			printf("Can't find source grid index %d in kbnd[*]\n",k);
			exit(-1);
		}
		if(abs(V.kbnd[im])>k){
			ih=im;
		}else{
			il=im;
		};
		im=int((il+ih)*0.5);

	}
	//printf("loc=%d, k=%d,%d\n",im, V.kbnd[im],k);
	return(im);
};


//-----------------------------------------------------
//
//		TIME-STEPPING OPERATION 1
//		  ( Velocity --> Stress)
//
//void PML :: v2s(double amu, double almb, double dt){
void CNTRL:: v2s(int it){
	int i,j,ib,jb;
	double dmpX,dmpY,alph,beta,xx,yy;
	double dv1,dv2;
	double dtr=1.0/dt;
	double a2m=almb+2.0*amu;
	double s110=0.0, s220=0.0;

//		Normal Stress: s11, s22
	//for(i=0; i<nx[0]-1 ;i++){
	for(i=0; i<Ndiv[0] ;i++){
		//dmpX=dmpx[i];
		xx=Xa[0]+dx[0]*(i+0.5);
		dmpX=dm.PML_dcy(0,xx); 
		alph=(dtr-dmpX)/(dtr+dmpX);
	//for(j=0; j<nx[1]-1 ;j++){
	for(j=0; j<Ndiv[1] ;j++){
		//dmpY=dmpy[j];
		yy=Xa[1]+dx[1]*(j+0.5);
		dmpY=dm.PML_dcy(1,yy);
		beta=(dtr-dmpY)/(dtr+dmpY);

		//dv1=(v1[i+1][j]-v1[i][j])/(dtr+dmpX);
		dv1=(v1.F[i+1][j]-v1.F[i][j])/(dtr+dmpX);
		//dv2=(v2[i][j+1]-v2[i][j])/(dtr+dmpY);
		dv2=(v2.F[i][j+1]-v2.F[i][j])/(dtr+dmpY);

		// p --> x
		// v --> y

		//s11p[i][j]*=alph;
		s11.Fx[i][j]*=alph;
		//s11v[i][j]*=beta;
		s11.Fy[i][j]*=beta;
		//s11p[i][j]+=(dv1*a2m/dx[0]);
		s11.Fx[i][j]+=(dv1*a2m/dx[0]);
		//s11p[i][j]+=(s110*2.*dmpX/(dtr+dmpX));	// <--------------
		s11.Fx[i][j]+=(s110*2.*dmpX/(dtr+dmpX));	// <--------------
		//s11v[i][j]+=(dv2*almb/dx[1]);
		s11.Fy[i][j]+=(dv2*almb/dx[1]);
		//s11[i][j]=s11p[i][j]+s11v[i][j];
		s11.F[i][j]=s11.Fx[i][j]+s11.Fy[i][j];

		//s22p[i][j]*=alph;
		s22.Fx[i][j]*=alph;
		//s22v[i][j]*=beta;
		s22.Fy[i][j]*=beta;
		//s22p[i][j]+=(dv1*almb/dx[0]);
		s22.Fx[i][j]+=(dv1*almb/dx[0]);
		//s22p[i][j]+=(s220*2.*dmpX/(dtr+dmpX)); // <--------------
		s22.Fx[i][j]+=(s220*2.*dmpX/(dtr+dmpX)); // <--------------
		//s22v[i][j]+=(dv2*a2m/dx[1]);
		s22.Fy[i][j]+=(dv2*a2m/dx[1]);
		//s22[i][j]=s22p[i][j]+s22v[i][j];
		s22.F[i][j]=s22.Fx[i][j]+s22.Fy[i][j];
	}
	}

//		Shear Stress: s12

	//for(i=1; i<nx[0]-1 ;i++){
	for(i=1; i<Ndiv[0]; i++){
		//dmpX=.5*(dmpx[i]+dmpx[i-1]);
		xx=Xa[0]+dx[0]*i;
		dmpX=dm.PML_dcy(0,xx); 
		alph=(dtr-dmpX)/(dtr+dmpX);
	//for(j=1; j<nx[1]-1 ;j++){
	for(j=1; j<Ndiv[1]; j++){
		yy=Xa[1]+dx[1]*j;
		//dmpY=.5*(dmpy[j]+dmpy[j-1]);
		dmpY=dm.PML_dcy(1,yy); 
		beta=(dtr-dmpY)/(dtr+dmpY);
		
		//dv2=(v2[i][j]-v2[i-1][j])*amu/(dtr+dmpX);
		dv2=(v2.F[i][j]-v2.F[i-1][j])*amu/(dtr+dmpX);
		//dv1=(v1[i][j]-v1[i][j-1])*amu/(dtr+dmpY);
		dv1=(v1.F[i][j]-v1.F[i][j-1])*amu/(dtr+dmpY);

		//s12p[i][j]*=alph;
		s12.Fx[i][j]*=alph;
		//s12v[i][j]*=beta;
		s12.Fy[i][j]*=beta;
		//s12p[i][j]+=(dv2/dx[0]);
		s12.Fx[i][j]+=(dv2/dx[0]);
		//s12v[i][j]+=(dv1/dx[1]);
		s12.Fy[i][j]+=(dv1/dx[1]);
		//s12[i][j]=s12p[i][j]+s12v[i][j];
		s12.F[i][j]=s12.Fx[i][j]+s12.Fy[i][j];
	}
	}
//		Stress Free B.C.	
	int k,l;
	for(k=0; k<s12.Nbnd; k++){
		l=s12.kbnd[k];
		s12.l2ij(abs(l),&i,&j);
		s12.F[i][j]=0.0;
		s12.Fx[i][j]=0.0;
		s12.Fy[i][j]=0.0;
	};
/*

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
*/
};

//-----------------------------------------------------
//
//		TIME-STEPPING OPERATION 2
//		  ( Stress --> Velocity )
//
//void PML :: s2v(double rho, double dt){
void CNTRL :: s2v(int itime){
	int i,j,ib,jb;
	double dmpX,dmpY,alph,beta,dS11,dS22,dS12;
	double dtr=1.0/dt,sgn;
	double rdx=dx[0]*rho, rdy=dx[1]*rho;
	double xx,yy;

//	################ Velocity v1 ##################
	int k,l;
	for(k=0;k<v1.Nin; k++){
		l=v1.kint[k];
		v1.l2ij(abs(l),&i,&j);
	//for(i=1; i<nx[0]-1 ;i++){
		xx=Xa[0]+dx[0]*i;
		yy=Xa[1]+dx[1]*(j+0.5);

		//dmpX=.5*(dmpx[i]+dmpx[i-1]);
		dmpX=dm.PML_dcy(0,xx); 
		alph=(dtr-dmpX)/(dtr+dmpX);
	//for(j=0; j<nx[1]-1; j++){
		//dmpY=dmpy[j];
		dmpY=dm.PML_dcy(1,yy); 
		beta=(dtr-dmpY)/(dtr+dmpY);

		//dS11=(s11[i][j]-s11[i-1][j])/((dtr+dmpX)*rdx);
		dS11=(s11.F[i][j]-s11.F[i-1][j])/((dtr+dmpX)*rdx);
		//dS12=(s12[i][j+1]-s12[i][j])/((dtr+dmpY)*rdy);
		dS12=(s12.F[i][j+1]-s12.F[i][j])/((dtr+dmpY)*rdy);

		//v1p[i][j]*=alph;
		v1.Fx[i][j]*=alph;
		//v1v[i][j]*=beta;
		v1.Fy[i][j]*=beta;
		//v1p[i][j]+=dS11;
		v1.Fx[i][j]+=dS11;
		//v1v[i][j]+=dS12;
		v1.Fy[i][j]+=dS12;
		//v1[i][j]=v1p[i][j]+v1v[i][j];
		v1.F[i][j]=v1.Fx[i][j]+v1.Fy[i][j];
	//}
	}
//		Boundary Nodes (stress-free B.C.)

	//i=0;
	//sgn=1.0;
	//dmpX=dmpx[i];
	
	double s110=0.0;
	for(k=0; k<v1.Nbnd; k++){
		if(v1.ksrc[k]) continue;
		l=v1.kbnd[k];
		v1.l2ij(abs(l),&i,&j);
		xx=Xa[0]+dx[0]*i;
		yy=Xa[1]+dx[1]*(j+0.5);
		sgn=1.0;
		if(l<0) sgn=-1.0;
		dmpX=dm.PML_dcy(0,xx);
		dmpY=dm.PML_dcy(1,yy);
//	for(ib=0;ib<2;ib++){
		//if(ib==1){
		//	i=nx[0]-1;
		//	dmpX=dmpx[i-1];
		//	sgn=-1.0;
		//}
		alph=(dtr-dmpX)/(dtr+dmpX);
//	for(j=0;j<nx[1]-1;j++){
		//dmpY=dmpy[j];
		beta=(dtr-dmpY)/(dtr+dmpY);
		//dS11=sgn*2.*(s11[i-ib][j]-s110)/((dtr+dmpX)*rdx);
		dS11=sgn*2.*(s11.F[i-ib][j]-s110)/((dtr+dmpX)*rdx);
		//dS12=(s12[i][j+1]-s12[i][j])/((dtr+dmpY)*rdy);
		dS12=(s12.F[i][j+1]-s12.F[i][j])/((dtr+dmpY)*rdy);
		//v1p[i][j]*=alph;
		v1.Fx[i][j]*=alph;
		//v1v[i][j]*=beta;
		v1.Fy[i][j]*=beta;
		//v1p[i][j]+=dS11;
		v1.Fx[i][j]+=dS11;
		//v1v[i][j]+=dS12;
		v1.Fy[i][j]+=dS12;
		//v1[i][j]=v1p[i][j]+v1v[i][j];
		v1.F[i][j]=v1.Fx[i][j]+v1.Fy[i][j];
//	}	
	}

//	################ Velocity v2 ##################
	for(k=0; k<v2.Nin; k++){
		l=v2.kint[k];
		v2.l2ij(abs(l),&i,&j);
		xx=Xa[0]+dx[0]*(i+0.5);
		yy=Xa[1]+dx[1]*j;
	//for(i=0; i<nx[0]-1; i++){
		//dmpX=dmpx[i];
		dmpX=dm.PML_dcy(0,xx); 
		alph=(dtr-dmpX)/(dtr+dmpX);
	//for(j=1; j<nx[1]-1; j++){
		//dmpY=.5*(dmpy[j]+dmpy[j-1]);
		dmpY=dm.PML_dcy(1,yy); 
		beta=(dtr-dmpY)/(dtr+dmpY);

		dS12=(s12.F[i+1][j]-s12.F[i][j])/((dtr+dmpX)*rdx);
		dS22=(s22.F[i][j]-s22.F[i][j-1])/((dtr+dmpY)*rdy);

		//v2p[i][j]*=alph;
		v2.Fx[i][j]*=alph;
		//v2v[i][j]*=beta;
		v2.Fy[i][j]*=beta;
		//v2p[i][j]+=dS12;
		v2.Fx[i][j]+=dS12;
		//v2v[i][j]+=dS22;
		v2.Fy[i][j]+=dS22;
		//v2[i][j]=v2p[i][j]+v2v[i][j];
		v2.F[i][j]=v2.Fx[i][j]+v2.Fy[i][j];
	//}
	}

//		Boundary Nodes (stress-free B.C.)
	//j=0;
	//sgn=1.0;
	//dmpY=dmpy[j];

	double s220=0.0;
	for(k=0; k<v2.Nbnd; k++){
		if(v2.ksrc[k]) continue;
		l=v2.kbnd[k];
		v2.l2ij(abs(l),&i,&j);
		sgn=1.0;
		if(l<0) sgn=-1.0;
		xx=Xa[0]+dx[0]*(i+0.5);
		yy=Xa[1]+dx[1]*j;
		dmpX=dm.PML_dcy(0,xx);
		dmpY=dm.PML_dcy(1,yy);
//	for(jb=0; jb<2; jb++){
	//	if(jb==1){
	//		j=nx[1]-1;
	//		sgn=-1.0;
//			dmpY=dmpy[j-1];
//		}	
		beta=(dtr-dmpY)/(dtr+dmpY);
//	for(i=0; i<nx[0]-1; i++){
//		dmpX=dmpx[i];
		alph=(dtr-dmpX)/(dtr+dmpX);

		//dS12=(s12[i+1][j]-s12[i][j])/((dtr+dmpX)*rdx);
		dS12=(s12.F[i+1][j]-s12.F[i][j])/((dtr+dmpX)*rdx);
		//dS22=2.0*sgn*(s22[i][j-jb]-s220)/((dtr+dmpY)*rdy);
		dS22=2.0*sgn*(s22.F[i][j-jb]-s220)/((dtr+dmpY)*rdy);

		//v2p[i][j]*=alph;
		v2.Fx[i][j]*=alph;
		//v2v[i][j]*=beta;
		v2.Fy[i][j]*=beta;
		//v2p[i][j]+=dS12;
		v2.Fx[i][j]+=dS12;
		//v2v[i][j]+=dS22;
		v2.Fy[i][j]+=dS22;
		//v2[i][j]=v2p[i][j]+v2v[i][j];
		v2.F[i][j]=v2.Fx[i][j]+v2.Fy[i][j];
	//}
	}

	TRNSDCR src;
	double bvl,tdly,tdly0,a0;
	int sft,iwv,i0,isgn;
	i0=nsrc*round;
	for(int isrc=0;isrc<nsrc;isrc++){
		if(ary.actv[isrc+i0]==0) continue;
		src=srcs[isrc]; // SOURCE
		iwv=src.iwv;

		//a0=ary.a0[isrc+i0];
		tdly0=ary.tdly[isrc+i0];
		//bvl=wvs[iwv].amp[itime];
		sft=1;
		isgn=src.nml;
		if(isgn ==-1) sft=0;

		switch(src.type){
		case 1:	// v1 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				if(src.Ctd<0.0) tdly=(k-src.ng+1)*src.Ctd+tdly0;	// negative incident angle
				//bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				bvl=wvs[iwv].val(itime*dt-tdly);
				//q1.F[i][j]=bvl*isgn;
			}
			break;
		case 2: // v2 source
			for(k=0; k<src.ng; k++){
				i=src.isrc[k];
				j=src.jsrc[k];
				tdly=src.Ctd*k+tdly0;
				if(src.Ctd<0.0) tdly=(k-src.ng+1)*src.Ctd+tdly0;	// negative incident angle
				//bvl=wvs[iwv].val(itime*dt-tdly)*a0;
				bvl=wvs[iwv].val(itime*dt-tdly);
				//q2.F[i][j]=bvl*isgn;
			}
			break;
		};

	}
};

//----------------------------------------------------
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

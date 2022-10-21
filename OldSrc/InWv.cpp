#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fdm2d.h"

using namespace std;

//-----------------------------------------------------------
InWv::InWv(int N){	// empty constructor
	t1=0.0; t2=1.0;
	dt=0.0; T0=0.1;
	if(N>1) dt=(t2-t1)/(N-1);
	wvtyp=1; 
	nbrst=3;
	Nt=N;
	mem_alloc();	
	
};
InWv::InWv(char *fname){
	int i;
	FILE *fp;
	char ch[128];
	fp=fopen(fname,"r");
	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&wvtyp);	

	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&Nt);	

	fgets(ch,128,fp);
	fscanf(fp,"%lf %d\n",&T0,&nbrst);	

	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&nsig); // narrowness factor used in Gaussian amp. modulation	

	fgets(ch,128,fp);
	fscanf(fp,"%lf %lf\n",&t1,&t2);	

	dt=(t2-t1)/(Nt-1);

	mem_alloc();
	gen_wvfm();

	double tb=nbrst*T0*0.5;
	if(nsig != 0) Amod(tb,nsig);
	if(wvtyp <0){	// read numerical waveform if wvtyp is negative
		fgets(ch,128,fp);
		for(i=0;i<Nt;i++){
			 fscanf(fp,"%lf\n",amp+i);	
		}
	}
	fclose(fp);
};

void InWv :: disp(){
	cout<<"**** waveform parameters ****\n";
	cout<<"wtyp="<<wvtyp<<'\n';
	cout<<"t1="<<t1<<" "<<"t2="<<t2<<'\n';
	cout<<"dt="<<dt<<'\n';
	cout<<"Nt="<<Nt<<'\n';
	cout<<"nbrst="<<nbrst<<'\n';
	cout<<"T0="<<T0<<'\n';
}

void InWv::mem_alloc(){
	amp=(double *)malloc(sizeof(double)*Nt);
}

void InWv:: gen_wvfm(){
	switch(wvtyp){
	case 1: // sine function
		gen_sin();
		break;
	case 2:	// cosine function
		gen_cos();
		break;
	};
}

void InWv:: gen_sin(){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
		if(tt-t1<=nbrst*T0) amp[i]=sin(omg*tt);
	}
}
void InWv:: gen_cos(){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
//		if(tt-t1<=nbrst*T0) amp[i]=1-cos(omg*tt);
		if(tt-t1<=nbrst*T0) amp[i]=cos(omg*tt);
	}
}
void InWv::out(char *fout){
	FILE *fp=fopen(fout,"w");
	int i;
	for(i=0;i<Nt;i++) fprintf(fp,"%lf %lf\n",t1+dt*i, amp[i]);
//	for(i=0;i<Nt;i++) fprintf(fp,"%lf\n",amp[i]);
	fclose(fp);
}
void InWv::DFTout(char *fout){
	FILE *fp=fopen(fout,"w");
	int i;
	for(i=0;i<Nt;i++) fprintf(fp,"%lf %lf %lf %lf \n",df*i, Amp[i].re,Amp[i].im,Amp[i].Abs());
	fclose(fp);
}

void InWv::Amod(
	double tb,	// mean 
	int nsig)	// narrowness 
{
	int i;
	double arg,sig=tb/nsig;
	for(i=0;i<Nt;i++){
		arg=(i*dt+t1-tb)/sig;
		arg*=arg;
		amp[i]*=exp(-arg*0.5);
	
	}
};
//---------------------------------------------------
void InWv::DFT(){
	int i,j;
	double omg,pi=4.0*atan(1.0),tt;
	fmax=1./dt;
	df=fmax/Nt;	
	Cmplx zi(0.0,1.0);
	double dw=2.0*pi*df;

	Amp=(Cmplx *)malloc(sizeof(Cmplx)*Nt);
	for(i=0;i<Nt;i++){
		Amp[i]=0.0;
		tt=t1+dt*i;
	for(j=0;j<Nt;j++){
		omg=dw*j;
		Amp[i]=Amp[i]+exp(zi*omg*tt)*amp[j];
	}	
	}

}

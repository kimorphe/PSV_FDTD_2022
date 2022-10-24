#define DB 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "wave1d.h"

using namespace std;

//------------------------------------------------------------
Wv1D::Wv1D(){
	//amp=0;
	//time=0;
	mllc=false;
	fft_stat=0;
	//Nt=0;
	//dt=0.0;
	//t1=0.0;t2=1.0;
	//sprintf(data_file,"not_specified");
};
Wv1D::Wv1D(int ndat){
	Nt=ndat;
	amp=(double *)malloc(sizeof(double)*Nt);
	time=(double *)malloc(sizeof(double)*Nt);
	mllc=true;
	fft_stat=0;
	strcpy(data_file,"");
	t1=0.0;
	t2=0.0;
	dt=0.0;
};
Wv1D::Wv1D(char *fname){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
	Wv1D::load(fname);
};
void Wv1D::gen_wv(char *fname){
	awv.gen(fname);
	Nt=awv.Nt;
	dt=awv.dt;
	amp=awv.amp;
	time=(double *)malloc(sizeof(double)*Nt);
	mllc=true;
	for(int i=0;i<Nt;i++) time[i]=dt*i;
	t1=time[0];
	t2=time[Nt-1];
};
double Wv1D::val(double tau){

	if(tau < t1) return(0.0);
	if(tau > t2) return(0.0);

	int i1=int((tau-t1)/dt);
	if(i1==Nt-1) return(amp[i1]);
	int i2=i1+1;
	double xi=(tau-t1)/dt-i1;
	return(amp[i1]*(1.-xi)+amp[i2]*xi);
};
void Wv1D::print_info(){
	printf("------wave data parameters-------\n");
	printf("File: %s\n",data_file);
	printf("(t1,t2)=%lf,%lf\n",t1,t2);
	printf("dt=%lf\n",dt);
	printf("Nt=%d\n",Nt);
	printf("---------------------------------\n");
}
int Wv1D::load(char fname[128]){

	FILE *fp=fopen(fname,"r");
	char cbff[128];

	if(fp==NULL){	// check if sucessfully opened
		printf("File %s not found!\n",fname);
		printf(" --> process terminated.");
		exit(-1);
	}
	Nt=0;
	while(fgets(cbff,128,fp)!=NULL){
		Nt++;	// count number of lines
	}	
	printf("Nt=%d\n",Nt);
	fclose(fp);

	if(!mllc){	// allocate memory
		amp=(double *)malloc(sizeof(double)*Nt);
		time=(double *)malloc(sizeof(double)*Nt);
		mllc=true;
	}

	fp=fopen(fname,"r");
	strcpy(data_file,fname);
	double sum=0.0;
	for(int i=0;i<Nt;i++){
		fscanf(fp,"%lf, %lf\n",time+i,amp+i);
		time[i]*=1.e06;	// to micro sec
		sum+=amp[i];
	};
	fclose(fp);
	sum/=Nt;
	for(int i=0;i<Nt;i++) amp[i]-=sum;
	t1=time[0];
	t2=time[Nt-1];
	dt=time[1]-time[0];
	return(0);
};
int Wv1D::FFT(int isgn){
	int p=ceil(log2(Nt));
	Np=pow(2,p);

	if(fft_stat==0){
		Amp=(complex<double> *)malloc(sizeof(complex<double>)*Np);
		phi=(double *)malloc(sizeof(double)*Np);
	}
	if(isgn==1){
		for(int i=0;i<Nt;i++) Amp[i]=complex<double>(amp[i],0.0);
		for(int i=Nt;i<Np;i++) Amp[i]=complex<double>(0.0,0.0);
	}

	fft(Amp,Np,isgn);
	fft_stat=isgn;

	int i;
	for(i=0;i<Np;i++) phi[i]=arg(Amp[i]); 
	df=1./dt/Np;
	return(Np);
};
void Wv1D::renew_phase(){
	int i;
	for(i=0;i<Np;i++) phi[i]=arg(Amp[i]);
};
void Wv1D::unwrap(double f0){	// unwraping from f=f0
	int i,nf=int(f0/df),nwrap;
	double dphi;
	double phi2,phi1;
	double PI2=8.0*atan(1.0);
	double aPI2=0.9*PI2;

	phi1=phi[nf];
	nwrap=0;
	for(i=nf+1;i<Np;i++){
		phi2=phi[i];
		dphi=phi2-phi1;
		if(dphi>aPI2) nwrap--;
		if(dphi<-aPI2) nwrap++;
		phi[i]+=(PI2*nwrap);
		phi1=phi2;
	};

	phi2=phi[nf];
	nwrap=0;
	for(i=nf-1;i>=0;i--){
		phi1=phi[i];
		dphi=phi2-phi1;
		if(dphi>aPI2) nwrap++;
		if(dphi<-aPI2) nwrap--;
		phi[i]+=(PI2*nwrap);
		phi2=phi1;
	}
};
void Wv1D::out_Amp(char *fn){
	FILE *fp=fopen(fn,"w");
	df=1./dt/Np;
	double freq;
	for(int i=0;i<Np;i++){
		freq=i*df;
		fprintf(fp,"%le %le %le\n",freq,Amp[i].real(),Amp[i].imag());
	}
	fclose(fp);
};
void Wv1D::out_Amp_polar(char *fn){
	FILE *fp=fopen(fn,"w");
	df=1./dt/Np;
	double freq;
	for(int i=0;i<Np;i++){
		freq=i*df;
		//fprintf(fp,"%le %le %le\n",freq,abs(Amp[i]),arg(Amp[i]));
		fprintf(fp,"%le %le %le\n",freq,abs(Amp[i]),phi[i]);
	}
	fclose(fp);
};

double Wv1D::L2(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double sum=0.0;
	for(i=i1;i<=i2;i++) sum+=(amp[i]*amp[i]);
	return(sqrt(sum)/(i2-i1+1));
};
double Wv1D::max(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double amax=fabs(amp[i2]);
	for(i=i1;i<i2;i++){
	       	if(amax < fabs(amp[i])) amax=fabs(amp[i]);
	}
	return(amax);
};
void Wv1D::out_Amp(char *fn,int ofst){
	FILE *fp=fopen(fn,"w");
	int j;
	double xx,x1,dx;
	double df;
	if(fft_stat==1){
		df=1./dt/Np;
		dx=df;
		x1=0.0;
	}else{
		x1=t1;
		dx=dt;
	}
	for(int i=0;i<Np;i++){
		j=(ofst+i)%Np;
		xx=(i-ofst)*dx;
		fprintf(fp,"%le %le %le %le\n",xx,Amp[j].real(),Amp[j].imag(),abs(Amp[j]));
	}
	fclose(fp);
};
void Wv1D::out_amp(char *fn,char dlm){
	FILE *fp=fopen(fn,"w");
	for(int i=0;i<Nt;i++) fprintf(fp,"%le%c %le\n",i*dt+t1,dlm,amp[i]);

	fclose(fp);
};
void Wv1D::Gauss(double tb, double sig){
	double arg;

	int i;
	for(i=0;i<Nt;i++){
		arg=time[i]-tb;
		arg*=arg;
		amp[i]*=exp(-arg*0.5);
	}
};
void Wv1D::Butterworth(double tb, double Tw_6dB){

	int p=4;
	double tt;
	double t0=Tw_6dB*0.5;
	double arg;
	for(int i=0; i<Nt;i++){
		tt=t1+dt*i;
		arg=(tt-tb)/t0;
		arg=pow(arg,p);
		amp[i]/=(1.+arg);
	};
};
double Wv1D::gdelay(){
	if(fft_stat==0) FFT(1);
	tg=(double *)malloc(sizeof(double)*(Nt-1));
	double *phi=(double *)malloc(sizeof(double)*Nt);
	double *tmp=(double *)malloc(sizeof(double)*Nt);

	int i;
	double phi0=arg(Amp[0]),phi1;
	double PI2=8.0*atan(1.0);
	double tol=0.60;
	tol*=PI2;
	FILE *fp=fopen("tg.out","w");

	// phase spectrum
	double cost,sint,zlen;
	for(i=0;i<Nt;i++){
		zlen=abs(Amp[i]);
		cost=real(Amp[i])/zlen;	
		sint=imag(Amp[i])/zlen;	
		tmp[i]=acos(cost);
		if( sint < 0.0) tmp[i]=PI2-tmp[i];
	};

	// unwrapping
	double dphi,ofst=0.0;
	phi[0]=tmp[0];
	for(i=1;i<Nt;i++){
		dphi=tmp[i]-tmp[i-1];
		phi[i]=tmp[i];
		if(dphi>tol){
			ofst-=PI2;
		}
		if(dphi<-tol){
			ofst+=PI2;
		}
		phi[i]+=ofst;
	};

	// group delay 
	double df=1./dt/Np;
	double Cf=PI2*df;
	double f1=0.05, f2=1.0;
	int nf1=int(f1/df), nf2=int(f2/df);
	for(i=0;i<Nt-1;i++){
		tg[i]=(phi[i+1]-phi[i])/Cf;
	};
	fclose(fp);
	double tgb=0.0;
	for(i=nf1;i<=nf2;i++) tgb+=tg[i];
	return(tgb/(nf2-nf1+1));
};
//------------------------------------------------------------
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax){
	wv1.FFT(1);
	wv2.FFT(1);

	Wv1D wv3(wv1.Np);
	wv3.FFT(1);
	double A1=0.0,A2=0.0;
	complex <double > Z1,Z2;
	Z1=complex<double>(0.0,0.0);
	Z2=complex<double>(0.0,0.0);
	wv3=wv1;
	for(int i=0;i<wv1.Np;i++){
		Z1=wv1.Amp[i];
		Z2=conj(wv2.Amp[i]);
		wv3.Amp[i]=Z1*Z2;
		A1+=abs(Z1*Z1);
		A2+=abs(Z2*Z2);
	};
	wv3.FFT(-1);
	A1=sqrt(A1*A2);
	int imax=0;
	double A0=abs(wv3.Amp[0].real()/A1);
	for(int i=0;i<wv1.Np;i++){
		wv3.Amp[i]/=A1;
		if(A0 < abs(wv3.Amp[i].real())){
			A0=abs(wv3.Amp[i].real());
			imax=i;
		}
	}
	double t0=imax*wv1.dt;
	if(imax>wv1.Np/2) t0=-(wv1.Np-imax)*wv1.dt;
	(*tmax)=t0;
	(*Amax)=A0;
	//printf("%d %le %le\n",imax,t0,A0);
	return(wv3);
};
//---------------------------------------------------------
int dft(complex<double> *Amp, int N, int isgn){
	
	int i,j;
	double PI=4.0*atan(1.0);
	complex<double> zi(0.0,1.0);
	complex<double> Wij;
	complex<double> *dat=(complex<double> *)malloc(sizeof(complex<double>)*N);

	for(i=0;i<N;i++) dat[i]=Amp[i];

	zi*=isgn;

//		Discrete Fourier Transform
	for(i=0;i<N;i++){
		Amp[i]=complex<double>(0.0,0.0);
	for(j=0;j<N;j++){
		Wij=exp((2.*PI/N*i*j)*zi);
		Amp[i]+=(dat[j]*Wij);
	}
		if(isgn==1) Amp[i]/=N;
	}

	return( ceil(log2(N)) );	
}; 

int fft(complex<double> *Amp, int N, int isgn){

	int p=ceil(log2(N));	
	if(N-pow(2,p)!=0){
		puts(" Error: Signal length must be a power of 2.");
		printf("N=%d, p=%lf\n",N,log2(N));
		puts(" --> process terminated.. ");
		exit(-1);
	}

	int i,k0,k_loc,nlen,irev;
	complex<double> tmp;

//	 Bit-reversal Scramble for FFT 
	for(i=0;i<N;i++){
		k_loc=i;
		k0=0;
		nlen=N;
		while(nlen>2){
			if(k_loc%2 != 0) k0+=(nlen/2); 
			k_loc/=2;
			nlen/=2;
		}
		irev=k0+k_loc;
		//printf("i=%d --> %d\n",i,irev);
		if(i < irev){
			tmp=Amp[i];
			Amp[i]=Amp[irev];
			Amp[irev]=tmp;
		};
	}


//		Butterfly Operation 
	int nseg;
	int j,i1,i2;
	double PI=4.0*atan(1.0);
	complex<double> zi(0.0,1.0);
	complex<double> X1,X2,W1,W2;

	zi*=isgn;

	nseg=N/2;
	nlen=1;
	while(nseg>0){	
		for(i=0;i<nseg;i++){
		for(j=0;j<nlen;j++){
			i1=i*(nlen*2)+j;
			i2=i1+nlen;
			X1=Amp[i1];
			X2=Amp[i2];
			W1=exp(zi*(PI/(double)nlen*j));
			W2=exp(zi*(PI/(double)nlen*(j+nlen)));
			Amp[i1]=X1+W1*X2;
			Amp[i2]=X1+W2*X2;
		}
		}
		nseg/=2;
		nlen*=2;
	}

	if(isgn==1) for(i=0;i<N;i++) Amp[i]/=N;

	return(p);
};
//------------------------------------------------------------
double DFT_prms::t(int i){ 
	return(ts+dt*i);
};
double DFT_prms::f(int i){
	return(df*i);
};

void DFT_prms::sigma(complex<double> *Amp){
	int i;
	double sig,PI=4.0*atan(1.0);
	double arg;

	for(i=1;i<Nt/2;i++){
		arg=(PI*i)/(100);
		sig=sin(arg)/arg;
		sig=pow(sig,6);
		Amp[i]*=sig;
		Amp[Nt-i]*=sig;
	};	
};
void DFT_prms::diff(complex<double> *Amp){
	complex<double> zi(0.0,1.0);
	complex<double> zw;
	int i;

	for(i=1;i<Nt/2;i++){
		zw=-zi*(dw*i);
		Amp[i]*=zw;
		Amp[Nt-i]*=(-zw);
	};	
};
void DFT_prms::integ(complex<double> *Amp){
	complex<double> zi(0.0,1.0);
	complex<double> zw;
	int i;

	for(i=1;i<Nt/2;i++){
		zw=-zi*(dw*i);
		Amp[i]/=zw;
		Amp[Nt-i]/=(-zw);
	};	
};

void DFT_prms::set_time(double t1, double t2, int n){
	double PI=4.0*atan(1.0);
	ts=t1;
	te=t2;
	Td=te-ts;
	dt=(te-ts)/n;

	df=1./Td;
	dw=2.*PI*df;
	
};
void DFT_prms::set_time(double dtau, int n){
	double PI=4.0*atan(1.0);
	Nt=n;
	dt=dtau;	
	te=Nt*dt;
	Td=te;
	df=1./Td;
	dw=2.*PI*df;
	ts=0.0;
};
//---------------------------------------------------------
InWv::InWv(){};	// empty constructor
InWv::InWv(int N){
	t1=0.0; t2=1.0;
	dt=0.0; T0=0.1;
	if(N>1) dt=(t2-t1)/(N-1);
	wvtyp=1; 
	nbrst=3;
	Nt=N;
	mem_alloc();	
	
};
void InWv::gen(char *fname){
//void InWv::gen(char fname[128]){
	int i;
	FILE *fp;
	char ch[128];
	fp=fopen(fname,"r");
	if(fp==NULL){
		printf("%s not found !!\n",fname);
		printf(" ---> abort process\n");
		exit(-1);
	};
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
	printf("**** waveform parameters ****\n");
	printf("wtyp=%d\n",wvtyp);
	printf("t1=%lf, t2=%lf\n",t1,t2);
	printf("dt=%lf\n",dt);
	printf("Nt=%d\n",Nt);
	printf("nbrst=%d\n",nbrst);
	printf("T0=%lf\n",T0);
}

void InWv::mem_alloc(){
	amp=(double *)malloc(sizeof(double)*Nt);
}

void InWv:: gen_wvfm(){
	double sgn=1.0;
	if(wvtyp<0) sgn=-1;
	switch(abs(wvtyp)){
	case 1: // sine function
		gen_sin(sgn);
		break;
	case 2:	// cosine function
		gen_cos(sgn);
		break;
	};
}

void InWv:: gen_sin(double sgn){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
		if(tt-t1<=nbrst*T0) amp[i]=sin(omg*tt)*sgn;
	}
}
void InWv:: gen_cos(double sgn){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
//		if(tt-t1<=nbrst*T0) amp[i]=1-cos(omg*tt);
		if(tt-t1<=nbrst*T0) amp[i]=cos(omg*tt)*sgn;
	}
}
void InWv::out(char *fout){
	FILE *fp=fopen(fout,"w");
	int i;
	for(i=0;i<Nt;i++) fprintf(fp,"%lf %lf\n",t1+dt*i, amp[i]);
//	for(i=0;i<Nt;i++) fprintf(fp,"%lf\n",amp[i]);
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

#if DB == 11 
int main(){
	char fname[128]="inwv0.dat";
	char fnout[128]="wv.out";

	Wv1D wv;
	wv.gen_wv(fname);
	wv.out_amp(fnout,',');
	wv.FFT(1);

	sprintf(fname,"awvt.out");
	wv.out_amp(fname,',');
	//awv.Butterworth(16.0,10.0);
	sprintf(fname,"awvw.out");
	wv.out_Amp(fname,0);

	return(0);
};
#endif

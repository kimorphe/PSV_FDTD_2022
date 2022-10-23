#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#ifndef __FFT__
	#define __FFT__
//	#include "fft.h"
#endif

using namespace std;

class InWv{
	public :
		double t1,t2,dt,T0;
		double *amp;
		int Nt,wvtyp,nbrst,nwv;
		int nsig; //narrowness factor used in Gaussian amplitude modulation
		InWv();
		//InWv(char *);	// constructor (data read from file char*)
		void gen(char *);	// constructor (data read from file char*)
		InWv(int Nt);	//
		void disp();
		void gen_sin(double sgn);
		void gen_cos(double sgn);
		void gen_wvfm();
		void out(char *);
		void Amod(double tb, int nsig);
	private:
		void mem_alloc();
};
int dft(complex<double> *Amp, int N, int isgn);
int fft(complex<double> *Amp, int N, int isgn);
class DFT_prms{
	public:
		int Nt;
		int p;	// p**2=Nt
		double dt;
		double dw;
		double df;
		double ts,te;
		double Td;
		void set_time(double t1, double t2, int n);
		void set_time(double dtau, int n);
		double t(int i); // return time_i
		double f(int i); // return freq_i
		double w(int i); // return omega_i
		void diff(complex<double> *Amp);
		void integ(complex<double> *Amp);
		void sigma(complex<double> *Amp);
	private:
};

class Wv1D{
	public:
		int Nt; // original data length
		int Np;	// FFT data length
		double *amp,*time,*tg;
		complex<double> *Amp;
		double *phi;
		double t1,t2,dt;
		double df;

		Wv1D();
		Wv1D(int Nt);
		Wv1D(char *fname);

		double val(double tt);
		int load(char fname[128]);
		void print_info();
		bool mllc;
		char data_file[128];
		int FFT(int sign);
		void unwrap(double f0);
		void renew_phase();
		void out_Amp(char *fn);
		void out_Amp_polar(char *fn);
		void out_Amp(char *fn,int ofst);
		void out_amp(char *fn,char dlm);
		int fft_stat;
		double L2(double t1,double t2);
		double max(double t1,double t2);

		void Butterworth(double tb, double Tw_6dB);
		void Gauss(double tb, double sig);
		double gdelay();

		InWv awv;
		//void gen_wv(char fname[128]);
		void gen_wv(char *fname);
	private:
	protected:
};
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax);
//------------------------------------------------------------

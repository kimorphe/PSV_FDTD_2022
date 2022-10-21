#include <stdio.h>
#include <math.h>
#include "fdm2d.h" 

Cmplx :: Cmplx(){};
Cmplx :: Cmplx(double a, double b){
	re=a;
	im=b;
};
double Cmplx :: Abs(){
	return(sqrt(re*re+im*im));
};
Cmplx Cmplx :: Cnjg(){
	Cmplx Z(re,-im);
	return(Z);
};
void Cmplx::set(double a, double b){
	re=a; im=b;
};

Cmplx Cmplx::operator+(Cmplx z){
	Cmplx w(re+z.re,im+z.im);
	return(w);
};
Cmplx Cmplx::operator-(Cmplx z){
	Cmplx w(re-z.re,im-z.im);
	return(w);
}
Cmplx Cmplx::operator*(Cmplx z){
	Cmplx w(re*z.re-im*z.im, re*z.im+im*z.re);
	return(w);
}
Cmplx Cmplx::operator*(double a){
	Cmplx w(a*re, a*im);
	return(w);
}
void Cmplx::operator=(double a){
	re=a;
	im=0.0;
}
Cmplx Cmplx::operator/(Cmplx z){
	Cmplx w(re*z.re+im*z.im, im*z.re-re*z.im);
	double r=z.re*z.re+z.im*z.im;
	w.re/=r;
	w.im/=r;
	return(w);
}

void Cmplx::disp(){
	printf("(%lf, %lf)\n",re,im);
};
//---------------------------------------------
Cmplx exp(Cmplx z){
	Cmplx w(exp(z.re)*cos(z.im),exp(z.re)*sin(z.im));
	return(w);
}
//--------------------------------------------
/*
int main(){
	Cmplx z1(1.0,2.0); 
	Cmplx z2(-1.0,2.0); 
	Cmplx zi(0.0,1.0); 
	Cmplx z;

	z=z1+z2;
	z.disp();

	z=z1-z2;
	z.disp();

	z=z1*z2;
	z.disp();

	z=z2*(z1/z2);
	z.disp();

	Cmplx w;
	w.set(100,100);
	w.disp();
	z=w*2.0;
	z.disp();

	int Nt=100;
	double T=4.0*atan(1.0)*2;
	double dt=T/(Nt-1);
	FILE *fp=fopen("cmplx.out","w");
	for(int i=0;i<Nt;i++){
		z.set(-0.1*i*dt,dt*i);
		w=exp(z);
		fprintf(fp,"%lf %lf %lf\n",dt*i,w.re,w.im);
	}

	w=123;
	w.disp();

	return 0;
}
*/

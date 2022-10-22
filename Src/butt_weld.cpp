#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//----------------------------------------------------------
class Bump{
	public:
		double a;	// half width
		double dh;	// height
		double h;	// elevation
		double Rd;	// radius
		double b;
		double sgn;
		Bump(double a, double dh);
		double height(double x);
		double xc;
		void shift(double xc, double h);
	private:
};

Bump::Bump(double aw, double ht){
	a=aw;
	dh=abs(ht);
	sgn=1.0;
	if( ht < 0.0) sgn=-1.0;
		
	Rd=(a*a+dh*dh)/dh*0.5;
	b=sqrt(Rd*Rd-a*a);

	xc=0.0;
	h=0.0;
};
void Bump::shift(double x_shift, double y_shift){
	xc=x_shift;
	h=y_shift;
};
double  Bump::height(double x){
	double y=h;
	x-=xc;
	if(abs(x)>a) return(y);
	y=sgn*(sqrt(Rd*Rd-x*x)-b)+h;
	return(y);
};
//----------------------------------------------------------
// BEAD
// a_top, dh_top, h_top
// a_btm, dh_btm, h_btm
// 
int main(){
	double x,dx,x1,x2;
	int i,Nx=201;
	double a=5.0, dh=2.5;
	Bump bp_t(a,dh);
	Bump bp_b(a*0.5, -dh*0.5);

	double xc=5.0;
	double h=12.0;

	bp_t.shift(xc,h);
	bp_b.shift(xc,0.0);

	x1=-20.0; 
	x2= 20.0;
	dx=(x2-x1)/(Nx-1);
	for(i=0;i<Nx;i++){
		x=x1+dx*i;
		printf("%lf, %lf, %lf\n",x,bp_t.height(x),bp_b.height(x));
	};
	return(0);
};

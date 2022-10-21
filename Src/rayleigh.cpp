#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double Rs(double sL, double sT, double s){

	double s2=s*s;
	double sL2=sL*sL, sT2=sT*sT;
	double R1,RL,RT;

	
	R1=(2.0*s2-sT2);
	R1*=R1;
	RL=sqrt(s2-sL2);
	RT=sqrt(s2-sT2);

	return(R1-4.0*s2*RL*RT);
};

double Rwave(double cL, double cT){
	double s0,cR;
	double sR,sL=1.0/cL, sT=1.0/cT; // Slownesses
	double tol=1.e-06;
	int it=0,it_max=100;
	double alph=1.5;
	double sa=sT, sb=alph*sT;
	double err,err0=fabs(Rs(sL,sT,sa));

	do{
		it++;
		s0=0.5*(sa+sb);
		err=Rs(sL,sT,s0);

		if( err < 0.0){
			sb=s0;
		}else{
			sa=s0;
		}
		if(it>= it_max) break;
	}while(fabs(err)/err0 > tol);

	sR=s0;
	cR=1.0/sR;
	printf("iteration=%d, relative error=%lf\n",it,fabs(err)/err0);

	return(cR);
};

int main(){

	double cL=2.0;	// L-wave velocity
	double cT=1.0;	// T-wave velocity
	double cR;	// R-wave velocity
	double sL=1.0/cL, sT=1.0/cT; // Slownesses
	double sR, double s,Rval;
	FILE *fp=fopen("rfun.out","w");

	double sa=sT, sb=alph*sT;
	int i,N=100;
	double ds=(sb-sa)/(N-1);

	for(i=0;i<N;i++){
		s=sa+ds*i;
		Rval=Rs(sL,sT,s);
		fprintf(fp,"%lf %lf\n",s,Rval);
	} 
	fclose(fp);	


	cR=Rwave(cL,cT);
	sR=1.0/cR;
	printf("cR=%lf, sR=%lf (cL=%lf, cT=%lf)\n",cR,sR, cL,cT);

	double beta=cL/cT; 
	beta*=beta;
	double nu=(beta-2.0)/(beta-1.0)*0.5;
	double cR_apprx=(0.862+1.14*nu)/(1+nu)*cT;
	printf("Rational approximation cR=%lf (Poisson ratio nu=%lf)\n",cR_apprx,nu);
	
	return(0);
};

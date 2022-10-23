#include <stdio.h>
#include <stdlib.h>
#include "fdm2d.h"

class CNTRL{
	public:
		int nele;
		int n_meas;
		void load(char *fn);
		int *actv;
		double *delay;
		int idx_a;
		int idx_d;
		bool pop_a();
		double pop_d();
	private:
};
bool CNTRL::pop_a(){
	bool active=true;
	if(actv[idx_a]==0) active=false;
	idx_a++;
	return(active);
};
double CNTRL::pop_d(){
	return(delay[idx_d++]);
};
void CNTRL::load(char *fn){
	char cbff[128];
	FILE *fp=fopen(fn,"r");
	if(fp==NULL){
		printf("Can't open %s\n",fn);
		exit(-1);
	}

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nele);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&n_meas);

	actv=(int *)malloc(sizeof(int)*nele*n_meas);
	delay=(double *)malloc(sizeof(double)*nele*n_meas);

	int i,j,k=0;
	for(i=0;i<n_meas;i++){
		fgets(cbff,128,fp);
	for(j=0;j<nele;j++){
		fscanf(fp,"%d %lf\n",actv+k,delay+k);
		printf("actv=%d delay=%lf\n",actv[k],delay[k]);
		k++;
	}
	}
};
int main(int argc, char *argv[]){
	CNTRL ctr;
	char fname[128]="array.inp";
	ctr.load(fname);

	int m,el;
	double dly;
	for(m=0; m < ctr.n_meas; m++) {
	for(el=0; el<ctr.nele; el++){
		dly=ctr.pop_d();
		if(ctr.pop_a()) printf("delay=%lf\n",dly);
		//if(printf("%d, %lf\n",ctr.pop_a(),ctr.pop_d());
	}
	}
	return(0);	
};

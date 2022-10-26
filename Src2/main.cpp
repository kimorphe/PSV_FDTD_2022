#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"fdm2d.h"

int main(int argc, char *argv[]){

	CNTRL ctr;
	char fgeom[128]="geom.inp";
	char ftset[128]="tset.inp";
	char fsrce[128]="src.inp";
	char farry[128]="array.inp";

	printf("********* log files *********\n");
	ctr.setup_domain(fgeom);
	ctr.time_setting(ftset);
	ctr.src_setting(fsrce);
	ctr.wvfm_setting();
	ctr.array_setting(farry);
	printf("*****************************\n");

	
	int it,isum,m;

	ctr.mark_src_grid();
/*	
	for(m=0;m<ctr.ary.nmeas;m++){
		printf("m=%d\n",m);
		isum=0;
		ctr.snapshot(m,isum++,0);
	for(it=0; it<ctr.Nt; it++){
		ctr.v2q(it); 
		ctr.q2v(it);
		ctr.capture(it);
		if(ctr.out_time(it)) ctr.snapshot(m,isum++,it);
	};
		//for(int i=0; i<ctr.nrec; i++) ctr.recs[i].fwrite();
		ctr.fwrite_ary();
		ctr.round++;
		ctr.clear();
	}
*/
	return(0);
};

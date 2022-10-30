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
	ctr.setup_domain(fgeom);		// define domain geometry & material properties
	ctr.time_setting(ftset);		// set time range and increment
	ctr.src_setting(fsrce);			// setup transducer elements
	ctr.wvfm_setting();			// define excitation signal(waveform)
	ctr.array_setting(farry);		// define array control 
	printf("*****************************\n");

	
	int it,isum,m;

	for(m=0;m<ctr.ary.nmeas;m++){
		ctr.mark_src_grid();	// ???????????????? coordination with array
		printf("m=%d\n",m);
		isum=0;
		//ctr.snapshot(m,isum++,0);
	for(it=0; it<ctr.Nt; it++){
		printf("it=%d\n",it);
		ctr.s2v(it);		// ????????????????  source term management
		ctr.v2s(it);
		//ctr.capture(it);
		//if(ctr.out_time(it)) ctr.snapshot(m,isum++,it);
	};
		//for(int i=0; i<ctr.nrec; i++) ctr.recs[i].fwrite();
		//ctr.fwrite_ary();
		ctr.round++;
		//ctr.clear();			
	}
	return(0);
};

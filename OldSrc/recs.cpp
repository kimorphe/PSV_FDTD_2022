#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fdm2d.h"
//#include "crack.h"

//using namespace std; 

int set_recs(char *fname, Dom2D dom, Recs **rec,int Nt, double dt){
	FILE *fp=fopen(fname,"r");
	int i,nrec,idir,ityp,npnt;
	char cbff[128];
	double Y1[2],Y2[2];

	while(fgets(cbff,7,fp)!=NULL){
		if(strcmp(cbff,"##RECS")==0){
			fgets(cbff,2,fp);
			fgets(cbff,128,fp);
			fscanf(fp,"%d\n",&nrec);
				(*rec)=(Recs *)malloc(sizeof(Recs)*nrec);
			fgets(cbff,128,fp);
			for(i=0;i<nrec;i++){
				fscanf(fp,"%d %d\n",&idir,&ityp);
				fscanf(fp,"%lf %lf\n",Y1,Y1+1);
				fscanf(fp,"%lf %lf\n",Y2,Y2+1);
				fscanf(fp,"%d\n",&npnt);
				(*rec)[i].setup(idir,ityp,Y1,Y2,npnt,dom,Nt,dt);
				(*rec)[i].ID=i;
				(*rec)[i].print(dom);
			}
			fgets(cbff,6,fp);
			if(strcmp(cbff,"##END")!=0){
				puts("Incorrect data format detected");
				puts("when reading receiver information.");
			}
		}
	};
	fclose(fp);
	return nrec;
}

void Recs :: setup(
	int jdir, // alignment 0:horizontal, 1: vertical
	int jtyp, // grid type v1=0, v2=1, s11=2, s12=3, s22=4 
	double *Y1, // End point 1
	double *Y2,  // End point 2
	int npnt, 	// number of points
	Dom2D dom, 	// domain data
	int nt, 	// number of time steps
	double dtau	// time increment
){
	
	int i,ntmp;
	int indx[2];
	int *ktmp=(int *)malloc(sizeof(int)*npnt);
	double dX=0.0;	// pitch 
	double xcod[2],Len;

	idir=jdir;	// alignment (0:horiz. 1:vert.)
	ityp=jtyp;	// grid type
	Nt=nt;
	dt=dtau;
	for(i=0;i<2;i++) {
		X1[i]=Y1[i];
		X2[i]=Y2[i];
	}

	ntmp=npnt; 
	npnt=0;
	Ng=0;
	Len=dist2d(X1,X2);
	if( ntmp > 1) dX=Len/(ntmp-1); 
	dom.gridNum(ityp);	// get number of grid  
	cod2indx(X1,ityp,dom.Xa,dom.dx,indx);
	switch(idir){
		case 0:	//horizontally aligned points
			if(indx[1] < 0) break;
			if(indx[1] >= dom.Nx[1]) break;
			xcod[1]=X1[1];
			npnt=ntmp;
			break;
		case 1:	// vertically aligned points
			if(indx[0] < 0) break;
			if(indx[0] >= dom.Nx[0]) break;
			xcod[0]=X1[0];
			npnt=ntmp;
			printf("npnt=%d\n",npnt);
			break;	
	};
	for(i=0;i<npnt;i++){
		xcod[idir]=X1[idir]+i*dX;
		cod2indx(xcod,ityp,dom.Xa,dom.dx,indx);
		if(indx[idir] <  0) continue;
		if(indx[idir] >= dom.Nx[idir]) continue;
//		if(dom.kcell[indx[0]][indx[1]] <0) continue;
		ktmp[Ng]=ij2l(indx[0],indx[1],ityp,dom.Ndiv);	
		Ng++;
	}

	krec=(int *)malloc(sizeof(int)*Ng);
	for(i=0;i<Ng;i++){
		krec[i]=ktmp[i];
	}

	mem_alloc2D(Ng,Nt,&val);
};
void Recs :: print(Dom2D dom){
	int i,j,l;
	double xcod[2];
	for(l=0;l<Ng;l++){
		l2ij(krec[l],ityp,dom.Ndiv,&i,&j);
		indx2cod(i,j,ityp,dom.Xa,dom.dx,xcod);
//		printf("xcod=%lf %lf\n",xcod[0],xcod[1]);
	}
};

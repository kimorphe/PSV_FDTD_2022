#include<iostream> 
#include<stdio.h>
#include<stdlib.h>
#include<math.h> 
#include<string.h>
#include "fdm2d.h"
using namespace std;


//	-----------------------------------------
//		CONSTRUCTOR
Fld2D :: Fld2D(Dom2D dom, char *ficon){
	
	int i,ndim=2;

	Ndiv[0]=dom.Ndiv[0];
	Ndiv[1]=dom.Ndiv[1];

	mem_alloc();	// memory allocation

	set_IC(dom,ficon);  // Initial Condition

	npml=0;
	for(i=0;i<ndim;i++){
		if( dom.nwa[i]> 0) npml++;
		if( dom.nwb[i]> 0) npml++;
	}

	pml=(PML *)malloc(sizeof(PML)*npml);

	i=0;
	if(dom.nwa[0]>0){
		pml[i].Xa[0]=dom.Xa[0]; 
		pml[i].Xa[1]=dom.Xa[1];
		pml[i].Xb[0]=dom.xa[0]; 
		pml[i].Xb[1]=dom.Xb[1];
		pml[i].setup(dom,ficon);
		pml[i].ID=0;
		i++;
	}
	if(dom.nwb[0]>0){
		pml[i].Xa[0]=dom.xb[0]; 
		pml[i].Xa[1]=dom.Xa[1];
		pml[i].Xb[0]=dom.Xb[0]; 
		pml[i].Xb[1]=dom.Xb[1];
		pml[i].setup(dom,ficon);
		pml[i].ID=1;
		i++;
	}
	if(dom.nwa[1]>0){
		pml[i].Xa[0]=dom.Xa[0]; 
		pml[i].Xa[1]=dom.Xa[1];
		pml[i].Xb[0]=dom.Xb[0]; 
		pml[i].Xb[1]=dom.xa[1];
		pml[i].setup(dom,ficon);
		pml[i].ID=2;
		i++;
	}
	if(dom.nwb[1]>0){
		pml[i].Xa[0]=dom.Xa[0]; 
		pml[i].Xa[1]=dom.xb[1];
		pml[i].Xb[0]=dom.Xb[0]; 
		pml[i].Xb[1]=dom.Xb[1];
		pml[i].setup(dom,ficon);
		pml[i].ID=3;
		i++;
	}
};
//------------- GRID NUMBER  ------------------

//	SET NUMBER OF GRIDS OF TYPE itype to Nx[2]
//      (identical with the function Dom2D:: gridNum) 
void Fld2D ::gridNum(int ityp){
	switch(ityp){
	case 0:	// v1-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1];
		break;
	case 1:	// v2-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1]+1;
		break;
	case 2:	// s11-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	case 3:	// s12-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1]+1;
		break;
	case 4:	// s22-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	};
	Ng=Nx[0]*Nx[1];
};
//----------------------------------------------

//		MEMORY ALLOCATION 

void Fld2D::mem_alloc(){

//		v1-grid 
	Fld2D::gridNum(0);
	mem_alloc2D(Nx[0],Nx[1],&v1);

//		v2-grid 
	Fld2D::gridNum(1);
	mem_alloc2D(Nx[0],Nx[1],&v2);

//		s11,s22-grid 
	Fld2D::gridNum(2);
	mem_alloc2D(Nx[0],Nx[1],&s11);
	mem_alloc2D(Nx[0],Nx[1],&s22);

//		s12-grid 
	Fld2D::gridNum(3);
	mem_alloc2D(Nx[0],Nx[1],&s12);

};

int Fld2D::set_IC(Dom2D dom, char *ficon){

	FILE *fp=fopen(ficon,"r");
	char cbff[128],fname[128];
	int i,j,l,it0=0;
	int init,inwv_typ,inwv_num;
	double ain,xp0[2];
	double gin,ginx,giny;
	
	
	if(fp == NULL){
		s110=0.0;
		s220=0.0;
		puts("Can't open icon.inp.");
		puts(" --> initial field set to identically zero.");
		return(-1);
	}	

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&init);
	printf("init=%d",init);
	std::fflush(stdout);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&ain);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&inwv_typ);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf \n",xp0,xp0+1);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&inwv_num);

	printf("inwvnum=%d",inwv_num);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",&s110,&s220);


	fclose(fp);
	
	sprintf(fname,"inwv%d.dat",inwv_num);
	InWv inwv=InWv(fname);

	for(l=0;l<5;l++){
		gridNum(l);
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		gin=0.0; 
		if(init==1) gin=plane_wave(i,j,l,it0,inwv_typ,ain,xp0,inwv,dom,&ginx,&giny);
		switch(l){
		case 0:
			v1[i][j]=gin;
			break;
		case 1: 
			v2[i][j]=gin;
			break;
		case 2: 
			s11[i][j]=gin+s110;
			break;
		case 3:
			s12[i][j]=gin;
			break;
		case 4:
			s22[i][j]=gin+s220;
			break;
		};
	}
	}
	}
	return(0);	
};

void mem_alloc2D(int nx, int ny, double ***pt){
	
	int i,ng=nx*ny;
	double *p1;

	p1=(double *)malloc(sizeof(double)*nx*ny);
	for(i=0;i<ng;i++) p1[i]=0.0;

	(*pt)=(double **)malloc(sizeof(double *)*nx);	
	for(i=0;i<nx;i++) (*pt)[i]=p1+i*ny;
}
//----------------------------------------------

//		GENERATE 1D INDEX FOR v1 
void Fld2D::gen_indx1(Dom2D dom){
	int i,j;
	int Nin,Nbnd,Nex;
	int il,ir;

	Fld2D :: gridNum(0);

	for(int m=0;m<2;m++){	//count grids at m=0, results stored at m=1; 
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		il=1; ir=1;
		if(i>0) il=dom.kcell[i-1][j];	
		if(i<Ndiv[0]) ir=dom.kcell[i][j];	
		switch(il+ir){
		case 0:	// interior grid 
			if(m==1) kint1[Nin]=i*Nx[1]+j;
			Nin++;
			break;
		case 1: // boundary grid
			if(m==1){
				kbnd1[Nbnd]=i*Nx[1]+j;
				if(il==1) kbnd1[Nbnd]*=-1;
			}
			Nbnd++; 
			break;
		case 2: // exterior grid
			Nex++;
			break;
		};
	}
	}

//		printf("Nin=%d Nbnd=%d Nex=%d\n",Nin,Nbnd,Nex);
//		printf("Ng=%d\n",Ng);
		if(m==0){
			kint1=(int *)malloc(sizeof(int)*Nin);
			kbnd1=(int *)malloc(sizeof(int)*Nbnd);
			Nint1=Nin;
			Nbnd1=Nbnd;
		}
	}
};
//		GENERATE 1D INDEX FOR v2 
void Fld2D::gen_indx2(Dom2D dom){
	int i,j;
	int Nin,Nbnd,Nex;
	int iu,id;

	Fld2D :: gridNum(1);

	for(int m=0;m<2;m++){
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		id=1; iu=1;
		if(j>0) id=dom.kcell[i][j-1];	
		if(j<Ndiv[1]) iu=dom.kcell[i][j];	
		switch(iu+id){
		case 0: // interior grid
			if(m==1) kint2[Nin]=i*Nx[1]+j;
			Nin++;
			break;
		case 1: // boundary grid
			if(m==1){
				kbnd2[Nbnd]=i*Nx[1]+j;
				if(id==1) kbnd2[Nbnd]*=-1; 
			}
			Nbnd++; 
			break;
		case 2: //exterior grid
			Nex++;
			break;
		};
	}
	}

//		printf("Nin=%d Nbnd=%d Nex=%d\n",Nin,Nbnd,Nex);
//		printf("Ng=%d\n",Ng);
		if(m==0){
			kint2=(int *)malloc(sizeof(int)*Nin);
			kbnd2=(int *)malloc(sizeof(int)*Nbnd);
			Nint2=Nin;
			Nbnd2=Nbnd;
		}
	}
};

//		GENERATE 1D INDEX FOR s12 
void Fld2D::gen_indx12(Dom2D dom){
	int i,j,k,m,alph;
	int k0,k1,k2,k3;
	int Nin,Nbnd,Nex;
	bool filled[4]; 

	Fld2D :: gridNum(3);

	for(m =0;m<2;m++){
		Nin=0;
		Nbnd=0;
		Nex=0;
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		k0=1; k1=1; k2=1; k3=1;
		for(k=0;k<4;k++) filled[k]=1;
		if(i==0){ 
			filled[1]=0; filled[2]=0;
		}
		if(i==Nx[0]-1){ 
			filled[0]=0; filled[3]=0;
		}
		if(j==0){ 
			filled[2]=0; filled[3]=0;
		}
		if(j==Nx[1]-1){ 
			filled[0]=0; filled[1]=0;
		}

		if(filled[0]==1) k0=dom.kcell[i][j];
		if(filled[1]==1) k1=dom.kcell[i-1][j];
		if(filled[2]==1) k2=dom.kcell[i-1][j-1];
		if(filled[3]==1) k3=dom.kcell[i][j-1];

		k=k0+k1+k2+k3;

		switch(k){
		case 0: // interior grid
			if(m==1) kint12[Nin]=i*Nx[1]+j;
			Nin++;
			break;
		case 1: // boundary grid
		case 3:
			if(m==1) kbnd12[Nbnd]=i*Nx[1]+j;
			Nbnd++; 
			break;
		case 2: 
			alph= k0*k1 + k1*k2 + k2*k3 + k3*k0;	
			if(alph==1){	// boundary grid
				 if(m==1) kbnd12[Nbnd]=i*Nx[1]+j;
				Nbnd++; 
			}
			if(alph==0){	// interior grid
				if(m==1) kint12[Nin]=i*Nx[1]+j;
				Nin++;
			}			
			break;	
		case 4: //exterior grid
			Nex++;
			break;
		};
	}
	}
		if(m==0){
			kint12=(int *)malloc(sizeof(int)*Nin);
			kbnd12=(int *)malloc(sizeof(int)*Nbnd);
			Nint12=Nin;
			Nbnd12=Nbnd;
		}
	}
};

int Fld2D:: del_src_grid(Src src){
	int i,j;
	int jbnd,jsld=0,i0=0,idel;
	int ngrd=src.ngrd;

	switch(src.ityp){
	case 0:	// v1-grid

		for(j=0; j<Nbnd1; j++){  // walking through kbnd1[Nbnd1]
			jbnd=kbnd1[j];
			idel=0;
			for(i=i0;i<ngrd;i++){ // scan src.ksrc[ngrd]
				if(src.ksrc[i]==jbnd){	// found !
					jsld++;
					i0=i;
					idel=1;
					break;
				}
			}
			if(idel==0) kbnd1[j-jsld]=kbnd1[j];
		}
		Nbnd1-=jsld;
		printf("Number of boundary grid deleted is %d\n",jsld);

		break;
	case 1: // v2-grid
		for(j=0; j<Nbnd2; j++){  // walking through kbnd2[Nbnd2]
			jbnd=kbnd2[j];
			idel=0;
			for(i=i0;i<ngrd;i++){ // scan src.ksrc[ngrd]
				if(src.ksrc[i]==jbnd){	// found !
					jsld++;
					i0=i;
					idel=1;
					break;
				}
			}
			if(idel==0) kbnd2[j-jsld]=kbnd2[j];
		}
		Nbnd2-=jsld;
		break;
	case 3: // s12-grid
		for(j=0; j<Nbnd12; j++){  // walking through kbnd12[Nbnd12]
			jbnd=kbnd12[j];
			idel=0;
			for(i=i0;i<ngrd;i++){ // scan src.ksrc[ngrd]
				if(src.ksrc[i]==jbnd){	// found !
					jsld++;
					i0=i;
					idel=1;
					break;
				}
			}
			if(idel==0) kbnd12[j-jsld]=kbnd12[j];
		}
		Nbnd12-=jsld;
		break;
	};


	return(jsld);

};

//----------------------------------------------
//
//	TIME-STEPPING : (Velocity --> Stress )
//
//
void Fld2D::v2s(Dom2D dom, double dt){
	int i,j,k,l,isum;
	int i0,j0,ngx,ngy;
	double almb,amu;
	double dx[2],dv1,dv2;

	double ds11,ds22,ds12;
	double Ax0,Ax1,Ax2;
	double Ay0,Ay1,Ay2;

	dx[0]=dom.dx[0];
	dx[1]=dom.dx[1];
		
//	rho=dom.rho;	//density
	almb=dom.almb;	// Lame constant
	amu=dom.amu;	// shear rigidity

	Ax0=(almb+2.*amu)*dt/dx[0];
	Ay0=almb*dt/dx[1];

	Ax1=almb*dt/dx[0];
	Ay1=(almb+2.*amu)*dt/dx[1];

	Ax2=amu*dt/dx[0];
	Ay2=amu*dt/dx[1];

//		--- Internal Grids ---
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		if(dom.kcell[i][j]!=0) continue;	

		dv1=v1[i+1][j]-v1[i][j];
		dv2=v2[i][j+1]-v2[i][j];

		ds11=Ax0*dv1+Ay0*dv2;
		ds22=Ax1*dv1+Ay1*dv2;
		s11[i][j]+=ds11;
		s22[i][j]+=ds22;

	}
	}

	gridNum(3);
	l=0;
	isum=0;	
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		if(kint12[l]!=isum){
			isum++;
			 continue;
		}
		l++;
		dv1=v1[i][j]-v1[i][j-1];
		dv2=v2[i][j]-v2[i-1][j];

		ds12=Ax2*dv2+Ay2*dv1;
		s12[i][j]+=ds12;
		isum++;
	}
	}

	for(l=0;l<Nbnd12;l++){
		l2ij(kbnd12[l],3,Ndiv,&i,&j);
		s12[i][j]=0.0;
	}
	

//		PML DOMAIN 
//	Import PML Field Data  ( PML --> Field )
	for(k=0;k<npml;k++){
		pml[k].v2s(amu,almb,dt);

		i0=pml[k].ia[0];
		j0=pml[k].ia[1];
		ngx=pml[k].nx[0]-1;
		ngy=pml[k].nx[1]-1;

		for(i=0;i<ngx;i++){
		for(j=0;j<ngy;j++){
			s11[i+i0][j+j0]=pml[k].s11[i][j];
			s22[i+i0][j+j0]=pml[k].s22[i][j];
		}
		}

		ngx=pml[k].nx[0];
		ngy=pml[k].nx[1];
		for(i=1;i<ngx-1;i++){
		for(j=1;j<ngy-1;j++){
			s12[i+i0][j+j0]=pml[k].s12[i][j];
		}
		}
	}

//	 Export Data to PML ( Field --> PML )
	for(k=0;k<npml;k++){
		i0=pml[k].ia[0];
		j0=pml[k].ia[1];

		ngx=pml[k].nx[0];
		ngy=pml[k].nx[1];

/*
		for(i=0;i<ngx;i++){
			pml[k].s12[i][0]=s12[i+i0][j0];
			pml[k].s12[i][ngy-1]=s12[i+i0][ngy-1+j0];
		}	
		for(j=0;j<ngy;j++){
			pml[k].s12[0][j]=s12[i0][j+j0];
			pml[k].s12[ngx-1][j]=s12[ngx-1+i0][j+j0];
		}	
*/
		switch(pml[k].ID){
		case 0:
			for(j=0;j<ngy;j++) pml[k].s12[ngx-1][j]=s12[ngx-1+i0][j+j0];
			break;
		case 1:
			for(j=0;j<ngy;j++) pml[k].s12[0][j]=s12[i0][j+j0];
			break;
		case 2:
			for(i=0;i<ngx;i++) pml[k].s12[i][ngy-1]=s12[i+i0][ngy-1+j0];
			break;
		case 3:
			for(i=0;i<ngx;i++) pml[k].s12[i][0]=s12[i+i0][j0];
			break;
		};
	}

};

//----------------------------------------------
//
//	TIME-STEPPING : ( Stress --> Velocity )
//
//
void Fld2D::s2v(Dom2D dom,double dt){
	int i,j,k,l,isum,id;
	int ngx,ngy,i0,j0;
	double dx[2];
	double ds11,ds22,ds12;
	double dv1,dv2;
	double rho;
	double dtr0,dtr1,sgn;
	double bvl;

	dx[0]=dom.dx[0];
	dx[1]=dom.dx[1];
		
	rho=dom.rho;
	dtr0=dt/(rho*dx[0]);
	dtr1=dt/(rho*dx[1]);

//		--- v1-grid --
	gridNum(0);
	l=0;
	isum=0;	
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		if(kint1[l]!=isum){
			isum++;
			 continue;
		}
		l++;
		ds11=s11[i][j]-s11[i-1][j];
		ds12=s12[i][j+1]-s12[i][j];
		dv1=dtr0*ds11+dtr1*ds12;
		v1[i][j]+=dv1;
		isum++;
	}
	}

	for(l=0;l<Nbnd1;l++){
		l2ij(abs(kbnd1[l]),0,Ndiv,&i,&j);
		id=0; sgn=-1.0;
		bvl=0.0;	// boundary value
		if(i==0 || i==Ndiv[0]) bvl=s110;
//		if(i==0 && dom.nwa[0]>0) bvl=s110;
//		if(i==Ndiv[0] && dom.nwb[0]>0) bvl=s110;
		if(kbnd1[l] > 0){
			 id=1; sgn=1.0;
		} 
//		ds11=sgn*(bvl-s11[i-id][j])*2.0;
		ds11=sgn*(bvl-s11[i-id][j])*1.5;
		ds12=s12[i][j+1]-s12[i][j];
		dv1=dtr0*ds11+dtr1*ds12;
		v1[i][j]+=dv1;
	}


//		--- v2-grid --
	gridNum(1);
	l=0;
	isum=0;	
	for(i=0;i<Nx[0];i++){
	for(j=0;j<Nx[1];j++){
		if(kint2[l]!=isum){
			isum++;
			 continue;
		}
		l++;
		ds12=s12[i+1][j]-s12[i][j];
		ds22=s22[i][j]-s22[i][j-1];
		dv2=dtr0*ds12+dtr1*ds22;
		v2[i][j]+=dv2;
		isum++;
	}
	}

	for(l=0;l<Nbnd2;l++){
		l2ij(abs(kbnd2[l]),1,Ndiv,&i,&j);

		id=0; sgn=-1.0;
		bvl=0.0;
		if(j==0 || j==Ndiv[1]) bvl=s220;
//		if(j==0 && dom.nwa[1]>0) bvl=s220;
//		if(j==Ndiv[1] && dom.nwb[1]>0) bvl=s220;
		if( kbnd2[l] >0){
			id=1; sgn=1.0;
		}
		ds12=s12[i+1][j]-s12[i][j];
//		ds22=sgn*(bvl-s22[i][j-id])*2.0;
		ds22=sgn*(bvl-s22[i][j-id])*1.5;
		dv2=dtr0*ds12+dtr1*ds22;
		v2[i][j]+=dv2;
	}

//		PML DOMAIN 
//		Import PML Field Data (PML --> Field) 
	for(k=0;k<npml;k++){ 
		pml[k].s2v(rho,dt);
		ngx=pml[k].nx[0];
		ngy=pml[k].nx[1];
		i0=pml[k].ia[0];
		j0=pml[k].ia[1];	

		for(i=1;i<ngx-1;i++){
		for(j=0;j<ngy-1;j++){
			v1[i+i0][j+j0]=pml[k].v1[i][j];
		}
		}
		for(i=0;i<ngx-1;i++){
		for(j=1;j<ngy-1;j++){
			v2[i+i0][j+j0]=pml[k].v2[i][j];
		}
		}
	}

//		Export Data to PML (Field --> PML) 
	for(k=0;k<npml;k++){ 

		ngx=pml[k].nx[0];
		ngy=pml[k].nx[1];
		i0=pml[k].ia[0];
		j0=pml[k].ia[1];	

/*
		for(j=0;j<ngy-1;j++){
			pml[k].v1[0][j]=v1[i0][j+j0];
			pml[k].v1[ngx-1][j]=v1[ngx-1+i0][j+j0];
		}	
		if(k==0){
		for(j=0;j<ngy-1;j++){
			pml[k].v1[0][j]=0.0;
		}
		}
		for(i=0;i<ngx-1;i++){
			pml[k].v2[i][0]=v2[i+i0][j0];
			pml[k].v2[i][ngy-1]=v2[i+i0][ngy-1+j0];
		}
*/

		switch(pml[k].ID){
		case 0 :
			for(j=0;j<ngy-1;j++) pml[k].v1[ngx-1][j]=v1[ngx-1+i0][j+j0];
			break;
		case 1 :
			for(j=0;j<ngy-1;j++) pml[k].v1[0][j]=v1[i0][j+j0]; 
			break;
		case 2 :
			for(i=0;i<ngx-1;i++) pml[k].v2[i][ngy-1]=v2[i+i0][ngy-1+j0];
			break;
		case 3 :
			for(i=0;i<ngx-1;i++) pml[k].v2[i][0]=v2[i+i0][j0];
			break;
		};
	}

};
void Fld2D :: apply_src(int it, Src src,Dom2D dom){
	int l,m,jwv,i,j,jt,id;
	double amp;
	int ityp=src.ityp;
	double rho=dom.rho;
	double dtr0=src.dt/(rho*dom.dx[0]);
	double dtr1=src.dt/(rho*dom.dx[1]);
	double sgn,ds11,ds22,ds12,dv1,dv2;

	for(m=0;m<src.ngrd;m++){
		l=src.ksrc[m];
		l2ij(abs(l),abs(ityp),Ndiv,&i,&j);
//		printf("%d,%d\n",i,j);
		jwv=src.wvID[m];
		amp=0.0;
		jt=it-src.idly[m];
		if(jt >= 0) amp=src.inwvs[jwv].amp[jt];

		switch (ityp){
			case  0: // v1-grid
			if(src.ifld==0){
					amp+=s110;
				id=1; sgn=1.0;
				if(l < 0){
					id=0; sgn=-1.0;
				} 
//				ds11=sgn*(amp-s11[i-id][j])*2.0;
				ds11=sgn*(amp-s11[i-id][j])*1.5;
				ds12=s12[i][j+1]-s12[i][j];
				dv1=dtr0*ds11+dtr1*ds12;
				v1[i][j]+=dv1;
			}else if(src.ifld==1){
				v1[i][j]=amp;
			}
				break;
			case  1: // v2-grid
			if(src.ifld==0){
				amp+=s220;
				id=1; sgn=1.0;
				if( l < 0){
					id=0; sgn=-1.0;
				}
				ds12=s12[i+1][j]-s12[i][j];
//				ds22=sgn*(amp-s22[i][j-id])*2.0;
				ds22=sgn*(amp-s22[i][j-id])*1.5;
				dv2=dtr0*ds12+dtr1*ds22;
				v2[i][j]+=dv2;
			}else if(src.ifld==1){
				v2[i][j]=amp;
			}
				break;
			case  2: // s11-grid
				s11[i][j]=amp;
				break;
			case  3: // s12-grid
				s12[i][j]=amp;
				break;
			case  4: // s22-grid
				s22[i][j]=amp;
				break;
		};
	}	
};
//	div and curl field output (2022/02/03)
void Fld2D :: snap_out_del(char *fname, double tout, Dom2D dom){
	int i,j,ndat[2];
	double x1[2],x2[2],sxy,vx,vy;
	FILE *fp=fopen(fname,"w");

	ndat[0]=dom.iYb[0]-dom.iYa[0]+1;
	ndat[1]=dom.iYb[1]-dom.iYa[1]+1;
	indx2cod(dom.iYa[0],dom.iYa[1],2,dom.Xa,dom.dx,x1);
	indx2cod(dom.iYb[0],dom.iYb[1],2,dom.Xa,dom.dx,x2);

	fprintf(fp,"#time=%lf\n",tout);
	fprintf(fp,"# Xa[0], Xb[0]\n");
		fprintf(fp,"%lf %lf\n",x1[0],x2[0]);
	fprintf(fp,"# Xa[1], Xb[1]\n");
		fprintf(fp,"%lf %lf\n",x1[1],x2[1]);
	fprintf(fp,"# Nx[0], Nx[1]\n");
		fprintf(fp,"%d %d\n",ndat[0],ndat[1]);
	fprintf(fp,"## ------ Divergence ------\n");
	for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
	for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
		vx=v1[i+1][j]-v1[i][j];
		vy=v2[i][j+1]-v2[i][j];
		fprintf(fp,"%lf %lf \n",vx,vy);	// div (v) components
	}
	}

	fprintf(fp,"## ------ Curl ------\n");
	for(i=dom.iYa[0]; i<dom.iYb[0]  ; i++){
	for(j=dom.iYa[1]; j<dom.iYb[1]  ; j++){
		vx=v1[i][j+1]-v1[i][j];
		vy=v2[i+1][j]-v2[i][j];
		fprintf(fp,"%lf %lf \n",vx,vy);	// curl (v) components
	}
	}
}
void Fld2D :: snap_out(int ityp, char *fname, double tout, Dom2D dom){

	int i,j,ndat[2];
	double x1[2],x2[2],sxy,vx,vy;
	FILE *fp=fopen(fname,"w");

	ndat[0]=dom.iYb[0]-dom.iYa[0]+1;
	ndat[1]=dom.iYb[1]-dom.iYa[1]+1;
	indx2cod(dom.iYa[0],dom.iYa[1],2,dom.Xa,dom.dx,x1);
	indx2cod(dom.iYb[0],dom.iYb[1],2,dom.Xa,dom.dx,x2);

	fprintf(fp,"#time=%lf\n",tout);
	fprintf(fp,"# Xa[0], Xb[0]\n");
		fprintf(fp,"%lf %lf\n",x1[0],x2[0]);
	fprintf(fp,"# Xa[1], Xb[1]\n");
		fprintf(fp,"%lf %lf\n",x1[1],x2[1]);
	fprintf(fp,"# Nx[0], Nx[1]\n");
		fprintf(fp,"%d %d\n",ndat[0],ndat[1]);
	fprintf(fp,"# amplitude\n");

	switch(ityp){
		case 0: 
		case 1: 
			for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
				vx=(v1[i][j]+v1[i+1][j])*0.5;
				vy=(v2[i][j]+v2[i][j+1])*0.5;
				fprintf(fp,"%lf %lf \n",vx,vy);
			}
			}
			break;
		case 2:
		case 3:
			for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
				sxy=s12[i][j]+s12[i+1][j]+s12[i][j+1]+s12[i+1][j+1];
				sxy*=0.25;
				fprintf(fp,"%lf %lf %lf\n",s11[i][j],s22[i][j],sxy);
			}
			}
			break;
		default :
			puts(" Invalid output field type specified!");
			puts(" --> process terminated.");
			exit(-1);
	}
	fclose(fp);

};
void Fld2D :: out(int ityp, char *fname, double tout, Dom2D dom){
	int i,j,ndat[2];
	double x1[2],x2[2];
	FILE *fp;
	fp=fopen(fname,"w");
//	gridNum(ityp);

	indx2cod(dom.iYa[0],dom.iYa[1],ityp,dom.Xa,dom.dx,x1);
	indx2cod(dom.iYb[0],dom.iYb[1],ityp,dom.Xa,dom.dx,x2);
	ndat[0]=dom.iYb[0]-dom.iYa[0]+1;
	ndat[1]=dom.iYb[1]-dom.iYa[1]+1;
	if(ityp==0){	// v1 field
		x1[0]-=0.5*dom.dx[0];
		x2[0]+=0.5*dom.dx[0];
		ndat[0]+=1;
	} 
	if(ityp==1){	// v2 field
		x1[1]-=0.5*dom.dx[1];
		x2[1]+=0.5*dom.dx[1];
		ndat[1]+=1;
	} 
	if(ityp==3){	// s12 field
		x1[0]-=0.5*dom.dx[0];
		x1[1]-=0.5*dom.dx[1];

		x2[0]+=0.5*dom.dx[0];
		x2[1]+=0.5*dom.dx[1];
		ndat[0]+=1;
		ndat[1]+=1;
	} 

	fprintf(fp,"#time=%lf\n",tout);
	fprintf(fp,"# Xa[0], Xb[0]\n");
	fprintf(fp,"%lf %lf\n",x1[0],x2[0]);
	fprintf(fp,"# Xa[1], Xb[1]\n");
	fprintf(fp,"%lf %lf\n",x1[1],x2[1]);
	fprintf(fp,"# Nx[0], Nx[1]\n");
	fprintf(fp,"%d %d\n",ndat[0],ndat[1]);
	fprintf(fp,"# amplitude\n");
	switch(ityp){
		case 0 :	// v1 field
			for(i=dom.iYa[0]; i<=dom.iYb[0]+1; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
				fprintf(fp,"%lf\n",v1[i][j]);
			}
			}
			break;
		case 1 : 	// v2 field
			for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]+1; j++){
				fprintf(fp,"%lf\n",v2[i][j]);
			}
			}
			break;
		case 2 : 	// s11,s22 field
		if(ityp ==2){
			for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
				fprintf(fp,"%lf\n",s11[i][j]);
			}
			}
		}else{
			for(i=dom.iYa[0]; i<=dom.iYb[0]  ; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]  ; j++){
				fprintf(fp,"%lf\n",s22[i][j]);
			}
			}
		}
			break;
		case 3 : 	// s12-field
			for(i=dom.iYa[0]; i<=dom.iYb[0]+1; i++){
			for(j=dom.iYa[1]; j<=dom.iYb[1]+1; j++){
				fprintf(fp,"%lf\n",s12[i][j]);
			}
			}
	};	
	fclose(fp);

}

//			Store and Export A-Scan waveforms 
void Fld2D :: out(int it, Recs rec, Dom2D dom){
	int i,j,k;	
	int ityp=rec.ityp;
	FILE *fp;
	char fname[128];
	double xcod[2];

	switch(ityp){
		case 0:	// v1
			for(k=0;k<rec.Ng;k++){
				l2ij(rec.krec[k],ityp,dom.Ndiv,&i,&j);
				rec.val[k][it]=v1[i][j];
			}
			break;
		case 1: // v2 
			for(k=0;k<rec.Ng;k++){
				l2ij(rec.krec[k],ityp,dom.Ndiv,&i,&j);
				rec.val[k][it]=v2[i][j];
			}
			break;
		case 2 : // s11
			for(k=0;k<rec.Ng;k++){
				l2ij(rec.krec[k],ityp,dom.Ndiv,&i,&j);
				rec.val[k][it]=s11[i][j];
			}
			break;
		case 3 : // s12
			for(k=0;k<rec.Ng;k++){
				l2ij(rec.krec[k],ityp,dom.Ndiv,&i,&j);
				rec.val[k][it]=s12[i][j];
			}
			break;
		case 4 : // s22
			for(k=0;k<rec.Ng;k++){
				l2ij(rec.krec[k],ityp,dom.Ndiv,&i,&j);
				rec.val[k][it]=s22[i][j];
			}
			break;
	};
	if(it==0){

		sprintf(fname,"rec%d.out",rec.ID);
		fp=fopen(fname,"w");

		fprintf(fp,"# dt, Nt\n");
		fprintf(fp,"%lf %d\n",rec.dt,rec.Nt);
		fprintf(fp,"# Ng, ityp(v1=0,v2=1,s11=2,s12=3,s22=3)\n");	
		fprintf(fp,"%d %d\n",rec.Ng,rec.ityp);
		fprintf(fp,"# xcod\n");	
		for(k=0;k<rec.Ng;k++){
			l2ij(rec.krec[k],abs(ityp),dom.Ndiv,&i,&j);
			indx2cod(i,j,abs(ityp),dom.Xa,dom.dx,xcod);
			fprintf(fp,"%lf %lf\n",xcod[0],xcod[1]);
		}	
		fprintf(fp,"# amp \n");
		fclose(fp);
	}
	if(it==rec.Nt-1){
		
		sprintf(fname,"rec%d.out",rec.ID);
		fp=fopen(fname,"a");
		for(k=0;k<rec.Ng;k++){
		for(i=0;i<rec.Nt;i++){
				fprintf(fp,"%lf\n",rec.val[k][i]);
		}
		}
		fclose(fp);
	}
}

void Fld2D :: outV1grid(Dom2D dom){

	int i,j,l;
	double xcod[2];
	FILE *fout=fopen("g1.out","w");

	gridNum(0);
	for(l=0;l<Nbnd1;l++){
		l2ij(abs(kbnd1[l]),Nx,&i,&j);
		indx2cod(i,j,0,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fprintf(fout,"\n");
	for(l=0;l<Nint1;l++){
		l2ij(kint1[l],Nx,&i,&j);
		indx2cod(i,j,0,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fclose(fout);
}
void Fld2D :: outV2grid(Dom2D dom){
	int i,j,l;
	double xcod[2];
	FILE *fout=fopen("g2.out","w");

	gridNum(1);
	for(l=0;l<Nbnd2;l++){
		l2ij(abs(kbnd2[l]),Nx,&i,&j);
		indx2cod(i,j,1,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fprintf(fout,"\n");
	for(l=0;l<Nint2;l++){
		l2ij(kint2[l],Nx,&i,&j);
		indx2cod(i,j,1,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fclose(fout);
}
void Fld2D :: outS12grid(Dom2D dom){
	int i,j,l;
	double xcod[2];
	FILE *fout=fopen("g12.out","w");

	gridNum(3);
	for(l=0;l<Nbnd12;l++){
		l2ij(abs(kbnd12[l]),Nx,&i,&j);
		indx2cod(i,j,3,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fprintf(fout,"\n");
	for(l=0;l<Nint12;l++){
		l2ij(kint12[l],Nx,&i,&j);
		indx2cod(i,j,3,dom.Xa,dom.dx,xcod);
		fprintf(fout,"%lf %lf\n",xcod[0],xcod[1]);	
	}
	fclose(fout);
}

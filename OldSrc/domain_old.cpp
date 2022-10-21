#include<iostream> 
#include<stdio.h>
#include<math.h> 
#include<string.h>
#include<stdlib.h>
#include "fdm2d.h"
using namespace std;


//  ---------CONSTRUCTOR -----------
Dom2D::Dom2D(char *fname){ //Contructor

	int i,j,k,ndim=2;
	double Ya[2],Yb[2];
	FILE *fp;
	char cbff[128],fkcell[128];
	fp=fopen(fname,"r");	

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fkcell);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf",xa,xa+1);
	fscanf(fp,"%lf %lf\n",xb,xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	int nchar_min=3;
	if(strlen(fkcell) > nchar_min){
		printf("read geometry from %s\n",fkcell);
		// read xa,xb,Ndiv & kcell from fkcell
		Dom2D::load_kcell(fkcell);
	};	

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Wa,Wa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Wb,Wb+1);

	int ndiv[2];
	ndiv[0]=Ndiv[0];
	ndiv[1]=Ndiv[1];
	for(i=0;i<ndim;i++){
		dx[i]=(xb[i]-xa[i])/Ndiv[i];
		Nx[i]=Ndiv[i]+1;
		nwa[i]=floor(Wa[i]/dx[i]+0.5);
		nwb[i]=floor(Wb[i]/dx[i]+0.5);
		Wa[i]=dx[i]*nwa[i];
		Wb[i]=dx[i]*nwb[i];

		Xa[i]=xa[i];
		Xb[i]=xb[i];
		Xa[i]=xa[i]-Wa[i];
		Xb[i]=xb[i]+Wb[i];
		Nx[i]+=(nwa[i]+nwb[i]);
		Ndiv[i]+=(nwa[i]+nwb[i]);
	}
	printf("nwa=%d %d\n",nwa[0],nwa[1]);
	printf("nwb=%d %d\n",nwb[0],nwb[1]);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf\n",&rho,&cL,&cT);
	cR=Rwave(cL,cT);
	amu=rho*cT*cT;
	almb=rho*(cL*cL-2.*cT*cT);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Ya,Ya+1);
	fscanf(fp,"%lf %lf\n",Yb,Yb+1);

	for(i=0;i<2;i++){
		iYa[i]=floor((Ya[i]-Xa[i])/dx[i]);
		iYb[i]=floor((Yb[i]-Xa[i])/dx[i]);
		if( iYa[i] <0) iYa[i]=0;
		if( iYb[i] >= Ndiv[i]) iYb[i]=Ndiv[i]-1;
	}

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %d\n",&tout_s,&tout_e,&Nout);

	fclose(fp);

	//mem_alloc();	// allocate & initialize kcell[i][j]
	ptmp=(long int *)malloc(sizeof(long int)*Ndiv[0]*Ndiv[1]);
	kcell=(long int **)malloc(sizeof(long int*)*Ndiv[0]);
	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}

	if(strlen(fkcell) > nchar_min){
		fp=fopen(fkcell,"r");
		for(i=0;i<5;i++){
		       	fgets(cbff,128,fp);
		}
		for(i=0;i<ndiv[0];i++){
		for(j=0;j<ndiv[1];j++){
			fscanf(fp,"%d\n",&k);
			kcell[i+nwa[0]][j+nwa[1]]=k;
		}
		}
		fclose(fp);
	}
};

void Dom2D::load_kcell(char *fname){
	FILE *fp;
	char cbff[128];

	fp=fopen(fname,"r");
	if(fp==NULL){
		printf("Can't find file %s\n",fname);
		printf("--> abort procss..");
		exit(-1);
	};

	xa[0]=0.0; 
	xa[1]=0.0;
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",xb,xb+1);

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Ndiv,Ndiv+1);

	printf("xa=%lf %lf\n",xa[0],xa[1]);
	printf("xb=%lf %lf\n",xb[0],xb[1]);
	printf("Ndiv=%d %d\n",Ndiv[0],Ndiv[1]);

	fgets(cbff,128,fp);
/*
	int i,j,k;
	ptmp=(long int *)malloc(sizeof(long int)*Ndiv[0]*Ndiv[1]);
	kcell=(long int **)malloc(sizeof(long int*)*Ndiv[0]);
	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fscanf(fp,"%d\n",&k);
		kcell[i][j]=k;
	}
	}
*/
	fclose(fp);
	
};
//  ----------- MEMORY ALLOCATION --------------
void Dom2D::mem_alloc(){
	int i,j;
	//long int *ptmp;
	ptmp=(long int *)malloc(sizeof(long int)*Ndiv[0]*Ndiv[1]);
	kcell=(long int **)malloc(sizeof(long int*)*Ndiv[0]);

	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kcell[i][j]=0;
	}}
};

//  ----------- DOMAIN PERFORATION  --------------
void Dom2D::perfo_ellip(char *fname){
	int i,j;
	double xcod[2],xc[2],a,b;
	double xx,yy;

	FILE *fp;
	char cbff[8];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo_ellip !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,8,fp) !=NULL){
		if(strcmp(cbff,"##Ellip")==0){
			fscanf(fp,"%lf %lf %lf %lf\n",xc,xc+1,&a,&b);
			printf("%lf %lf %lf %lf\n",xc[0],xc[1],a,b);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				xx=(xcod[0]-xc[0])/a;
				yy=(xcod[1]-xc[1])/b;
				if(xx*xx+yy*yy <1.0) kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
};
//  ----------- DOMAIN PERFORATION  --------------
void Dom2D::Cut(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[7];
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo !");
		puts(fname);
		exit(-1);
	}
	double xs[2];	// start point
	double XL;	// slit length
	double wd;	// half width 
	double alph;	// angle in deg
	double xe[2];	// end point
	double PI=atan(1.0)*4.0;
	double rd;

	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Wire")==0){
			fscanf(fp,"%lf %lf \n",xs,xs+1);
			fscanf(fp,"%lf %lf %lf \n",&XL, &wd, &alph);
			rd=wd*0.5;
			//printf("xs=%lf %lf\n",xs[0],xs[1]);
			//printf("Xw=%lf alph=%lf\n",Xw,alph);
			//exit(-1);
			alph=alph/180.*PI;
			xe[0]=xs[0]+XL*cos(alph);
			xe[1]=xs[1]+XL*sin(alph);

			Dom2D::WireCut(xs,xe,rd);
		}
	};
};
void Dom2D::WireCut(double xc1[2], double xc2[2],double rd){
	//double xc1[2]={11.0, 12.0};
	//double xc2[2]={13.0,  8.0};
	printf("xc1=%lf %lf\n",xc1[0],xc1[1]);
	printf("xc2=%lf %lf\n",xc2[0],xc2[1]);
	printf("rd=%lf\n",rd);
	//double rd=0.1;
	int Nstep;
	int i,j,k;

	double xc[2],dxc[2];

	double rx=xc2[0]-xc1[0];
	double ry=xc2[1]-xc1[1];
	double rr=sqrt(rx*rx+ry*ry);
	Nstep=int(rr/dx[0]);
	dxc[0]=rx/(Nstep-1);
	dxc[1]=ry/(Nstep-1);


	for(i=0;i<Nstep;i++){
		xc[0]=xc1[0]+dxc[0]*i;
		xc[1]=xc1[1]+dxc[1]*i;
		Dom2D::perfo_tight(xc,rd);
	};

};
void Dom2D::perfo_tight(double xc[2],double rd){
	int i,j;
	double xcod[2];
	bool io;
	Circ cdat;
	cdat.xc[0]=xc[0];
	cdat.xc[1]=xc[1];
	cdat.radi=rd;


	int i1=int((xc[0]-rd-Xa[0])/dx[0])-1;
	int j1=int((xc[1]-rd-Xa[1])/dx[1])-1;
	int i2=int((xc[0]+rd-Xa[0])/dx[0])+1;
	int j2=int((xc[1]+rd-Xa[1])/dx[1])+1;

	if(i1<0) i1=0;
	if(j1<0) j1=0;
	if(i2>=Ndiv[0]-1) i2=Ndiv[0]-1;
	if(j2>=Ndiv[1]-1) j2=Ndiv[1]-1;

	for(i=i1; i<=i2; i++){
		xcod[0]=Xa[0]+dx[0]*(i+0.5);	
	for(j=j1; j<=j2; j++){
		xcod[1]=Xa[1]+dx[1]*(j+0.5);	
		io=cdat.isin(xcod);	
		if(io==true) kcell[i][j]=1;
	}
	}

};
//  ----------- DOMAIN PERFORATION  --------------
void Dom2D::perfo(char *fname){
	int i,j;
	double xcod[2];
	bool io;
	Circ cdat;

	FILE *fp;
	char cbff[7];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::perfo !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Circ")==0){
			fscanf(fp,"%lf %lf %lf\n",cdat.xc,cdat.xc+1, &cdat.radi);

			for(i=0;i<Ndiv[0];i++){
				xcod[0]=Xa[0]+dx[0]*(i+0.5);	
			for(j=0;j<Ndiv[1];j++){
				xcod[1]=Xa[1]+dx[1]*(j+0.5);	
				io=cdat.isin(xcod);	
				if(io==true){
				 kcell[i][j]=1;
				}
			}
			}
		}
	};
	fclose(fp);
};
void Dom2D::polygon(char *fname){
	FILE *fp;
	char cbff[7];
	int i,j,k,np,iin;
	double *x,*y;
	double x1[2],x2[2],tx[2],nx[1],xf[2];
	double xmin,ymin,xmax,ymax,r12,h;
	double tmp;
	int imin,imax,jmin,jmax;


	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::polygon !");
		puts(fname);
		puts("--> process terminated");
		exit(-1);
	}
			x=(double *)malloc(sizeof(double)*10);
			y=(double *)malloc(sizeof(double)*10);


	i=0;
	printf("Ndiv=%d %d\n",Ndiv[0],Ndiv[1]);
	while(fgets(cbff,7,fp)!=NULL){
		if(strcmp(cbff,"##Poly")==0){
			printf("Flag ##Poly found !!\n");
			fscanf(fp,"%d\n",&np);	// number vertices
			printf("np=%d\n",np);
			//x=(double *)malloc(sizeof(double)*np);
			//y=(double *)malloc(sizeof(double)*np);

			for(int jp=0;jp<np;jp++){
				fscanf(fp,"%lf %lf\n",x+jp,y+jp);
				if(jp==0){
				       	xmin=x[0]; ymin=y[0];
				       	xmax=x[0]; ymax=y[0];
				}else{
					if(x[jp]< xmin) xmin=x[jp];
					if(x[jp]> xmax) xmax=x[jp];
					if(y[jp]< ymin) ymin=y[jp];
					if(y[jp]> ymax) ymax=y[jp];
				};
			}
			printf("xmin=%lf, xmax=%lf\n",xmin,xmax);
			printf("ymin=%lf, ymax=%lf\n",ymin,ymax);

			imin=int(floor((xmin-Xa[0])/dx[0]));
			imax=floor((xmax-Xa[0])/dx[0]);
			jmin=int(floor((ymin-Xa[1])/dx[1]));
			jmax=floor((ymax-Xa[1])/dx[1]);

			if(imin <0) imin=0;
			if(jmin <0) jmin=0;
			if(imin>=Ndiv[0]) imin=Ndiv[0]-1;
			if(jmin>=Ndiv[1]) jmin=Ndiv[1]-1;

			if(imax>=Ndiv[0]) imax=Ndiv[0]-1;
			if(jmax>=Ndiv[1]) jmax=Ndiv[1]-1;
			printf("imin,jmin=%d %d\n",imin,jmin);
			printf("imax,jmax=%d %d\n",imax,jmax);
			printf("Ndiv=%d %d\n",Ndiv[0],Ndiv[1]);

			//for(i=imin; i<=imax; i++){
			//for(i=imin; i<=imax; i++){
			for(i=0; i<Ndiv[0]; i++){
				xf[0]=Xa[0]+dx[0]*(i+0.5);
				//fflush(stdout);
			//for(j=jmin; j<=jmax; j++){
			//for(j=jmin; j<=jmax; j++){
			for(j=0; j<Ndiv[1]; j++){
				xf[1]=Xa[1]+dx[1]*(j+0.5);
				iin=0;
				for(k=0;k<np;k++){
					x1[0]=x[k%np];
					x2[0]=x[(k+1)%np];
					x1[1]=y[k%np];
					x2[1]=y[(k+1)%np];

					tx[0]=x2[0]-x1[0];
					tx[1]=x2[1]-x1[1];
					r12=sqrt(tx[0]*tx[0]+tx[1]*tx[1]);

					tx[0]/=r12;
					tx[1]/=r12;

					nx[0]=-tx[1];
					nx[1]= tx[0];

					h=(xf[0]-x1[0])*nx[0]+(xf[1]-x1[1])*nx[1];
					if(h>0.e0) iin++;
				}

				//if(j<0) puts("j!");
				//if(i>=Ndiv[0]) puts("i!!");
				//if(j>=Ndiv[1]) puts("j!!");
				if(iin==np){
					Dom2D::kcell[i][j]=1;
			//		puts("Interior point !");
				};
			}
			//printf("isum=%d\n",isum);
			}
		}
	};

};
void Dom2D::perfo_para(char *fname){
	double ew[2],eh[2];  // unit vectors (basis)
	double ewh[2],ehh[2];  // unit vectors (basis)
	double xf[2],x0[2];	// lower left corner point
	double W,H,alph;
	char cbff[7];
	double xw,xh;
	FILE *fp;
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::slit !");
		puts(fname);
		exit(-1);
	}


	double det;
	int i,j;
	double pi=4.*atan(1.0);
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Para")==0){
			fscanf(fp,"%lf %lf\n",x0,x0+1);
			fscanf(fp,"%lf %lf %lf\n",&W,&H,&alph);
			ew[0]=W;
			ew[1]=0.0;
			alph=alph/180.*pi;
			eh[0]=H*sin(alph);
			eh[1]=H*cos(alph);
			det=ew[0]*eh[1]-ew[1]+eh[0];
			ewh[0]=eh[1]/det;
			ewh[1]=-eh[0]/det;

			ehh[0]=-ew[1]/det;
			ehh[1]= ew[0]/det;

			for(i=0;i<Ndiv[0];i++){
				xf[0]=Xa[0]+dx[0]*(i+0.5)-x0[0];
			for(j=0;j<Ndiv[1];j++){
				xf[1]=Xa[1]+dx[1]*(j+0.5)-x0[1];

				xw=xf[0]*ewh[0]+xf[1]*ewh[1];
				xh=xf[0]*ehh[0]+xf[1]*ehh[1];
				if(xw<0.0) continue;
				if(xh<0.0) continue;
				if(xw>1.0) continue;
				if(xh>1.0) continue;
				kcell[i][j]=1;
			}
			}
		}
	};
}
void Dom2D::angled_slit(char *fname){
	double xc[2],e1[2],e2[2],xf[2];
	double a,b,alph,x1,x2;
	double xs[4],ys[4];
	double xmin,xmax,ymin,ymax;
	double pi=4.0*atan(1.0);
	int i,j,i1,i2;
	int imin,imax,jmin,jmax;
	FILE *fp;
	char cbff[7];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::slit !");
		puts(fname);
		exit(-1);
	}

	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Slit")==0){
			fscanf(fp,"%lf %lf\n",xc,xc+1);
			fscanf(fp,"%lf %lf %lf\n",&a,&b,&alph);
			alph=alph/180.*pi;

			e1[0]=cos(alph);
			e1[1]=sin(alph);

			e2[0]=-e1[1];
			e2[1]=e1[0];

			xs[0]=xc[0]+a*e1[0]+b*e2[0];
			ys[0]=xc[1]+a*e1[1]+b*e2[1];

			xs[1]=xc[0]-a*e1[0]+b*e2[0];
			ys[1]=xc[1]-a*e1[1]+b*e2[1];

			xs[2]=xc[0]-a*e1[0]-b*e2[0];
			ys[2]=xc[1]-a*e1[1]-b*e2[1];

			xs[3]=xc[0]+a*e1[0]-b*e2[0];
			ys[3]=xc[1]+a*e1[1]-b*e2[1];

			xmin=xs[0]; xmax=xs[0]; 
			ymin=ys[0]; ymax=ys[0]; 
			for(i=1;i<4;i++){
				if(xs[i]<xmin) xmin=xs[i];
				if(xs[i]>xmax) xmax=xs[i];
				if(ys[i]<ymin) ymin=ys[i];
				if(ys[i]>ymax) ymax=ys[i];
			}

			
			imin=floor((xmin-Xa[0])/dx[0]);
			imax=floor((xmax-Xa[0])/dx[0]);
			jmin=floor((ymin-Xa[1])/dx[1]);
			jmax=floor((ymax-Xa[1])/dx[1]);

			if(imin <0) imin=0;
			if(jmin <0) jmin=0;
			if(imax>=Ndiv[0]) imax=Ndiv[0]-1;
			if(jmax>=Ndiv[1]) jmax=Ndiv[1]-1;

			for(i=imin; i<=imax; i++){
				xf[0]=Xa[0]+dx[0]*i-xc[0];
			for(j=jmin; j<=jmax; j++){
				xf[1]=Xa[1]+dx[1]*j-xc[1];


				x1=xf[0]*e1[0]+xf[1]*e1[1];
				x2=xf[0]*e2[0]+xf[1]*e2[1];
				i1=0; i2=0;
				if(fabs(x1)<a) i1=1;
				if(fabs(x2)<b) i2=1;
				if(i1+i2==2) kcell[i][j]=1;
			}
			}
/*
			for(i=0;i<2;i++){
				ixs[i]=floor((xs[i]-Xa[i])/dx[i]);
				ixe[i]=floor((xe[i]-Xa[i])/dx[i]);
				if( ixs[i] <0) ixs[i]=0;
				if( ixe[i] >= Ndiv[i]) ixe[i]=Ndiv[i]-1;
			}

			for(i=ixs[0]; i<=ixe[0]; i++){
			for(j=ixs[1]; j<=ixe[1]; j++){
				 kcell[i][j]=1;
			}
			}
*/
		}
	};
	

}
void Dom2D::slit(char *fname){
	int i,j;
	double xs[2],xe[2];
	int ixs[2],ixe[2];

	FILE *fp;
	char cbff[7];
	fp=fopen(fname,"r");
	if(fp==NULL){
		puts("Can't open file from Dom2D::slit !");
		puts(fname);
		exit(-1);
	}
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Rect")==0){
			printf("Flag ##Rect found !!\n");
			fscanf(fp,"%lf %lf\n",xs,xs+1);
			fscanf(fp,"%lf %lf\n",xe,xe+1);

			for(i=0;i<2;i++){
				ixs[i]=floor((xs[i]-Xa[i])/dx[i]);
				ixe[i]=floor((xe[i]-Xa[i])/dx[i]);
				if( ixs[i] <0) ixs[i]=0;
				if( ixe[i] >= Ndiv[i]) ixe[i]=Ndiv[i]-1;
			}

			for(i=ixs[0]; i<=ixe[0]; i++){
			for(j=ixs[1]; j<=ixe[1]; j++){
				 kcell[i][j]=1;
			}
			}
		}
	};
	fclose(fp);
}
//-------------EXPORT GEOMETRY DATA  ------------------
void Dom2D :: out_kcell(){

	FILE *fp=fopen("kcell.dat","w");
	int i,j;
	double Ya[2],Yb[2];


	fprintf(fp,"# xa[2], xb[2] (computational domian)\n");
	fprintf(fp," %lf %lf\n %lf %lf\n",xa[0],xa[1],xb[0],xb[1]);

	fprintf(fp,"# Xa[2], Xb[2] (physical domain) \n");
	fprintf(fp," %lf %lf\n %lf %lf\n",Xa[0],Xa[1],Xb[0],Xb[1]);

	fprintf(fp,"# Ya[2], Yb[2] (imaging area) \n");
	for(i=0;i<2;i++){
		Ya[i]=Xa[i]+dx[i]*(iYa[i]+0.5);
		Yb[i]=Xa[i]+dx[i]*(iYb[i]+0.5);
	}
	fprintf(fp,"%lf %lf\n %lf %lf\n",Ya[0],Ya[1],Yb[0],Yb[1]);

	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# kcell[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp, "%ld\n",kcell[i][j]);
	}
	}

	fflush(fp);
};

//-------------COURANT NUMBER  ------------------
void Dom2D :: CFL(double dt){

	cfl0=cL*dt/dx[0]*sqrt(2.0);
	cfl1=cT*dt/dx[0]*sqrt(2.0);
	printf("CFL(L-wave)=%lf\n",cfl0);
	printf("   (T-wave)=%lf\n",cfl1);
	printf("cL=%lf, cT=%lf\n",cL,cT);
	printf("dt=%f dx=(%f, %f)\n",dt,dx[0],dx[1]);
}

//------------- GRID NUMBER  ------------------

//	SET NUMBER OF GRIDS OF TYPE itype to Nx[2]
//      (identical with the function Fld2D:: gridNum) 
void Dom2D ::gridNum(int ityp){
	switch(ityp){
	case 0:	// v1-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1];
		break;
	case 1:	// v2-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1]+1;
		break;
	case 2:	// s11
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	case 3:	// s12-grid
		Nx[0]=Ndiv[0]+1;
		Nx[1]=Ndiv[1]+1;
		break;
	case 4: //s22-grid
		Nx[0]=Ndiv[0];
		Nx[1]=Ndiv[1];
		break;
	};
	Ng=Nx[0]*Nx[1];
};

//------------- Circ CLASS MEMBER FUNCTION  ------------------
bool Circ :: isin(double *x){
	double dist;

	dist=(x[0]-xc[0])*(x[0]-xc[0]);
	dist+=(x[1]-xc[1])*(x[1]-xc[1]);
	dist=sqrt(dist);

	if(dist <= radi){
		return true;
	}{
		return false;
	}
};

//------------- Rayleigh Wave Phase Velocity --------------
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


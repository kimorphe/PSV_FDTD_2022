#include<stdio.h>
#include<stdlib.h>


class Kcell{
	public:
		int **kcell;
		int Nx,Ny;
		double Wd;
		double Ht;
		void load(char fname[128]);
		void load2(char fname[128]);
		int is_Int(int i, int j);
		int is_Bnd(int i, int j);
		int is_Bnd4(int i, int j);
		int Walk(int i, int j);
		int mdir;
		double xa[2],xb[2];
		double Xa[2],Xb[2];
		double Ya[2],Yb[2];
		double dx[2];
		double xcod(int i);
		double ycod(int j);
	private:
};

double Kcell::xcod(int i){
	return(Xa[0]+dx[0]*(i+0.5));
};
double Kcell::ycod(int j){
	return(Xa[1]+dx[1]*(j+0.5));
};
void Kcell::load2(char fname[128]){
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("Can't find file %s\n",fname);
		exit(-1);
	};
	int Ndiv[2];
	fgets(cbff,128,fp);
	puts(cbff);
	fscanf(fp,"%lf %lf\n",xa,xa+1);
	fscanf(fp,"%lf %lf\n",xb,xb+1);
	fgets(cbff,128,fp);
	puts(cbff);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	puts(cbff);
	fscanf(fp,"%lf %lf\n",Ya,Ya+1);
	fscanf(fp,"%lf %lf\n",Yb,Yb+1);
	fgets(cbff,128,fp);
	puts(cbff);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);
	fgets(cbff,128,fp);
	puts(cbff);

	dx[0]=(Xb[0]-Xa[0])/Ndiv[0];
	dx[1]=(Xb[1]-Xa[1])/Ndiv[1];
	printf("dx,dy=%lf %lf\n",dx[0],dx[1]);

	Nx=Ndiv[0];
	Ny=Ndiv[1];
	kcell=(int **)malloc(sizeof(int *)*Nx);
	int *pt=(int *)malloc(sizeof(int)*Nx*Ny);
	int i,j,itmp,n;
	for(i=0;i<Nx;i++) kcell[i]=pt+i*Ny;

	n=0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		fscanf(fp,"%d\n",&itmp);
       		kcell[i][j]=itmp+1;	// 1:interior, 2: exterior
		n++;
	}
	}
	printf("n=%d\n",n);
	mdir=0;
	fclose(fp);
};

void Kcell::load(char fname[128]){
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("Can't find file %s\n",fname);
		exit(-1);
	};
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",&Wd,&Ht);
	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",&Nx,&Ny);
	fgets(cbff,128,fp);

	kcell=(int **)malloc(sizeof(int *)*Nx);
	int *pt=(int *)malloc(sizeof(int)*Nx*Ny);
	int i,j,itmp,n;
	for(i=0;i<Nx;i++) kcell[i]=pt+i*Ny;

	n=0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		fscanf(fp,"%d\n",&itmp);
       		kcell[i][j]=itmp+1;	// 1:interior, 2: exterior
		n++;
	}
	}
	mdir=0;
	fclose(fp);
};
int Kcell::is_Int(int i, int j){
	int typ=0;	// set to false (saying "No., this is not interior cell") 
	if(i<0) return(typ);
	if(j<0) return(typ);
	if(i>=Nx) return(typ);
	if(j>=Ny) return(typ);
	if(abs(kcell[i][j])==2) return(typ); // 1:interior, 2:exterior cell
	typ=1;		// set to true (saying "Yes, this is an interior cell")
	return(typ);
};
int Kcell::is_Bnd4(int i, int j){
	int I,J;
	int ii[4]={ 0, 1, 0,-1};
	int jj[4]={-1, 0, 1, 0};
	int m,md;
	int typ=0;	// non boundary cell
	if(is_Int(i,j)==0) return(typ);
	for(m=0;m<4;m++){
		I=ii[m]+i;
		J=jj[m]+j;
		md=Kcell::is_Int(I,J);
		if(md==0){
			typ=1;
			return(typ);
		};
	};
	return(typ);
};
int Kcell::is_Bnd(int i, int j){
	int I,J;
	int ii[8]={ 0,  1, 1, 1, 0, -1, -1,-1};
	int jj[8]={-1, -1, 0, 1, 1,  1,  0,-1};
	int m,md;
	int typ=0;	// non boundary cell
	if(is_Int(i,j)==0) return(typ);
	for(m=0;m<8;m++){
		I=ii[m]+i;
		J=jj[m]+j;
		md=Kcell::is_Int(I,J);
		if(md==0){
			typ=1;
			return(typ);
		};
	};
	return(typ);
};
int Kcell::Walk(int i, int j){
	// typ= 0(interior) , -1(out), 1,2,3,4 (bnd)
	// typ
	//
	bool intr=true;
	if(i<0) intr=false;
	if(j<0) intr=false;
	if(i>=Nx) intr=false;
	if(j>=Ny) intr=false;
	
	if(intr != true){
		puts("exterior point given to Walk");
		exit(-1);
	};

	int k,l,m;
	int I,J;
	int id=0;
	int ii[4]={ 0, 1, 0,-1};
	int jj[4]={-1, 0, 1, 0};
	for(m=0;m<4;m++){
		I=ii[m]+i;
		J=jj[m]+j;
		l=Kcell::is_Int(I,J);
		if(l==0) return(m);
	};
	return(-1);

};

int main(){
	//char fname[128]="bead_kcell.dat";
	char fname[128]="kcell.dat";

	FILE *fp=fopen("bnd.out","w");
	FILE *fp2=fopen("cell.out","w");
	Kcell KCL;

	KCL.load2(fname);
	//KCL.load(fname);

	int ii[8]={ 0,  1, 1, 1, 0, -1, -1,-1};
	int jj[8]={-1, -1, 0, 1, 1,  1,  0,-1};

	int m,i1,j1,i;
	i1=0; j1=10;
	int Nstep=KCL.Nx*20;
	int itry,jtry,ibnd;
	double xcod,ycod;
	m=0;
	for(i=0; i<Nstep; i++){

		itry=i1+ii[m];
		jtry=j1+jj[m];
		if(KCL.is_Int(itry,jtry)==0){
			m++;
			m=m%8;
			continue;
		}
		if(KCL.kcell[itry][jtry]<=0){
			m++;
			m=m%8;
			continue;
		};

		ibnd=KCL.is_Bnd4(itry,jtry);
		if(ibnd==0){
			m++;
			m=m%8;
			continue;
		}

		KCL.kcell[i1][j1]=-abs(KCL.kcell[i1][j1]);
		i1=itry;
		j1=jtry;
		//printf("%d, %d \n",i1,j1);
		xcod=KCL.xcod(i1);
		ycod=KCL.ycod(j1);
		fprintf(fp,"%lf, %lf \n",xcod,ycod);
	};
	/*
	int j,kval;
	for(i=0;i<KCL.Nx;i++){
	for(j=0;j<KCL.Ny;j++){
		kval=KCL.kcell[i][j];
		if(kval <0) fprintf(fp2,"%d, %d\n",i,j);
	}
	}
	*/
	return(0);

};

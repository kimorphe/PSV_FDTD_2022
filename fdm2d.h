void mem_alloc2D(int nx, int ny, double ***pt);

//	********************************************
//
//		DOMAIN CLASS
//
//	********************************************
//-------------Circ Class--------------------
class Circ{
	public:
		double xc[2];	// center 
		double radi;	// radius
		bool isin(double* x);
	private:
};

//-------------Dom2D Class ------------------
class Dom2D{
	public:
		double Xa[2],Xb[2],dx[2];// Physical Domain
		int iYa[2],iYb[2];	 // Snapshot Window(2d-index)
		double tout_s,tout_e; int Nout;  // timestep for output
		double Wa[2], Wb[2];	 // PML Size
		double xa[2],xb[2];	 // Computational Domain
		double rho,cT,cL,cR,almb,amu; //density & stiffness constants
		double cfl0,cfl1;	//Counrant Numbers
		int Nx[2],Ndiv[2],nwa[2],nwb[2];
		int Ng;	// number of grids of a specified type
		long int *ptmp;
		long int **kcell;
		void perfo(char *fn);
		void perfo_ellip(char *fn);
		void perfo_tight(double xc[2],double rd);
		void perfo_para(char *fn);
		void slit(char *fn);
		void Cut(char *fn);
		void WireCut(double xs[2], double xe[2], double rd);
		void angled_slit(char *fn);
		void polygon(char *fn);
		void out_kcell();
		Dom2D(char *fname);
		void CFL(double dt);
		void gridNum(int ityp);
		void load_kcell(char *fn);
	private:
		void mem_alloc();
};

//	********************************************
//
//		SOURCE and INCIDENT WAVE CLASS
//
//	********************************************


class Cmplx{	// Complex Number Class
	public:
		double re,im;
		double Abs();
		Cmplx(double a, double b);
		Cmplx();
		Cmplx Cnjg();
		void disp();
		void set(double a, double b);
		Cmplx operator+(Cmplx z);
		Cmplx operator-(Cmplx z);
		Cmplx operator*(Cmplx z);
		Cmplx operator/(Cmplx z);
		Cmplx operator*(double z);
		void operator=(double z);
	private:
};
Cmplx exp(Cmplx z);
class InWv{
	public :
		double t1,t2,dt,T0;
		double *amp;
		int Nt,wvtyp,nbrst,nwv;
		int nsig; //narrowness factor used in Gaussian amplitude modulation
		InWv(char *);	// constructor (data read from file char*) 
		InWv(int Nt);	// 
		void disp();
		void gen_sin();		
		void gen_cos();		
		void gen_wvfm();
		void out(char *);
		double df,fmax;
		Cmplx *Amp;
		void DFT();
		void DFTout(char *);
		void Amod(double tb, int nsig);
	private:
		void mem_alloc();
};		
/*
class InWv{
	public :
		double t1,t2,dt,T0;
		double *amp;
		int Nt,wvtyp,nbrst,nwv;
		InWv(char *);
		void disp();
		void gen_sin();		
		void gen_cos();		
		void gen_wvfm();
	private:
		void mem_alloc();
		
};
*/
class Src{
	public:
		double Xsa[2],Xsb[2];	// corner point
		InWv *inwvs;		// waveform class array 
		int *ksrc, *wvID; 	// 1D grid index, waveform IDs --> inwv**.dat
		int ngrd,nwv,inmb;	// number of grids, waveforms. Waveform asigning sequence
		int ityp;		// source grid type (0=v1, 1=v2, 3=s12)
		int ifld;		// field type (0:stress, 1:velocity)
		int Nt; double dt;	// total time steps, time increment
		double ain;		// incident angle measured couter-clockwise from x1-axis
		double xf[2];		// focal point
		int dly_seq;		// 0:steering, 1: focusing, -1: uniform delay 
		int is_sync;		// synchrhonize array 0: Yes, 1: No
		double delay;		// delay time used when array elements are synchronized
		int *idly;		// delay steps of excitation signal 
		int wvtyp;		// incident wave type (0:P, 1:SV, 2:R)
		double cin;		// incident wave vel.(0: P, 1:SV-wave, 2:R-wave)
		Src(int ii);		// constructor(source type)
		Src();			// empty default constructor
		void gen_Src(char *fname); // set source parameters reading data from *fname
		void gen_Grid(Dom2D dom); // generate source grids 
		double set_Wvfm();	// generate numerical waveform data 
		void isType();		// show source type
		void out(int isrc, int *Nx, Dom2D dom); // print source parameters to "source.out" 
		void set_Delay(int *Nx, Dom2D dom); // set phase delay
	private:
};
int set_srcs(
	char *fname, 	// data file "src.inp"
	Dom2D dom, 	// instance of Dom2D class 
	Src **src	// pointer to source class array
);
int sync_src(
	char *fname,  // data file "src.inp"
	Dom2D dom,Src // Domain data
	**src	//pointer to source class array
);

//	********************************************
//
//		PML CLASS
//
//	********************************************


class PML{
	public:
		double Xa[2],Xb[2];
		double dx[2]; //signed mesh size 
		int ia[2],ib[2],nx[2];
		int ID;		// PML ID {0(left),1(right),2(bottom),3(top)}
		int Nx[2];	// use 'gridNum(ityp)' to fill 
		double *dmpx;	//damping coefficient
		double *dmpy;	//damping coefficient
		double **v1p, **v1v, **v1;
		double **v2p, **v2v, **v2;
		double **s11p, **s11v, **s11;	
		double **s22p, **s22v, **s22;	
		double **s12p, **s12v, **s12;	
		double s110,s220;	// initial stress
//		PML( Dom2D dom, double a, double b);
		void setup(Dom2D dom, char *ficon);
		void mem_alloc();
		void s2v(double rho, double dt);
		void v2s(double amu, double almb, double dt);
		double gmm;
		int mexp;
		int set_IC(Dom2D dom, char *ficon);
		void gridNum(int ityp);
	private:
};

//	********************************************
//
//		RECEIVER CLASS	
//
//	********************************************

class Recs{
	public:
		int idir;	// alignment (0: horiz. 1:vert)
		int ityp;	// grid type(0:v3, 1: s31, 2: s32)
		int Ng;		// number of receiver grids 
		int ID;		// receiver group ID
		int *krec;	// 1D index set
		double **val;	// time-history data val[Ng][Nt]
		double X1[2],X2[2]; // end points
		int Nt;	//number of time steps
		double dt;	//time increment
		void setup(
			int jdir, // alignment 0:horizontal, 1: vertical
			int jtyp, // grid type 0:v3, 1:s31, 2;s32
			double *Y1, // End point 1
			double *Y2,  // End point 2
			int npnt, 	// number of points
			Dom2D dom, 	// domain data
			int Nt, 	// number of time steps
			double dt	// time increment
		);
		void print(Dom2D dom);
	private:
};

int set_recs(char *fname, Dom2D dom, Recs **rec, int Nt, double dt);
//--------------------------------------------------------------
double plane_wave(
	int ig, int jg,	// 2D grid index 	
	int igrd,	// grid type (0: v1,1:v2, 2:s11,3:s12,4:s22)
	int it,		// time step
	int ityp,	// wave type (L-wave=0, T-wave=1)
	double ain,	// incident angle [deg]
	double *xp0,	// zero phase point
	InWv inwv,	// wave data(class)
	Dom2D dom,	// domain data
	double *ginx, 	// incident field || (para)
	double *giny	// incident field _|_ (vert)
);

//----------- Fld2D Class---------------------
class Fld2D{
	public:
		double **v1,**v2;	 // velocity fields
		double **s11,**s22,**s12;// stress fields
		double s110,s220;	// pre-stress
		int *kint1, *kbnd1; // 1D index for v1-grid walk through
		int *kint2, *kbnd2; // 1D index for v2-grid walk through
		int *kint12, *kbnd12; // 1D index for s12-grid walk through
		int Nint1, Nbnd1; // number of v31-grids
		int Nint2, Nbnd2; // number of v32-grids
		int Nint12, Nbnd12; // number of s12-grids
		int npml;	//number of PML patches
		PML *pml;	// PML Class Array
		Fld2D(Dom2D dom, char *ficon);    // Constructor ( N2[2]:number of cells)
		int Ndiv[2];	// Number of cells on x and y axes
		int Nx[2];	// Number of grids on x and y axes
		int Ng;		//total number of FD cells
		void gridNum(int ii); //Set grid number of type 'ii' to Nx[2] 
		void gen_indx1(Dom2D dom); // generate 1D index sets: kint1 & kbnd1 
		void gen_indx2(Dom2D dom); // generate 1D index sets: kint2 & kbnd2
		void gen_indx12(Dom2D dom); // generate 1D index sets: kint12 & kbnd12
		void s2v(Dom2D dom, double dt );
		void v2s(Dom2D dom, double dt);
		void apply_src(int it,Src src,Dom2D dom);
		void out(int ityp, char *fname, double tout, Dom2D dom);
		void out(int it, Recs rec, Dom2D dom);
		void snap_out(int ityp, char *fname, double tout, Dom2D dom);
		void snap_out_del(char *fname, double tout, Dom2D dom);
		void outV1grid(Dom2D dom);
		void outV2grid(Dom2D dom);
		void outS12grid(Dom2D dom);
		int set_IC(Dom2D dom, char *ficon);
		int del_src_grid(Src src);
	private:
		void mem_alloc();
};

//	********************************************
//
//		VARIOUS INDEX OPERATOINS
//
//	********************************************
int ij2l(
	int i, int j, 
	int *Nx
);
int ij2l(
	int i, int j, 
	int ityp,	//grid type 
	int *Ndiv	//number of cells
);

void l2ij(
	int l, 
	int *Nx,
	int *i, int *j //return value
);
void l2ij(
	int l, 
	int ityp, 	//grid type
	int *Ndiv,	// number of cells
	int *i, int *j //return value
);

void Ndiv2Nx(
	int *Ndiv,
	int ityp, 
	int *Nx		//return value
);


void indx2cod(
	int i, int j, 	// 2D index 
	int ityp, 	// grid type
	double *Xa, 	// Lower left point
	double *dx, 	// grid spacing
	double *xcod	// coordinate (return value)
);

void indx_ofst(
	int ityp,
	int *i0
);
void cod2indx(
	double *xcod, 
	int ityp, 
	double *Xa, 
	double *dx, 
	int *indx
);
double dist2d(double *x1, double *x2);

double Rs(double sL, double sT, double s);
double Rwave(double cL, double cT);

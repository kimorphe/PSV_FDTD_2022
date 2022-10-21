
class Crack{
	public:
		double X1[2], X2[2];	// crack-tip coordinate
		double Len;	// crack length
		int Ng;		// number of grids
		double *c0;	// coheasive stress (static)	
		double *c1;	// coheasive stress (dynamic)
		double *v1p, *v1m;	// velocity 
		double *v2p, *v2m;	// velocity
		double *s11p, *s11m;	// normal stress
		double *s22p, *s22m;	// normal stress
		double *s11b, *s22b;	// normal stress (mean) 
		double *s12;		// shear stress
		double *u1p,*u1m;	// crack-face velocity 
		double *u2p,*u2m;	// crack-face displacement
		double cod0;		// initial COD (set a huge number for open(linear) crack )
		int *iop;		// open/close status (0:open, 1:close)
		int *stat;		// state (0:stick, 1:slip, -1: open) 
		int *kv1,*kv2; // 1D index sets 
		int idir;	// orientation 0: horizontal, 1: vertical
		Crack();
		void gen_Crack(double *Y1, double *Y2, Dom2D dom);
		void print_grid(int num, Dom2D dom); 
		void s2v(Fld2D fld, Dom2D dom, double dt);
		void v2s(Fld2D fld, Dom2D dom, double dt);
		double IsSlip(int l, double str_tmp);
		void set_c01(int nseg, double *xi, double *c_0, double *c_1);
		void COD_hrz(int l, double dt, double rho, double *dx, double s22b_now);
		void COD_vrt(int l, double dt, double rho, double *dx, double s22b_now);
		void out_cod(FILE *fp, int it, int Nt, double dt, Dom2D dom);
		void out_csd(FILE *fp, int it, int Nt, double dt, Dom2D dom, Fld2D fld);
		void setIC(Dom2D dom, Fld2D fld);
	private:
};

int make_cracks(char *fname, Dom2D dom, Crack **crk);

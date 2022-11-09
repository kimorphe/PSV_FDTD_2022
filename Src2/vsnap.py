import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import bnd
import matplotlib.ticker

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.colors import ListedColormap
import copy

class Mask:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()   # comment line
        dat=fp.readline()
        dat=dat.strip().split(" ");
        xlim=[ float(dat[0]),float(dat[1])];

        fp.readline()   # comment line
        dat=fp.readline()
        dat=dat.strip().split(" ");
        ylim=[float(dat[0]),float(dat[1])];

        fp.readline()   # comment line
        dat=fp.readline()
        dat=dat.strip().split(" ");
        Ndiv=[int(dat[0]),int(dat[1])]
        fp.readline()

        mk=fp.readlines()
        mk=np.array(mk)
        mk=mk.astype(int)
        mk=np.reshape(mk,Ndiv)

        self.mk=np.transpose(mk)
        print(np.shape(mk))

        fp.close()
class Vfld:
	def __init__(self,fname):
		fp=open(fname,"r");
		tstmp=fp.readline();
		#txt=self.tstmp.split("=")
		#self.time=float(txt[1])
		self.time=float(fp.readline());

		fp.readline();
		tmp=fp.readline().lstrip().split(",");
		self.Xa=list(map(float,tmp));

		fp.readline();
		tmp=fp.readline().lstrip().split(",");
		self.Xb=list(map(float,tmp));
		self.xlim=[self.Xa[0],self.Xb[0]]
		self.ylim=[self.Xa[1],self.Xb[1]]

		fp.readline();
		tmp=fp.readline().lstrip().split(",");
		self.Ng=list(map(int,tmp));
		print("xlim=",self.xlim)
		print("ylim=",self.ylim)

		fp.readline();

		dat=fp.readlines();

		ndat=self.Ng[0]*self.Ng[1];
		self.v1=np.zeros(ndat);
		self.v2=np.zeros(ndat);
		ndat=0
		for row in dat:
			item=row.lstrip().split(" ");
			self.v1[ndat]=float(item[0])
			self.v2[ndat]=float(item[1])
			ndat+=1

		fp.close()
		self.v1=np.reshape(self.v1,[self.Ng[0],self.Ng[1]])
		self.v2=np.reshape(self.v2,[self.Ng[0],self.Ng[1]])
		self.v1=np.transpose(self.v1)
		self.v2=np.transpose(self.v2)
		self.v=np.sqrt(self.v1*self.v1+self.v2*self.v2);
	def draw1(self,ax,vmin=-1.e-08,vmax=0.01,cmap="nipy_spectral",Fac=1):
		rng=[self.xlim[0],self.xlim[1],self.ylim[0], self.ylim[1]];
		V=np.sqrt(self.v1*self.v1+self.v2*self.v2)+1.e-08;
		my_cmap=plt.cm.jet
		my_cmap=copy.copy(plt.cm.get_cmap("jet"))
		#my_cmap.set_over("w")
		my_cmap.set_under("k")
		#img=ax.imshow(V*Fac,extent=rng,vmin=vmin,vmax=vmax,cmap=my_cmap,origin="lower",rasterized=True);
		img=ax.imshow(V*Fac,extent=rng,vmin=vmin,vmax=vmax,cmap=my_cmap,origin="lower",rasterized=False);
		ax.set_xlabel("x[mm]",fontsize=12)
		ax.set_ylabel("y[mm]",fontsize=12)
		return(img)


if __name__=="__main__":

    #vmsk=Mask()
    #vmsk.load("vmask.dat")

    dir_name="./"
    #geom=bnd.Crv(dir_name+"bnd.out")

    nfile=101;
    inc=1;
    n1=0

    fig=plt.figure();
    ax=fig.add_subplot(1,1,1)
    ax_divider=make_axes_locatable(ax)
    cax=ax_divider.append_axes("right",size="5%",pad="3%")	

    isum=0
    fmt=matplotlib.ticker.ScalarFormatter(useMathText=True)
    Fac=100
    Fac=20
    for k in range(n1,n1+nfile,inc):
        fname="v"+str(k)+".out";
        print(fname)
        vf=Vfld(fname);
        #vf.v-=10*vmsk.mk
        if isum==0:
            img=vf.draw1(ax,cmap="jet",Fac=Fac,vmax=1.0);
            #plt.colorbar(cax1,orientation='horizontal')
            fig.colorbar(img,cax=cax,orientation="vertical",format=fmt)
        else:
            img.set_data(vf.v*Fac)
        isum+=1
        #ax.text(20,30,"t={:.2f}".format(vf.time)+"$\mu$sec",color="w")
        ax.set_title("t={:.2f}".format(vf.time)+"$\mu$sec",color="k",fontsize=14)#,x=20,y=30)
        #geom.draw(ax,"w",lw=1.0)
        #ax.set_xlim([-25,20])

        fig.savefig(str(k)+".png",bbox_inches="tight",dpi=300)
        #fig.savefig(str(k)+".pdf",bbox_inches="tight",dpi=300)
        #pdf.savefig(fig,dpi=400)
        #ax.clear()
        #plt.pause(0.1)
        #raw_input("press enter to continue");
    #pdf.close()

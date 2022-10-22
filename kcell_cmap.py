import numpy as np
import matplotlib.pyplot as plt

import matplotlib.collections as mcll
from matplotlib import patches

from matplotlib.colors import ListedColormap


class KCELL:
    def __init__(self):
        colors=["gray","k"]
        cmap=ListedColormap(colors,name="my_cmap")
        #cmap.set_under("k")
        #cmap.set_over("pink")
        #cmap.set_bad("w")
        self.cmap=cmap
    def load(self,fname):
        fp=open(fname,"r")

        cmt=fp.readline()
        dat=fp.readline().strip().split(" ");
        x=float(dat[0]); y=float(dat[1])
        Xa=np.array([x,y])

        dat=fp.readline().strip().split(" ");
        x=float(dat[0]); y=float(dat[1])
        Xb=np.array([x,y])
        Wd=Xb-Xa;

        for k in range(6):
            fp.readline()

        cmt=fp.readline()
        dat=fp.readline().strip().split(" ");
        Nx=int(dat[0])
        Ny=int(dat[1])
        Z=[]
        cmt=fp.readline()
        for row in fp:
            Z.append(int(row))

        Z=np.array(Z)
        Z=np.reshape(Z,[Nx,Ny])

        self.Nx=Nx
        self.Ny=Ny
        self.Ndiv=np.array([Nx,Ny])
        self.Z=np.transpose(Z)
        self.Xa=Xa
        self.Xb=Xb
        self.Wd=Wd
        self.dx=self.Wd/self.Ndiv
    def show(self,ax):
        dx=self.dx;
        Xa=self.Xa;
        Xb=self.Xb;
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        #ax.imshow(self.Z,aspect=1.0,extent=ext,cmap="gray",origin="lower",vmin=-1.0,vmax=1)
        ax.imshow(self.Z,aspect=1.0,extent=ext,cmap=self.cmap,origin="lower")
        ax.set_ylim([Xa[1]-10,Xb[1]+10])
        #ax.set_xlim([Xa[0]+10,Xb[0]-10])
        ax.grid(True)
        #ax.grid(True)
    def pix2patch(self):
        n=0
        ptchs=[]
        ptc=patches.Rectangle([0.0,0.0],width=20*self.dx[0],height=20*self.dx[1])
        for ix in range(self.Nx):
            for iy in range(self.Ny):
                val=self.Z[iy,ix]
                if val==1:
                    xx=self.Xa[0]+ix*self.dx[0];
                    yy=self.Xa[1]+iy*self.dx[1];
                    #ptchs.append(patches.Rectangle([xx,yy],width=self.dx[0],height=self.dx[1]))
                    #patches.Rectangle([xx,yy],width=self.dx[0],height=self.dx[1])
                    ptc.set_xy([xx,yy])
                    ptchs.append(ptc)
                    

        #self.xx=xx*self.dx[0]+self.Xa[0]
        #self.yy=yy*self.dx[1]+self.Xa[1]
        ptc.set_color("k")
        self.ptc=ptc
        self.ptchs=ptchs


if __name__=="__main__":
    K=KCELL()

    fname="kcell.dat"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    K.load(fname)
    ax.set_facecolor("k")
    #K.pix2patch()
    K.show(ax)
    #fname=fname.replace("out","png")
    #fig.savefig(fname,bbox_inches="tight")

    #ax.plot(K.xx,K.yy,".")
    #for m in range(5):
    #    ax.add_patch(K.ptchs[m])


    plt.show()


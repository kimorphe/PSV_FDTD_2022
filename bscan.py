#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys


class REC:
    def __init__(self,fname):
        fp=open(fname,"r");
        fp.readline();
        tmp=fp.readline().lstrip().split(" ");
        self.dt=float(tmp[0]);
        self.Nt=int(tmp[1]);
        self.time=self.dt*np.array(range(self.Nt))
        fp.readline();
        tmp=fp.readline().lstrip().split(" ");
        self.Ng=int(tmp[0]);
        self.ityp=int(tmp[1]);
        fp.readline();
        self.xsrc=[];
        self.ysrc=[];
        for k in range(self.Ng):
            self.xsrc.append(0.);
            self.ysrc.append(0.);
            (self.xsrc[k],self.ysrc[k])=fp.readline().lstrip().split(" ");
        self.xsrc=list(map(float,self.xsrc));
        self.ysrc=list(map(float,self.ysrc));
        self.xsrc=np.array(self.xsrc);
        self.ysrc=np.array(self.ysrc);
        fp.readline();
        dat=fp.readlines();
        dat=np.array(list(map(float,dat)))
        print(np.shape(dat))
        print(self.Nt*self.Ng);
        self.dat=np.reshape(dat,(self.Ng,self.Nt))
    def bscan(self,v1=-0.001,v2=0.001,cmap="gray"):
        fig=plt.figure();
        ax=fig.add_subplot(1,1,1);
        x1=self.ysrc[0];
        x2=self.ysrc[-1];
        x1=self.xsrc[0];
        x2=self.xsrc[-1];
        t1=0.;
        t2=self.dt*self.Nt;
        im=ax.imshow(self.dat,extent=[t1,t2,x1,x2],cmap=cmap,aspect="auto",vmin=v1,vmax=v2,interpolation="bicubic",origin="lower");
        return fig,ax;
    def bscan2(self,ax,v1=-0.001,v2=0.001,cmap="gray"):
        x1=self.ysrc[0];
        x2=self.ysrc[-1];
        x1=self.xsrc[0];
        x2=self.xsrc[-1];
        t1=0.;
        t2=self.dt*self.Nt;
        im=ax.imshow(self.dat,extent=[t1,t2,x1,x2],cmap=cmap,aspect="auto",vmin=v1,vmax=v2,interpolation="bicubic",origin="lower");
    def get_amp(self,rnum,time):
        it=int(time/self.dt)
        err=False
        if rnum <0:
            err=True
        if rnum >self.Ng:
            err=True
        if err:
            print("Invalid receiver number given to get_amp ",rnum)
            return(0)

        if it <0:
            return(0)
        if it >= self.Nt:
            return(0)

        return(self.dat[rnum,it])
    def delay_sum(self,cR=2.8): 
        ydat=np.zeros(self.Nt)
        for k in range(self.Ng):
            Tdiff=(self.xsrc[k]-self.xsrc[0])/cR
            for m in range(self.Nt):
                tt=m*self.dt+Tdiff
                ydat[m]+=self.get_amp(k,tt)
        self.Ysum=ydat/self.Ng;

if __name__=="__main__":


    fname="rec0.out"
    dir_name="."

    narg=len(sys.argv)
    if narg>1:
        dir_name=sys.argv[1]
    if narg>2:
        fname=sys.argv[2]

    rec=REC(dir_name+"/"+fname)
    fig,ax=rec.bscan();

    fig2=plt.figure()
    ax=fig2.add_subplot(111)

    vels=np.linspace(2.7,3.1,9)
    dA=0.005
    A0=0
    for cR in vels:
        rec.delay_sum(cR=cR)
        ax.plot(rec.time,rec.Ysum+A0,"-",label=str(cR)+"km/s")
        A0+=dA
    
    #ax.plot(rec.time,rec.dat[0,:],"-b")
    ax.grid(True)
    ax.legend()
    plt.show()


		
		

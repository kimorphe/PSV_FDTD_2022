#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys


class REC:
    def __init__(self,fname):
        #fp=open(fname,"r");
        with open(fname,"r") as fp: # open file ( B-scan  rec*.out)
            fp.readline(); # comment line
            tmp=fp.readline().lstrip().split(" ");
            self.dt=float(tmp[0]);  # time increment dt
            self.Nt=int(tmp[1]);    # total time steps Nt
            self.time=self.dt*np.array(range(self.Nt))  # make time axis data

            fp.readline(); # comment line
            tmp=fp.readline().lstrip().split(" ");
            self.Ng=int(tmp[0]);    # number of grid (=waveforms)
            self.ityp=int(tmp[1]);  # grid type (0: v1, 1: v2)

            fp.readline(); # comment line
            self.xsrc=[]; # x-grid coordinate  
            self.ysrc=[]; # y-grid coordinate
            for k in range(self.Ng):
                self.xsrc.append(0.);
                self.ysrc.append(0.);
                (self.xsrc[k],self.ysrc[k])=fp.readline().lstrip().split(" ");
            self.xsrc=list(map(float,self.xsrc));
            self.ysrc=list(map(float,self.ysrc));
            self.xsrc=np.array(self.xsrc);
            self.ysrc=np.array(self.ysrc);

            fp.readline(); # comment line
            dat=fp.readlines();
            dat=np.array(list(map(float,dat)))
        self.dat=np.reshape(dat,(self.Ng,self.Nt))  # B-scan Data (Ng x Nt)
        self.tlim=[0.0,self.time[-1]]

    def bscan(self,v1=-0.05,v2=0.05,cmap="jet"):
        fig=plt.figure();
        ax=fig.add_subplot(111);
        x1=self.ysrc[0];
        x2=self.ysrc[-1];
        x1=self.xsrc[0];
        x2=self.xsrc[-1];
        t1=0.;
        t2=self.dt*self.Nt;
        im=ax.imshow(self.dat,extent=[t1,t2,x1,x2],cmap=cmap,aspect="auto",vmin=v1,vmax=v2,interpolation="bicubic",origin="lower");
        return fig,ax;
    def get_amp(self,rnum,time):

        err=False
        if rnum <0:
            err=True
        if rnum >self.Ng:
            err=True

        if err:
            print("Invalid receiver number given to get_amp ",rnum)
            return(0)

        it=int(time/self.dt)
        if it <0:
            return(0.)
        if it >= self.Nt:
            return(0.)
        if it == self.Nt-1:
            return(self.dat[rnum,-1])

        xi=(time-it*self.dt)/self.dt
        eta=1.-xi;
        amp=eta*self.dat[rnum,it]+xi*self.dat[rnum,it+1]

        return(amp)

    def delay_and_sum(self,th_deg,cp):
        th=np.deg2rad(th_deg)
        X0=self.xsrc[0];
        ydat=np.zeros(self.Nt)
        for k in range(self.Ng):
            Xk=self.xsrc[k];
            dly=(Xk-X0)*np.sin(th)/cp
            for m in range(self.Nt):
                tt=m*self.dt-dly
                ydat[m]+=self.get_amp(k,tt)
        self.Ysum=ydat/self.Ng;


if __name__=="__main__":


    dir_name="."; fname="rec0.out";
    #-----------------Parse console input ----------------
    narg=len(sys.argv)
    if narg>1: # data directory provided
        dir_name=sys.argv[1]
    if narg>2: # data directory & file name provided
        fname=sys.argv[2]       
    #-----------------------------------------------------

    rec=REC(dir_name+"/"+fname)  # create REC class instance
    fig,ax=rec.bscan(); # show B-scan plot 

    fig2=plt.figure(figsize=[8,4])
    ax=fig2.add_subplot(111)
    cps=[3.0]# phase velocity
    ths=[50,60.0,65.];     # transducer angle
    A0=0


    for cp in cps:
        for th in ths:
            rec.delay_and_sum(th,cp)
            ax.plot(rec.time,rec.Ysum+A0,"-",label=str(th)+"[deg]")

    ax.set_xlim(rec.tlim)
    ax.grid(True)
    plt.show()


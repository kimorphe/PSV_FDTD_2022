import numpy as np
import matplotlib.pyplot as plt

class dom2d:	# Domain Class
	def __init__(self,fname):
		fp=open(fname,"r");
		fp.readline();	
		tmp=fp.readline().lstrip().split(" "); self.xa=np.array(list(map(float,tmp)))
		tmp=fp.readline().lstrip().split(" "); self.xb=np.array(list(map(float,tmp)))
		fp.readline();	
		tmp=fp.readline().lstrip().split(" "); self.Xa=np.array(list(map(float,tmp)))
		tmp=fp.readline().lstrip().split(" "); self.Xb=np.array(list(map(float,tmp)))
		fp.readline();	
		tmp=fp.readline().lstrip().split(" "); self.Ya=np.array(list(map(float,tmp)))
		tmp=fp.readline().lstrip().split(" "); self.Yb=np.array(list(map(float,tmp)))
		fp.readline();	
		tmp=fp.readline().lstrip().split(" "); self.Ndiv=np.array(list(map(int,tmp)))
		fp.readline();
		tmp=fp.readlines(); self.kcell=np.array(list(map(int,tmp)));
		self.Ndat=self.Ndiv[0]*self.Ndiv[1];

		dx=((self.Xb-self.Xa)/self.Ndiv);
		self.xcod=(np.arange(self.Ndiv[0])+0.5)*dx[0]+self.Xa[0];
		self.ycod=(np.arange(self.Ndiv[1])+0.5)*dx[1]+self.Xa[1];

	def show(self):
		print("xa  =",self.xa[0],self.xa[1]);
		print("xb  =",self.xb[0],self.xb[1]);
		print("Xa  =",self.Xa[0],self.Xa[1]);
		print("Xb  =",self.Xb[0],self.Xb[1]);
		print("Ya  =",self.Ya[0],self.Ya[1]);
		print("Yb  =",self.Yb[0],self.Yb[1]);
		print("Ndiv=",self.Ndiv[0],self.Ndiv[1]);
		print("Ndat=",self.Ndat);
		print("Shape of kcell=",np.shape(self.kcell));

	def draw_dom(self):
		Y,X=np.meshgrid(self.ycod,self.xcod)
		print(np.shape(X))
		Kcell=np.reshape(self.kcell,(self.Ndiv[0],self.Ndiv[1]))	
		indx=np.arange(self.Ndiv[1],0,-1)-1;
		for k in range(self.Ndiv[0]):
			Kcell[k]=Kcell[k][indx];

		Kcell=np.transpose(Kcell);
		print(np.shape(Kcell))
		fig=plt.figure();
		ax=fig.add_subplot(1,1,1)
		rng=[self.Xa[0],self.Xb[0],self.Xa[1],self.Xb[1]]
		ax.imshow(Kcell,extent=rng,vmin=-1,vmax=1);
		#ax.pcolor(X,Y,Kcell)	# <-- very slow

		Xa=self.Xa;
		Xb=self.Xb;
		XX=[Xa[0],Xb[0],Xb[0],Xa[0],Xa[0]];
		YY=[Xa[1],Xa[1],Xb[1],Xb[1],Xa[1]];

		xa=self.xa;
		xb=self.xb;
		xx=[xa[0],xb[0],xb[0],xa[0],xa[0]];
		yy=[xa[1],xa[1],xb[1],xb[1],xa[1]];

		Ya=self.Ya;
		Yb=self.Yb;
		UU=[Ya[0],Yb[0],Yb[0],Ya[0],Ya[0]];
		VV=[Ya[1],Ya[1],Yb[1],Yb[1],Ya[1]];

		ax.plot(XX,YY,'k',lw=2)
		ax.plot(xx,yy,'b',lw=2)
		ax.plot(UU,VV,'w',lw=2)

		return ax

if __name__=="__main__":
    dm=dom2d("kcell.dat")
    dm.show()
    dm.draw_dom()

    plt.show()

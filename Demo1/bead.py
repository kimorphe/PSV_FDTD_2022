import numpy as np
import matplotlib.pyplot as plt

class Bump:
    def __init__(self,a,dh):
        self.a=a
        self.dh=dh
        self.R=(a*a+dh*dh)/dh*0.5
        self.b=np.sqrt(self.R*self.R-a*a)
    def Yval(self,X):
        Y=[]
        for x in X:
            if abs(x)>self.a:
                Y.append(0.0)
            else:
                R=self.R
                b=self.b
                Y.append(np.sqrt(R*R-x*x)-b)
        return(np.array(Y))

if __name__=="__main__":
    a=5.0
    dh=2.0
    h=12

    bp=Bump(a,dh)
    X=np.linspace(-20,20,201)
    Y=bp.Yval(X)

    bpl=Bump(2.0,1)
    Yl=bpl.Yval(X)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)
    ax.plot(X,Y+h)
    ax.plot(X,-Yl)
    ax.set_aspect(1.0)

    plt.show()

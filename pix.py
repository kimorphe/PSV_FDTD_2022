import numpy as np
import matplotlib.pyplot as plt

fig=plt.figure()
ax=fig.add_subplot()
ax.set_aspect(1.0)

A=np.random.random([20,20])
B=np.random.random([10,10])

ax.imshow(A,extent=[0,1,0,1])

from matplotlib import patches
R1=patches.Rectangle(xy=(0.25,0.25), width=0.2, height=0.4,facecolor="w",edgecolor="k")
R2=patches.Rectangle(xy=(0.5,0.2), width=0.3, height=0.05,facecolor="w",edgecolor="k")

#ax.add_patch(R1)
#ax.add_patch(R2)

import matplotlib.collections as mcll

PC=mcll.PatchCollection([R1,R2])
PC.set(color="w",alpha=0.5,edgecolor="k")
ax.add_collection(PC)
PC.set_label("rectangles")
ax.legend()

plt.show()

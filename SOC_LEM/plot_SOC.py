import matplotlib.pyplot as plt
import numpy as np

SOC = np.load('final_SOC.npy')
z = np.linspace(-.995,.995,200)
x = 35
y = 73
plt.figure(1)
im = plt.imshow(SOC[:,:,75])
plt.scatter(y,x,color='r')
plt.colorbar(im)
plt.figure(2)
plt.plot(SOC[x,y,:],z)
plt.show()

import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('postProcessing/probes/0/T', comments='#')

t = data[:,0]
p0 = data[:,1]
p1 = data[:,2]

fig, axs = plt.subplots(1,2,figsize=(9,3))
axs[0].plot(t, p0, 'r')
axs[1].plot(t, p1, 'b')
plt.show()
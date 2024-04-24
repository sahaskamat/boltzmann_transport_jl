import numpy as np
import matplotlib.pyplot as plt



ax = plt.figure().add_subplot(projection='3d')

orbit_python = np.loadtxt("orbitpoints_python.dat",delimiter=" ")
orbit_julia = np.loadtxt("orbitpoints_julia.dat",delimiter=" ")

ax.scatter(orbit_python[:,0],orbit_python[:,1],orbit_python[:,2],s=1)
ax.scatter(orbit_julia[:,0],orbit_julia[:,1],orbit_julia[:,2],s=1)

plt.show()
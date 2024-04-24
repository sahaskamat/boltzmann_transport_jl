import matplotlib.pyplot as plt
import numpy as np


#boltzmann_transport data, high res
thetalist,rhoxylist = np.genfromtxt("rhoZZvsT_python.dat",delimiter=" ",unpack=True)
plt.plot(thetalist,rhoxylist,ls="-",marker="o",ms=2,label="Python")

#marching cube
thetalist,rhoxylist = np.genfromtxt("rhoZZvsT_julia.dat",delimiter=" ",unpack=True,usecols=[0,1])
plt.plot(thetalist,(rhoxylist),ls="-",marker="o",ms=2,label="Julia")

plt.xlim(0,90)
plt.ylabel(r"$\rho_{zz}$ ($m\Omega$ cm )")
plt.xlabel(r'$\theta^\circ$')
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

fig,axes = plt.subplots()

for H in [10,20,50,100]:
    Tlist,rhoxylist = np.genfromtxt("rhoXYvsT_julia"+str(H)+".dat",delimiter=" ",unpack=True)
    axes.plot(Tlist,rhoxylist/10000,label="B = "+str(H)+" T")

axes.set_ylabel(r"$\rho_{xy}$ (m$\Omega$-cm)")
axes.set_xlabel(r"$T$ (K)")
axes.legend()
axes.text(230,-2e-7,r"NdLSCO x=0.24"+ "\n" r"$B$//$c$")

fig2,axes2 = plt.subplots()

for H in [10,20,50,100]:
    Tlist,rhoxylist = np.genfromtxt("rhoXYvsT_julia"+str(H)+".dat",delimiter=" ",unpack=True)
    axes2.plot(Tlist,(rhoxylist - rhoxylist[0])/10000,label="B = "+str(H)+" T")

axes2.set_ylabel(r"$\rho_{xy} - \rho_{xy0}$ (m$\Omega$-cm)")
axes2.set_xlabel(r"$T$ (K)")
axes2.legend()
axes2.text(230,0,r"NdLSCO x=0.24"+ "\n" r"$B$//$c$")

plt.show(block=True)

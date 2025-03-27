import numpy as np
import matplotlib.pyplot as plt
import hamiltonian as hamiltonian
from constants import CaF
from numpy.linalg import eigh
import calculate as calc
import scipy.constants 


#%%
muN = scipy.constants.physical_constants['nuclear magneton'][0]

Nmax=5
H_0,H_B,H_dc,H_ac = \
    hamiltonian.build(Nmax,CaF,zeeman=True,Edc=True,ac=True)


I = 0 #W/m^2
E = np.linspace(0, 5e6, int(40)) #V/m
B = 1 #T

H = H_0[..., None]+\
    H_B[..., None]*B+\
    H_dc[..., None]*E+\
    H_ac[..., None]*I 
H = H.transpose(2,0,1)

energies, states = calc.solve(H, Nmax, CaF, False)

d = calc.electric_moment(states, Nmax, CaF)

#%%Plot
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(E/1e5, d/CaF['d0'],ls = '-')#,color = '#fd411e')
# ax.set_xlim(0, 5)
ax.set_ylabel(r"$\rm{Electrice~Dipole~Moment ~(d_0)}$")
ax.set_xlabel(r"${\rm Electric~Field~ (kV/m)}$")
ax.tick_params (direction = 'in')
plt.tight_layout()
plt.show()



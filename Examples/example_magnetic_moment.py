import numpy as np
import matplotlib.pyplot as plt
import  hamiltonian as hamiltonian
import calculate as calc
from constants import CaF
from scipy import constants
h = constants.h
from numpy.linalg import eigh


#%%
muN = constants.physical_constants['nuclear magneton'][0]

Nmax=2
H_0,H_B,H_dc,H_ac = \
    hamiltonian.build(Nmax,CaF,zeeman=True,Edc=True,ac=True)


I = 0 #W/m^2
E = 0 #V/m
B = np.linspace(1, 60, int(600))*1e-4 #T

H = H_0[..., None]+\
    H_B[..., None]*B+\
    H_dc[..., None]*E+\
    H_ac[..., None]*I 
H = H.transpose(2,0,1)

energies, states, label_list = calc.solve(H, Nmax, CaF, True, B)

mu = calc.magnetic_moment(states, Nmax, CaF)

#%%N=0
fig = plt.figure()
ax = fig.add_subplot()

for i in range(0,6):
    ax.plot(B*10**4,mu[:,i]/muN,ls = '-',label = label_list[i])#,color = '#c20078')
# ax.set_ylim(0,250)
ax.set_xlim(0, 30)
ax.set_ylabel(r"$\rm{Magnetic~moment (\mu/ \mu_N)}$")
ax.set_xlabel(r"${\rm Magnetic~field~(G)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$\rm{N=0}$')
plt.legend(title= '(N, F, mF)',frameon=False, loc =  'center left', bbox_to_anchor = (1,0.5))
plt.tight_layout()
plt.show()
#%% N=1
fig = plt.figure()
ax = fig.add_subplot()

for i in range(7,16):
    ax.plot(B*10**4,mu[:,i]/muN,ls = '-',label = label_list[i])#,color = '#c20078')
# ax.set_ylim(0,250)
ax.set_xlim(0, 30)
ax.set_ylabel(r"$\rm{Magnetic~moment (\mu/ \mu_N)}$")
ax.set_xlabel(r"${\rm Magnetic~field~(G)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$\rm{N=1}$')
plt.legend(title= '(N, F, mF)',frameon=False, loc =  'center left', bbox_to_anchor = (1,0.5))
plt.tight_layout()
plt.show()



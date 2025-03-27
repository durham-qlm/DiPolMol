'''Generates a simple Zeeman plot for SrF as an example'''

import numpy as np
import matplotlib.pyplot as plt
import  hamiltonian as hamiltonian
import calculate as calc
from constants import CaF, SrF
from scipy.constants import h


color_list = ['#fedf08','#fcb001','#fd411e','#c20078','#c20078','#c20078','#c20078','#c20078']

#%%
Nmax=2
H_0,H_B,H_dc, H_ac = \
    hamiltonian.build(Nmax,SrF,zeeman=True,Edc=False,ac=False)
  

B = np.linspace(0,100,5000)*1e-4 #G

H = H_0[..., None]+\
    H_B[..., None]*B
H = H.transpose(2,0,1)

energies, states, label_list = calc.solve(H, Nmax, SrF, True, B)


#%%
plt.figure()
for x in range(0,len(label_list)):
     plt.plot(B,np.transpose(energies)[x]*1e-9/h,ls = '-',label = label_list[x])

plt.ylim(20,21)

plt.legend(title= '(N, F, mF)',frameon=False, bbox_to_anchor=(1.,1.05))
plt.ylabel("Energy/h (GHz)")
plt.xlabel("Magnetic Field (T)")
plt.show()

#%% Plot N=0
fig = plt.figure()
ax = fig.add_subplot()

for x in range(4):
    ax.plot(B*10**4,np.transpose(energies)[x]*1e-6/h-energies[0][0]*1e-6/h,ls = '-',label = label_list[x])
ax.set_xlim(0, 100)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Magnetic~Field~ (G)}$")
plt.legend(title= '(N, F, mF)',frameon=False, bbox_to_anchor=(1.,1.05))

ax.tick_params (direction = 'in')
ax.set_title(r'$N=0$')
plt.tight_layout()
plt.show()

#%% Plot N=1
fig = plt.figure()
ax = fig.add_subplot()

for x in range(4,16):
    ax.plot(B*10**4,np.transpose(energies)[x]*1e-6/h-energies[0][4]*1e-6/h,ls = '-',label = label_list[x])
    print(np.transpose(energies)[x][0]*1e-6/h-energies[0][0]*1e-6/h)
ax.set_ylim(-80,250)
ax.set_xlim(0, 100)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Magnetic~Field~ (G)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$N=1$')
plt.legend(title= '(N, F, mF)',frameon=False, bbox_to_anchor=(1.,1.05))
plt.tight_layout()
plt.show()

#%% Plot N=2
fig = plt.figure()
ax = fig.add_subplot()

for x in range(16,36):
    ax.plot(B*10**4,np.transpose(energies)[x]*1e-6/h-energies[0][17]*1e-6/h,ls = '-',label = label_list[x])
ax.set_ylim(-150,250)
ax.set_xlim(0, 100)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Magnetic~Field~ (G)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$N=2$')
plt.legend(title= '(N, F, mF)',frameon=False, bbox_to_anchor=(1.,1.05))
plt.tight_layout()
plt.show()
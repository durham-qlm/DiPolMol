'''Generates a simple dc Stark plot for CaF as an example
         !Takes a while to run due to high Nmax!
'''
import numpy as np
import matplotlib.pyplot as plt
import  hamiltonian as hamiltonian
import calculate as calc
from constants import BaF138
from scipy.constants import h
from numpy.linalg import eigh


Nmax=2
H_0,H_B,H_dc,H_ac = \
    hamiltonian.build(Nmax,BaF138,zeeman=True,Edc=True,ac=False)

I = 0.0# #W/cm^2
E = np.linspace(0, 3, int(2000))*1e5 #V/m
B = 0.0 #T

H = H_0[..., None]+\
    H_B[..., None]*B+\
    H_dc[..., None]*E
   # H_ac[..., None]*I 
H = H.transpose(2,0,1)

energies, states = calc.solve(H, Nmax, BaF138, False)

#%% Plot N=0
fig = plt.figure()
ax = fig.add_subplot()

for x in range(4):
    ax.plot(E*10**-3,np.transpose(energies)[x]*1e-6/h-energies[0][0]*1e-6/h,ls = '-',color = '#fd411e')#, label=label_list[x])
ax.set_xlim(0, 300)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Electric~Field~ (kV/m)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$N=0$')
# ax.legend(title= '(N, F)')
plt.tight_layout()
plt.show()

#%% Plot N=1
fig = plt.figure()
ax = fig.add_subplot()

for x in range(4,16):
    ax.plot(E*10**-3,np.transpose(energies)[x]*1e-6/h-energies[0][4]*1e-6/h,ls = '-',color = '#c20078')#,label=label_list[x])
ax.set_ylim(0,250)
ax.set_xlim(0, 300)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Electric~Field~ (kV/m)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$N=1$')
# ax.legend(title= '(N, F)')

plt.tight_layout()
plt.show()

#%% Plot N=2
fig = plt.figure()
ax = fig.add_subplot()

for x in range(17,36):
    ax.plot(E*10**-3,np.transpose(energies)[x]*1e-6/h-energies[0][17]*1e-6/h,ls = '-',color = '#fcb001')#,label=label_list[x])
ax.set_ylim(0,250)
ax.set_xlim(0, 300)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"${\rm Electric~Field~ (kV/m)}$")
ax.tick_params (direction = 'in')
ax.set_title(r'$N=2$')
# ax.legend(title= '(N, F)')

plt.tight_layout()
plt.show()

#%%
# np.save('results//dc_E_list.npy',E)
# np.save('results//dc_energies.npy',energies)
# np.save('results//dc_states.npy',states)
# np.save('results//dc_labels.npy',label_list)
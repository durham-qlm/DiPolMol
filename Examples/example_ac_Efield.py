import numpy as np
import matplotlib.pyplot as plt
import  hamiltonian as hamiltonian
import calculate as calc
from constants import CaF
from numpy.linalg import eigh
from scipy import constants
h = constants.h

color_list = ['#fedf08','#fedf08','#fedf08','#fcb001','#fd411e','#fd411e','#fd411e','#c20078','#c20078','#c20078','#c20078','#c20078']

#%% Import and calculate required constants

Consts = CaF
Nmax = 2

# [alpha_0,alpha_1,alpha_2] = calc.alpha_012(780, CaF)
# #Overwrite polarisability component values in constants file
# Consts['alpha_0'] = 0 #Looking at transition energy so set isotropic component to zero
# Consts['alpha_1'] = alpha_1
# Consts['alpha_2'] = alpha_2


#%% change beta
# n = 30
# beta_list = np.linspace(0,np.pi/2,n)   
# I = 30e9 #W m^-2
# E = 0 #V/m
# B = 300e-4 #T

# energy_list = []
# for i in range(n):
#     Consts['Beta'] = beta_list[i]
#     H_0,H_B,H_dc,H_ac = hamiltonian.build(Nmax,Consts,zeeman=True,Edc=True,ac=True)

#     H = H_0[..., None]+\
#         H_B[..., None]*B+\
#         H_dc[..., None]*E+\
#         H_ac[..., None]*I 
#     H = H.transpose(2,0,1)
#     energies, states = calc.solve(H, Nmax, CaF, False)
#     energy_list.append(energies)


#%%Plot beta change

# fig = plt.figure()
# ax = fig.add_subplot()

# for x in range(4,10):
#     ax.plot(beta_list*180/np.pi,np.transpose(energy_list)[13-x][0]*1e-6/h-energy_list[0][0][4]*1e-6/h,ls = '-')
# ax.set_ylim(-0,120)
# ax.set_xlim(0, 90)
# ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
# ax.set_xlabel(r"$\beta (^\circ)$")
# ax.tick_params (direction = 'in')
# plt.tight_layout()
# plt.show()



#%%Calculate zero light field energies
I=0
H_0,H_B,H_dc,H_ac = hamiltonian.build(Nmax,Consts,zeeman=False,Edc=False,ac=False)

H = H_0[..., None]
H = H.transpose(2,0,1)
energies_0, states_0 = calc.solve(H, Nmax, CaF, False)

#%% change lambda
n = 20
lambda_list = np.linspace(490,800,n)   
I = 30e9 #W m^-2
E = 0 #V/m
B = 0 #T
Consts['Beta'] = 0

energy_list = []
for i in range(n):
    [alpha_0,alpha_1,alpha_2] = calc.alpha_012(lambda_list[i], CaF)
    #Overwrite polarisability component values in constants file
    Consts['alpha_0'] = alpha_0
    Consts['alpha_1'] = alpha_1
    Consts['alpha_2'] = alpha_2
    H_0,H_B,H_dc,H_ac = hamiltonian.build(Nmax,Consts,zeeman=False,Edc=False,ac=True)

    H = H_0[..., None]+\
        H_ac[..., None]*I 
    H = H.transpose(2,0,1)
    energies, states = calc.solve(H, Nmax, CaF, False, B)
    energy_list.append(energies)

print(np.transpose(energy_list)[13-4][0]*1e-6/h-energy_list[0][0][4]*1e-6/h)
#%%Plot lambda change

fig = plt.figure()
ax = fig.add_subplot()
x_list = range(4,16)
for x in x_list:
    ax.plot(lambda_list,np.transpose(energy_list)[x][0]*1e-6/h-energies_0[0][x]*1e-6/h,ls = '',marker = '.')#-energy_list[0][0][4]*1e-6/h
ax.hlines(0,500,650)
ax.set_ylim(-175,-60)
ax.set_xlim(640, 650)
ax.set_ylabel(r"$\rm{Energy~shift/h~(MHz)}$")
ax.set_xlabel(r"$\lambda (nm)$")
ax.tick_params (direction = 'in')
plt.tight_layout()
plt.show()
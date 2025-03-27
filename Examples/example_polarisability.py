# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:08:48 2023

@author: xnkr87
"""

import matplotlib.pyplot as plt
from constants import CaF
import numpy as np
import calculate as calc
from scipy import constants
h = constants.h
c=constants.c
e_0 = constants.epsilon_0
params = {'font.size'       : 35,
          'font.family'     : 'Helvetica',
          'font.monospace'  : 'Computer Modern',
          'axes.labelsize'  : 35,
          'backend'         : 'ps',
          'legend.fontsize' : 35,
          'xtick.labelsize' : 35,
          'ytick.labelsize' : 35,
          'text.usetex'     : True,
          'axes.linewidth'  : 2.5,
          'xtick.major.size': 8,
          'xtick.major.width': 2.5,          
          'ytick.major.width': 2.5,          
          'ytick.major.size': 8,
          'xtick.minor.size': 4,
          'xtick.minor.width': 1.8,          
          'ytick.minor.width': 1.8,          
          'ytick.minor.size': 4,
          'figure.figsize'  : (12,7)}

plt.rcParams.update(params)
color_list = ['#c20078','#069af3','#2baf6a','#0d75f8','#0d75f8','#0d75f8','#c79fef']
#%% Hannah test - delete later

# a_list = np.asarray(calc.alpha_012(545.08, CaF))
# print(a_list/(h)*10**4)#cm^2
#%% Calculate
label_list = [0,1,2]
n = 5000
l_list = np.linspace(640,650,n)
shape = (n,3)
a_array = np.zeros(shape)
for i in range(n):
    a_array[i] = calc.alpha_012(l_list[i],CaF)
  

#%% Plotting

fig = plt.figure()
ax = fig.add_subplot()
ax.axhline(y=0,color = 'Black')
for x in range(3):
    ax.plot(l_list,np.transpose(a_array)[x]/h*10**4,ls = '',marker = '.',label = r'$'+str(label_list[x])+'$',color=color_list[x])

ax.set_ylim(-200,200)
ax.set_xlim(500,660)
plt.legend(title= r'$k$',frameon=False, bbox_to_anchor=(1.,0.9))
ax.set_ylabel(r"$\alpha_{(k)}^{'}~\rm{(Hz/(W/cm^2))}$")
ax.set_xlabel(r"$\lambda (nm)$")
ax.tick_params (direction = 'in')
plt.tight_layout()
plt.show()



import numpy as np
from numpy.linalg import eigh
import hamiltonian as hamiltonian
from constants import SrF
import calculate as calc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.constants import h
from pylab import *

#%%
Nmax=2
H_0,H_B,H_dc, H_ac = \
      hamiltonian.build(Nmax,SrF,zeeman=True,Edc=False,ac=False)  
    
I = 0 #W/m^2
E = 0 #V/m
B = 60e-6#T

H = H_0[..., None]+\
    H_B[..., None]*B
    
H = H.transpose(2,0,1)

energies, states, label_list = calc.solve(H, Nmax, SrF, True, B)


#%% Create matrix of tdm
def tdmatrix(gs):
    dm = np.round(calc.transition_dipole_moment(Nmax,SrF,+1,states,gs),6)
    dz = np.round(calc.transition_dipole_moment(Nmax,SrF,0,states,gs),6)
    dp = np.round(calc.transition_dipole_moment(Nmax,SrF,-1,states,gs),6)
    shape = (len(label_list),3)
    tdm_list = np.zeros(shape)#tdm, polarisation, energy
    for i in range(len(label_list)):
        if dm[i]+dz[i]+dp[i] ==0:
            tdm_list[i,0]=0
            tdm_list[i,1]=0
            tdm_list[i,2]= energies[0][i]*1e-6/h 
        else:    
            tdm_list[i,0]=dm[i]+dz[i]+dp[i]
            tdm_list[i,1]=(dm[i]*1+dz[i]*2+dp[i]*3)/((tdm_list[i,0]))
            tdm_list[i,2]= energies[0][i]*1e-6/h   

    return (tdm_list)



#%% Table with colour gradient
gs_list  = [6,7,10,11,15] #random selection of ground states
cmap_dm = cm.get_cmap('Blues', 1000)    # PiYG
dm_colors = []
for i in range(cmap_dm.N):
    rgba = cmap_dm(i)
    # rgb2hex accepts rgb or rgba
    dm_colors.append(colors.rgb2hex(rgba))
    
cmap_dz = cm.get_cmap('Purples', 1000)    # PiYG
dz_colors = []
for i in range(cmap_dz.N):
    rgba = cmap_dz(i)
    # rgb2hex accepts rgb or rgba
    dz_colors.append(colors.rgb2hex(rgba))
    
cmap_dp = cm.get_cmap('Reds', 1000)    # PiYG
dp_colors = []
for i in range(cmap_dp.N):
    rgba = cmap_dp(i)
    # rgb2hex accepts rgb or rgba
    dp_colors.append(colors.rgb2hex(rgba))
    
    
cell_text = []
colour_array = []
row_labels = []
column_labels = []
for i in range(len(label_list)): 
    column_labels.append((r'(' + str(label_list[i][0]) + ',' + str(label_list[i][1]) + ',' + str(label_list[i][2]) + ')'))
    
for x in range(len(gs_list)):
    gs = gs_list[x]
    text_list = []
    colour_list = []
    ground_state = label_list[gs]
    row_labels.append(r'(' + str(ground_state[0]) + ',' + str(ground_state[1]) + ',' + str(ground_state[2]) + ')')
    td_matrix = tdmatrix(gs)
    for i in range(len(label_list)):
        text_list.append(str(np.abs(round(td_matrix[i,0],3))))
        colour_scale = int(td_matrix[i,1])
        if colour_scale ==1:
            colour_index = abs(int(float(td_matrix[i,0])*1000))
            # print(colour_index)
            colour_list.append(dm_colors[colour_index])
        elif colour_scale ==2:
            colour_index = abs(int(float(td_matrix[i,0])*1000))
            colour_list.append(dz_colors[colour_index])
        elif colour_scale ==3:
            colour_index = abs(int(float(td_matrix[i,0])*1000))
            colour_list.append(dp_colors[colour_index])
        else:
            colour_list.append('#ffffff')
    cell_text.append(text_list)
    colour_array.append(colour_list)


    #%% Table

fig = plt.figure(figsize=(5,2))
ax = fig.add_subplot(111, frameon=False, xticks=[], yticks=[])
ax.axis('tight')
ax.axis('off')
the_table=plt.table(cellText=np.transpose(cell_text), rowLabels=column_labels, colLabels=row_labels ,# colWidths = [0.03]*vals.shape[1],
                      loc='center',cellColours=np.transpose(colour_array),cellLoc='center'
                      )
the_table.auto_set_font_size(False)
the_table.set_fontsize(10)
plt.show()


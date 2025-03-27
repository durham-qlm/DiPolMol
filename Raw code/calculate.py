import hamiltonian as hamiltonian
import numpy as np
import scipy.constants
from sympy.physics.wigner import wigner_3j as sympy_wig_3j
from sympy.physics.wigner import wigner_6j as sympy_wig_6j
import matplotlib.pyplot as plt
from constants import CaF
###############################################################################
# Start by definining a bunch of constants that are needed for the code       #
###############################################################################

'''
    Important note!

    All units in this code are SI i.e. elements in the Hamiltonian have units
    of Joules. Outputs will be on the order of 1e-30

'''

h = scipy.constants.h
muN = scipy.constants.physical_constants['nuclear magneton'][0]
bohr = scipy.constants.physical_constants['Bohr radius'][0]
epsilon_0 = scipy.constants.epsilon_0
c = scipy.constants.c

#Defining Wigner symbols to ensure on float errors
def wigner_3j(a, b, c, d, e, f):
    a = float(a)
    b = float(b)
    c = float(c)
    d = float(d)
    e = float(e)
    f = float(f)
    return sympy_wig_3j(a,b,c,d,e,f)

def wigner_6j(a, b, c, d, e, f):
    a = float(a)
    b = float(b)
    c = float(c)
    d = float(d)
    e = float(e)
    f = float(f)
    return sympy_wig_6j(a,b,c,d,e,f)

#Frequencies of each vibrational transition needed for polarisability calculation
def XXfreq(v1,v2,consts):
    f = 100*c*(consts['XT'][v2]+2*consts['XB'][v2] - consts['XT'][v1]) 
    return f
def XAfreq(Xv,Av,Omega,consts):
    f = 100*c*(consts['AT'][Av]+consts['AA'][Av]*(Omega-1) - consts['XT'][Xv])
    return f
def XBfreq(Xv,Bv,consts):
    f = 100*c*consts['BT'][Bv] - consts['XT'][Xv]
    return f




def alpha_012(l, consts):
    '''Calculates the polarisability of the X(v=0) state for a given wavelength l.
    
    Args: 
        l (float): Wavelength at which to calculate the polarisability at. (m)
        consts (dict): Dictionary of constants for the molecule to be calculated.
        
    Returns:
        alpha (list): List containing alpha_0, alpha_1 and alpha_2, the scalar, vector 
        and tensor components of molecular polarisability respectively.'''
        
    f=c/(l*10**-9)
    XXmu = consts['d0']
    XXvibmu = consts['XXvibmu']
    XAmu = consts['XAmu']
    XBmu = consts['XBmu']
    XAfc = consts['XAfc']
    XBfc = consts['XBfc']
    #There are three alpha terms: alpha parallel (Sigma overlap) and alpha perpendicular (Pi overlap for both Omega values)
    parallel_terms = [[XXfreq(0,0,consts),XXmu],[XXfreq(0,1,consts),XXvibmu],[XBfreq(0, 0,consts),XBmu*(XBfc[0][0]**0.5)]]
    perp_terms1 = [[XAfreq(0,0,0.5,consts), XAmu*(XAfc[0][0]**0.5)],[XAfreq(0,1,0.5,consts), XAmu*(XAfc[0][1]**0.5)]]
    perp_terms3 = [[XAfreq(0,0,1.5,consts), XAmu*(XAfc[0][0]**0.5)],[XAfreq(0,1,1.5,consts), XAmu*(XAfc[0][1]**0.5)]]
    a_par = 0
    a_perp1 = 0
    a_perp3 = 0
    for i in range(len(parallel_terms)):
        a_par += 1/h*(1/(parallel_terms[i][0]+f)+1/(parallel_terms[i][0]-f))*parallel_terms[i][1]**2
    for i in range(len(perp_terms1)):
        a_perp1 += 1/h*(1/(perp_terms1[i][0]+f)+1/(perp_terms1[i][0]-f))*perp_terms1[i][1]**2
        a_perp3 += 1/h*(1/(perp_terms3[i][0]+f)+1/(perp_terms3[i][0]-f))*perp_terms3[i][1]**2
    alpha_0 = 1/3*(a_par+a_perp1+a_perp3)*10**4/(2*c*epsilon_0)*10**-4
    alpha_1 = 1/2*(f/XAfreq(0,0,0.5,consts)*a_perp1-f/XAfreq(0,0,1.5,consts)*a_perp3)*10**4/(2*c*epsilon_0)*10**-4
    alpha_2 = 1/3*(2*a_par - a_perp1 - a_perp3)*10**4/(2*c*epsilon_0)*10**-4
    return [alpha_0,alpha_1, alpha_2]


def dipole(Nmax,Consts,M):
    ''' Generates the induced dipole moment operator for a Rigid rotor.
    Expanded to cover state vectors in the uncoupled hyperfine basis.

    Args:
        Nmax (int): maximum rotational states.
        Consts (dict): Dictionary of constants for the molecule to be calculated.
        M (float): index indicating the helicity of the dipole field. -1 = S+, 0 = Pi, +1 = S-

    Returns:
        Dmat (np.ndarray): - dipole matrix
    '''


    I = Consts['I']
    S = Consts['S']
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    dmat = np.zeros((shape,shape),dtype= complex)
    dmat = (np.kron(dmat,np.kron(np.identity(Ishape),
                                                    np.identity(Sshape))))
    i=0
    j=0
    
    for N1 in range(0,Nmax+1):
        for J1 in np.arange(np.abs(N1-S),(N1+S+1),1):
            for F1 in np.arange(J1-I,(J1+I+1),1):
                for mF1 in np.arange(-F1,(F1+1),1):
                    for N2 in range(0,Nmax+1):
                        for J2 in np.arange(np.abs(N2-S),(N2+S+1),1):
                            for F2 in np.arange(J2-I,(J2+I+1),1):
                                for mF2 in np.arange(-F2,+(F2+1),1):
                                    dmat[i,j] = np.round((-1+0j)**(F2-mF2+J1+J2+F1+N2+1))*\
                                        ((2*F1+1)*(2*F2+1)*(2*J1+1)*(2*J2+1))**0.5*((2*N1+1)*(2*N2+1))**0.5*(-1)**N2*wigner_3j(N2,1,N1,0,0,0)*\
                                           wigner_3j(F2,1,F1,-mF2,M,mF1)*wigner_6j(J1,F1,I,F2,J2,1)*wigner_6j(N1,J1,0.5,J2,N2,1)
                                    
                                    i+=1
                    i=0
                    j+=1

    
    return dmat

    
def transition_dipole_moment(Nmax,Consts,M,states,gs,locs=None):
    ''' Function to calculate the Transition Dipole Moment between a state gs
    and a range of states. Returns the TDM in units of the permanent dipole
    moment (d0).

    Args:
        Nmax (int): Maximum rotational quantum number in original calculations
        Consts (dict): Dictionary of constants for the molecule to be calculated
        M (float): Helicity of Transition, -1 = S+, 0 = Pi, +1 = S-
        States (np.ndarray): matrix for eigenstates of problem output from np.linalg.eig
        gs (int): index of ground state.

    kwargs:
        locs (list of ints): optional argument to calculate for subset of States, should be an
                array-like.

    Returns:
        TDM (list of floats): transition dipole moment between gs and States.
    
    '''
    states = states[0]
    dipole_op = dipole(Nmax,Consts,M)

    gs = np.conj(states[:,gs])
    if locs != None :
        states =  states[:,locs]

    tdm =  np.einsum('i,ij,jk->k',gs,dipole_op,states).real

    return tdm


def magnetic_moment(States, Nmax, Consts):
    '''Returns the magnetic moments of each eigenstate.
    
    Args:
        States (np.ndarray): matrix for eigenstates of problem output from np.linalg.eig
        Nmax (int): Maximum rotational quantum number in original calculations
        Consts (dict): Dictionary of constants for the molecule to be calculated
        
    Returns:
        mu (list of floats): magnetic moment for each eigenstate in States.
        
    '''
    
    muz = -1*hamiltonian.H_zee(Nmax,Consts,Consts['I'],Consts['S'])
    
    mu =np.einsum('ijk,jl,ilk->ik',
            np.conjugate(States),muz,
            States)
    return mu


def electric_moment(States, Nmax, Consts):
    '''Returns the electric dipole moments of each eigenstate
    
    Args:
        States (np.ndarray): matrix for eigenstates of problem output from np.linalg.eig
        Nmax (int): Maximum rotational quantum number in original calculations
        Consts (dict): Dictionary of constants for the molecule to be calculated
        
    Returns:
        d (list of floats): electric dipole moment for each eigenstate in States
    '''
    
    dz = -1*hamiltonian.H_dc(Nmax,Consts,Consts['I'],Consts['S'])
    
    d =np.einsum('ijk,jl,ilk->ik',
            np.conjugate(States),dz,
            States)
    return d



def sort_smooth(energy, states):
    ''' Sort states to remove false avoided crossings.
    This is a function to ensure that all eigenstates plotted change
    adiabatically, it does this by assuming that step to step the eigenstates
    should vary by only a small amount (i.e. that the  step size is fine) and
    arranging states to maximise the overlap one step to the next.
    
    Args:
        Energy (np.ndarray) : np.ndarray containing the eigenergies, as from np.linalg.eig
        States (np.ndarray): np.ndarray containing the states, in the same order as Energy
    Returns:
        Energy (np.ndarray) : np.ndarray containing the eigenergies, as from np.linalg.eig
        States (np.ndarray): np.ndarray containing the states, in the same order as Energy E[x,i] -> States[x,:,i]
    '''
    ls = np.arange(states.shape[2],dtype="int") #from 0 to number of states
    
    number_iterations = len(energy[:,0]) #ie the number of B field values
    
    for i in range(2, number_iterations): 
        '''
        This loop sorts the eigenstates such that they maintain some
        continuity. Each eigenstate should be chosen to maximise the overlap
        with the previous.
        '''
        #calculate the overlap of the ith and jth eigenstates
        overlaps = np.einsum('ij,ik->jk',
                                np.conjugate(states[i-1,:,:]),states[i,:,:])     
        orig2 = states[i,:,:].copy() 
        orig1 = energy[i,:].copy()
        #insert location of maximums into array ls
        np.argmax(np.abs(overlaps),axis=1,out=ls)
        
        for k in range(states.shape[2]): 
            l = ls[k] 
            if l!=k:
                energy[i,k] = orig1[l].copy()
                states[i,:,k] = orig2[:,l].copy()
    return energy, states





def label_FmF_states(States, Nmax, consts, B):
    """
    Function that takes the array of eigenstates as generated from the hamiltonian 
    (and processed with sort_smooth to avoid false avoided crossings) and assigns 
    each an (N, F, mF) label.
    
    Inputs:
        States (np.ndarray) : np.ndarray containing the eigenstates, as from np.linalg.eig
        Nmax (int) : Maximum N state to consider in the calculations
        B (int/float/list/array) : Magnetic field values used to calculate the eigenstates
        consts (dict): Dictionary of constants for the molecule to be calculated
        
    Returns:
        FmF_labels or F_labels (list): list of [N, F, mF] or [N, F] labels, one for each of the eigenstates in the input
    """
    
    S = consts['S']
    I= consts['I']
    
    #initialise error counters to 0
    zeromagfielderror = 0
    F_error = 0
    mF_error = 0
    
    startingstateindex = 0
    
    #check the B fields input by the user, select a non-zero B field to label at
    if type(B) == float or type(B) == int:
        States = States[0]
        if B == 0:
            zeromagfielderror += 1
    else:
        if B[0] != 0:
            States = States[0]  
        else:
            if len(States) > 1:
                States = States[1] 
                startingstateindex += 1
            else:
                States = States[0]
                zeromagfielderror += 1

    

    #generate a list of the possible (N, J, F, mF) states 
    statelist = []
    for N1 in range(0,Nmax+1):
        for J1 in np.arange(np.abs(N1-S),(N1+S+1),1):
            for F1 in np.arange(J1-I,(J1+I+1),1):
                for mF1 in np.arange(-F1,(F1+1),1):
                    state = [N1, J1, F1, mF1]
                    statelist.append(state)
    
    FmF_labels = []
    F_labels = []
    
    #now match the (N,J,F,mF) labels with the eigenstates
    for i in range(len(States)): #For each of the eigenstates
        statebreakdown = []
        for n in range(len(States)): #loop through each value in the eigenstate
            if np.round(States[:,i][n],1) != 0: #consider the non-zero coeff
                statebreakdown.append(statelist[n]) 
                #following can be useful to check states
                #print(f'{States[:,i][n]} coeff of {statelist[n]}' )
                
        #Confirming consistency in the F,mF labelling of the states
        if len(statebreakdown) > 1: #only do this if there are multiple (F,mF) states to consider
            avgF = np.mean(statebreakdown, axis=0)[2]
            avgmF = np.mean(statebreakdown, axis=0)[3]
            
            if avgF != statebreakdown[0][2] or statebreakdown[0][2] != statebreakdown[1][2]:
                F_error += 1
                
            if avgmF != statebreakdown[0][3] or statebreakdown[0][3] != statebreakdown[1][3]:
                mF_error += 1
                F_labels.append([int(statebreakdown[0][0]), int(statebreakdown[0][2])])
                
            if avgF == statebreakdown[0][2] and statebreakdown[0][2] == statebreakdown[1][2] and avgmF == statebreakdown[0][3] and statebreakdown[0][3] == statebreakdown[1][3]:
                FmF_labels.append([int(statebreakdown[0][0]), int(statebreakdown[0][2]), int(statebreakdown[0][3])])
                F_labels.append([int(statebreakdown[0][0]), int(statebreakdown[0][2])])
                
        else: #if we are only considering one state
            FmF_labels.append([int(statebreakdown[0][0]), int(statebreakdown[0][2]), int(statebreakdown[0][3])])
            F_labels.append([int(statebreakdown[0][0]), int(statebreakdown[0][2])])
    
    
    #Various error messages 
    if zeromagfielderror != 0:
        if mF_error != 0 and F_error == 0:
            print('mF labelling failed as attempted to label at zero magnetic field: please input a non-zero magnetic field. Returning N,F labels instead.')
        if F_error != 0:
            print('mF, F labelling failed as attempted to label at too high a magnetic field: please input a smaller starting magnetic field value.')
    else:
        if F_error != 0:
            if startingstateindex == 0:
                print('mF, F labelling failed as attempted to label at too high a magnetic field: please input a smaller starting magnetic field value.')
            else:
                print('mF, F labelling failed as attempted to label at too high a magnetic field: please reduce the spacing between your input magnetic field values. Alternatively, start plotting at a small, non-zero magnetic field.')
        if mF_error != 0 and F_error ==0:
            if startingstateindex == 0:
                print('mF labelling failed as attempted to label at too high a magnetic field: please input a smaller starting magnetic field value. Returning N,F labels instead.')
            else:
                print('mF labelling failed as attempted to label at too high a magnetic field: please reduce the spacing between your input magnetic field values. Alternatively, start plotting at a small, non-zero magnetic field. Returning N,F labels instead.')
       
    
    if F_error == 0 and mF_error == 0:
        return FmF_labels
    if F_error == 0 and mF_error != 0:
        return F_labels
    if F_error != 0:
        return None




def solve(H, Nmax, consts, label, B=None):
    """Function that combines the diagonalisation of the Hamiltonian generated by the
    hamiltonian.build method, and the sorting of the states performed by the 
    calculate.sort_smooth function. Generates the sorted set of eigenstates and 
    eigenenergies for the given Hamiltonian.
    
    Inputs:
        H (np.ndarray): Hamiltonian to generate the eigenstates and eigenenergies for,
                        as generated by hamiltonian.build method.
        Nmax (int): Maximum rotational state to consider in the calculations.
        consts (dict): Dictionary of constants for the molecule to be calculated
        B (list/array): Magnetic field values used to calculate the eigenstates
        label (boolean): If True, return the F, mF state labels for the eigenstates
    """
    #generate the eigenenergies and eigenstates
    energies, states = np.linalg.eigh(H)
    #apply sort_smooth to the eigenstates to maintain consistency with the labelling
    energies, states = sort_smooth(energies, states)	
    
    if label:
        if type(B) == float or type(B) == int or type(B) == list or type(B) == np.ndarray:
            #ie if a suitable type of B value is provided
            label_list = label_FmF_states(states, Nmax, consts, B)
        
            return energies, states, label_list  
        
        else:
            print('Please provide magnetic field value(s) to label at.')
            return energies, states, None
            
        
    else:
        return energies, states
import numpy as np
from sympy.physics.wigner import wigner_3j as sympy_wig_3j
from sympy.physics.wigner import wigner_6j as sympy_wig_6j

from sympy import KroneckerDelta
from scipy.linalg import block_diag
import scipy.constants
from scipy.special import sph_harm
import matplotlib.pyplot as plt
'''
This module contains the main code to calculate the hyperfine structure of
singlet -sigma molecules. In usual circumstances most of the functions within
are not user-oriented.

Example:
    Basic usage of this module is for accessing the eigenstates and
    eigenvalues of the molecule in question. This is most easily done
    by combining this module with the user's favourite linear algebra module.
    For instance to find the zero-field hyperfine states of Molecule::

        $ from diatom import Hamiltonian
        $ from np import linalg as la
        $ H0,Hz,HDC,HAC = Hamiltonian.Build_Hamiltonians(5,Molecule)
        $ ev,es = la.eigh(H0)
'''


###############################################################################
# Start by definining constants that are needed for the code                  #
###############################################################################

'''
    Important note!

    All units in this code are SI i.e. elements in the Hamiltonian have units
    of Joules. Outputs will be on the order of 1e-30

'''

h = scipy.constants.h
muN = scipy.constants.physical_constants['nuclear magneton'][0]
bohr = scipy.constants.physical_constants['Bohr radius'][0]
eps0 = scipy.constants.epsilon_0
c = scipy.constants.c
pi = np.pi

DebyeSI = 3.33564e-30
""" Conversion factor from debyes to J/V/m """

###############################################################################
# Functions for the calculations to use                                       #
###############################################################################

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

def wigner_D(l, m, alpha, beta, gamma):
    ''' The Wigner D matrix with labels l and m.

    Calculates the Wigner D Matrix for the given Alpha,beta,gamma in radians.
    The wigner-D matrices represent rotations of angular momentum operators.
    The indices l and m determine the value of the matrix.
    The second index (m') is always zero.

    The input angles are the x-z-x euler angles

    Args:
        l (int) : order of wigner Matrix
        m (float): first index of Wigner Matrix
        alpha,beta,gamma (float) : x,z,x Euler angles in radians
    Returns:
        D (float) : Value of the wigner-D matrix
    '''
    prefactor = np.sqrt((4*np.pi)/(2*l+1))
    function = np.conj(sph_harm(m,l,alpha,beta))
    return prefactor*function


def H_rot(Nmax,consts,I,S):
    ''' Calculates <N,J,F,mF| H_rot |N',J',F', mF'>. The rotational hamiltonian H_rot
        is diagonal in this basis, and this function iterates over N, J, F, mF, N',
        J', F', mF' to build a matrix and evaluate the diagonal elements.
        This function is based off Jesus Aldegunde's FORTRAN 77 code. 
    
        Args:
            Nmax (int): Maximum rotational quantum number to calculate.
            I (float): Nuclear spin. We only consider molecules where one consistuent
                       atom has no nuclear spin, so this is the non zero I of the other.
            S (float): Electronic spin.
            consts (dict): Dictionary of constants for the molecule to be calculated.
    
        Returns:
            H (np.ndarray): Hamiltonian in joules
     '''
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    H_ROT = np.zeros((shape,shape), dtype= complex)
    H_ROT = (np.kron(H_ROT,np.kron(np.identity(Ishape),
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
                                    
                                    
                                    H_ROT[i,j]= consts['Brot']*(N1*(N1+1))*KroneckerDelta(N1,N2)*\
                                        KroneckerDelta(J1,J2)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)+\
                                            consts['Drot']*(N1*(N1+1))**2*KroneckerDelta(N1,N2)*\
                                                KroneckerDelta(J1,J2)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)
                                                
                                            
                                    i+=1
                    i=0
                    j+=1
                    
    #final check for NaN errors, mostly this is due to division by zero or
    # multiplication by a small prefactor. it is safe to set these terms to 0
    H_ROT[np.isnan(H_ROT)] =0

    #return the matrix, in the full coupled basis.

    return H_ROT


def H_hf(Nmax,consts,I,S):
    ''' Calculates <N,J,F,mF| H_hf |N',J',F', mF'> by iterating over N, J, F, mF, N',
        J', F', mF' to build a matrix and evaluate the elements.
        This function is based off Jesus Aldegunde's FORTRAN 77 code. 
    
        Args:
    
            Nmax (int) - maximum rotational quantum number to calculate
            I (float): Nuclear spin. We only consider molecules where one consistuent
                       atom has no nuclear spin, so this is the non zero I of the other.
            S (float): Electronic spin.
            consts (dict): Dictionary of constants for the molecule to be calculated.
    
        Returns:
            H (np.ndarray): Hamiltonian in joules
     '''
     
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    H_HF = np.zeros((shape,shape),dtype= complex)
    H_HF = (np.kron(H_HF,np.kron(np.identity(Ishape),
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
                                    H_HF[i,j]= consts['gamma']/2*(J1*(J1+1)-N1*(N1+1)-S*(S+1))*KroneckerDelta(N1,N2)*\
                                        KroneckerDelta(J1,J2)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)+\
                                            (consts['b']+consts['c']/3)*(np.round((-1+0j)**(J1+J2+F1+N1),0)*\
                                                KroneckerDelta(N1,N2)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)*\
                                                    3/2*((2*J2+1)*(2*J1+1))**0.5*wigner_6j(F1,0.5,J2,1,J1,0.5)*\
                                                        wigner_6j(0.5, J2, N1, J1, 0.5, 1)) +\
                                                10**0.5*consts['c']*np.round((-1+0j)**(J1+F1+N2+0.5),0)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)*\
                                                    ((2*N2+1)*(2*N1+1)*(2*J2+1)*(2*J1+1))**0.5*wigner_6j(F1,0.5,J2,1,J1,0.5)*\
                                                        wigner_6j(J2,J1,1,0.5,1.5,N1)*wigner_6j(N1,N2,2,0.5,1.5,J2)*wigner_3j(N2,2,N1,0,0,0)+\
                                                            consts['CC']*KroneckerDelta(N1,N2)*KroneckerDelta(F1,F2)*KroneckerDelta(mF1,mF2)*\
                                                                np.round((-1+0j)**(2*J1+N2+F2),0)*(3/2)**0.5*((2*J2+1)*(2*J1+1)*N1*(N1+1)*(2*N1+1))**0.5*\
                                                                    wigner_6j(F1,0.5,J2,1,J1,0.5)*wigner_6j(N2,J2,0.5,J1,N1,1)
                                                
                                            
                                    i+=1
                    i=0
                    j+=1
                    
    #final check for NaN errors, mostly this is due to division by zero or
    # multiplication by a small prefactor. it is safe to set these terms to 0
    H_HF[np.isnan(H_HF)] =0

    #return the matrix, in the full coupled basis.
    return H_HF

def H_dc(Nmax,consts,I,S):
    ''' Calculates the effect of the anisotropic DC light shift for a rigid-rotor
        like molecule, , <N,J,F,mF| H_dc |N',J',F', mF'>.
    
        This function  based on Jesus Aldegunde's FORTRAN 77 code and iterates 
        over N, J, F, mF, N', J', F', mF' to build a matrix and evaluate the elements. 
    
        Args:
            Nmax (int) - maximum rotational quantum number to calculate
            I (float): Nuclear spin. We only consider molecules where one consistuent
                       atom has no nuclear spin, so this is the non zero I of the other.
            S (float): Electronic spin.
            constst (dict): Dictionary of constants for the molecule to be calculated.
    
        Returns:
            H (np.ndarray): Hamiltonian in joules
     '''
     
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    HDC = np.zeros((shape,shape),dtype= complex)
    HDC = (np.kron(HDC,np.kron(np.identity(Ishape),
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
                                    HDC[i,j] = -1*consts['d0']*np.round((-1+0j)**(F2-mF2+J1+J2+F1+N2+1))*\
                                        ((2*F1+1)*(2*F2+1)*(2*J1+1)*(2*J2+1))**0.5*((2*N1+1)*(2*N2+1))**0.5*(-1)**N2*wigner_3j(N2,1,N1,0,0,0)*\
                                           wigner_3j(F2,1,F1,-mF2,0,mF1)*wigner_6j(J1,F1,I,F2,J2,1)*wigner_6j(N1,J1,0.5,J2,N2,1)
                                    
                                    i+=1
                    i=0
                    j+=1
    #final check for NaN errors, mostly this is due to division by zero or 
    # multiplication by a small prefactor. it is safe to set these terms to 0
    HDC[np.isnan(HDC)] =0

    #return the matrix, in the full coupled basis.
    return HDC

def H_ac(Nmax,consts,I,S):
    ''' Calculates the effect of the anisotropic light shift for a rigid-rotor
        like molecule, , <N,J,F,mF| H_ac |N',J',F', mF'>.
    
        This function  based on Jesus Aldegunde's FORTRAN 77 code and iterates 
        over N, J, F, mF, N', J', F', mF' to build a matrix and evaluate the elements. 
    
        Args:
    
            Nmax (int) - maximum rotational quantum number to calculate
            I (float): Nuclear spin. We only consider molecules where one consistuent
                       atom has no nuclear spin, so this is the non zero I of the other.
            S (float): Electronic spin.
            constst (dict): Dictionary of constants for the molecule to be calculated.
    
        Returns:
            H (np.ndarray): Hamiltonian in joules
     '''
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    HAC = np.zeros((shape,shape),dtype= complex)
    HAC = (np.kron(HAC,np.kron(np.identity(Ishape),
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
                                    M=mF2-mF1  
                                    HAC[i,j]= -1*consts['alpha_2']*\
                                        (wigner_D(2,M,0,consts['Beta'],0)*\
                                         ((-1)**(N1+N2)+1)*np.round((-1+0j)**(F2-mF2+F1-J2+J1+I+0.5),0)*((2*F1+1)*(2*F2+1))**0.5*\
                                                ((2*N1+1)*(2*N2+1)*(2*J1+1)*(2*J2+1))**0.5*wigner_6j(J2,F2,I,F1,J1,2)*\
                                                wigner_3j(F2,2,F1,-mF2,M,mF1)*wigner_3j(J1,0.5,N1,-0.5,0.5,0)*\
                                                wigner_3j(J2,0.5,N2,-0.5,0.5,0)*wigner_3j(J2,2,J1,-0.5,0,0.5)) -\
                                        consts['alpha_1']*(wigner_D(2,M,0,consts['Beta'],0)*((-1)**(N1+N2)+1)*np.round((-1+0j)**(F2-mF2+F1-J2+J1+I+0.5),0)*\
                                                        ((2*F1+1)*(2*F2+1))**0.5*((2*N1+1)*(2*N2+1)*(2*J1+1)*(2*J2+1))**0.5*\
                                                        wigner_6j(J2,F2,I,F1,J1,1)*wigner_3j(F2,1,F1,-mF2,M,mF1)*wigner_3j(J1,0.5,N1,-0.5,0.5,0)*\
                                                    wigner_3j(J2,0.5,N2,-0.5,0.5,0)*wigner_3j(J2,1,J1,-0.5,0,0.5)) -\
                                            consts['alpha_0']*(((-1)**(N1+N2)+1)*np.round((-1+0j)**(F2-mF2+F1-J2+J1+I+0.5),0)*\
                                                    ((2*F1+1)*(2*F2+1))**0.5*((2*N1+1)*(2*N2+1)*(2*J1+1)*(2*J2+1))**0.5*\
                                                    wigner_6j(J2,F2,I,F1,J1,0)*wigner_3j(F2,0,F1,-mF2,M,mF1)*wigner_3j(J1,0.5,N1,-0.5,0.5,0)*\
                                                wigner_3j(J2,0.5,N2,-0.5,0.5,0)*wigner_3j(J2,0,J1,-0.5,0,0.5))                                            
                                    i+=1
                    i=0
                    j+=1
                    
    #final check for NaN errors, mostly this is due to division by zero or
    # multiplication by a small prefactor. it is safe to set these terms to 0
    HAC[np.isnan(HAC)] =0

    #return the matrix, in the full coupled basis.
    return HAC

def H_B(Nmax,consts,I,S):
    ''' Calculates the Zeeman shift for a molecule in a magnetic field, <N,J,F,mF| H_zee |N',J',F', mF'>
    
        This function  based on Jesus Aldegunde's FORTRAN 77 code and iterates 
        over N, J, F, mF, N', J', F', mF' to build a matrix and evaluate the elements. 
    
        Args:
    
            Nmax (int) - maximum rotational quantum number to calculate
            I (float): Nuclear spin. We only consider molecules where one consistuent
                       atom has no nuclear spin, so this is the non zero I of the other.
            S (float): Electronic spin.
            constst (dict): Dictionary of constants for the molecule to be calculated.
    
        Returns:
            H (np.ndarray): Hamiltonian in joules
     '''
    Ishape = int(2*I+1)
    Sshape = int(2*S+1)
    shape = np.sum(np.array([2*x+1 for x in range(0,Nmax+1)]))
    HB = np.zeros((shape,shape),dtype= complex)
    HB = (np.kron(HB,np.kron(np.identity(Ishape),
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
                                    
                                    
                                    HB[i,j]= (consts['MuS']+consts['MuL'])*((-1+0j)**(F2-mF2+2*J2+F1+N1+1)*KroneckerDelta(N1,N2)*\
                                                    (3/2*(2*F1+1)*(2*F2+1)*(2*J1+1)*(2*J2+1))**0.5*\
                                                        wigner_6j(J1,F1,0.5,F2,J2,1)*wigner_6j(0.5,J2,N1,J1,0.5,1)*\
                                                            wigner_3j(F2,1,F1,-mF2,0,mF1)) +\
                                        consts['MuL']*(np.round((-1+0j)**(F1+F2-mF2+2*J2+N1+N2+0.5),0)*\
                                            ((2*F1+1)*(2*F2+1)*(2*J1+1)*(2*J2+1)*(2*N1+1)*(2*N2+1))**0.5*\
                                                wigner_3j(F2,1,F1,-mF2,0,mF1)*wigner_6j(J1,F1,I,F2,J2,1)*\
                                                    np.sum((-1)**(O)*\
                                                              wigner_3j(J2,1,J1,-O,0,O)*wigner_3j(J2,S,N2,O,-O,0)*\
                                                                  wigner_3j(J1,S,N1,O,-O,0)*O for O in [-0.5,0.5]))+\
                                            consts['MuR']*((-1+0j)**(F1-mF1+J1+J2+F1+N1+3)*\
                                                            ((2*F1+1)*(2*F2+1)*(2*J1+1)*(2*J2+1))**0.5*\
                                                                wigner_6j(J1,F1,0.5,F2,J2,1)*wigner_6j(N1,J1,0.5,J2,N2,1)*\
                                                                    wigner_3j(F1,1,F2,-mF1,0,mF2)) +\
                                            consts['MuN']*((-1+0j)**(2*F1-mF1+J1+3/2)*\
                                                            (3/2*(2*F1+1)*(2*F2+1))**0.5*KroneckerDelta(J1,J2)\
                                                                *wigner_6j(0.5,F2,J2,F1,0.5,1)*\
                                                                    wigner_3j(F1,1,F2,-mF1,0,mF2))
                                            
                                    i+=1
                    i=0
                    j+=1
    #final check for NaN errors, mostly this is due to division by zero or
    # multiplication by a small prefactor. it is safe to set these terms to 0
    HB[np.isnan(HB)] =0

    #return the matrix, in the full coupled basis.

    return HB



# This is the main build function and one that the user will actually have to use.
def build(Nmax,constants,zeeman=False,Edc=False,ac=False):
    ''' Return the hyperfine hamiltonian.

        This function builds the hamiltonian matrices for evaluation so that
        the user doesn't have to rebuild them every time and we can benefit from
        np's ability to do distributed multiplication.

        Args:
            Nmax (int) - Maximum rotational level to include
            Constants (Dictionary) - Dictionary of constants for the molecule to be calculated.
            B,EDC,AC (Boolean) - Switches for turning off parts of the total Hamiltonian.
                                      This can save significant time on calculations where DC 
                                      and AC fields are not required due to nested for loops

        Returns:
            H0,HB,HDC,HAC (np.ndarray): Each of the terms in the Hamiltonian.
    '''
    
    I = constants['I']
    S = constants['S']


    H0 = H_rot(Nmax, constants, I, S) + H_hf(Nmax,constants,I,S) 
    if zeeman:
        HB = H_B(Nmax,constants,I,S)
    else:
        HB =0.
    if Edc:
        Hdc = H_dc(Nmax,constants,I,S)
    else:
        Hdc =0.
    if ac:
        Hac = H_ac(Nmax,constants,I,S) 
    else:
        Hac =0.
    return H0,HB,Hdc,Hac



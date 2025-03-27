import scipy.constants
import numpy as np

###############################################################################
# Molecular Constants
###############################################################################
#Here are some starting dictionaries for various molecules. 
#References given, but check up to date if precision needed!
e = scipy.constants.e
h = scipy.constants.h
muN = scipy.constants.physical_constants['nuclear magneton'][0]
muB = scipy.constants.physical_constants['Bohr magneton'][0]
bohr = scipy.constants.physical_constants['Bohr radius'][0]
eps0 = scipy.constants.epsilon_0
a0 = scipy.constants.physical_constants['atomic unit of length'][0]
c = scipy.constants.c
DebyeSI = 3.33564e-30
inv_cm_to_Hz = 100*c

CaF = {"I":0.5,
       "S":0.5,
       "d0": 1.02041*10**-29,#3.07*DebyeSI,
       "Brot":10267.5387*10**6*h,
       "Drot":0.01406*10**6*h,
       "gamma": 39.65891*10**6*h,
       "b" : 109.1839*10**6*h,
       "c": 40.119*10**6*h,
       "CC": 28.76*10**3*h,
       "MuS": 2.002*muB,
       "MuL": -1.86*10**-3*muB, #10.1103/PhysRevLett.124.063001
       "MuN": 5.585*muN,
       "MuR": -5.13*10**-5*muB, #10.1103/PhysRevLett.124.063001
       "alpha_0":1.4*10**-3*h, #h*Hz/(W/m^2) at 780 nm calculated in Caldwell2020 https://arxiv.org/abs/1910.10689 NB these are alpha/2ce_0
       "alpha_1":3*10**-5*h,
       "alpha_2":-8*10**-4*h,
       "Beta":0,#angle of polarisation for polarisability function
       #Transition dipole moments
       "XXvibmu" : 0.0999*e*a0, #
       "XAmu" : 2.03*1e-29,#2.3439*e*a0, #P. J. Dagdigian, H. W. Cruse, and R. N. Zare, J. Chem. Phys. 60, 2330 (1974).
       "XBmu" : 1.71*e*a0, #   P. J. Dagdigian, H. W. Cruse, and R. N. Zare, J. Chem. Phys. 60, 2330 (1974).
       #FC factors from M. Pelegrini, C. S. Vivacqua, O. Roberto-Neto, F. R. Ornellas, and F. B. Machado.  Braz. J. Phys. 35.4 A (2005), pp. 950–956. X-A
       # X-B M. Dulick, P. F. Bernath, and R. W. Field. Can. J. Phys. 58.5 (1980), pp. 703–712.
       "XAfc" : [[0.964,0.036,0],[0.035,0.895,0.07],[0.001,0.065,0.83]],
       "XBfc" : [[9.992*10**-1,7.27*10**-4,3.809*10**-5],[7.396*10**-4,9.973*10**-1,1.814*10**-3],[2.473*10**-5,1.873*10**-3,9.945*10**-1]],
       #Terms for vibrational transitions [v=0,1,2] Journal of Molecular Spectroscopy 197, 289–296 (1999)
       #Only need for v=0 in X and B, and then V=0,1 in A
       #In cm^-1
       "XT" : np.array([0,582.8478,1159.9473]), #electronic splitting
       "XB" : np.array([0.34248818,0.34005359, 0.33762897]),#rotational constants
       "AT" : np.array([16529.653,17118.103,17700.49]), #electronic splitting
       "AA" : np.array([71.491, 71.614, 71.737]), #Lambda doubling
       "BT" : np.array([18832.031, 19398.254, 19958.276]), #electronic splitting
       "BB" : np.array([0.341058, 0.338469, 0.335903]), # rotational constants
       }

BaF = {"I":0.5,
       "S":0.5,
       "d0": 1.05739e-29,#3.17*DebyeSI,
       "Brot":6473.9586572e6*h,# from PRA 94, 063415 (2016) https://journals.aps.org/pra/pdf/10.1103/PhysRevA.94.063415
       "Drot":5.5296816e3*h,#1085*10**7*29979.2458*h,
       "gamma": 80.95547199999999e6*h,
       "b" : 63.509e6*h,
       "c": 8.224e6*h,
       "CC": 0,
       "MuS": 2.002*muB,
       "MuL": -0.028*muB,#Where did I get this?
       "MuN": 5.585*muN,
       "MuR" : 0,
       #"alpha_0":,
       #"alpha_1":,
       #"alpha_2":,
       "Beta":0 #angle of polarisation for polarisability function
       }


SrF = {"I": 0.5,
       "S": 0.5,
       "d0": 3.4963*DebyeSI,
       "Brot": 7487.6*10**6*h,
       "Drot":0.0075*10**6*h,
       "gamma": 74.795*10**6*h,
       "b" : 97.0827*10**6*h,
       "c":30.2675*10**6*h,
       "CC": 0.0023*10**6*h,
       "MuS": 2.002*muB, # g_s * muB
       "MuN": 5.585*muN,
       "MuL": -4.97*10**-3*muB, #10.1103/PhysRevLett.124.063001
       "MuR": -4.77*10**-5*muB, #10.1103/PhysRevLett.124.063001
       "Beta":0, #angle of polarisation for polarisability function
       #Transition dipole moments
       "XXvibmu" : 0.0999*e*a0, #MISSING
       "XAmu" : 2.45*e*a0, #P. J. Dagdigian, H. W. Cruse, and R. N. Zare, J. Chem. Phys. 60, 2330 (1974).
       "XBmu" : 1.94*e*a0, #   Berg, L. E. et al. Chem. Phys. Lett. 248, (1996)
       #FC factors from Hao et al. J. Chem. Phys. 151, 034302 (2019)
       "XAfc" : [[0.9789,0.02054,4*10**-4],[0.02102,0.9377,0.03969],[2.72*10**-5,4.158*10**-2,8.978*10**-1]],
       "XBfc" : [[9.961*10**-1,3.866*10**-3,3.604*10**-6],[3.856*10**-3,9.881*10**-1,8*10**-3],[1.343*10**-5,7.959*10**-3,9.796*10**-1]],
       #Terms for vibrational transitions [v=0,1,2] - for v=0 find in J. Barry thesis table 2.6
       "XT" : np.array([0,0,0])*inv_cm_to_Hz, #electronic splitting
       "XB" : np.array([7.4876,7.44123,7.395])*10**9,#rotational constants J. Barry thesis table 2.8
       
       "AT" : np.array([15072.09,0,0])*inv_cm_to_Hz, #electronic splitting
       "AA" : np.array([281.46138,0,0])*inv_cm_to_Hz, #Lambda doubling
       
       "BT" : np.array([17267.41,0,0])*inv_cm_to_Hz, #electronic splitting
       "BB" : np.array([0.249396,0,0])*inv_cm_to_Hz, # rotational constants
       }
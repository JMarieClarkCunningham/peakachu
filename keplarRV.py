from __future__ import print_function, division
import numpy as np
from matplotlib import rc
from PyAstronomy import pyasl
from astropy import constants as const
from astropy import units as u
'''
This function calculates synthetic radial velocity curves using the form:
RV = gamma + K * (e * cos(w) + cos(theta(t) + w)), Keplarian MarkleyKESolver

THIS IS A GREAT DOCSTRING BUT IT TRADITIONALLY GOES INSIDE THE FUNCTION :)
__________________________________________________________________________

                               Parameters:
__________________________________________________________________________
phases :  `np.array`               Array of phases, typically from 0 to 1
 omega : `int`                     Argument of periapsis     [degrees] ... REALLY NOW
     K : `int`                     Velocity semi-amplitude   [km/s]
   ecc : `int`                     Eccentricity of the system
 gamma : `int`                     Systemic Velocity         [km/s]
__________________________________________________________________________

                                Returns:
__________________________________________________________________________
KepRVCurve  :  a 1D array of RVs corresponding to phases
'''
def getkepRV(phases, ecc, omega, gamma, K):
    ks = pyasl.MarkleyKESolver()                  # Initializes the Keplarian Solver
    #theta = []                                    # Makes empty list for theta  NOT NEEDED
    KepRVCurve = []                               # Makes empty list for RV curve data

##  Run Keplarian solver MarkleyKESolver() to find the mean anomaly as a function
##  of phase [M_ph], the eccentric anomaly E_ph, and the True anomaly, Tru, or
##  theta, the missing piece of the RV curve calculation
    for phase in phases:                          # Loop over phases
        MeAn = 2. * np.pi * phase                 # Solve for mean anomaly, M(phase)
        EcAn = ks.getE(MeAn, ecc)                   # Solve for eccentric anomaly, E

        # Compute the true anomaly
        cosTru = (np.cos(EcAn) - ecc) / (1 - ecc*np.cos(EcAn))
        sinTru = (np.sqrt(1 - ecc**2) * np.sin(EcAn)) / (1 - ecc*np.cos(EcAn))
        Tru = np.arctan2(sinTru, cosTru)

        RV = gamma + K * (ecc * np.cos(omega) + np.cos(Tru + omega))

##   How exciting, now we've got all the pieces needed to solve for the Keplarian
##   projected RV curve!
        #result = np.add(system,costhetaplusomega)  # NO THIS WON'T WORK
        KepRVCurve.append(RV)               # add RV curve to the empty list

    return KepRVCurve

def getK(Porb, M1, M2, incl, ecc):
    '''
    Original version by Joni is in comments,
    Slightly condensed edits by Meredith appear as-is
    '''
#     M1_app = M1.to(u.kg) + M2.to(u.kg)
#     P_orb = Porb.to(u.second)
#     twopiG = 2. * cons.pi * cons.G
#     twopiGoverP = np.divide(twopiG, P_orb)
#     part_1 = twopiGoverP**(1/3)
#     M2sin_i = M2.to(u.kg) * np.sin(incl.to(u.deg))
#     M1pow2third = (M1.to(u.kg) + M2.to(u.kg))**(2/3)
#     part_2 = M2sin_i / M1pow2third
#     part_a = part_1 * part_2
#     part_n = np.sqrt(1 - ecc**2)
#     part_3 = 1 / part_n
# 
#     K = np.multiply(part_a, part_3)
#     return(K)

    twopiG = 2. * np.pi * const.G
    twopiGoverP = np.divide(twopiG, Porb)
    part_1 = twopiGoverP**(1/3)

    M2sin_i = M2 * np.sin(incl)
    M1pow2third = (M1 + M2)**(2/3)
    part_2 = M2sin_i / M1pow2third
    
    part_3 = 1. / np.sqrt(1 - ecc**2)

    K = np.multiply(part_1*part_2, part_3)
    
    # ALL IN ONE LINE MWAHAHA -MR
    K1_calc = ((2*np.pi*const.G / Porb)**(1/3)) * ((M2*np.sin(incl)) / (M1+M2)**(2/3)) * (1/np.sqrt(1 - ecc**2))

    return K
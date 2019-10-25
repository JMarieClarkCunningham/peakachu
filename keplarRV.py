from __future__ import print_function, division
import numpy as np
from matplotlib import rc
from PyAstronomy import pyasl
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
        #system = gamma + K * (ecc * np.cos(omega))  # THIS NEEDS TO HAPPEN LAST

        MeAn = 2. * np.pi * phase                 # Solve for mean anomaly, M(phase)
        EcAn = ks.getE(MeAn, ecc)                   # Solve for eccentric anomaly, E
        TruTop = (np.cos( EcAn - ecc ))             # Find the denominator of arctan2
        TruBot = ((np.sqrt(1 - ecc*ecc)) * np.sin(EcAn))    # Find the numerator of arctan2
        Tru = np.arctan2(TruTop, TruBot)          # Solve for true anomaly, Tru
        costhetaplusomega = np.cos(Tru + omega)
        #theta.append(Tru)                         # theta is the true anomaly  NOT NEEDED

        RV = gamma + K * (ecc * np.cos(omega) + np.cos(Tru + omega))

##   How exciting, now we've got all the pieces needed to solve for the Keplarian
##   projected RV curve!
        #result = np.add(system,costhetaplusomega)  # NO THIS WON'T WORK
        KepRVCurve.append(RV)               # add RV curve to the empty list

    return KepRVCurve

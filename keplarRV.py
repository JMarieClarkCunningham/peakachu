from __future__ import print_function, division
import numpy as np
from matplotlib import rc
from PyAstronomy import pyasl
'''
This function calculates synthetic radial velocity curves using the form:
RV = gamma + K * (e * cos(w) + cos(theta(t) + w), Keplarian MarkleyKESolver
__________________________________________________________________________

                               Parameters:
__________________________________________________________________________
phases :  np.arange(0, 1, 0.001)   Array from 0 to 1 in 0.001 step size
 omega : `int`                     Argument of periapsis     [degrees]
     K : `int`                     Velocity semi-amplitude   [km/s]
   ecc : `int`                     Eccentricity of the system
 gamma : `int`                     Systemic Velocity         [km/s]
__________________________________________________________________________

                                Returns:
__________________________________________________________________________
KepRVCurve  :  a 2D array of phases and RV: RV(phases)
'''
def getkepRV(phases, ecc, omega, gamma, K):
    ks = pyasl.MarkleyKESolver()                  # Initializes the Keplarian Solver
    theta = []                                    # Makes empty list for theta
    KepRVCurve = []                               # Makes empty list for RV curve data

##  Run Keplarian solver MarkleyKESolver() to find the mean anomaly as a function
##  of phase [M_ph], the eccentric anomaly E_ph, and the True anomaly, Tru, or
##  theta, the missing piece of the RV curve calculation
    for phase in phases:                          # Loop over phases
        system = gamma + K * (ecc * np.cos(omega))
        MeAn = 2. * np.pi * phase                 # Solve for mean anomaly, M(phase)
        EcAn = ks.getE(MeAn, ecc)                   # Solve for eccentric anomaly, E
        TruTop = (np.cos( EcAn - ecc ))             # Find the denominator of arctan2
        TruBot = ((np.sqrt(1 - ecc*ecc)) * np.sin(EcAn))    # Find the numerator of arctan2
        Tru = np.arctan2(TruTop, TruBot)          # Solve for true anomaly, Tru
        costhetaplusomega = np.cos(Tru + omega)
        theta.append(Tru)                         # theta is the true anomaly
##   How exciting, now we've got all the pieces needed to solve for the Keplarian
##   projected RV curve!
        result = np.add(system,costhetaplusomega)
        KepRVCurve.append(result)               # add RV curve to the empty list

    return KepRVCurve

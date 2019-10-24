from __future__ import print_function, division
import numpy as np
from matplotlib import rc
from PyAstronomy import pyasl
'''
This function calculates synthetic radial velocity curves using the form:
RV = gamma + K * (e * cos(w) + cos(theta(t) + w), Keplarian MarkleyKESolver
                               Parameters:
__________________________________________________________________________
omega : `int`
      Argument of periapsis [degrees]
K :  `int`
e : `int`
      eccentricity of the system
phases : `list`
      2D list with phases of the system
'''
def getE(phases, e, K):
    ks = pyasl.MarkleyKESolver()
####################### THIS IS WHERE YOU CHANGE THINGS ########################
    bfinfile = '6864859/6864859Outfile.txt'     # Where the BF_RV phase data lives
    bfrvdata = np.loadtxt(bfinfile, comments='#', unpack=True)
    phases = bfrvdata[1]                        # len(BF_ph) = 19 [BF phases]
    gamma=6.54; K=93945.5                       # Keblat (Windemuth et al. 2019)
    e = 0.63462; omega=3.101696                 # (Cunningham et al. 2019)
    theta = []; M_ph = []; E_ph = [];           # theta, the true anomaly
    KepRVCurve = []


##  Run Keplarian solver MarkleyKESolver() to find the mean anomaly as a function
##  of phase [M_ph], the eccentric anomaly E_ph, and the True anomaly, Tru, or
##  theta, the missing piece of the RV curve calculation

                                            # Using phase and e, solve for M
    for phase in enumerate(phases):         # USING BF DATA: FIND M, E, Tru
        M = 6.283 * phases                  # Solve for mean anomaly, M(phase)
        E = ks.getE(M, e)                   # Solve for eccentric anomaly, E
        TruTop = (np.cos( E - e ))          # Find the denominator of arctan2
        TruBot = (np.sqrt(1-e^2)*np.sin(E)) # Find the numerator of arctan2
        Tru = np.arctan2(TruTop, TruBot)    # Solve for true anomaly, Tru
        M_ph.append(M)                      # Append empty lists
        E_ph.append(E)                      # to include what we used the keplarian
        theta.append(Tru)                   # solver to determine

##   How exciting, now we've got all the pieces needed to solve for the Keplarian
##   projected RV curve!

        result = gamma + K * (e * np.cos( omega ) + np.cos( theta + omega ))
        kepRV.append(result)
        kepcurve = kepRV(theta[i], gamma=6.54, K=93945.5, e = 0.63462, omega=3.101696)
        KepRVCurve.append(kepcurve)
    print(kepcurve)


    return KepRVCurve

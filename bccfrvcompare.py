from __future__ import print_function, division
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
from PyAstronomy import pyasl


'''
Yay! If you've made it this far, you've got some peaks fitted with some gaussians
maybe even have errors on those gaussians, and now, you're looking to plot some
RVs. Well, here's lookin at you kid, cause this here program, it does just that!
'''

##  Let's read in fit estimates (RV estimates) from peakfit.py [CCF RVs]      ##
##     _____________________________________________________________          ##
##     [0]:Date (JD), [1]:RV 1, [2]:RV 1 err, [3]:RV 1, [4]:RV 1 err          ##
ccfrvdata = np.loadtxt('6864859/6864859ccfpeakout.txt', comments='#', unpack=True)
jd = ccfrvdata[0]
ccf_rv1 = ccfrvdata[1]; ccf_rv1_err = ccfrvdata[2]
ccf_rv2 = ccfrvdata[3]; ccf_rv2_err = ccfrvdata[4]

##  Next, let's need to skip any [CCF RVs] that have rv_err = 0, as these are ##
##  where the Gaussfitter in peakfit.py [LevMarLSQFitter] failed to fit the   ##
##  peaks in the CCF ApStar spectra, 6 cases of this in KIC 6864859.          ##

for idx, err in enumerate(ccf_rv1_err):
    if err == 0:
        ccf_rv1[idx] = None
        ccf_rv1_err[idx] = None
for idx, err in enumerate(ccf_rv2_err):
    if err == 0:
        ccf_rv2[idx] = None
        ccf_rv2_err[idx] = None


##  Let's  read in fit estimates (RV estimates) from BF_python.py [BF RVs]
##  _____________________________________________________________________________________
##  |[0]:Date (BJD)[1]:phase[3]:BF RV 1[4]:BF RV 1 err[5]:BF RV 2[6]:BF RV 2 err|
bfinfile = '6864859/6864859Outfile.txt'
bfrvdata = np.loadtxt(bfinfile, comments='#', unpack=True)
bjd = bfrvdata[0]
phase = bfrvdata[1]
bf_rv1 = bfrvdata[3]; bf_rv1_err = bfrvdata[4]
bf_rv2 = bfrvdata[5]; bf_rv2_err = bfrvdata[6]

##  Next, let's need to skip any [BF RVs] that have rv_err = 0, as these are  ##
##  where the Gaussparty in BF_python.py failed to fit the peaks in the BF    ##
##  ApVisit spectra. Also, just an FYI, not to tout on about how BF >> CCF    ##
##  in RV extraction... but ...                                               ##
##  KIC 6864859: NO RVS HAD TO BE SKIPPED BECAUSE THEIR FITTING FAILED.       ##

for idx, err in enumerate(bf_rv1_err):
    if err == 0:
        bf_rv1[idx] = None
        bf_rv1_err[idx] = None
for idx, err in enumerate(bf_rv2_err):
    if err == 0:
        bf_rv2[idx] = None
        bf_rv2_err[idx] = None

### print(ccf_rv1, ccf_rv1_err, ccf_rv2, ccf_rv2_err)  # Sanity √ print-party
### print(bf_rv1, bf_rv1_err, bf_rv2, bf_rv2_err)      # Sanity √ print-party

''' Systematic Effects, Time, APOGEE RV offset, Barycentric Velocities, etc '''

## Get the systemic RV of the star according to APOGEE
apSystemic = 82  ## TODO: FIND ACTUAL SYSTEMIC RV PARAMETER !

##  Recall: The Barycentric Julian Date (BJD) is the Julian Date (JD) corrected
##  for differences in the Earth's position with respect to the barycentre of
##  the Solar System
dateoffset = 2454833.                                  # BJD - dateoffset = JD
                                                       # for Unfolded RV vs time
rv1line = np.isfinite(bf_rv1)                          # for the dotted RV lines
rv2line = np.isfinite(bf_rv2)

##  Doubling the arrays to allow this folded phase functionality [if the arrays
##  aren't doubled, we don't have RV information to plot all the way to phase 2]

phase_double = np.concatenate((np.array(phase),np.array(phase)+1.0), axis=0)
bf_rv1_double = np.concatenate((bf_rv1,bf_rv1), axis=0)
bf_rv1_err_double = np.concatenate((bf_rv1_err,bf_rv1_err), axis=0)
bf_rv2_double = np.concatenate((bf_rv2,bf_rv2), axis=0)
bf_rv2_err_double = np.concatenate((bf_rv2_err,bf_rv2_err), axis=0)

ccf_rv1_double = np.concatenate((ccf_rv1,ccf_rv1), axis=0)
ccf_rv1_err_double = np.concatenate((ccf_rv1_err,ccf_rv1_err), axis=0)
ccf_rv2_double = np.concatenate((ccf_rv2,ccf_rv2), axis=0)
ccf_rv2_err_double = np.concatenate((ccf_rv2_err,ccf_rv2_err), axis=0)

##  Let's compute the RV curve (modeled in Cunningham et al. 2019) with the   ##
##  Keplarian solver keblat, and add it to the second subplot!                ##
##  RV = gamma + K*(e * cos(w) + cos(omega(t) + w)                            ##

##  We'll need to solve Kepler's equation to turn our timestamps into theta   ##
##  and omega using pyastronomy:MarkleyKESolver()
##  First, let us define some constants we will need for RV curve calculation ##
KIC = '6864859'
ks = pyasl.MarkleyKESolver()             # Instantiates the keplarian solver
M = (2*np.pi*phase).any                  # Mean Anomaly [.any troubleshooting MarkleyKESolver]
e = 0.63462                              # eccentricity (Cunningham et al. 2019)

##  Solves Kepler's Equation for a set of mean anomaly and eccentricity. Uses ##
##  the algorithm presented by Markley 1995.                                  ##
print("Eccentric anomaly: ", ks.getE(M, e))

ecosomega = -0.634115                    # ecos(ω) (Cunningham et al. 2019)
omega = 3.101696                         # Yay! Inverse Trig! √ Working [J=MC^2]
#K = # velocity semi amplitude

###            Now, let us set up this figure you came for...                ###
fig = plt.figure(1, figsize=(8, 6))
##  Plot formatting, axis formatting [makin' it look good tho']
timestart = 1710; timeend = 1990
RVmin = -30; RVmax = 60
phasemin = 0.1; phasemax = 0.9
plt.axis([timestart, timeend, RVmin, RVmax])
plt.ylabel("Radial Velocity (km s$^{-1}$)", size='large')

##  Titles for the subplots
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center',
    rotation='vertical', size='large')
fig.text(0.14, 0.13, 'Folded')
fig.text(0.14, 0.55, 'Unfolded')
fig.text(0.45, 0.9, 'KIC 6864859', size='large')

##########################  TOP PLOT subplot(2,1,1): ##########################
#####################  Unfolded RV VS time (BJD-2454833)  #####################
ax2 = plt.subplot(2,1,1)
plt.axis([timestart, timeend, RVmin, RVmax])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset))

##  Dotted lines for both CCF and BF RV curve
plt.plot(bjd[rv1line]-dateoffset, bf_rv1[rv1line] - apSystemic, color="orange",
    mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2line]-dateoffset, bf_rv2[rv2line] - apSystemic, color="orange",
    mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(jd-dateoffset, ccf_rv1, color="skyblue", mfc=None, mec=None, lw=1.5,
    ls=':')
plt.plot(jd-dateoffset, ccf_rv2, color="skyblue", mfc=None, mec=None, lw=1.5,
    ls=':')

##  Plot the BF RV estimates and their errors
for idx, date in enumerate(bjd):                         # BJD - dateoffset = JD
    plt.errorbar(date-dateoffset, bf_rv1[idx] - apSystemic, yerr=bf_rv1_err[idx],
        fmt='ko', color="C0", ecolor="C0", mfc="C0", mec="C0", ms=8, lw=1.5,
        label='BF RV$_{1}$')
    plt.errorbar(date-dateoffset, bf_rv2[idx] - apSystemic, yerr=bf_rv2_err[idx],
        fmt='ko', color="C0", ecolor="C0", mfc="white", mec="C0", ms=8, lw=1.5,
        label='BF RV$_{2}$')

## Plot the CCF RV estimates and their errors
for idx, date in enumerate(jd):
    plt.errorbar(date-dateoffset, ccf_rv1[idx], yerr=ccf_rv1_err[idx], fmt='ko',
        color="C1", ecolor="C1", mfc="C1", mec="C1", ms=8, lw=1.5,
        label='CCF RV$_{1}$')
    plt.errorbar(date-dateoffset, ccf_rv2[idx], yerr=ccf_rv2_err[idx], fmt='ko',
        color="C1", ecolor="C1", mfc="white", mec="C1", ms=8, lw=1.5,
        label='CCF RV$_{2}$')

########################  BOTTOM PLOT subplot(2,1,2): ##########################
#############################  Folded RV VS Phase   ############################
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')

##  Loop through the enumerated doubled phase array and plot BF and CCF RVs
##  This is why the arrays were doubled, lines 69 - 81
for idx, ph in enumerate(phase_double):
    ##  Plot the Folded BF RV estimates and their errors
    plt.errorbar(phase_double[idx], bf_rv1_double[idx] - apSystemic, yerr=bf_rv1_err_double[idx],
        marker='o', color="C1", mec="C1", ecolor="C1", mfc="C1", ms=8, ls='None',
        lw=1.5, label='BF RV$_{1}$')
    plt.errorbar(phase_double[idx], bf_rv2_double[idx] - apSystemic, yerr=bf_rv1_err_double[idx],
        marker='o', color="C1", mec="C1", ecolor="C1", mfc="white", ms=8,
        ls='None', lw=1.5, label='BF RV$_{2}$')

    ##  Plot the Folded CCF RV estimates and their errors
    plt.errorbar(phase_double[idx], ccf_rv1_double[idx], yerr=ccf_rv1_err_double[idx],
        marker='o', color="C0", mec="C0", ecolor="C0", mfc="C0", ms=8, ls='None',
        lw=1.5, label='CCF RV$_{1}$')
    plt.errorbar(phase_double[idx], ccf_rv2_double[idx], yerr=ccf_rv2_err_double[idx],
        marker='o', color="C0", mec="C0", ecolor="C0", mfc="white", ms=8,
        ls='None', lw=1.5, label='CCF RV$_{2}$')
    ##  MNake a legend
    if idx == 0:
        plt.legend(ncol=2, loc=1, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35),
            columnspacing=0.7)
plt.xlabel("Orbital Phase", size='large')

##  Legend formatting [makin' it look good tho']

plt.show()
#plt.savefig('6864859BFCCFRV.png')

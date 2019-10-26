from __future__ import print_function, division
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
#sys.path.append('/peakachu/keplarRV.py')
from keplarRV import getkepRV, getK
from astropy import constants as const
from astropy import units as u

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
#phases = bfrvdata[1]
bf_rv1 = bfrvdata[3]; bf_rv1_err = bfrvdata[4]
bf_rv2 = bfrvdata[5]; bf_rv2_err = bfrvdata[6]

# Compute orbital phases
dateoffset = 2454833  # Kepler time offset
BJD0 = 158.3189733 + dateoffset  # time of primary eclipse from paper
#BJD0 = 2454955.5563  # original BJD0 used in RV modeling
Porb = 40.8778425 * u.day
phases = []
for i in range(0, len(bjd)):
    fracP = (bjd[i] - BJD0) / Porb.value
    if fracP < 0:
        phases.append(1 + (fracP % 1))
        cycle = int(fracP) - 1
    else:
        phases.append((fracP % 1))

phases = [phase - 0.065 if phase > 0.065 else phase - 0.065 + 1 for phase in phases]
# per the paper, periastron passage is 0.065 off from deepest eclipse in phase

## Get the systemic RV of the star according to APOGEE
apSystemic = 82  # empirical RV offset observed between BF and CCF due introduced by APOGEE

synphases = np.arange(0, 1, 0.01)

# Data from Clark Cunningham et al. paper
gamma = 93.9455
M1 = 1.411 * u.Msun
M2 = 1.354 * u.Msun
incl = 1.5415  * u.rad
ecc = 0.63462

# Calculate K properly
K1_calc = getK(Porb, M1, M2, incl, ecc)
print('K1 is', K1_calc.to(u.km / u.s))
K1 = K1_calc.to(u.km / u.s).value
K2 = -1*K1*M1/M2
print('K2 is', K2)

# Calculate omega properly
esinw = -0.0254
ecosw = -0.634115
sinw = esinw / ecc
cosw = ecosw / ecc
#print(np.arccos(cosw), np.arcsin(sinw))  # trig is annoying

# Try a few different omegas to see what happens
# Fun exercise: try values between 0 and 2*pi! it's illustrative!
#for omega, color in zip([3.1, 3.2, np.arccos(cosw)], ['C2', 'C3', 'C4']):
#    keplercurve1 = getkepRV(synphases, ecc, omega, gamma, K1)
#    keplercurve2 = getkepRV(synphases, ecc, omega, gamma, K2)
#    plt.plot(synphases, keplercurve1, marker=None, ls='-', color=color, label='omega='+str(omega))
#    plt.plot(synphases, keplercurve2, marker=None, ls='-', color=color)
#plt.errorbar(np.array(phases), bf_rv1, yerr=bf_rv1_err, marker='o', ls='None', lw=1.5, label='BF RV$_{1}$')
#plt.errorbar(np.array(phases), bf_rv2, yerr=bf_rv2_err, marker='o', ls='None', lw=1.5, label='BF RV$_{2}$')
#plt.legend()
#plt.show()

omega = np.arccos(cosw)  # sure why not

synRV1 = getkepRV(synphases, ecc, omega=omega, gamma=gamma, K=K1)
synRV2 = getkepRV(synphases, ecc, omega=omega, gamma=gamma, K=K2)


###            Now, let us set up this figure you came for...                ###
fig = plt.figure(1, figsize=(10, 8))
##  Plot formatting, axis formatting [makin' it look good tho']
timestart = 1710; timeend = 2005
#RVmin = -30; RVmax = 60
RVmin = 50; RVmax = 150
phasemin = 0.1; phasemax = 0.99
plt.axis([timestart, timeend, RVmin, RVmax])
plt.ylabel("Radial Velocity (km s$^{-1}$)", size='x-large')

##  Titles for the subplots
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center',
    rotation='vertical', size='x-large')
fig.text(0.14, 0.13, 'Folded')
fig.text(0.14, 0.55, 'Unfolded')
fig.text(0.45, 0.9, 'KIC 6864859', size='x-large')

##########################  TOP PLOT subplot(2,1,1): ##########################
#####################  Unfolded RV VS time (BJD-2454833)  #####################
ax2 = plt.subplot(2,1,1)
plt.axis([timestart, timeend, RVmin, RVmax])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset), size='x-large')

plt.errorbar(bjd-dateoffset, bf_rv1, yerr=bf_rv1_err, 
    marker='o', ls=':', ms=8, ecolor="C1", mfc="C1", mec="C1", color="C1")
plt.errorbar(bjd-dateoffset, bf_rv2, yerr=bf_rv2_err,
    marker='o', ls=':', ms=8, ecolor="C1", mfc="white", mec="C1", color="C1")

plt.errorbar(jd-dateoffset, ccf_rv1 + apSystemic, yerr=ccf_rv1_err,
    marker='o', ls=':', ms=8, ecolor="C0", mfc="C0", mec="C0", color="C0")
plt.errorbar(jd-dateoffset, ccf_rv2 + apSystemic, yerr=ccf_rv2_err,
    marker='o', ls=':', ms=8, ecolor="C0", mfc="white", mec="C0", color="C0")
    

########################  BOTTOM PLOT subplot(2,1,2): ##########################
#############################  Folded RV VS Phase   ############################
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')

# Plot the synthetic RV curves
plt.plot(synphases, synRV1, marker=None, color='C4', ls=':')
plt.plot(synphases, synRV2, marker=None, color='C4', ls=':')

##  Loop through the enumerated doubled phase array and plot BF and CCF RVs
##  This is why the arrays were doubled, lines 69 - 81
for idx, ph in enumerate(np.array(phases)):
    ##  Plot the Folded BF RV estimates and their errors [[color='C1' ORANGE FOR BF]]
    plt.errorbar(ph, bf_rv1[idx], yerr=bf_rv1_err[idx],
        marker='o', color="C1", mec="C1", ecolor="C1", mfc="C1", ms=8, ls='None',
        lw=1.5, label='BF RV$_{1}$')
    plt.errorbar(ph, bf_rv2[idx], yerr=bf_rv2_err[idx],
        marker='o', color="C1", mec="C1", ecolor="C1", mfc="white", ms=8, ls='None',
        lw=1.5, label='BF RV$_{2}$')

    ##  Plot the Folded CCF RV estimates and their errors [[color='C0' BLUE FOR CCF]]
    plt.errorbar(ph, ccf_rv1[idx] + apSystemic, yerr=ccf_rv1_err[idx],
        marker='o', color="C0", mec="C0", ecolor="C0", mfc="C0", ms=8, ls='None',
        lw=1.5, label='CCF RV$_{1}$')
    plt.errorbar(ph, ccf_rv2[idx] + apSystemic, yerr=ccf_rv2_err[idx],
        marker='o', color="C0", mec="C0", ecolor="C0", mfc="white", ms=8, ls='None',
        lw=1.5, label='CCF RV$_{2}$')
    ##  Make a legend
    if idx == 0:
        plt.legend(ncol=2, loc=1, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35),
            columnspacing=0.7)

plt.xlabel("Orbital Phase", size='x-large')
plt.show()


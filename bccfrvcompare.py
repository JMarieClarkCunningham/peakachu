import numpy as np
import matplotlib.pyplot as plt
'''
Yay! If you've made it this far, you've got some peaks fitted with some gaussians
maybe even have errors on those gaussians, and now, you're looking to plot from
RVs
'''

##  Let's read in fit estimates (RV estimates) from peakfit.py [CCF RVs]
##  _____________________________________________________________
##  [0]:Date (JD), [1]:RV 1, [2]:RV 1 err, [3]:RV 1, [4]:RV 1 err
ccfrvdata = np.loadtxt('6864859/6864859ccfpeakout.txt', comments='#', unpack=True)
jd = ccfrvdata[0]
ccf_rv1 = ccfrvdata[1]; ccf_rv1_err = ccfrvdata[2]
ccf_rv2 = ccfrvdata[3]; ccf_rv2_err = ccfrvdata[4]

##  Next, let's need to skip any [CCF RVs] that have rv_err = 0, as these are where the
##  Gaussfitter in peakfit.py [LevMarLSQFitter] failed to fit the peaks in the CCF
##  ApStar spectra.

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
##  [0]:Date (BJD), [1]:phase, [3]:BF RV 1, [4]:BF RV 1 err, [5]:BF RV 2, [6]:BF RV 2 err
bfinfile = '6864859/6864859Outfile.txt'
bfrvdata = np.loadtxt(bfinfile, comments='#', unpack=True)
bjd = bfrvdata[0]
phase = bfrvdata[1]
bf_rv1 = bfrvdata[3]; bf_rv1_err = bfrvdata[4]
bf_rv2 = bfrvdata[5]; bf_rv2_err = bfrvdata[6]

##  Next, let's need to skip any [BF RVs] that have rv_err = 0, as these are where the
##  Gaussparty in BF_python.py failed to fit the peaks in the BF ApVisit spectra.
##  Also, just an FYI, not to tout on about how BF >> CCF in RV extraction, but in
##  KIC 6864859: NO RVS HAD TO BE SKIPPED BECAUSE THEIR FITTING FAILED.

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

## Now, let's set up this figure you came for...
fig = plt.figure(1, figsize=(13, 9))

## Plot the BF RV estimates and their errors
for idx, date in enumerate(bjd):                         # BJD - dateoffset = JD
    plt.errorbar(date-dateoffset, bf_rv1[idx], yerr=bf_rv1_err[idx], fmt='ko', color="darkorange", mfc="darkorange", mec="darkorange", ms=10, lw=1.5, label='BF RV$_{1}$')
    plt.errorbar(date-dateoffset, bf_rv2[idx], yerr=bf_rv2_err[idx], fmt='ko', color="darkorange", mfc="white", mec="darkorange", ms=10, lw=1.5, label='BF RV$_{2}$')

## Plot the CCF RV estimates and their errors
for idx, date in enumerate(jd):
    plt.errorbar(date-dateoffset, ccf_rv1[idx], yerr=ccf_rv1_err[idx], fmt='ko', color="dodgerblue", mfc="dodgerblue", mec="dodgerblue", ms=10, lw=1.5, label='CCF RV$_{1}$')
    plt.errorbar(date-dateoffset, ccf_rv2[idx], yerr=ccf_rv2_err[idx], fmt='ko', color="dodgerblue", mfc="white", mec="dodgerblue", ms=10, lw=1.5, label='CCF RV$_{2}$')

#plt.legend(ncol=2, loc=1, numpoints=1, frameon=True)
rv1line = np.isfinite(bf_rv1)
rv2line = np.isfinite(bf_rv2)

## Stuff to make the plot look nice
timestart = 1920; timeend = 1990
RVmin = -50; RVmax = 150
plt.axis([timestart, timeend, RVmin, RVmax])
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset))
plt.ylabel("Radial Velocity (km s$^{-1}$)")
## Add some dotted lines, see if any primary/secondary
## switchy witchy stuff is going on [J=MC^2]
plt.plot(bjd[rv1line]-dateoffset, bf_rv1[rv1line], color="orange", mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2line]-dateoffset, bf_rv2[rv2line], color="orange", mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(jd-dateoffset, ccf_rv1, color="skyblue", mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(jd-dateoffset, ccf_rv2, color="skyblue", mfc=None, mec=None, lw=1.5, ls=':')

plt.show()
plt.savefig('6864859BFvsCCFRVs.png')

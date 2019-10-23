import numpy as np
import matplotlib.pyplot as plt
'''
Yay! If you've made it this far, you've got some peaks fitted with some gaussians
maybe even have errors on those gaussians, and now, you're looking to plot from
RVs
'''

# Let's read in fit estimates (RV estimates) from peakfit.py [CCF RVs]
# [0]:Date (JD), [1]:RV 1, [2]:RV 1 err, [3]:RV 1, [4]:RV 1 err
ccfrvdata = np.loadtxt('6864859/6864859ccfpeakout.txt', comments='#', unpack=True)
jd = ccfrvdata[0]
ccf_rv1 = ccfrvdata[1]; ccf_rv1_err = ccfrvdata[2]
ccf_rv2 = ccfrvdata[3]; ccf_rv2_err = ccfrvdata[4]

## Recall: The Barycentric Julian Date (BJD) is the Julian Date (JD) corrected
## for differences in the Earth's position with respect to the barycentre of
## the Solar System
dateoffset = 2454833.                     # BJD - dateoffset = JD
                                          # for Unfolded RV vs time

## Let's  read in fit estimates (RV estimates) from BF_python.py [BF RVs]
bfinfile = '6864859/6864859Outfile.txt'
bfrvdata = np.loadtxt(bfinfile, comments='#', unpack=True)

bjd = bfrvdata[0]
phase = bfrvdata[1]
bf_rv1 = bfrvdata[3]; bf_rv1_err = bfrvdata[4]
bf_rv2 = bfrvdata[5]; bf_rv2_err = bfrvdata[6]

#print(bf_rv1, bf_rv1_err, bf_rv2, bf_rv2_err)

## Recall: Most of the figures with RV estimates in Cunningham et al (2019)
## are unfolded phase, so we will begin with that, but first: the calculations
## i.e. doubling the arrays to allow this folded phase functionality

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
#    plt.errorbar(date-dateoffset, bf_rv1[idx]*bf_stretch, yerr=bf_rv1_err[idx], fmt='ko', color="darkorange", mfc="darkorange", mec="darkorange", ms=10, lw=1.5)
#    plt.errorbar(date-dateoffset, bf_rv2[idx]*bf_stretch, yerr=bf_rv2_err[idx], fmt='ko', color="darkorange", mfc="white", mec="darkorange", ms=10, lw=1.5)
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

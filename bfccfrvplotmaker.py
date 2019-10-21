import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import RegularPolyCollection
from BF_functions import user_rc
from astropy.io import fits
from astropy.time import Time
'''
Radial velocity plotter!
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase

You need to have these columns in your input file: TIME, PHASE, RV1, RV1_ERR, RV2, RV2_ERR

Update September 2015:
Has two flag options
1. apply some shift to the RVs before plotting them (doShift)
2. read in another set of calculated RVs from cols 8,9,10,11 and plot RV1 = RV_col8-RV_col3
   and RV2 = RV_col10-RV_col5, plus a line at RV = 0 (compareRVs)

Update June 2016:
Simplified some options; no longer manually sets point shape as a function of "source" string.
If you want that functionality, use an older version of this code... it was messy.
**NOTE that any RV value with an error bar = 0 is not plotted!**
'''
###CCF

KIC = '6131659'; windowrows = 4; windowcols = 7
dir = os.path.join('data', KIC)
ccffile = 'apStar-r8-2M19370697+4126128.fits' #6131659
bffile = '6131659BFOutALL.txt'
BCVfile = '6131659bjdinfile.txt'
ccfinfile = os.path.join('data', ccffile)
bfinfile = os.path.join(dir, bffile)
bjdinfile = os.path.join(dir, BCVfile)

# Read in relevant CCF info from apStar file
hdu = fits.open(ccfinfile)
# The CCF information lives in the 9th HDU, because of course it does
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0][2:]
CCFerrors = hdu9['CCFERR'][0][2:]
CCFxaxis = hdu9['CCFLAG'][0]  # pixels, needs to be turned into RVs
CCF_delta = hdu9['CCFDW'][0]  # Delta (log_10 lambda)
#print(hdu9['VHELIO'][0] - hdu9['VREL'][0])
# (the data model website says this log lambda spacing per lag step of
# 6.d-6 corresponds to 4.145 km/s)
# Get the systemic RV of the star according to APOGEE
apSystemic = 82  ## TODO: FIND ACTUAL SYSTEMIC RV PARAMETER !!!
# Get the barycentric velocities of each visit according to APOGEE
# This is called VHELIO for super obvious reasons, Jen Sobeck private communication, for reals
apBCVs = hdu9['VHELIO'][0]
CCF_rvaxis = [CCFxaxis * 4.145 + bcv for bcv in apBCVs]

# Get the timestamp for each CCF visit from the 0th HDU
hdu0 = hdu[0].header
ccftimes = []
for idx in range(1, len(CCFvalues)+1):
    headercard = 'HJD' + str(idx)
    ccftimes.append(hdu0[headercard] + 2400000.)

ccftimesAstropy = []
for ccftime in ccftimes:
    ccftimesAstropy.append(Time(ccftime, scale='utc', format='jd'))

dateoffset = 2454833. # this value will be subtracted from bjds in pane vs. time

#sysname = '5285607'; filename = '5285607Outfile_take2.txt'
#timestart = 980; timeend = 1020
#phasemin = 0.5; phasemax = 1.5
#RVmin = -45; RVmax = 180

#sysname = '6449358'; filename = 'data/6449358/6449358Outfile.txt'
#timestart = 1700; timeend = 2000
#phasemin = 0.5; phasemax = 1.5
#RVmin = 0; RVmax = 140

colors = user_rc()

sysname = '6131659'; filename = 'data/6131659/6131659outfileALL.txt'
timestart = 1520; timeend = 2000
phasemin = -0.1; phasemax = 1.0
RVmin = 0; RVmax = 200

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2


###BF RV DATA###
# usecols=(0,1,3,4,5,6) # this is the default, with RVs in 3,4,5,6 not 8,9,10,11
rvdata = np.loadtxt(filename, comments='#', unpack=True)

bjd = rvdata[0]
phase = rvdata[1]
rv1 = rvdata[3]; rverr1 = rvdata[4]
rv2 = rvdata[5]; rverr2 = rvdata[6]
rv3 = rvdata[7]; rverr3 = rvdata[8]

# Skip any RV values that have 0 for error bars or are already None
for idx, err in enumerate(rverr1):
    if err == 0:
        rv1[idx] = None
        rverr1[idx] = None
for idx, err in enumerate(rverr2):
    if err == 0:
        rv2[idx] = None
        rverr2[idx] = None
rv1mask = np.isfinite(rv1)
rv2mask = np.isfinite(rv2)
rv3mask = np.isfinite(rv3)

# Double the BF arrays so we can plot any phase from 0 to phase 2... assuming phase is in range (0,1)
rv1_double = np.concatenate((rv1,rv1), axis=0)
rv2_double = np.concatenate((rv2,rv2), axis=0)
rv3_double = np.concatenate((rv3,rv3), axis=0)
phase_double = np.concatenate((np.array(phase),np.array(phase)+1.0), axis=0)
rverr1_double = np.concatenate((rverr1,rverr1), axis=0)
rverr2_double = np.concatenate((rverr2,rverr2), axis=0)
rverr3_double = np.concatenate((rverr3, rverr3), axis=0)

#Double the CCF arrays so we can plot any phase from 0 to phase 2...
'''
rv1_double = np.concatenate((rv1,rv1), axis=0)
rv2_double = np.concatenate((rv2,rv2), axis=0)
rv3_double = np.concatenate((rv3,rv3), axis=0)
phase_double = np.concatenate((np.array(phase),np.array(phase)+1.0), axis=0)
rverr1_double = np.concatenate((rverr1,rverr1), axis=0)
rverr2_double = np.concatenate((rverr2,rverr2), axis=0)
rverr3_double = np.concatenate((rverr3, rverr3), axis=0)
'''
# Set up the figure
fig = plt.figure(1, figsize=(13,9))
'''
# Unfolded RV vs time (BJD-2454833)
ax2 = plt.subplot(2,1,1)
plt.axis([timestart, timeend, RVmin, RVmax])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
# dotted lines to guide the eye
plt.plot(bjd[rv1mask]-dateoffset, rv1[rv1mask], color=colors[15], mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2mask]-dateoffset, rv2[rv2mask], color=colors[15], mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv3mask]-dateoffset, rv3[rv3mask], color=colors[15], mfc=None, mec=None, lw=1.5, ls=':')
for idx, date in enumerate(bjd):
    plt.errorbar(date-dateoffset, rv1[idx], yerr=rverr1[idx], fmt='ko', color=colors[15], mfc=colors[6], mec=colors[14], ms=10, lw=1.5)
    plt.errorbar(date-dateoffset, rv2[idx], yerr=rverr2[idx], fmt='ko', color=colors[15], mfc=colors[2], mec=colors[14], ms=10, lw=1.5)
    plt.errorbar(date-dateoffset, rv3[idx], yerr=rverr3[idx], fmt='ko', color=colors[15], mfc=colors[8], mec=colors[14], ms=10, lw=1.5)
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset))
'''
axlims = [-140, 140, -0.06, 0.56]
ccfyoffset = -0.12
bfyamp = 9 # arbitrary factor to stretch BF values for clarity

#BF_vs_CCF
fig = plt.figure(1, figsize=(15,10))
ax = fig.add_subplot(111)
ax.tick_params(top='off', bottom='off', left='off', right='off')
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xlabel('Radial Velocity (km s$^{-1}$)', labelpad=20, size='x-large')
ax.set_ylabel('Arbitrary CCF or BF amplitude', labelpad=20, size='x-large')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
plt.title('KIC ' + KIC, size='x-large')

# Loop over and plot CCF data
for idx, CCFdata in enumerate(CCFvalues):
    ax0 = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    plt.plot(CCF_rvaxis[idx], CCFvalues[idx] + ccfyoffset)
    plt.text(23, 0.5, ccftimesAstropy[idx].iso[0:10], size=10)
    ## TURN ON FOR TIMESTAMP TROUBLESHOOTING
    #plt.text(25, 0.5, '{0:.3f}'.format(ccftimes[idx] - 2450000.), size=10, color='C0')
    #plt.text(25, 0.3, idx, size=10, color='C0')
    plt.subplots_adjust(wspace=0, hspace=0)

# Read in relevant BF info from the BF infile
bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
BFrvaxis = bfdata[0]
BFvalues = bfdata[1]
BFerrors = bfdata[2]

# Get the timestamp for each BF
with open(bfinfile) as bfinfo:
    bftimes = [float(line[13:]) for line in bfinfo if 'timestamp' in line]

# Save the indices that separate each BF's visit in the input file
visitidx = 0
BFindices = []
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
        visitidx = visitidx + 1
        BFindices.append(idx)

# Read in barycentric velocity correction info from the BJD infile
BCVdata = np.loadtxt(bjdinfile, comments='#', usecols=(2,), unpack=True)
BCVdata = BCVdata[1:]

# Loop over and plot BF data
for idx in range(0, len(bftimes)):
    ax1 = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if (idx < 20):  # set to 18 for 686 and 20 for 613
        ax1.set_xticklabels(())
    if (idx!=0 and idx!=7 and idx!=14 and idx!=21):
        ax1.set_yticklabels(())
    ## TURN ON FOR TIMESTAMP TROUBLESHOOTING
    #plt.text(25, 0.4, '{0:.3f}'.format(bftimes[idx] - 2450000.), size=10, color='C1')
    #plt.text(25, 0.2, idx, size=10, color='C1')
    try:
        plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] - apSystemic + BCVdata[idx],
                 bfyamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
    except: # handle the final case where there is no idx+1
        try:
            plt.plot(BFrvaxis[BFindices[idx]::] - apSystemic + BCVdata[idx],
                    bfyamp*BFvalues[BFindices[idx]::])
        except:
            print('You\'re missing a BF where an APOGEE CCF exists')
            continue

plt.show()

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
for idx, ph in enumerate(phase_double):
    plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=colors[6], mec=colors[14], ecolor=colors[6], ms=10, ls='None', lw=1.5)
    plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=colors[2], mec=colors[14], ecolor=colors[2], ms=10, ls='None', lw=1.5)
    plt.errorbar(phase_double[idx], rv3_double[idx], yerr=rverr3_double[idx], marker='o', color=colors[8], mec=colors[14], ecolor=colors[8], ms=10, ls='None', lw=1.5)
plt.xlabel("Orbital Phase")

# Draw vertical lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Option for a legend and labels (note: for a legend you will need to add a label to the plt.errorbar commands)
#plt.legend(ncol=2, loc=1, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35), columnspacing=0.7)
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', rotation='vertical', size='large')
fig.text(0.14, 0.13, 'Folded')
fig.text(0.14, 0.55, 'Unfolded')
fig.text(0.14, 0.9, sysname, size='large')

plt.show()

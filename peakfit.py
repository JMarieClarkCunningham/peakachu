from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl
from scipy import ndimage
import pandas as pd
import gaussfitter as gf


##############################      GAUSSPARTY     #############################
def gaussparty(gausspars, nspec, CCFvaluelist, CCFxaxis, error_array, amplimits, threshold, widlimits):
    '''
    Fits 2 gaussians to some data
    '''
    param = []
    with open(gausspars) as f1:
        for line in f1:
            if line[0] != '#':
                param.append( line.rstrip() )
    #param = np.loadtxt(gausspars, comments='#')
    bffitlist = []
    gauss1 = [[] for i in range(nspec)]
    gauss2 = [[] for i in range(nspec)]
    #error_array = np.ones(len(CCFvaluelist[0]))*0.01 # dummy array with 0.01 error values
    print(' ')
    print('Gaussian fit results: peak amplitude, width, rvraw, rvraw_err')
    print ('-------------------------------------------------------------')
    for i in range(0, nspec):
        if '#' in param[i]:
            commentbegin = param[i].find('#')
            partest = param[i][0:commentbegin].split()
        else:
            partest = param[i].split()
        partest = [float(item) for item in partest]
        ngauss = 2
        #minpars = [amplimits[0], float(partest[1])-threshold, widlimits[0]]
        #maxpars = [amplimits[1], float(partest[1])+threshold, widlimits[1]]
        #minpars.extend([amplimits[2], float(partest[4])-threshold, widlimits[2]])
        #maxpars.extend([amplimits[3], float(partest[4])+threshold, widlimits[3]])
        bffit = gf.multigaussfit(CCFxaxis[i], CCFvaluelist[i], ngauss=ngauss,
                params=partest, err=error_array[i], quiet=True, shh=True, veryverbose=False)
        bffitlist.append(bffit)
        print('{0:.3f} {1:.2f} {2:.4f} {3:.4f} \t {4:.3f} {5:.2f} {6:.4f} {7:.4f}'.format(
            bffit[0][0], bffit[0][2], bffit[0][1], bffit[2][1],
            bffit[0][3], bffit[0][5], bffit[0][4], bffit[2][4]))
    return bffitlist

'''
This is a truncated version of the BF_python.py software by Meredith Rawls
to fit Gaussians via the method of least squares to BF and CCF peaks.
This results of this iterative fitting are then written to an out
file for use in the RV plotting software bfccfrvplotmaker.py

INPUT
gausspars:  your best initial guesses for fitting gaussians to the CC(B)F peaks
            the parameters are [amp1, offset1, width1, amp2, offset2, width2]
            the top line is ignored (template), but must have six values
            one line per observation
            IMPORTANT: If the system being fitted is a suspected triple system,
            gausspars had better have 9 values per observation, i.e. include
            amp3, offset3, width3.

OUTPUT
ccfpeakout: a output .txt file created with # columns: [CCFRV1, CCFRV1_err,
            CCFRV2, CCFRV2_err, CCFRV3, CCFRV3_err]
'''

#################      READ IN DATA FOR TARGET OF INTEREST     #################

#################                   6131659                    #################
KIC = '6131659'; windowrows = 4; windowcols = 7        # KIC num; rows; columns
nspec = 27                                             # KIC 6131659 has 27 visits
gausspars = '6131659/6131659gaussparsCCF.txt'          # gausspars
outfile = '6131659/6131659ccfpeakout.txt'               # outfile
ccfinfile = 'ApStar/apStar-r8-2M19370697+4126128.fits' # CCF infile
#bffile = '6131659BFOutALL.txt'                         # BF infile
#BCVfile = '6131659bjdinfile.txt'                       # Barycentric correction
#infiles =   'data/6131659/6131659infiles.txt'          # BF infile
#bjdinfile = 'data/6131659/6131659bjdinfile.txt'        # bjd infile

#################                    6864859                  ##################
#KIC = '6864859'; windowrows = 4; windowcols = 7         # KIC num; rows; columns
#nspec = 24                                              # KIC 6864859 has 25 visits
#gausspars = '6864859/6864859gaussparsCCF.txt'           # gausspars
#outfile = '6864859/6864859ccfpeakout.txt'               # outfile
#ccfinfile = 'ApStar/apStar-r6-2M19292405+4223363.fits'  # 6864859  #r8 is missing one visit
#bffile = '6131659BFOutALL.txt'                          # BF infile
#BCVfile = '6131659bjdinfile.txt'                        # Barycentric correction
#infiles =   'data/6131659/6131659infiles.txt'           # BF infile
#bjdinfile = 'data/6131659/6131659bjdinfile.txt'         # bjd infile

#################    FOR CCF PEAKS FROM ApStar  [CCF INPUT]    #################
hdu = fits.open(ccfinfile)
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0][2:]                  # CCF peak information
CCFerrors = hdu9['CCFERR'][0][2:]               # CCF peak errors
CCFxaxis = hdu9['CCFLAG'][0]                    # pixels, needs to be turned into RVs
CCF_delta = hdu9['CCFDW'][0]                    # Delta (log_10 lambda)
apBCVs = hdu9['VHELIO'][0]                      # Barycentric velocities of each visit
CCF_rvaxis = [CCFxaxis * 4.145 + bcv for bcv in apBCVs]
#print(CCFvalues[0])
#print(len(CCFvalues), len(CCF_rvaxis))
                                                # BARYCENTRIC CORRECTION TO CCF
hdu0 = hdu[0].header                            # 56 - 60 get the timestamp for
ccftimes = []                                   # each visit which lines 62 - 66
for idx in range(1, len(CCFvalues)+1):          # calculate the dateoffset, which
    headercard = 'HJD' + str(idx)               # will be subtracted from bjds
    ccftimes.append(hdu0[headercard] + 2400000.)

ccftimesAstropy = []
for ccftime in ccftimes:
    ccftimesAstropy.append(Time(ccftime, scale='utc', format='jd'))

dateoffset = 2454833.

#################     FOR BF PEAKS FROM ApVisit  [BF INPUT]    #################
# This same interative Gaussian peak-fitting process was completed on ApVisit  #
# spectra for seven SEBs using the BF_python software suite by Meredith Rawls  #
# see Cunningham et al. (2019). The results of this analysis are only read in  #
# here, but functionality exists for peakfit.py to fit Gaussians via MLS to BF #
# and CCF peaks. This allows for comparison between BF and CCF RV extractions. #

#bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
#BFrvaxis = bfdata[0]
#BFvalues = bfdata[1]
#BFerrors = bfdata[2]
'''
with open(bfinfile) as bfinfo:                  # retrieve timstamps from BF infile
    bftimes = [float(line[13:]) for line in bfinfo if 'timestamp' in line]

visitidx = 0                                    # save the indices that separate
BFindices = []                                  # each BF's visit in the input file
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
        visitidx = visitidx + 1
        BFindices.append(idx)

BCVdata = np.loadtxt(bjdinfile, comments='#', usecols=(2,), unpack=True)
BCVdata = BCVdata[1:]                           # read in barycentric velocity
                                                # correction info from BJD infile
'''
#################             IMPORTANT DEFINITIONS            #################
amplimits = [0,1.2, 0,1.2, 0,1.2]               # limits, gaussian normalized amp
                                                # [min1,max1,min2,max2]
threshold = 50                                  # margin for gaussian position
                                                # [raw RV in km/s]
#widlimits = [0,9, 0,9, 0,11]; rvneg = 0; rvpos = 199; ymin = -0.15; ymax = 1.18 # 6131659
widlimits = [20,40, 20,40, 20,40]               # limits for Gaussian width (km/s)
                                                # [min1, max1, min2, max2]
rvneg = -300; rvpos = 300

#rvstd = 0; bcvstd = 0                           # RV and BCV model template (km/s)
                                                 # [0, 0] if using a model
#smoothstd = 1.5                                 # stdev of Gaussian to smooth peaks by
                                                # (~slit width in pixels)
#w00 = 15170; n = 6000; stepV = 4.0              # (lower res, apStar)
#timestart = 1520; timeend = 2000                # observation times (APOGEE)
#phasemin = -0.1; phasemax = 1.0                 # phase limits
#RVmin = 0; RVmax = 200                          # RV limits
#colors = bff.user_rc()                          # color consistancy in project
#bf_ind = svd.getRVAxis(r, 1) + rvstd - bcvstd   # obtains indices in RV space that
#                                                # correspond to the BF
#svd = pyasl.SVD()                               # single value decomposition
#svd.decompose(newspeclist[0], m)
#singularvals = svd.getSingularValues()          # don't believe I need to decompose
#CCF_yaxis = CCFvalues[idx] + ccfyoffset

#################             MAIN PLOT FORMATTING             #################

fig = plt.figure(1, figsize=(12,6))
ax = fig.add_subplot(111)
ax.tick_params(top='off', bottom='off', left='off', right='off')
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xlabel('Radial Velocity (km s$^{-1}$)', labelpad=20, size='x-large')
ax.set_ylabel('Arbitrary CCF amplitude', labelpad=20, size='x-large')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')

#################    PRELIMINARY VIEW OF PEAKS TO BE FITTED    #################
plt.axis([rvneg, rvpos, -0.2, float(nspec)+1])
axlims = [-140, 140, -0.06, 0.56]
ccfyoffset = 0                              # offset of CCF amplitude
#bfyamp = 9                                      # arbitrary factor to stretch BF
yoffset = 0                                     # values for clarity

# LOOP OVER AND PLOT CCFs FROM ApStar
#for idx, CCFdata in enumerate(CCFvalues):
#    ax0 = fig.add_subplot(windowrows, windowcols, idx+1)
#    plt.axis(axlims)
#    plt.plot(CCF_rvaxis[idx], CCFvalues[idx] + ccfyoffset)
#    plt.text(23, 0.5, ccftimesAstropy[idx].iso[0:10], size=10)
#    #plt.axhline(y=yoffset, color=colors[15], ls=':')
#    plt.subplots_adjust(wspace=0, hspace=0)
#    yoffset = yoffset + 1.0
#plt.show()

#################     FIT PEAKS WITH TWO OR MORE GAUSSIANS     #################
# IMPORTANT: You must have fairly decent guesses in the gausspars file for the #
# following code to work. Otherwise, this will fail and resultant outputs will #
# look terrible. Remember: if the system being fitted is a suspected trinary,  #
# gausspars had better have 9 values per observation.                          #
#peakfitlist = bff.gaussparty(gausspars, nspec, filenamelist, CCF_rvaxis, CCFvalues, amplimits, threshold, widlimits)
#peakfitlist = bff.gaussparty(nspec, gausspars, filenamelist, CCFvalues, CCF_rvaxis, amplimits, threshold, widlimits)
peakfitlist = gaussparty(gausspars, nspec, CCFvalues, CCF_rvaxis, CCFerrors, amplimits, threshold, widlimits)
#print(peakfitlist) #prints the peak fitting results
rvraw1 = []; rvraw2 = []; rvraw1_err = []; rvraw2_err = []; rvraw3 = []; rvraw3_err = []
amp1 = []; amp2 = []; width1 = []; width2 = []
for peakfit in peakfitlist:
    rvraw1.append(peakfit[0][1]) # indices are [parameter, BF, error array][amp,rv,width x N]
    rvraw2.append(peakfit[0][4]) # [0,1,2] is amp,rv,width for star1; [3,4,5] is same for star2, etc.
    rvraw1_err.append(peakfit[2][1])
    rvraw2_err.append(peakfit[2][4])
    amp1.append(peakfit[0][0])
    amp2.append(peakfit[0][3])
    width1.append(peakfit[0][2])
    width2.append(peakfit[0][5])

#print(rvraw1, amp1, width1)

# handy little gaussian function maker
def gaussian(x, amp, mu, sig): # i.e., (xarray, amp, rv, width)
    return amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

xmin = -120
xmax = 120
ymin = -0.2
ymax = 0.8

for i in range(0, nspec):
    ax = fig.add_subplot(windowrows, windowcols, i+1) # out of range if windowcols x windowrows < nspec
    ax.yaxis.set_major_locator(MultipleLocator(0.4)) #increments of y axis tic marks
    if (i!=0 and i!=7 and i!=14 and i!=21):
        ax.set_yticklabels(())
    if i < nspec-windowcols:
        ax.set_xticklabels(())
    plt.subplots_adjust(wspace=0, hspace=0.0, bottom=0.2) #6131659
#    plt.subplots_adjust(wspace=0, hspace=0.0, bottom=0.2) #6449358
#    plt.plot_adjust(wspace=0, hspace=0)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(axis='both', which='major')
    #plt.text(xmax - 0.19*(np.abs(xmax-xmin)), 0.60*ymax, i)
    #plt.text(xmax - 0.19*(np.abs(xmax-xmin)), 0.60*ymax, '%.3f $\phi$' % (phase[i]), size='small')
    #plt.text(xmax - 0.26*(np.abs(xmax-xmin)), 0.35*ymax, '%s' % (datetimelist[i].iso[0:10]), size='small')
    #plt.plot(bf_ind, bfsmoothlist[i], color=colors[14], lw=1.5, ls='-', label='Smoothed BF')

    plt.plot(CCF_rvaxis[i], CCFvalues[i], color='0.5', lw=2, ls='-', label='ApStar CCFs')
    plt.plot(CCF_rvaxis[i], peakfitlist[i][1], color='C0', lw=2, ls='-', label='Two Gaussian fit')
    gauss1 = gaussian(CCF_rvaxis[i], amp1[i], rvraw1[i], width1[i])
    gauss2 = gaussian(CCF_rvaxis[i], amp2[i], rvraw2[i], width2[i])
    plt.plot(rvraw1[i], 0.1, color='C3', marker='|', ms=15)#, label='RV 1')
    plt.plot(rvraw2[i], 0.1, color='C1', marker='|', ms=15)#, label='RV 2')
    plt.plot(CCF_rvaxis[i], gauss1, color='C3', lw=3, ls='--')#, label='Gaussian fit 1')
    plt.plot(CCF_rvaxis[i], gauss2, color='C1', lw=3, ls='--')#, label='Gaussian fit 2')
    # OPTION TO PLOT VERTICAL LINE AT ZERO
    #plt.axvline(x=0, color=colors[15])

    # in this situation, the legend is printed to the right of the final subplot
    if i == nspec-1:
        ax.legend(bbox_to_anchor=(2.6,0.6), loc=1, borderaxespad=0.,
                  frameon=False, handlelength=2, prop={'size':12})

plt.show()
plt.savefig('6131659prelim.jpg')

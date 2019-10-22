import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from scipy import optimize
import gaussfitter as gf
from astropy.modeling import models, fitting
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
            CCFRV2, CCFRV2_err]
'''


#def gaussbetter(model, CCFxaxis, CCFvalues, CCFerrors):
def gaussbetter(pars, CCFvalues, CCFxaxis, CCFerrors):
    '''
    Fit 2 gaussians to some data, astropy style.

    Unlike gaussparty, this fits one visit at a time. I should finish this docstring.
    '''
    g1 = models.Gaussian1D(pars[0], pars[1], pars[2], fixed={'stddev': True})
    g2 = models.Gaussian1D(pars[3], pars[4], pars[5], fixed={'stddev': True})
    # the parameter order is amp, mean (aka position aka RV), stddev
    # as written, the widths are not fit, and are held fixed at the input values
    gg_init = g1 + g2
    fitter = fitting.LevMarLSQFitter()#http://docs.astropy.org/en/stable/api/astropy.modeling.fitting.LevMarLSQFitter.html
    gg_fit = fitter(gg_init, CCFxaxis, CCFvalues, weights=1./np.array(CCFerrors))     #This fits the combined Gaussians (g1 + g2)
    cov = fitter.fit_info['param_cov']#Should the cov calculation be here, instead??? [JMC^2]
    return cov
    return gg_fit


def gaussparty(gausspars, nspec, CCFvaluelist, CCFxaxis, error_array, ngauss=2,
               fixed=[False, False, False],
               limitedmin=[False, False, True],
               limitedmax=[False, False, False],
               minpars=[0, 0, 0], maxpars=[0, 0, 0]):
    '''
    Fit 2 gaussians to some data.

    Parameters
    ----------
    gausspars : `str`
        Filename containing initial Gaussian parameter guesses
    nspec : `int`
        Number of spectra (visits)
    CCFvaluelist : `list`
        2D list with CCF values, one per visit
    CCFxaxis : `list`
        2D list with CCF RV x-axis values corresponding to CCFvaluelist, one per visit
    error_array: `list`
        2D list with errors corresponding to CCFvaluelist
    fixed : `list` of bools
        Is parameter fixed? (amplitude, position, width)
    limitedmin/minpars : `list`
        Set lower limits on each parameter (default: width > 0)
    limitedmax/maxpars : `list`
        Set upper limits on each parameter
    minpars : `list` (amplitude, position, width)
    maxpars : `list` (amplitude, position, width)

    Returns
    -------
    gaussfitlist : `list`
        2D list of results from calling multigaussfit
    '''
    assert len(CCFvaluelist) == nspec
    assert len(CCFxaxis) == nspec
    assert len(error_array) == nspec

    param = []
    with open(gausspars) as f1:
        for line in f1:
            if line[0] != '#':
                param.append(line.rstrip())

    assert len(param) == nspec

    gaussfitlist = []
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
        gaussfit = gf.multigaussfit(CCFxaxis[i], CCFvaluelist[i], ngauss=ngauss,
                                    params=partest, err=error_array[i], fixed=fixed,
                                    limitedmin=limitedmin, limitedmax=limitedmax,
                                    minpars=minpars, maxpars=maxpars,
                                    quiet=True, shh=True, veryverbose=False)
        gaussfitlist.append(gaussfit)
        print('{0:.3f} {1:.2f} {2:.4f} {3:.4f} \t {4:.3f} {5:.2f} {6:.4f} {7:.4f}'.format(
              gaussfit[0][0], gaussfit[0][2], gaussfit[0][1], gaussfit[2][1],
              gaussfit[0][3], gaussfit[0][5], gaussfit[0][4], gaussfit[2][4]))
    return gaussfitlist


def gaussian(x, amp, mu, sig):  # i.e., (xarray, amp, rv, width)
    '''
    Handy little gaussian function maker
    '''
    return amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# THIS IS WHERE YOU CHANGE THINGS

# KIC = '6131659'
KIC = '6864859'
fitter = 'astropy'
# fitter = 'gaussfitter'
gausspars = os.path.join(KIC, KIC+'gaussparsCCF.txt')
outfile = os.path.join(KIC, KIC+'ccfpeakout.txt')

if KIC == '6131659':
    windowrows = 4
    windowcols = 7
    nspec = 27
    ccfinfile = os.path.join(KIC, 'apStar-r8-2M19370697+4126128.fits')
    # BCVfile = '6131659bjdinfile.txt'
elif KIC == '6864859':
    windowrows = 4
    windowcols = 7
    nspec = 25  # KIC 6864859 has 25 visits...
    ccfinfile = os.path.join(KIC, 'apStar-r6-2M19292405+4223363.fits')
    # BCVfile = '6131659bjdinfile.txt'
else:
    raise NotImplementedError()

# Load the CCF data from the apStar file
hdu = fits.open(ccfinfile)
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0][2:]
CCFerrors = hdu9['CCFERR'][0][2:]
CCFxaxis = hdu9['CCFLAG'][0]  # pixels, needs to be turned into RVs
apBCVs = hdu9['VHELIO'][0]  # Barycentric velocities
CCF_rvaxis = [CCFxaxis * 4.145 + bcv for bcv in apBCVs]

# Get the timestamps for each observation
hdu0 = hdu[0].header
ccftimes = []
for idx in range(1, len(CCFvalues)+1):
    headercard = 'HJD' + str(idx)
    ccftimes.append(hdu0[headercard] + 2400000.)
ccftimesAstropy = []
for ccftime in ccftimes:
    ccftimesAstropy.append(Time(ccftime, scale='utc', format='jd'))

# Optional initial plot of CCFs
doInitialPlot = False

if doInitialPlot:
    rvneg = -300
    rvpos = 300
    fig = plt.figure(1, figsize=(12, 6))
    ax = fig.add_subplot(111)
    plt.axis([rvneg, rvpos, -0.2, float(nspec) + 1])
    axlims = [-140, 140, -0.06, 0.56]
    ccfyoffset = 0
    yoffset = 0
    for idx, CCFdata in enumerate(CCFvalues):
        ax0 = fig.add_subplot(windowrows, windowcols, idx + 1)
        plt.axis(axlims)
        plt.plot(CCF_rvaxis[idx], CCFvalues[idx] + ccfyoffset)
        plt.text(23, 0.5, ccftimesAstropy[idx].iso[0:10], size=10)
        # plt.axhline(y=yoffset, color=colors[15], ls=':')
        plt.subplots_adjust(wspace=0, hspace=0)
        yoffset = yoffset + 1.0
    plt.show()

# IMPORTANT: You must have fairly decent guesses in the gausspars file for the
# following code to work. Otherwise, this will fail and resultant outputs will
# look terrible.

# Time to fit the peaks

if fitter == 'gaussfitter':
    peakfitlist = gaussparty(gausspars, nspec, CCFvalues, CCF_rvaxis, CCFerrors)
    rvraw1 = []
    rvraw2 = []
    rvraw1_err = []
    rvraw2_err = []
    amp1 = []
    amp2 = []
    width1 = []
    width2 = []
    for peakfit in peakfitlist:
        rvraw1.append(peakfit[0][1])  # indices are [parameter, BF, error array][amp,rv,width x N]
        rvraw2.append(peakfit[0][4])  # [0,1,2] is amp,rv,width for star1; [3,4,5] is same for star2, etc.
        rvraw1_err.append(peakfit[2][1])
        rvraw2_err.append(peakfit[2][4])
        amp1.append(peakfit[0][0])
        amp2.append(peakfit[0][3])
        width1.append(peakfit[0][2])
        width2.append(peakfit[0][5])
    bestFitModelList = []
    for i in range(0, nspec):
        bestFitModelList.append(peakfitlist[i][1])

elif fitter == 'astropy':
    bestFitModelList = []
    param = []
    rvraw1 = []
    rvraw2 = []
    rvraw1_err = []
    rvraw2_err = []
    amp1 = []
    amp2 = []
    width1 = []
    width2 = []
    with open(gausspars) as f1:
        for line in f1:
            if line[0] != '#':              #Skips commented out lines
                param.append(line.rstrip())
    assert len(param) == nspec
    partestlist = []
    for par in param:                       #Skips commented out lines
        if '#' in par:
            commentbegin = par.find('#')
            partest = par[0:commentbegin].split()
        else:
            partest = par.split()
        partest = [float(item) for item in partest]
        partestlist.append(partest)
    for idx in range(idx, nspec): #used to be: for idx in range(0, nspec)
        pars = partestlist[idx]
        result = gaussbetter(pars, CCFvalues[idx],
                             CCF_rvaxis[idx], CCFerrors[idx])
        bestFitModel = result(CCF_rvaxis[idx])
        bestFitModelList.append(bestFitModel)
        rvraw1.append(result.mean_0.value)
        rvraw2.append(result.mean_1.value)
        #cov = LevMarLSQFitter.fit_info['param_cov']
        parnames = [n for n in result.param_names if n not in ['stddev_0', 'stddev_1']]
        parvals = [v for (n, v) in zip(result.param_names, result.parameters) if n not in ['stddev_0', 'stddev_1']]
        for i, (name, value) in enumerate(zip(parnames, parvals)):
            print('{}: {} +/- {}'.format(name, value, np.sqrt(cov[i][i])))

        # rvraw1_err.append(result.SOMETHING???)  # TODO: get errors on fit parameters
        # rvraw2_err.append(result.SOMETHING???)
        rvraw1_err.append(0)  # DELETE THIS LATER
        rvraw2_err.append(0)  # DELETE THIS LATER
        amp1.append(result.amplitude_0.value)
        amp2.append(result.amplitude_1.value)
        width1.append(result.stddev_0.value)
        width2.append(result.stddev_1.value)

else:
    raise NotImplementedError()

# Print fit results to the outfile
with open(outfile, 'w') as f2:
    print('# time [JD], RV1 [km/s], error1 [km/s], RV2 [km/s], error2 [km/s]', file=f2)
    print('#', file=f2)
    for i in range(1, nspec): #used to be: for i in range(1, nspec)
        print ('{0:.9f} {1:.5f} {2:.5f} {3:.5f} {4:.5f}'.format(ccftimesAstropy[i].jd,
                                                                rvraw1[i],
                                                                rvraw1_err[i],
                                                                rvraw2[i],
                                                                rvraw2_err[i]), file=f2)
print('Time and RVs written to %s.' % outfile)

# Plotting time
fig = plt.figure(1, figsize=(12, 6))
ax = fig.add_subplot(111)
ax.tick_params(top=False, bottom=False, left=False, right=False)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xlabel('Radial Velocity (km s$^{-1}$)', labelpad=20, size='x-large')
ax.set_ylabel('Arbitrary CCF amplitude', labelpad=20, size='x-large')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')

xmin = -120
xmax = 120
ymin = -0.2
ymax = 0.8

for i in range(0, nspec):
    ax = fig.add_subplot(windowrows, windowcols, i + 1)  # out of range if windowcols x windowrows < nspec
    ax.yaxis.set_major_locator(MultipleLocator(0.4))  # increments of tick marks
    ax.xaxis.set_major_locator(MultipleLocator(80))
    if (i != 0 and i != 7 and i != 14 and i != 21):  # this is a function of windowcols!
        ax.set_yticklabels(())
    if i < nspec-windowcols:
        ax.set_xticklabels(())
    plt.subplots_adjust(wspace=0, hspace=0.0, bottom=0.1)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(axis='both', which='major')
    plt.text(xmax - 0.19*(np.abs(xmax-xmin)), 0.60*ymax, i)
    plt.text(xmax - 0.6*(np.abs(xmax-xmin)), 0.8*ymax, '%s' % (ccftimesAstropy[i].iso[0:10]), size='small')
    plt.plot(CCF_rvaxis[i], CCFvalues[i], color='0.5', lw=2, ls='-', label='ApStar CCFs')
    plt.plot(CCF_rvaxis[i], bestFitModelList[i], color='C0', lw=2, ls='-', label='Two Gaussian fit')
    gauss1 = gaussian(CCF_rvaxis[i], amp1[i], rvraw1[i], width1[i])
    gauss2 = gaussian(CCF_rvaxis[i], amp2[i], rvraw2[i], width2[i])
    plt.plot(CCF_rvaxis[i], gauss1, color='C3', lw=3, ls='--')
    plt.plot(CCF_rvaxis[i], gauss2, color='C1', lw=3, ls='--')
    plt.plot(rvraw1[i], 0.1, color='C3', marker='|', ms=15)
    plt.plot(rvraw2[i], 0.1, color='C1', marker='|', ms=15)
    # plt.axvline(x=0, color=colors[15])  # optional vertical line at 0

    # in this situation, the legend is printed to the right of the final subplot
    if i == nspec-1:
        ax.legend(bbox_to_anchor=(2.6, 0.6), loc=1, borderaxespad=0.,
                  frameon=False, handlelength=2, prop={'size': 12})

plt.show()
plt.savefig(KIC + 'prelim.jpg')

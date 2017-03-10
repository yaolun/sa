# This is a python 3 script
def AlphaClassification(obj, phot_dir, spec=None, usephot=True, waverange=None, plotdir=None):
    import numpy as np
    from astropy.io import ascii
    import astropy as apy
    import os
    from Spectrophotometry import Spectrophotometry
    import matplotlib.pyplot as plt
    import astropy.constants as const
    c = const.c.cgs.value

    if waverange == None:
        waverange = [2.0, 25.0]

    # read in the photometry ASCII file
    if usephot:
        if not os.path.exists(phot_dir+obj.lower()+'.txt'):
            print('No photometry data found')
            if (spec is None) & (phot is None):
                return None
            # else:
            #     phot = apy.table.Table(names=['wavelength', 'flux(Jy)', 'error(Jy)'])
        else:
            phot = ascii.read(phot_dir+obj.lower()+'.txt')
            phot = phot[(phot['wavelength'] >= waverange[0]) & (phot['wavelength'] <= waverange[1])]
            phot = phot[(phot['Name'] == 'IRAC') + (phot['Name'] == 'MIPS') + (phot['Name'] == '2MASS')]

    if spec is not None:
        # The names of all filters overlapping with IRS spectra
        # fil_names = np.array(['IRAC Channel 1','IRAC Channel 2','IRAC Channel 3','IRAC Channel 4',
        #              'MIPS 24um'])
        fil_waves = np.array([3.6, 4.5, 5.8, 8.0, 24.0])

        # getting the min and max of the input spectrum
        waverange = [spec['Wavelength(um)'].min(), spec['Wavelength(um)'].max()]
        wave_filter = (fil_waves >= waverange[0]) & (fil_waves <= waverange[1])

        # empty array to store the calculated spectrophotometry
        specphot = np.empty_like(fil_waves[wave_filter])
        specphot_err = np.empty_like(fil_waves[wave_filter])

        # loop through all wavelengths
        for i in range(len(fil_waves[wave_filter])):
            specphot[i], specphot_err[i] = Spectrophotometry(spec, fil_waves[wave_filter][i])
            if np.isnan(specphot_err[i]):
                # apply a 30% error on Spitzer flux
                # need a reference for this number
                specphot_err[i] = 0.3*specphot[i]
        print(fil_waves[wave_filter])

        # If 'phot' is not defined yet, create a new table, otherwise update the old one
        try:
            phot
        except NameError:
            arr = {'wavelength': fil_waves[wave_filter],
                   'flux(Jy)': specphot,
                   'error(Jy)': specphot_err}
            phots = apy.table.Table(arr)
        else:
            print(np.hstack((phot['wavelength'].data, fil_waves[wave_filter])))
            arr = {'wavelength': np.hstack((phot['wavelength'].data, fil_waves[wave_filter])),
                   'flux(Jy)': np.hstack((phot['flux(Jy)'].data, specphot)),
                   'error(Jy)': np.hstack((phot['error(Jy)'], specphot_err))}
            phots = apy.table.Table(arr)
    elif (spec is None) & (usephot is False):
        return None
    else:
        phots = phot

    # need to fit alpha = dlog10(lambda*F_lambda) / dlog10(lambda)
    # first plot the photometry into log10(lambda*F_lambda) vs log10(lambda)

    sorter = np.argsort(phots['wavelength'])
    phots = phots[sorter]

    # exclude NaN value in flux, and negative flux
    phots = phots[(np.isnan(phots['flux(Jy)']) == False) & (phots['flux(Jy)'] >= 0)]
    #  ÷(np.isnan(phots['error(Jy)']) == False)]

    print(phots)

    # interrupt the script if the available fluxes are fewer than 3
    if len(phots) < 3:
        return None

    # unit: erg/s/cm2/cm
    phots['F_lambda'] = phots['flux(Jy)']*c/(phots['wavelength']/1e4)**2 *1e-23
    phots['F_lambda_err'] = phots['error(Jy)']*c/(phots['wavelength']/1e4)**2 *1e-23

    x = np.log10(phots['wavelength'].data*1e-4)
    y = np.log10(phots['wavelength'].data*1e-4*phots['F_lambda'].data)
    y_err_hi = np.log10(phots['F_lambda'].data+phots['F_lambda_err'].data)-np.log10(phots['F_lambda'].data)
    y_err_low = np.log10(phots['F_lambda'].data)-np.log10(phots['F_lambda'].data-phots['F_lambda_err'].data)

    print(y_err_hi, y_err_low)

    # fitting part
    # p, cov = np.polyfit(x, y, 1, w=1/y_err**2, cov=True)
    # use my own fitting code
    from leastsqfit import lin_leastsqfit
    [y_fit, y_fit_err, p, cov, s_min] = lin_leastsqfit(x, y, y_err_low)
    # y_fit = p[0]*x + p[1]
    print(phots)
    print(p)
    print(cov)

    if plotdir != None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # d, = ax.plot(x, y, 'o')
        # fit, = ax.plot(x, y_fit, color='k')

        d = ax.errorbar(x, y, yerr=[y_err_hi, y_err_low], marker='o')
        fit, = ax.plot(x, y_fit, color='k')
        ax.fill_between(np.squeeze(x), y_fit+y_fit_err, y_fit-y_fit_err, facecolor='Grey', edgecolor='None', alpha=0.7)

        ax.set_xlabel(r'$log(\lambda)$', fontsize=20)
        ax.set_ylabel(r'$log(\lambda F_{\lambda})$', fontsize=20)
        ax.legend([d, fit], ['Data', 'Fit'], loc='upper left')
        ax.text(0.58,0.1, r'$\alpha$={:5.2f}$\pm${:5.2f}'.format(p[1,0], cov[1,1]**0.5),
                transform=ax.transAxes, fontsize=18)

        fig.savefig(plotdir+obj+'_AlphaClassification.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    return p[1,0], cov[1,1]**0.5

obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301','TMR1',
            'TMC1A','TMC1','IRAS15398','RNO91','GSS30-IRS1',
            'VLA1623','WL12','RCrA-IRS5A','L483','B335',
            'DKCha']
spec = 1
phot_dir = '/Users/yaolun/data/herschel_phot/cops-spire_alpha/'
spec_dir = '/Volumes/SD-Mac/Dropbox/cops-spire/IRS_spec/reformatted/'
# plotdir = '/Volumes/SD-Mac/Dropbox/cops-spire/AlphaClassification/'
plotdir = '/Users/yaolun/research/cops-spire/AlphaClassification/'
from astropy.io import ascii
import os

for obj in obj_list:
    if obj != 'VLA1623':
        continue
    print(obj)
    if os.path.exists(spec_dir+obj.lower()+'.txt'):
        spec = ascii.read(spec_dir+obj.lower()+'.txt', comment='%')
    else:
        spec = None
    alpha = AlphaClassification(obj, phot_dir, usephot=True, plotdir=plotdir, spec=None)
    print(alpha)

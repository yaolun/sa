# This is a python 3 script
def AlphaClassification(obj, phot_dir, spec=None, waverange=None, plotdir=None):
    import numpy as np
    from astropy.io import ascii
    from Spectrophotometry import Spectrophotometry
    import matplotlib.pyplot as plt
    import astropy.constants as const
    c = const.c.cgs.value

    if waverange == None:
        waverange = [2.0, 25.0]

    # read in the photometry ASCII file
    phot = ascii.read(phot_dir+obj.lower()+'.txt')
    phot = phot[(phot['wavelength'] >= waverange[0]) & (phot['wavelength'] <= waverange[1])]

    # something needs to be put in here for including spectrophotometry with photometry

    # interrupt the script if the available fluxes are fewer than 2
    if len(phot) < 3:
        return None

    # need to fit alpha = dlog10(lambda*F_lambda) / dlog10(lambda)
    # first plot the photometry into log10(lambda*F_lambda) vs log10(lambda)

    sorter = np.argsort(phot['wavelength'])
    phot = phot[sorter]

    # unit: erg/s/cm2/cm
    phot['F_lambda'] = phot['flux(Jy)']*c/(phot['wavelength']/1e4)**2 *1e-23

    x = np.log10(phot['wavelength'].data*1e-4*phot['F_lambda'].data)
    y = np.log10(phot['wavelength'].data*1e-4)
    print(x,y)
    # fitting part
    p, cov = np.polyfit(x, y, 1, cov=True)
    y_fit = p[0]*x + p[1]

    if plotdir != None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        d, = ax.plot(x, y, 'o')
        fit, = ax.plot(x, y_fit, color='k')
        ax.set_xlabel(r'$log(\lambda F_{\lambda})$', fontsize=20)
        ax.set_ylabel(r'$log(\lambda)$', fontsize=20)
        ax.legend([d, fit], ['Data', 'Fit'])
        ax.text(0.58,0.1, r'$\alpha$={:5.2f}$\pm${:5.2f}'.format(p[0], cov[0,0]**0.5),
                transform=ax.transAxes, fontsize=18)

        fig.savefig(plotdir+obj+'_AlphaClassification.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    return p[0], cov[0,0]**0.5

obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301','TMR1',
            'TMC1A','TMC1','IRAS15398','RNO91','GSS30-IRS1',
            'VLA1623','WL12','RCrA-IRS5A','L483','B335',
            'DKCha']
spec = 1
phot_dir = '/Users/yaolun/data/herschel_phot/'
plotdir = '/Volumes/SD-Mac/Dropbox/cops-spire/AlphaClassification/'
for obj in obj_list:
    print(obj)
    alpha = AlphaClassification(obj, phot_dir, plotdir=plotdir)
    print(alpha)

def AlphaHyperion(rtout, aperfile, dstar, wave_center, lbollsmm=False):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from hyperion.model import ModelOutput, Model
    from hyperion.util.constants import pc, c, lsun, au
    from astropy.io import ascii

    def getAlpha(spec, wave_center, plot=False, plotname=None):
        """
        spec = {'Wavelength(um)': wave,
                'Flux_Density(Jy)': flux,
                'Uncertainty(Jy)': unc}
        """
        from astropy.modeling import models, fitting, powerlaws
        from scipy.interpolate import interp1d
        import astropy.constants as const
        import numpy as np
        c = const.c.cgs.value

        # initial guess
        x_ref = c/1e5/wave_center
        f_flux = interp1d(c/1e5/spec['Wavelength(um)'], spec['Flux_Density(Jy)'])
        amp = f_flux(x_ref)
        alpha = 0

        # unit conversions
        freq_dum = (c/1e5/spec['Wavelength(um)'])\
                   [(c/1e5/spec['Wavelength(um)']>= x_ref-100) & (c/1e5/spec['Wavelength(um)']<= x_ref+100)]
        flux_dum = spec['Flux_Density(Jy)']\
                   [(c/1e5/spec['Wavelength(um)']>= x_ref-100) & (c/1e5/spec['Wavelength(um)']<= x_ref+100)]

        pow_model = powerlaws.PowerLaw1D(amp, x_ref, alpha)
        fitter = fitting.LevMarLSQFitter()
        fit = fitter(pow_model, freq_dum, flux_dum)

        alpha = -fit.alpha.value
        if fitter.fit_info['param_cov'] is None:
            alpha_err = np.nan
        else:
            alpha_err = fitter.fit_info['param_cov'][2,2]**0.5

        if plot:
            # to avoid X server error
            import matplotlib as mpl
            mpl.use('Agg')
            #
            import matplotlib.pyplot as plt

            fig = plt.figure(figsize=(10,6))
            ax = fig.add_subplot(111)
            ax.plot(freq_dum, fit(freq_dum), '-', color='k')
            ax.text(0.35, 0.15, r'$\alpha_{250,350,500} = %3.2f, %3.2f, %3.2f$' % (fitted_alpha[0], fitted_alpha[1], fitted_alpha[2]),
                    transform=ax.transAxes, fontsize=18)
            ax.text(0.35, 0.05, r'$\sigma_{\alpha}\,(250,350,500) = %5.3f, %5.3f, %5.3f$' % (fitted_alpha_err[0], fitted_alpha_err[1], fitted_alpha_err[2]),
                    transform=ax.transAxes, fontsize=18)
            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on()
            ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
            ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
            ax.set_xlabel('Frequency [GHz]', fontsize=18)
            ax.set_ylabel('Flux Density [Jy]', fontsize=18)
            ax.set_xlim([400,2000])

            fig.savefig(plotname+'Alpha500Hyperion.pdf', format='pdf', dpi=300, bbox_inches='tight')
            fig.clf()

        return (alpha, alpha_err)

    def lsubmm(low_wave, spec, dist):
        """
        spec = 'Wavelength(um)' and 'Flux_Density(Jy)'
        dist: distance in parsec
        """
        import sys
        import os
        sys.path.append(os.path.expanduser('~')+'/programs/misc/hyperion/')
        from l_bol import l_bol

        l = l_bol(spec['Wavelength(um)'][spec['Wavelength(um)'] >= low_wave],
                  spec['Flux_Density(Jy)'][spec['Wavelength(um)'] >= low_wave], dist)

        return l

    # get aperture list
    aperture = ascii.read(aperfile)
    # create the non-repetitive aperture list and index array
    aper_reduced = list(set(aperture['aperture(arcsec)']))

    # read in Hyperion output
    m = ModelOutput(rtout)

    # list for storing outputs
    alpha = []
    alpha_err = []
    aperture_list = []

    # option for calculating Lbol/Lsumm at the same time
    if lbollsmm:
        lbol = []
        lsmm = []

    # fig = plt.figure()
    # ax = fig.add_subplot(111)

    # iterate through apertures
    for i in range(0, len(aper_reduced)):
        # getting simulated SED from Hyperion output. (have to match with the reduced index)
        sed_dum = m.get_sed(group=i, inclination=0,aperture=-1,distance=dstar*pc, uncertainties=True)
        # get aperture from the output SED
        aperture_list.append(sed_dum.ap_min/(dstar*pc)/(1/3600.*np.pi/180.)*2)

        # ax.plot(sed_dum.wav, sed_dum.val, label='{:4.1f}'.format(aperture_list[i]))

        # construct 'spec' dictionary
        sorter = np.argsort(sed_dum.wav)
        spec = {'Wavelength(um)': sed_dum.wav[sorter],
                'Flux_Density(Jy)': (sed_dum.val/(c/sed_dum.wav*1e4)*1e23)[sorter],
                'Uncertainty(Jy)': (sed_dum.unc/(c/sed_dum.wav*1e4)*1e23)[sorter]}

        # option for calculating Lbol/Lsumm at the same time
        if lbollsmm:
            # model_name = os.path.basename(rtout).split('.')[0]
            # data = ascii.read(lbollsmm+model_name+'_sed_w_aperture.txt')
            # specphot = {'Wavelength(um)': data['wave'],
            #             'Flux_Density(Jy)': data['vSv']/(c/np.array(data['wave'])*1e4)*1e23}
            # lbol.append(lsubmm(specphot['Wavelength(um)'].min(), specphot, dstar))
            # lsmm.append(lsubmm(350.0, specphot, dstar))

            lbol.append(lsubmm(spec['Wavelength(um)'].min(), spec, dstar))
            lsmm.append(lsubmm(350.0, spec, dstar))

        # get alpha
        plotname = '/home/bettyjo/yaolun/test/'+aperture_list[-1]+'_'
        print(plotname)
        alpha_dum, alpha_err_dum = getAlpha(spec, wave_center, plot=True, plotname=plotname)
        alpha.append(alpha_dum)
        alpha_err.append(alpha_err_dum)

    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.legend(ncol=2, loc='best')
    # fig.savefig('/Users/yaolun/test/sed_test.pdf', format='pdf', dpi=300, bbox_inches='tight')
    if not lbollsmm:
        return np.array(aperture_list), np.array(alpha), np.array(alpha_err)
    else:
        return np.array(aperture_list), np.array(alpha), np.array(alpha_err), np.array(lbol), np.array(lsmm)

# rtout = '/Users/yaolun/test/model4/model4.rtout'
# aperfile = '/Users/yaolun/research/bhr71/best_calibrated/aperture.txt'
# lbollsmm = '/Users/yaolun/research/bhr71/hyperion/controlled/'
# dstar = 200.0
# wave_center = 500.0
# results = AlphaHyperion(rtout, aperfile, dstar, wave_center, lbollsmm=lbollsmm)

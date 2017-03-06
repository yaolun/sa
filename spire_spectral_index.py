def spire_spectral_index(outdir, obsid, obj):

    # to avoid X server error
    import matplotlib as mpl
    mpl.use('Agg')
    #
    # Note that the spectral index works in frequency
    import matplotlib.pyplot as plt
    from astropy.modeling import models, fitting, powerlaws
    from scipy.interpolate import interp1d
    import astropy.constants as const
    from astropy.io import ascii
    c = const.c.cgs.value

    # Read in spectra
    spire_sect = ascii.read(outdir+obsid+'_spire_sect.txt', data_start=4)

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)

    spire, = ax.plot(c/1e5/spire_sect['wave_segm1_0'], spire_sect['flux_segm1_0'], linewidth=1, color='b')
    ax.plot(c/1e5/spire_sect['wave_segm2_0'], spire_sect['flux_segm2_0'], linewidth=1, color='r')

    f_ssw = interp1d(c/1e5/spire_sect['wave_segm2_0'], spire_sect['flux_segm2_0'])
    f_slw = interp1d(c/1e5/spire_sect['wave_segm1_0'], spire_sect['flux_segm1_0'])
    fitted_alpha = []
    fitted_alpha_err = []

    for band in [250, 350, 500]:
        x_ref = c/1e5/band
        if band != 250:
            amp = f_slw(x_ref)
            freq_dum = (c/1e5/spire_sect['wave_segm1_0'])\
                       [(c/1e5/spire_sect['wave_segm1_0'] >= x_ref-200) & (c/1e5/spire_sect['wave_segm1_0'] <= x_ref+200)]
            flux_dum = spire_sect['flux_segm1_0']\
                       [(c/1e5/spire_sect['wave_segm1_0'] >= x_ref-200) & (c/1e5/spire_sect['wave_segm1_0'] <= x_ref+200)]
        else:
            amp = f_ssw(x_ref)
            freq_dum = (c/1e5/spire_sect['wave_segm2_0'])\
                       [(c/1e5/spire_sect['wave_segm2_0'] >= x_ref-200) & (c/1e5/spire_sect['wave_segm2_0'] <= x_ref+200)]
            flux_dum = spire_sect['flux_segm2_0']\
                       [(c/1e5/spire_sect['wave_segm2_0'] >= x_ref-200) & (c/1e5/spire_sect['wave_segm2_0'] <= x_ref+200)]
        alpha = 0

        pow_model = powerlaws.PowerLaw1D(amp, x_ref, alpha)
        fitter = fitting.LevMarLSQFitter()
        fit = fitter(pow_model, freq_dum, flux_dum)

        ax.plot(freq_dum, fit(freq_dum), '-', color='k')
        # take negative sign because the frequency array is reversed
        fitted_alpha.append(-fit.alpha.value)
        fitted_alpha_err.append(fitter.fit_info['param_cov'][2,2]**0.5)
        print fit.alpha.value, '+/-', fitter.fit_info['param_cov'][2,2]**0.5

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

    fig.savefig(outdir+obj+'_spire_alpha.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

    # write out the alpha index
    foo = open(outdir+obj+'_alpha.txt','w')
    foo.write('250um \t 350um \t 500um \n')
    foo.write('%8.6f \t %8.6f \t %8.6f \n' % (fitted_alpha[0], fitted_alpha[1], fitted_alpha[2]))
    foo.write('%8.6f \t %8.6f \t %8.6f \n' % (fitted_alpha_err[0], fitted_alpha_err[1], fitted_alpha_err[2]))
    foo.close()
#
# # COPS-SPIRE objects selected from successfully SECT reduction
# # COPS-SPIRE source list
# obsid_spire = [1342242620,1342242621,1342245084,1342245094,1342245857,
#                1342247625,1342248246,1342248249,1342249053,1342249470,
#                1342249474,1342249475,1342249476,1342249477,1342250509,
#                1342250510,1342250512,1342250515,1342251285,1342251286,
#                1342251287,1342251290,1342253646,1342253649,1342253652,
#                1342254037,1342252897]
#
# # SECT cannot converage at L1489 1342249473, L1527 1342250511, HH100 1342252897
# # mapping observation IRS44/46 1342251289
#
# obj_list_spire = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
#                   'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
#                   'L1455-IRS3','B1-a','B1-c','IRAS03301','TMR1',
#                   'TMC1A','TMC1','IRAS15398','RNO91','GSS30-IRS1',
#                   'VLA1623','WL12','RCrA-IRS5A','L483','B335',
#                   'DKCha','HH100']
#
# # archive_dir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/'
# archive_dir = '/Volumes/SD-Mac/CDF_archive_v2/'
#
# #
# for o in obsid_spire:
#     obj = obj_list_spire[obsid_spire.index(o)]
#     # if obj != 'BHR71':
#     #     continue
#     print obj
#     spire_spectral_index(archive_dir+obj+'/spire/data/', str(o), obj)

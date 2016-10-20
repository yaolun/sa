def spire_spectral_index(outdir, obsid):

    # Note that the spectral index works in frequency
    from astropy.modeling import models, fitting, powerlaws
    from scipy.interpolate import interp1d
    import astropy.constants as const
    from astropy.io import ascii
    c = const.c.cgs.value

    # Read in spectra
    spire_sect = ascii.read(outdir+obsid+'spire_sect.txt', data_start=4)

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)

    spire, = ax.plot(c/1e5/spire_sect['wave_segm1_0'], spire_sect['flux_segm1_0'], linewidth=1, color='b')
    ax.plot(c/1e5/spire_sect['wave_segm2_0'], spire_sect['flux_segm2_0'], linewidth=1, color='r')

    f_ssw = interp1d(c/1e5/spire_sect['wave_segm2_0'], spire_sect['flux_segm2_0'])
    f_slw = interp1d(c/1e5/spire_sect['wave_segm1_0'], spire_sect['flux_segm1_0'])
    fitted_alpha = []

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
        print -fit.alpha

    ax.text(0.9, 0.1, r'$\rm{\alpha_{250,350,500} = %3.2f, %3.2f, %3.2f}$' % (fitted_alpha), transform=ax.transAxes)
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
    ax.set_xlabel('Frequency [GHz]', fontsize=18)
    ax.set_ylabel('Flux Density [Jy]', fontsize=18)
    ax.set_xlim([400,2000])

    fig.savefig(outdir+obsid+'_spire_alpha.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

    # write out the alpha index
    foo = open(outdir+obsid+'_alpha.txt','r')
    foo.write('250um \t 350um \t 500um \n')
    foo.write('%8.6f \t %8.6f \t %8.6f \n' % (fitted_alpha))
    foo.close()

# COPS-SPIRE objects selected from successfully SECT reduction
Obsid=[1342242620,1342242621,1342245084,1342245094,1342245857,
       1342247625,1342248246,1342248249,1342249053,1342249470,
       1342249474,1342249475,1342249476,1342249477,1342250509,
       1342250510,1342250512,1342250515,1342251285,1342251286,
       1342251287,1342251290,1342253646,1342253649,1342253652,
       1342254037]
outdir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/'

for o in Obsid:
    spire_spectral_index(outdir, o)

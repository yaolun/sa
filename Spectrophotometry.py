def Spectrophotometry(spec, phot_wave, filter_dir=None):

    import numpy as np
    from scipy.interpolate import interp1d
    from astropy.io import ascii
    from phot_filter import phot_filter  # from specrtal_analysis
    import os

    # preset uncertainty is NaN
    unc_aper = np.nan

    # mask the NaN flux
    spec = spec[np.isnan(spec['Flux_Density(Jy)']) == False]
    # sort
    sorter = np.argsort(spec['Wavelength(um)'].data)
    spec = spec[sorter]

    if filter_dir == None:
        filter_dir = os.path.expanduser('~')+'/programs/spectra_analysis/'

    # function for properly calculating uncertainty of spectrophotometry value
    def unc_spectrophoto(wl, unc, trans):
        # adopting smiliar procedure as Trapezoidal rule
        # (b-a) * [ f(a) + f(b) ] / 2
        #
        return ( np.sum( trans[:-1]**2 * unc[:-1]**2 * (wl[1:]-wl[:-1])**2 ) / np.trapz(trans, x=wl)**2 )**0.5

    # apply the filter function
    # decide the filter name
    if phot_wave == 70:
        fil_name = 'Herschel PACS 70um'
    elif phot_wave == 100:
        fil_name = 'Herschel PACS 100um'
    elif phot_wave == 160:
        fil_name = 'Herschel PACS 160um'
    elif phot_wave == 250:
        fil_name = 'Herschel SPIRE 250um'
    elif phot_wave == 350:
        fil_name = 'Herschel SPIRE 350um'
    elif phot_wave == 500:
        fil_name = 'Herschel SPIRE 500um'
    elif phot_wave == 3.6:
        fil_name = 'IRAC Channel 1'
    elif phot_wave == 4.5:
        fil_name = 'IRAC Channel 2'
    elif phot_wave == 5.8:
        fil_name = 'IRAC Channel 3'
    elif phot_wave == 8.0:
        fil_name = 'IRAC Channel 4'
    elif phot_wave == 24:
        fil_name = 'MIPS 24um'
    elif phot_wave == 850:
        fil_name = 'SCUBA 850WB'
    else:
        fil_name = None

    if fil_name != None:
        filter_func = phot_filter(fil_name, filter_dir)
        filter_func = filter_func[(filter_func['wave']/1e4 >= min(spec['Wavelength(um)']))*\
                                  ((filter_func['wave']/1e4 >= 54.8)+(filter_func['wave']/1e4 <= 36.0853))*\
                                  ((filter_func['wave']/1e4 <= 95.05)+(filter_func['wave']/1e4 >=103))*\
                                  ((filter_func['wave']/1e4 <= 190.31)+(filter_func['wave']/1e4 >= 195))*\
                                  (filter_func['wave']/1e4 <= max(spec['Wavelength(um)']))]

        f = interp1d(spec['Wavelength(um)'], spec['Flux_Density(Jy)'])
        flux_aper = np.trapz(f(filter_func['wave']/1e4)*\
                                  filter_func['transmission'],x=filter_func['wave']/1e4 )/\
                       np.trapz(filter_func['transmission'], x=filter_func['wave']/1e4)
        if 'Uncertainity(Jy)' in spec.colnames:
            f_unc = interp1d(spec['Wavelength(um)'], spec['Uncertainity(Jy)'])
            unc_aper = unc_spectrophoto(filter_func['wave']/1e4,
                                        f_unc(filter_func['wave']/1e4),
                                        filter_func['transmission'])
    else:
        # use a rectangle function the average the simulated SED
        # apply the spectral resolution
        if (phot_wave < 50.) & (phot_wave >= 5):
            res = 60.
        elif phot_wave < 5:
            res = 10.
        else:
            res = 1000.
        ind = np.where((spec['Wavelength(um)'] < phot_wave*(1+1./res)) & (spec['Wavelength(um)'] > phot_wave*(1-1./res)))
        if len(ind[0]) != 0:
            flux_aper = np.nanmean(spec['Flux_Density(Jy)'][ind])
            if 'Uncertainty(Jy)' in spec.colnames:
                unc_aper = np.nanmean(spec['Uncertainty(Jy)'][ind])
        else:
            f = interp1d(spec['Wavelength(um)'], spec['Flux_Density(Jy)'])
            flux_aper = f(phot_wave)
            if 'Uncertainty(Jy)' in spec.colnames:
                f_unc = interp1d(spec['Wavelength(um)'], spec['Uncertainty(Jy)'])
                unc_aper = f_unc(phot_wave)

    return flux_aper, unc_aper

# from astropy.io import ascii
# # spec = ascii.read('/Volumes/SD-Mac/CDF_archive_v2/BHR71/pacs/data/BHR71_pacs_weighted.txt')
# spec = ascii.read('/Users/yaolun/research/bhr71/best_calibrated/BHR71_spire_corrected_continuum.txt')
# phot_wave = 500.
# (flux_aper, unc_aper) = Spectrophotometry(spec, phot_wave)
# print(flux_aper, unc_aper)

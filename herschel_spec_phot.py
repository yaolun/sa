def herschel_spec_phot(wl, flux, pacs=True, spire=True, filter_func=True):
    import numpy as np
    import os
    from scipy.interpolate import interp1d
    import sys
    sys.path.append(os.path.expanduser('~')+'/programs/spectra_analysis/')
    from phot_filter import phot_filter
    phot_wl = []
    if pacs == True:
        phot_wl.extend([70,100,160])
    if spire == True:
        phot_wl.extend([250,350,500])
    phot_wl = np.array(phot_wl)
    phot_flux = np.empty_like(phot_wl)

    # clean up NaN value in observation
    wl = wl[np.isnan(flux) == False]
    flux = flux[np.isnan(flux) == False]

    # warning message for checking both PACS and SPIRE bands
    if (min(wl) >= 60) or (max(wl) <= 650):
        print 'Incomplete observed spectra detected!  Aborting...'
        return None

    for i in range(len(phot_wl)):
        if filter_func == False:
            res = 3     # temp. resolution to mimic the photometry filters
            ind = np.where((wl < phot_wl[i]*(1+1./res)) & (wl > phot_wl[i]*(1-1./res)))
            if len(ind[0]) != 0:
                phot_flux[i] = np.nanmean(flux[ind])
        else:
            # apply the filter function
            # decide the filter name
            if phot_wl[i] == 70:
                fil_name = 'Herschel PACS 70um'
            elif phot_wl[i] == 100:
                fil_name = 'Herschel PACS 100um'
            elif phot_wl[i] == 160:
                fil_name = 'Herschel PACS 160um'
            elif phot_wl[i] == 250:
                fil_name = 'Herschel SPIRE 250um'
            elif phot_wl[i] == 350:
                fil_name = 'Herschel SPIRE 350um'
            elif phot_wl[i] == 500:
                fil_name = 'Herschel SPIRE 500um'

            filter_func = phot_filter(fil_name)

            # trim the filter function
            # if phot_wl[i] in [70,100,160]:
            filter_func = filter_func[(filter_func['wave']/1e4 >= max(54.8,min(wl)))*\
                                      ((filter_func['wave']/1e4 <= 95.05)+(filter_func['wave']/1e4 >=103))*\
                                      ((filter_func['wave']/1e4 <= 190.31)+(filter_func['wave']/1e4 >= 195))]
            # elif phot_wl[i] in [250,350,500]:
            #   filter_func = filter_func[(filter_func['wave']/1e4 >= 195)]

            f = interp1d(wl, flux)
            # print fil_name
            # print filter_func['wave']/1e4
            # print min(wl), max(wl)
            phot_flux[i] = np.trapz(filter_func['wave']/1e4, f(filter_func['wave']/1e4)*filter_func['transmission'])/np.trapz(filter_func['wave']/1e4, filter_func['transmission'])

    return phot_wl, phot_flux
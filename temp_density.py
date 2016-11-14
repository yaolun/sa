def envmass_fit(dustpath, obsdir, obj, plotdir):

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    import astropy.constants as const
    c = const.c.cgs.value
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    pc = const.pc.cgs.value
    AU = const.au.cgs.value
    MS = const.M_sun.cgs.value
    mh = const.m_p.cgs.value + const.m_e.cgs.value

    # dust model
    # Weingartner & Draine
    # dustfile = ascii.read('/Users/yaolun/Google Drive/dust_model/Weingartner_Draine/kext_albedo_WD_MW_5.5A_30_D03_nocomment.all',
    #                       header_start=64, data_start=67)

    # Ossenkopf & Hennings
    # dustfile = ascii.read('/Users/yaolun/programs/misc/oh5_hyperion.txt', names=['nu', 'albedo','chi','g'])

    # Use Neal's original version of OH5 dust
    dustfile = np.genfromtxt(dustpath, skip_header=2, skip_footer=1).T

    # read in observations
    import sys
    sys.path.append('/Users/yaolun/programs/misc/hyperion/')
    from get_obs import get_obs
    obs = get_obs(obsdir, obj=obj, spec=False)  # flux in Jy

    def bbfunc(omegab, kappa):
        def bbfunc_dum(nu, *p):
            tdust = p[0]
            Nd = p[1]

            # kappa in cm2 of dust
            # nd = 1/cm2
            import astropy.constants as const
            import numpy as np

            h = const.h.cgs.value
            c = const.c.cgs.value
            k = const.k_B.cgs.value

            bb = 2.*h*nu**3/c**2 / (np.exp(h*nu/k/tdust)-1)

            return bb * kappa * Nd * omegab
        return bbfunc_dum

    def bbfunc2(nu, omegab, kappa, *p):
        tdust = p[0]
        Nd = p[1]

        # kappa in cm2/g of dust
        # nd = g/cm2
        import astropy.constants as const
        import numpy as np

        h = const.h.cgs.value
        c = const.c.cgs.value
        k = const.k_B.cgs.value

        bb = 2.*h*nu**3/c**2 / (np.exp(h*nu/k/tdust)-1)

        return bb * kappa * Nd * omegab

    # get the kappa_nu from dust model
    # the unit is cm2/g
    # this part depends on the format of input dust model
    # kappa = (1-dustfile['albedo']) * dustfile['chi']
    # kappa = dustfile['chi']
    # Neal's version of OH5 dust
    kappa = (dustfile[0,:]+dustfile[1,:])
    nu_kappa = c/dustfile[3,:]*1e4

    # for oh5, md = 1/(1.0395087779657002e14)
    md = 1/(1.0395087779657002e14)

    # take down to 160 um
    trimmer = (obs['phot'][0] >= 100) & (obs['phot'][0] < 1400)
    #
    nu = c/obs['phot'][0][trimmer]/1e-4
    f_phot = obs['phot'][1][trimmer]*1e-23
    f_phot_err = obs['phot'][2][trimmer]*1e-23

    # sort
    sorter = np.argsort(nu)
    nu = nu[sorter]
    f_phot = f_phot[sorter]
    f_phot_err = f_phot_err[sorter]

    # load aperture
    aperture = ascii.read(obsdir+obj.lower()+'.txt')
    aper = []
    for wl in obs['phot'][0][trimmer]:
        if wl in aperture['wavelength']:
            aper.append(aperture['aperture(arcsec)'][aperture['wavelength'] == wl][0])
        else:
            print 'No aperture found for %f um' % wl
            aper_dum = raw_input('What is the aperture size?')
            aper.append(aper_dum)
    aper = np.array(aper)
    omegab = np.pi*(aper/2)**2/4.25e10

    # interpolate kappa
    from scipy.interpolate import interp1d
    f = interp1d(nu_kappa, kappa)
    kappa_int = f(nu)

    from scipy.optimize import curve_fit
    p = curve_fit(bbfunc(omegab, kappa_int), nu, f_phot, sigma=f_phot_err, p0=[20, 1e22])
    print p[0]
    perr = np.sqrt(np.diag(p[1]))
    print perr

    area = (90.*200/2*AU)**2*np.pi
    print p[0][1]*area * md *100 /MS * 0.315*pc/(90*200./2*AU)
    print perr[1]*area * md *100 /MS * 0.315*pc/(90*200./2*AU)

    # evaluate RJ limit
    print 'hv/kT = ', h*nu/k/p[0][0]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    obs_trimmer = (obs['phot'][0] > 60) & (obs['phot'][0] < 1400)
    nu = c/obs['phot'][0][obs_trimmer]/1e-4
    f_phot = obs['phot'][1][obs_trimmer]*1e-23
    f_phot_err = obs['phot'][2][obs_trimmer]*1e-23
    kappa_int = f(nu)

    # re-get the aperture for observational data
    aper_obs = []
    for wl in obs['phot'][0][obs_trimmer]:
        if wl in aperture['wavelength']:
            aper_obs.append(aperture['aperture(arcsec)'][aperture['wavelength'] == wl][0])
        else:
            print 'No aperture found for %f um' % wl
            aper_dum = raw_input('What is the aperture size?')
            aper_obs.append(aper_dum)
    aper_obs = np.array(aper_obs)
    print aper_obs

    omegab_obs = np.pi*(aper_obs/2)**2/4.25e10

    obs_yerr_hi = np.log10(f_phot*nu+f_phot_err*nu)-np.log10(f_phot*nu)
    obs_yerr_low = np.log10(f_phot*nu)-np.log10(f_phot*nu-f_phot_err*nu)

    obs_plot = ax.errorbar(np.log10(c/nu*1e4), np.log10(f_phot*nu), yerr=[obs_yerr_hi, obs_yerr_low],
                           fmt='o', markersize=10, markeredgewidth=1.5,
                           mec='b', ecolor='b', color='None', elinewidth=1.2, capthick=1.2, barsabove=True)

    yerr_hi = np.log10(bbfunc2(nu, omegab_obs, kappa_int, *p[0]+perr)*nu) - np.log10(bbfunc2(nu, omegab_obs, kappa_int, *p[0])*nu)
    yerr_low = np.log10(bbfunc2(nu, omegab_obs, kappa_int, *p[0])*nu) - np.log10(bbfunc2(nu, omegab_obs, kappa_int, *p[0]-perr)*nu)

    fit = ax.errorbar(np.log10(c/nu*1e4), np.log10(bbfunc2(nu, omegab_obs, kappa_int, *p[0])*nu),
                      yerr=[yerr_low, yerr_hi], fmt='o', linestyle='-', markersize=10, markeredgewidth=1.5,
                      color='None', mec='g',ecolor='g', elinewidth=1.2, capthick=1.2, barsabove=True)

    ax.legend([obs_plot, fit], [r'$\rm{Photometry}$', r'$\rm{Fit\,with\,OH5\,dust}$'],
              loc='best', numpoints=1, fontsize=16)

    ax.text(0.1, 0.2, r'$\rm{T_{dust}= %3.1f \pm %3.1f \,K}$' % (p[0][0], perr[0]), transform=ax.transAxes, fontsize=24)
    # expo = np.floor(np.log10(p[0][1]))
    # ax.text(0.1, 0.1, r'$\rm{N_{dust}= %2.1f \pm %2.1f \times 10^{%d}\,cm^{-2}}$' % (p[0][1]/10**expo, perr[1]/10**expo, expo),
    #         transform=ax.transAxes, fontsize=24)
    ngas = p[0][1]*md*100/mh/2.37
    ngas_err = perr[1]*md*100/mh/2.37
    expo = np.floor(np.log10(ngas))
    ax.text(0.1, 0.1, r'$\rm{N_{gas}= %2.1f \pm %2.1f \times 10^{%d}\,cm^{-2}}$' % (ngas/10**expo, ngas_err/10**expo, expo),
            transform=ax.transAxes, fontsize=24)

    ax.set_xlim([np.log10(50), np.log10(1500)])

    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
    ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
    ax.set_xlabel(r'$\rm{log(Wavelength)\,[\mu m]}$', fontsize=18)
    ax.set_ylabel(r'$\rm{log(\nu S_{\nu})\,[erg\,s^{-1}\,cm^{-2}]}$', fontsize=18)

    fig.savefig(plotdir+obj+'_greybb_fit.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.clf()

    return p, perr

dustpath = '/Users/yaolun/Google Drive/dust_model/Dust_OH5_Evans_Shirley/sigma.oh5.ref'
obsdir = '/Volumes/SD-Mac/CDF_archive/'
plotdir = '/Users/yaolun/test/greybb_fit/'
obj_list = ['RCrA-IRS7B','RCrA-IRS7C','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301',
            'TMR1','TMC1A','TMC1','IRAS15398','RNO91',
            'GSS30-IRS1','VLA1623','WL12','RCrA-IRS5A','L483',
            'B335','DKCha']
for o in obj_list:
    envmass_fit(dustpath, obsdir, o, plotdir)

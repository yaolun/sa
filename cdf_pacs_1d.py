def cdf_pacs_1d(osbid, objname, outdir, fits_for_header, line_fitting_dir):
    """
    cubedir: the directory contains FITS rebinned cube
    outcubedir: the directory contains ASCII style spectra for each spaxel
    out1ddir: the directory contains ASCII style 1-D spectrum
    """

    from pacs_weight import pacs_weight
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii, fits
    import pidly
    idl = pidly.IDL('/Applications/exelis/idl83/bin/idl')

    # outdir = '/Users/yaolun/bhr71/best_calibrated/'
    # cubedir = '/Users/yaolun/bhr71/data/HSA/'
    # photpath = '/Users/yaolun/bhr71/best_calibrated/bhr71.txt'

    cubedir = outdir+'data/fits/'

    cubefile = [cubedir+'OBSID_'+obsid[0]+'_blue_finalcubes_slice00_os8_sf7.fits',
                cubedir+'OBSID_'+obsid[0]+'_red_finalcubes_slice00_os8_sf7.fits',
                cubedit+'OBSID_'+obsid[1]+'_blue_finalcubes_slice00_os8_sf7.fits',
                cubedir+'OBSID_'+obsid[1]+'_red_finalcubes_slice00_os8_sf7.fits',]

    idl('.r '+line_fitting_dir+'get_pacs.pro')
    idl.pro('get_pacs', outdir=outdir+'data/cube/', objname=objname, filename=cubefile, suffix='os8_sf7')

    wl, flux = pacs_weight(outdir+'data/cube/', objname, 31.8, outdir+'data/', fits_for_header, suffix='os8_sf7')

    idl('.r '+line_fitting_dir+'gauss.pro')
    idl('.r '+line_fitting_dir+'extract_pacs.pro')

    idl.pro('extract_pacs', indir=outdir+'data/', filename=objname+'_pacs_weighted',
            outdir=outdir+'advanced_products/',
            plotdir=outdir+'advanced_products/plots/', noiselevel=3, ra=0, dec=0, global_noise=20,
            localbaseline=10, opt_width=1, continuum=1, flat=1, object=objname, double_gauss=1, fixed_width=1)
            

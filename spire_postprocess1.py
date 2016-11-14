def spire_postprocess1(indir,outdir,obsid,obj):
    """
    Dedicated for 1-D spectra analysis.
    convert the ASCII format of the spectra to a better format for reading, and perform spectrum trimming.
    """
    from astropy.io import ascii, fits
    import numpy as np

    # read in the spectrum
    spire_spec = ascii.read(indir+obsid+'spire_sect.txt', data_start=4)
    # convert it to the usual format
    spire_wl = np.hstack((spire_spec['wave_segm1_0'][spire_spec['wave_segm1_0'] >= 310].data,
                spire_spec['wave_segm2_0'][(spire_spec['wave_segm2_0'] < 310) & (spire_spec['wave_segm2_0'] > 195)].data))
    spire_flux = np.hstack((spire_spec['flux_segm1_0'][spire_spec['wave_segm1_0'] >= 310].data,
                spire_spec['flux_segm2_0'][(spire_spec['wave_segm2_0'] < 310) & (spire_spec['wave_segm2_0'] > 195)].data))

    sorter = np.argsort(spire_wl)
    spire_wl = spire_wl[sorter].data
    spire_flux = spire_flux[sorter].data

    # Write to file
    foo = open(outdir+obj+'/spire/data/'+obj+'_spire_corrected.txt','w')
    foo.write('%s \t %s \n' % ('Wavelength(um)', 'Flux_Density(Jy)'))
    for i in range(len(spire_wl)):
        foo.write('%f \t %f \n' % (spire_wl[i], spire_flux[i]))
    foo.close()

    # perform line fitting
    import pidly
    # grad13yy
    # idl = pidly.IDL('/Applications/exelis/idl83/bin/idl')
    # bettyjo
    idl = pidly.IDL('/opt/local/exelis/idl83/bin/idl')
    idl('.r '+os.path.expanduser('~')+'/programs/line_fitting/extract_spire.pro')
    # read the coordinates from cube fits file
    hdu = fits.open(obsid+'_spectrum_extended_HR_aNB_15.fits')
    ra = hdu[0].header['RA']
    dec = hdu[0].header['DEC']
    idl.pro('extract_spire', indir=outdir+obj+'/spire/data/', filename=obj+'_spire_corrected',
            outdir=outdir+obj+'/spire/advanced_products/',
            plotdir=outdir+obj+'/spire/advanced_products/plots/', localbaseline=10,
            global_noise=20, ra=ra, dec=dec, noiselevel=3, fx=1, object=obj, flat=1, continuum=1,
            double_gauss=1, no_plot=0)

def spire_cubeprocess(indir, outdir, obsid, obj):
    """
    Dedicated for cube analysis
    1. extract the spectrum of each pixel from SPIRE cube.
    2. perform line fitting on each pixel.
    """
    from astropy.io import ascii, fits
    import numpy as np
    import pidly

    # grad13yy
    # idl = pidly.IDL('/Applications/exelis/idl83/bin/idl')
    # bettyjo
    idl = pidly.IDL('/opt/local/exelis/idl83/bin/idl')
    idl('.r '+os.path.expanduser('~')+'/programs/line_fitting/get_spire.pro')
    idl('.r '+os.path.expanduser('~')+'/programs/line_fitting/extract_spire.pro')
    # spectra extraction
    idl.pro('get_spire', outdir=outdir+obj+'/spire/data/cube/', object=obj,
            filename=indir+obsid+'_spectrum_extended_HR_aNB_15.fits', brightness=1)
    # line fitting
    # get coordinates
    # not sure the variable "ra_slw" etc will work via python
    idl('.r '+os.path.expanduser('~')+'/programs/line_fitting/get_radec_spire.pro')
    # not going to work in this way
    idl.pro('get_radec_spire', filename=indir+obsid+'_spectrum_extended_HR_aNB_15.fits',
            pix_slw=pix_slw, ra_slw=ra_slw, dec_slw=dec_slw, slw=1)
    idl.pro('get_radec_spire', filename=indir+obsid+'_spectrum_extended_HR_aNB_15.fits',
            pix_ssw=pix_ssw, ra_ssw=ra_ssw, dec_ssw=dec_ssw, ssw=1)

    idl.pro('extract_spire', indir=outdir+obj+'/spire/data/cube/',
            outdir=outdir+obj+'/spire/advanced_products/cube/',
            plotdir=outdir+obj+'/spire/advanced_products/cube/plots/',
            localbaseline=10, global_noise=20, ra=ra_slw, dec=dec_slw, coordpix=pix_slw,
            slw=1, noiselevel=3, brightness=1, object=obj, flat=1, continuum_sub=1,
            current_pix=1, double_gauss=1)
    idl.pro('extract_spire', indir=outdir+obj+'/spire/data/cube/',
            outdir=outdir+obj+'/spire/advanced_products/cube/',
            plotdir=outdir+obj+'/spire/advanced_products/cube/plots/',
            localbaseline=10, global_noise=20, ra=ra_ssw, dec=dec_ssw, coordpix=pix_ssw,
            ssw=1, noiselevel=3, brightness=1, object=obj, flat=1, continuum_sub=1,
            current_pix=1, double_gauss=1)


Obsid=[1342242620,1342242621,1342245084,1342245094,1342245857,
       1342247625,1342248246,1342248249,1342249053,1342249470,
       1342249474,1342249475,1342249476,1342249477,
       1342250509,1342250510,1342250512,1342250515,1342251285,
       1342251286,1342251287,1342251290,1342253646,1342253649,
       1342253652,1342254037]

obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301',
            'TMR1','TMC1A','TMC1','IRAS15398','RNO91',
            'GSS30-IRS1','VLA1623','WL12','RCrA-IRS5A','L483',
            'B335','DKCha']

indir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/'
outdir = '/home/bettyjo/yaolun/CDF_SPIRE_reduction/'
for obsid in Obsid:
    spire_postprocess1(indir, outdir, str(obsid), obj_list[Obsid.index(obsid)])
    spire_cubeprocess(indir, outdir, str(obsid), obj_list[Obsid.index(obsid)])

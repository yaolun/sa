def process_phot(pacsindir, spireindir, outdir, obj):
    """
    compile photometry from PACS and SPIRE
    """
    from astropy.io import ascii
    import numpy as np
    import os

    foo = open(outdir+obj.lower()+'.txt', 'w')
    foo.write('wavelength \t flux(Jy) \t error(Jy) \t aperture(arcsec) \n')

    # read in PACS photometry
    if os.path.exists(pacsindir+obj+'_phot_sect.txt'):
        pacsdata = ascii.read(pacsindir+obj+'_phot_sect.txt', data_start=4)

        # construct non-redundant wavelength points
        pacswave = np.sort(list(set(pacsdata['wavelength(um)'].data)))
        pacsphot = np.empty_like(pacswave)
        pacsphot_err = np.empty_like(pacswave)

        # calculate the mean flux and uncertainty for each wavelength point
        for i in range(len(pacswave)):
            pacsphot[i] = np.mean(pacsdata['flux(Jy)'][pacsdata['wavelength(um)'] == pacswave[i]])
            pacsphot_err[i] = np.sqrt(np.sum(pacsdata['uncertainty(Jy)'][pacsdata['wavelength(um)'] == pacswave[i]]**2))/\
                          len(pacsdata['wavelength(um)'][pacsdata['wavelength(um)'] == pacswave[i]])
        for ip in range(len(pacswave)):
            foo.write('%f \t %f \t %f \t 31.8 \n' % (pacswave[ip], pacsphot[ip], pacsphot_err[ip]))

    # read in SPIRE photometry
    if os.path.exists(spireindir+obj+'_spire_phot.txt'):
        spiredata = ascii.read(spireindir+obj+'_spire_phot.txt', data_start=4)

        # construct non-redundant wavelength points
        spirewave = np.sort(list(set(spiredata['wavelength(um)'].data)))
        spirephot = np.empty_like(spirewave)
        spirephot_err = np.empty_like(spirewave)
        spirephot_aper = np.empty_like(spirewave)

        # calculate the mean flux and uncertainty for each wavelength point
        for i in range(len(spirewave)):
            spirephot[i] = np.mean(spiredata['flux(Jy)'][spiredata['wavelength(um)'] == spirewave[i]])
            spirephot_err[i] = np.sqrt(np.sum(spiredata['uncertainty(Jy)'][spiredata['wavelength(um)'] == spirewave[i]]**2))/\
                          len(spiredata['wavelength(um)'][spiredata['wavelength(um)'] == spirewave[i]])
            spirephot_aper[i] = np.mean(spiredata['aperture(arcsec)'][spiredata['wavelength(um)'] == spirewave[i]])
        for i in range(len(spirewave)):
            foo.write('%f \t %f \t %f \t %f \n' % (spirewave[i], spirephot[i], spirephot_err[i], spirephot_aper[i]))
    foo.close()

pacsindir = '/Volumes/SD-Mac/pacsphot_cdf/'
spireindir = '/Volumes/SD-Mac/CDF_archive/photometry/'
obj_list = ['RCrA-IRS7B','RCrA-IRS7C','HH46','L723-MM','L1014',
            'L1157','Ced110','BHR71','IRAS03245','L1551-IRS5',
            'L1455-IRS3','B1-a','B1-c','IRAS03301',
            'TMR1','TMC1A','TMC1','IRAS15398','RNO91',
            'GSS30-IRS1','VLA1623','WL12','RCrA-IRS5A','L483',
            'B335','DKCha']
outdir = '/Volumes/SD-Mac/CDF_archive/'
for o in obj_list:
    process_phot(pacsindir, spireindir, outdir, o)

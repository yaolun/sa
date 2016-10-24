def spire_postprocess1(indir,outdir,obsid,obj):
    from astropy.io import ascii
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
    foo = open(outdir+obj+'_spire_corrected.txt','w')
    foo.write('%s \t %s \n' % ('Wavelength(um)', 'Flux_Density(Jy)'))
    for i in range(len(spire_wl)):
        foo.write('%f \t %f \n' % (spire_wl[i], spire_flux[i]))
    foo.close()

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

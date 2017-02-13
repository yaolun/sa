def pacs_weight(cubedir, obj, aper_size, outdir, fits_for_header,
                suffix='', plot=None, wavelimit=None):
    """
    cubedir: The directory contains the flux and coordinates of each spaxel in ASCII files
    plot: If a path is given, then it will make a plot of an overview of aperture_masked
          spaxel configuration.
    fits_for_header: The FITS file contains the required keywords -
                     CROTA2, CDELT1 or CDELT2
    """
    from astropy.io import ascii, fits
    import numpy as np

    def CircularAperture(x, y, r, pix=10000.):

        grid_x, grid_y = np.meshgrid(np.linspace(0,pix-1,pix), np.linspace(0,pix-1,pix))
        grid_x = grid_x - (pix-1)/2.
        grid_y = grid_y - (pix-1)/2.
        grid_dist = ((grid_x-x)**2+(grid_y-y)**2)**0.5
        aperture = np.where(grid_dist <= r, np.full((int(pix), int(pix)), 1, dtype=int), 0)

        return aperture

    def Total_MaskedAperture(aperture, mask):

        aperture_masked = np.where(mask != 0, aperture, 0)

        return aperture_masked, np.sum(aperture_masked)

    def Mask(x, y, edge, rot, pix=10000.):
        from scipy.ndimage.interpolation import rotate

        init_arr = np.full((pix, pix), 1, dtype=int)
        grid_x, grid_y = np.meshgrid(np.linspace(0,pix-1,pix), np.linspace(0,pix-1,pix))

        # RA offset (x) and Dec offset (y) need to be rotated, because the spaxels are rotated.
        rotMatrix = np.array([[np.cos(np.radians(rot)), -np.sin(np.radians(rot))],
                              [np.sin(np.radians(rot)),  np.cos(np.radians(rot))]])
        coord_rot = np.squeeze(np.dot(rotMatrix, np.array([x,y]).reshape(2,1)))

        dist_x = abs(grid_x - (pix-1)/2. - coord_rot[0])
        dist_y = abs(grid_y - (pix-1)/2. - coord_rot[1])

        aper_trim = np.where((dist_x <= edge/2.) & (dist_y <= edge/2.), init_arr, 0)

        return aper_trim

    mean_coord = {'RA':[], 'Dec':[]}
    coord = {'RA':[], 'Dec':[]}
    foo_cen = ascii.read(cubedir+obj+'_pacs_pixel13_'+suffix+'_coord.txt')
    cen_ra, cen_dec = foo_cen['RA(deg)'][0], foo_cen['Dec(deg)'][0]
    mean_cen_ra, mean_cen_dec = np.mean(foo_cen['RA(deg)']), np.mean(foo_cen['Dec(deg)'])

    for i in range(1,26):
        foo = ascii.read(cubedir+obj+'_pacs_pixel'+str(i)+'_'+suffix+'_coord.txt')
        mean_coord['RA'].append((np.mean(foo['RA(deg)']) - mean_cen_ra)*np.cos(np.radians(mean_cen_dec))*3600.)
        mean_coord['Dec'].append((np.mean(foo['Dec(deg)']) - mean_cen_dec)*3600.)

        coord['RA'].append((foo['RA(deg)'][0]-cen_ra)*np.cos(np.radians(cen_dec))*3600.)
        coord['Dec'].append((foo['Dec(deg)'][0]-cen_dec)*3600.)

    # do 60" x 60"
    pix = 1000.
    factor = pix/60.
    weight = []
    aperture = CircularAperture(0, 0, aper_size/2.*factor, pix=pix)
    pix_size = abs(fits.open(fits_for_header)[1].header['CDELT1']*3600)
    ideal = (pix_size*factor)**2

    rot_angle = fits.open(fits_for_header)[1].header['CROTA2']

    a = np.full((int(pix), int(pix)), 0)
    for i in range(0,25):
        ra = coord['RA'][i]
        dec = coord['Dec'][i]
        masked_arr, arr_sum = Total_MaskedAperture(aperture, Mask(ra*factor, dec*factor, pix_size*factor, rot_angle, pix=pix))
        a = a + masked_arr
        weight.append(arr_sum)

    # somehow there is about 1% error even if the whole spaxel should be within the circular aperture
    for i in range(len(weight)):
        if np.abs(weight[i]/ideal) >= 0.99:
            weight[i] = 1
        else:
            weight[i] = weight[i]/ideal

    # plot
    if plot != None:
        import matplotlib.pyplot as plt
        from scipy.ndimage.interpolation import rotate

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        factor_rot = len(rotate(a, -rot_angle)[:,0])/pix

        ax.imshow(rotate(a, -rot_angle), extent=[-30*factor_rot,30*factor_rot,-30*factor_rot,30*factor_rot],
                  cmap='viridis', vmax=2)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on()
        ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
        ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)
        ax.set_xlabel('RA offset [arcsec]', fontsize=18)
        ax.set_ylabel('Dec offset [arcsec]', fontsize=18)

        fig.savefig(plot+'pacs_masked.pdf', format='pdf', dpo=300, bbox_inches='tight')

    # get the flux of each spaxel and apply the weighting
    foo_cen = ascii.read(cubedir+obj+'_pacs_pixel13_'+suffix+'.txt')
    if wavelimit == None:
        trimmer = np.ones_like(foo_cen['Wavelength(um)'], dtype=bool)
    else:
        trimmer = (foo_cen['Wavelength(um)'] > wavelimit[0]) & (foo_cen['Wavelength(um)'] < wavelimit[1])

    wl = foo_cen['Wavelength(um)'][trimmer]
    flux = np.zeros_like(foo_cen['Flux_Density(Jy)'][trimmer])

    for i in range(1,26):
        foo = ascii.read(cubedir+obj+'_pacs_pixel'+str(i)+'_'+suffix+'.txt')
        # set NaN values to zero
        foo['Flux_Density(Jy)'][np.isnan(foo['Flux_Density(Jy)'])] = 0

        flux = flux + foo['Flux_Density(Jy)'][trimmer]*weight[i-1]

    # write the weighted spectrum into a file
    foo = open(outdir+obj+'_pacs_weighted_'+suffix+'.txt','w')
    foo.write('{} \t {}\n'.format('Wavelength(um)', 'Flux_Density(Jy)'))
    for i in range(len(wl)):
        if flux[i] != 0:
            foo.write('{} \t {}\n'.format(wl[i], flux[i]))
    foo.close()

    print 'Weighted spectrum saved at ', outdir+obj+'_pacs_weighted_'+suffix+'.txt'

    return wl, flux

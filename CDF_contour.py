def CDF_contour(linename, objname, fitting_table,
                plot=True, cont=False, pix_cen=None, output=None, print_obj=False,
                sigma_floor=0, nofits=False, ggd37=False):

    from scipy.interpolate import griddata
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from astropy.io import ascii, fits
    import numpy as np
    import matplotlib.pyplot as plt

    # get the line data
    fitting = ascii.read(fitting_table)
    data = fitting[(fitting['Object'] == objname) & (fitting['Line'] == linename) & (fitting['Pixel_No.'] != 'c')]

    # determine the name of the central spaxel
    if pix_cen == None:
        pix_cen_list = ['SLWC3', 'SSWD4','13']
        pix_cen = data['Pixel_No.'][(data['Pixel_No.'] == pix_cen_list[0])+\
                                    (data['Pixel_No.'] == pix_cen_list[1])+\
                                    (data['Pixel_No.'] == pix_cen_list[2])]

    # get the RA & Dec of the central spaxel
    ra_cen = data['RA(deg)'][data['Pixel_No.'] == pix_cen].data
    dec_cen = data['Dec(deg)'][data['Pixel_No.'] == pix_cen].data

    # calculate the RA and Dec separations
    plot_ra = (data['RA(deg)'].data-ra_cen)*np.cos(np.radians(dec_cen))*3600
    plot_dec = (data['Dec(deg)'].data-dec_cen)*3600

    # determine the size of the contour plot
    size = np.ceil(max(abs(plot_ra).max(), abs(plot_dec).max())/10)*10

    # round up the ranges of RA and Dec separations, and calculate the number of points in between
    ra_range = [np.ceil(plot_ra.max()/10)*10,
                np.ceil(plot_ra.min()/10)*10,
                np.ceil((plot_ra.max()-plot_ra.min())/10)*10]
    dec_range = [np.ceil(plot_dec.max()/10)*10,
                 np.ceil(plot_dec.min()/10)*10,
                 np.ceil((plot_dec.max()-plot_dec.min())/10)*10]

    # create the rebinned grid for RA and Dec.  Use oversample of 4.
    ra_rebin = np.linspace(ra_range[0], ra_range[1], ra_range[2]*4)
    dec_rebin = np.linspace(dec_range[1], dec_range[0], dec_range[2]*4)

    # use the rebinned RA and Dec to re-grid the contours, both line and continuum
    if ggd37:
        selector = (plot_dec <= 0.)
        plot_ra = plot_ra[selector]
        plot_dec = plot_dec[selector]
        data = data[selector]

    z = griddata((plot_ra, plot_dec), data['Str(W/cm2)'].data, (ra_rebin[None,:], dec_rebin[:,None]), method='cubic')
    z_cont = griddata((plot_ra, plot_dec), data['Base(W/cm2/um)'].data*data['FWHM(um)'].data*1.086,
                      (ra_rebin[None,:], dec_rebin[:,None]), method='cubic')

    # calculate the noise floor for the line emission
    sigma = np.nanmin(data['Str(W/cm2)']/data['SNR'].data)
    if sigma_floor != 0:
        z_floor = sigma*sigma_floor
    else:
        z_floor = np.nanmin(z)

    # create the figure and axis objects
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot the contour with color and lines
    levels = np.linspace(z_floor, np.nanmax(z), 10)[1:]
    ax.contour(ra_rebin, dec_rebin, z, levels, linewidths=1.5, cmap='Reds')

    # whether show the continuum as image or line emission as image
    if cont:
        im = ax.imshow(z_cont, cmap='Blues', origin='lower',
                       extent=[ra_range[0], ra_range[1],dec_range[1],dec_range[0]])
        im_label = 'F_{cont.}'
    else:
        im = ax.imshow(z, cmap='Blues', origin='lower', extent=[ra_range[0], ra_range[1],dec_range[1],dec_range[0]])
        im_label = 'F_{line}'

    # set the bad pixel to white
    im.cmap.set_bad('w', 1.)
    # setup ticks and tick labels
    ax.set_xticks(np.linspace(-np.ceil(-ra_range[1]/10)*10, np.ceil(ra_range[0]/10)*10, 5, dtype='int'))
    ax.set_xticklabels(np.linspace(-np.ceil(-ra_range[1]/10)*10, np.ceil(ra_range[0]/10)*10, 5, dtype='int'))
    ax.set_yticks(np.linspace(-np.ceil(-dec_range[1]/10)*10 , np.ceil(dec_range[0]/10)*10 , 5, dtype='int'))
    ax.set_yticklabels(np.linspace(-np.ceil(-dec_range[1]/10)*10, np.ceil(dec_range[0]/10)*10, 5, dtype='int'))

    # create the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)
    cb.solids.set_edgecolor("face")
    cb.ax.minorticks_on()
    cb.ax.set_ylabel(r'$\rm'+im_label+'\,[W\,cm^{-2}]$',fontsize=16)
    cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cb_obj,fontsize=12)

    # set the x, y labels
    ax.set_xlabel(r'$\rm{RA\,offset\,[arcsec]}$', fontsize=18)
    ax.set_ylabel(r'$\rm{Dec\,offset\,[arcsec]}$', fontsize=18)
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax.minorticks_on()
    ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=5,length=5)
    ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=5,length=2.5)
    ax.set_aspect('equal', 'datalim')

    # print the object name
    if print_obj:
        ax.text(0.9, 0.9, objname, transform=ax.transAxes, fontsize=18, ha='right')

    # output the interpolated 2D array and the RA/Dec arrays into FITS and ASCII files
    if output != None:
        if not nofits:
            hdulist = fits.HDUList([fits.PrimaryHDU(z), fits.ImageHDU(z_cont)])
            hdulist.writeto(output, overwrite=True)

        # write out the RA/Dec arrays
        foo = open(output.split('.')[0]+'_interpolated_RA.txt', 'w')
        # the coordinates of the reference pixel
        foo.write('# Pixel 1: {:<12.8f} / {:<12.8f}\n'.format(ra_cen[0], dec_cen[0]))
        foo.write('{:<12s}\n'.format('RA_offset'))
        for i, ra_dum in enumerate(ra_rebin):
            foo.write('{:<12.8f}\n'.format(ra_dum))
        foo.close()

        foo = open(output.split('.')[0]+'_interpolated_Dec.txt', 'w')
        # the coordinates of the reference pixel
        foo.write('# Pixel 1: {:<12.8f} / {:<12.8f}\n'.format(ra_cen[0], dec_cen[0]))
        foo.write('{:<12s}\n'.format('Dec_offset'))
        for i, dec_dum in enumerate(dec_rebin):
            foo.write('{:<12.8f}\n'.format(dec_dum))
        foo.close()

        fig.savefig(output.split('.')[0]+'_contour.pdf', format='pdf', dpi=300, bbox_inches='tight')
        print('Figure saved at ', output.split('.')[0]+'_contour.pdf')

    return (ra_rebin, dec_rebin), (z, z_cont)

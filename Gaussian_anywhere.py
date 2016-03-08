def gauss2d(x, y, sigmax, sigmay=None):
    import numpy as np
    # sigmax needs to be in pixel coordinates
    if sigmay == None:
        sigmay = sigmax
    gauss = 1/(2*np.pi*sigmax*sigmay) * np.exp( -(x**2/2./sigmax**2 + y**2/2./sigmay**2) )
    return gauss
def Gaussian_anywhere(ra_offset, dec_offset, pixel_aper, gauss_beam=31.8, size=1201., phys_size=120.):
    # Create a 2-D cartian grid for mapping 2-D Gaussian
    import numpy as np
    # convert from aperture (diameter) to radius
    radius = pixel_aper/2.

    grid_x, grid_y = np.meshgrid(np.linspace(0,size-1,size), np.linspace(0,size-1,size))
    grid_x = grid_x - (size-1)/2.
    grid_y = grid_y - (size-1)/2.
    grid_gauss2d = gauss2d(grid_x,grid_y, sigmax=(size-1)/phys_size*(gauss_beam/2.354))
    dA = ((1/((size-1)/2.))*phys_size)**2

    # convert from physcial coordinates to pixel coordinates
    x = (ra_offset-phys_size/2.) * (size-1)/2./(phys_size/2.) + (size-1)/2.
    y = (dec_offset-phys_size/2.) * (size-1)/2./(phys_size/2.) + (size-1)/2.
    r_pix = radius * (size-1)/phys_size
    grid_dist = ((grid_x-x)**2+(grid_y-y)**2)**0.5
    gauss2d_mask = np.where(grid_dist<=r_pix, grid_gauss2d,0)

    return np.sum(gauss2d_mask)

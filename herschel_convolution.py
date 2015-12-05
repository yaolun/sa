def herschel_convolution(image, outdir, wave=None, plot=False, scale=1., dist=1.0, noise=None):
    """
    image: 2d array contains the image at certain wavelength.
    scale: The conversion factor that converts from an element to physical size in cm.
           The easiest way to calculate this is size of the object divides the total number 
           of pixels across the image.
    dist:  The distance to the object in parsec.  If not given, the default value is 1 pc.
    noise: The detection limit of the observation in unit of erg/s/cm2/Hz/sr
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import astropy.constants as const
    from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel

    # constant setup
    AU = const.au.cgs.value
    pc = const.pc.cgs.value

    beam_size = 9.4   # in arcsec
    beam_size = 1.0
    # calculate the image angular size
    w = len(image[0,:]) * scale / (dist*pc) / 2. * 180/np.pi*3600.

    kernel = Gaussian2DKernel(beam_size * dist * AU / scale)
    syn_image = convolve(image, kernel)

    if plot:
        # plot the image
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        cmap = plt.cm.CMRmap
        if noise == None:
            floor = -22
        else:
            floor = np.log10(noise)
        im = ax.imshow(np.log10(syn_image), vmin= floor, vmax= -12,
                  cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.solids.set_edgecolor("face")
        cb.ax.minorticks_on()
        cb.ax.set_ylabel(r'$\rm{log(I_{\nu})\,[erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}]}$',fontsize=14)
        cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cb_obj,fontsize=12)

        # fix the tick label font
        ticks_font = font_manager.FontProperties(family='STIXGeneral',size=18)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in cb.ax.get_yticklabels():
            label.set_fontproperties(ticks_font)

        ax.set_xlabel(r'$\rm{RA\,Offset\,(arcsec)}$', fontsize=18)
        ax.set_ylabel(r'$\rm{Dec\,Offset\,(arcsec)}$', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)

        fig.savefig(outdir+'syn_image_'+str(wave)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
        fig.clf()

    return syn_image

filename = '/Users/yaolun/bhr71/hyperion/cycle9/model34.rtout'
outdir = '/Users/yaolun/test/'
dist = 178.
wave = 500.

from hyperion.model import ModelOutput
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from mpl_toolkits.axes_grid1 import make_axes_locatable

# constant setup
pc = const.pc.cgs.value

m = ModelOutput(filename)
image = m.get_image(group=22, inclination=0, distance=dist * pc, units='MJy/sr')
# Find the closest wavelength
iwav = np.argmin(np.abs(wave - image.wav))

# Calculate the image width in arcseconds given the distance used above
# get the max radius
rmax = max(m.get_quantities().r_wall)
w = np.degrees(rmax / image.distance) * 3600.

# Image in the unit of MJy/sr
# Change it into erg/s/cm2/Hz/sr
factor = 1e-23*1e6
# avoid zero in log
# flip the image, because the setup of inclination is upside down
val = image.val[::-1, :, iwav] * factor + 1e-30

# plot the image
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

cmap = plt.cm.CMRmap
im = ax.imshow(np.log10(val), vmin= -22, vmax= -12,
          cmap=cmap, origin='lower', extent=[-w, w, -w, w], aspect=1)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{log(I_{\nu})\,[erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}]}$',fontsize=14)
cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
plt.setp(cb_obj,fontsize=12)

# fix the tick label font
ticks_font = font_manager.FontProperties(family='STIXGeneral',size=18)
for label in ax.get_xticklabels():
    label.set_fontproperties(ticks_font)
for label in ax.get_yticklabels():
    label.set_fontproperties(ticks_font)
for label in cb.ax.get_yticklabels():
    label.set_fontproperties(ticks_font)

ax.set_xlabel(r'$\rm{RA\,Offset\,(arcsec)}$', fontsize=18)
ax.set_ylabel(r'$\rm{Dec\,Offset\,(arcsec)}$', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=18)

fig.savefig(outdir+'image_'+str(wave)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()


syn_im = herschel_convolution(val, outdir, wave=wave, plot=True, scale=2*rmax/len(val[0,:]), dist=dist, noise=1e-19)
print np.shape(syn_im)
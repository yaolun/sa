# OI velocity shift analysis for L1551-IRS5
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy.io import ascii
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

c = const.c.cgs.value

def read_fitting_oi(filepath,noiselevel):
    from astropy.io import ascii
    line_label = []
    data = ascii.read(filepath)
    header = data.colnames
    data = data[(data['SNR']>noiselevel) & (data['Str(W/cm2)']>0)]

    ind_line = []
    for i in range(0, len(data['Line'])):
        if 'OI3P1-3P2' in data['Line'][i]:
            #wl, wl_sig, flux, flux_sig, E_u, A, g
            ind_line.append(i)

    if len(ind_line) > 0:
    	line_data = data[ind_line[0]]
    	return line_data, line_data['Line']
    else:
    	return [],[]

indir = '/Users/yaolun/test/L1551-IRS5/cube/'
oi_ref_wl = 63.18367004
ra_cen = 67.89196135
dec_cen = 18.13468391
ra = []
dec = []
oi_wl = []

# read in [OI] line from fitting
for i in range(1,26):
	oi_dum, oi_name_dum = read_fitting_oi(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt', 3)
	if len(oi_name_dum) > 0:
		oi_wl.append(oi_dum[2])
		ra.append(oi_dum[14])
		dec.append(oi_dum[15])
	else:
		data = ascii.read(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt')
		data = data[data['Line'] == 'OI3P1-3P2']
		oi_wl.append(oi_ref_wl)
		ra.append(data['RA(deg)'].data[0])
		dec.append(data['Dec(deg)'].data[0])

oi_wl = np.array(oi_wl)
ra = np.array(ra)
dec = np.array(dec)

velo = (oi_wl-oi_ref_wl)/oi_ref_wl*c / 1e5 # unit km/s

# projection and convert to arcsec
plot_ra = (ra-ra_cen)*np.cos(np.radians(dec_cen))*3600
plot_dec = (dec-dec_cen)*3600

# for i in range(len(velo)):
# 	print plot_ra[i], plot_dec[i], velo[i]

# make the plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# ra_rebin = np.linspace(min(plot_ra), max(plot_ra), 50)
# dec_rebin = np.linspace(min(plot_dec), max(plot_dec), 50)
ra_rebin = np.linspace(20, -20, 50)
dec_rebin = np.linspace(-20, 20, 50)

z = griddata((plot_ra, plot_dec), velo, (ra_rebin[None:], dec_rebin[:,None]), method='cubic')

# plot the contour with color and lines
ax.contour(ra_rebin, dec_rebin, z, 15, linewidths=0.5,colors='k')
im = ax.imshow(z, cmap='jet', origin='lower', extent=[20,-20,-20,20])#,\
    # norm=LogNorm(vmin=velo.min(), vmax=velo.max()))
# Blues_r
ax.set_xticks(np.linspace(-20, 20, 5))
ax.set_xticklabels(np.linspace(-20, 20, 5))
ax.set_yticks(np.linspace(-19.5, 19.5, 5))
ax.set_yticklabels(np.linspace(-20, 20, 5))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{\delta}v\,[km\,s^{-1}]$',fontsize=16)
cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
plt.setp(cb_obj,fontsize=12)

# ori_data = ax.plot(plot_ra,plot_dec, 'bo', markersize=5)
# plt.gca().invert_xaxis()

ax.set_xlabel(r'$\rm{RA\,offset\,[arcsec]}$', fontsize=20)
ax.set_ylabel(r'$\rm{Dec\,offset\,[arcsec]}$', fontsize=20)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on() 
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

fig.savefig('/Users/yaolun/test/L1551-IRS5_velo_shift.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()
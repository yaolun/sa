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
def read_fitting_cii(filepath,noiselevel):
    from astropy.io import ascii
    line_label = []
    data = ascii.read(filepath)
    header = data.colnames
    data = data[(data['SNR']>noiselevel) & (data['Str(W/cm2)']>0)]

    ind_line = []
    for i in range(0, len(data['Line'])):
        if 'CII2P3_2-2P1_2' in data['Line'][i]:
            #wl, wl_sig, flux, flux_sig, E_u, A, g
            ind_line.append(i)

    if len(ind_line) > 0:
        line_data = data[ind_line[0]]
        return line_data, line_data['Line']
    else:
        return [],[]

def great_config(array_name, ax, rot=0, shift=(0,0), color='r', alpha=1, linestyle='--'):
    import numpy as np
    import matplotlib.pyplot as plt
    # configuration detail from Risacher+2016 (IEEE), Table VI
    if array_name == 'HFA':
        beam_dia = 6.3
        spacing = 13.8
    if array_name == 'LFA':
        beam_dia = 15.5 # at 1.9THz, 11.8 at 2.5THz
        spacing = 34.0

    array = [np.array([0,0]),
             np.array([spacing, 0]),
             np.array([-spacing,0]),
             np.array([spacing*np.cos(np.radians(60.0)), spacing*np.sin(np.radians(60.0))]),
             np.array([-spacing*np.cos(np.radians(60.0)), spacing*np.sin(np.radians(60.0))]),
             np.array([spacing*np.cos(np.radians(60.0)), -spacing*np.sin(np.radians(60.0))]),
             np.array([-spacing*np.cos(np.radians(60.0)), -spacing*np.sin(np.radians(60.0))])]

    rot_matrix = np.matrix([[np.cos(np.radians(rot)), -np.sin(np.radians(rot))],
                     [np.sin(np.radians(rot)), np.cos(np.radians(rot))]])

    for i in range(len(array)):
        array[i] = np.array([array[i][0]+shift[0], array[i][1]+shift[1]])
        array[i] = np.dot(array[i], rot_matrix)[0]

    for pix in array:
        ax.add_artist(plt.Circle((pix[0,0],pix[0,1]), beam_dia/2, color=color,
                                 fill=False, linestyle=linestyle, linewidth=2, alpha=alpha))

    return None


indir = '/Users/yaolun/test/L1551-IRS5/cube/'
oi_ref_wl = 63.18367004
cii_ref_wl = 157.6922760
ra_cen = 67.89196135
dec_cen = 18.13468391

# OI part
ra = []
dec = []
oi_wl = []
oi_width = []

# read in [OI] line from fitting
for i in range(1,26):
    oi_dum, oi_name_dum = read_fitting_oi(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt', 3)
    if len(oi_name_dum) > 0:
        oi_wl.append(oi_dum[2])
        ra.append(oi_dum[14])
        dec.append(oi_dum[15])
        oi_width.append(max([oi_dum[6]/oi_dum[2]*c/1e5, 186]))
    else:
        data = ascii.read(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt')
        data = data[data['Line'] == 'OI3P1-3P2']
        oi_wl.append(oi_ref_wl)
        ra.append(data['RA(deg)'].data[0])
        dec.append(data['Dec(deg)'].data[0])

print max(oi_width), min(oi_width), np.mean(oi_width)

oi_wl = np.array(oi_wl)
ra = np.array(ra)
dec = np.array(dec)

velo = (oi_wl-oi_ref_wl)/oi_ref_wl*c / 1e5 # unit km/s

# projection and convert to arcsec
plot_ra = (ra-ra_cen)*np.cos(np.radians(dec_cen))*3600
plot_dec = (dec-dec_cen)*3600

# for i in range(len(velo)):
#   print plot_ra[i], plot_dec[i], velo[i]

# make the plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# ra_rebin = np.linspace(min(plot_ra), max(plot_ra), 100)
# dec_rebin = np.linspace(min(plot_dec), max(plot_dec), 100)
ra_rebin = np.linspace(20, -20, 100)
dec_rebin = np.linspace(-20, 20, 100)

z = griddata((plot_ra, plot_dec), velo, (ra_rebin[None,:], dec_rebin[:,None]), method='cubic')

# masked array for masking NaN values
masked_z = np.ma.array (z, mask=np.isnan(z))

# add the GREAT pixel configuration
great_config('HFA', ax, rot=-33)
# shift along the outflow direction
great_config('HFA', ax, rot=-33, shift=(6.9,0), color='m', linestyle='-')
great_config('HFA', ax, rot=-33, shift=(-6.9,0), color='cyan')

# plot the contour with color and lines
ax.contour(ra_rebin, dec_rebin, z, 10, linewidths=0.5,colors='k')
im = ax.imshow(masked_z, cmap='jet', origin='lower', extent=[20,-20,-20,20])#,\
    # norm=LogNorm(vmin=velo.min(), vmax=velo.max()))
im.cmap.set_bad('w', 1.)
# Blues_r
ax.set_xticks(np.linspace(-20, 20, 5))
ax.set_xticklabels(np.linspace(-20, 20, 5))
ax.set_yticks(np.linspace(-20 , 20 , 5))
ax.set_yticklabels(np.linspace(-20, 20, 5))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{\Delta}v\,[km\,s^{-1}]$',fontsize=16)
cb_obj = plt.getp(cb.ax.axes, 'yticklabels')
plt.setp(cb_obj,fontsize=12)

ori_data = ax.plot(plot_ra,plot_dec, 'bo', markersize=5)
# plt.gca().invert_xaxis()

ax.set_xlabel(r'$\rm{RA\,offset\,[arcsec]}$', fontsize=20)
ax.set_ylabel(r'$\rm{Dec\,offset\,[arcsec]}$', fontsize=20)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on()
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

fig.savefig('/Users/yaolun/test/L1551-IRS5_oi_velo_shift.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

# CII part
ra = []
dec = []
cii_wl = []
cii_width = []

# read in [OI] line from fitting
for i in range(1,26):
    cii_dum, cii_name_dum = read_fitting_cii(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt', 3)
    if len(cii_name_dum) > 0:
        cii_wl.append(cii_dum[2])
        ra.append(cii_dum[14])
        dec.append(cii_dum[15])
        cii_width.append(cii_dum[6]/cii_dum[2] * c/1e5)
    else:
        data = ascii.read(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt')
        data = data[data['Line'] == 'CII2P3_2-2P1_2']
        cii_wl.append(cii_ref_wl)
        ra.append(data['RA(deg)'].data[0])
        dec.append(data['Dec(deg)'].data[0])

cii_wl = np.array(cii_wl)
ra = np.array(ra)
dec = np.array(dec)

print max(cii_width), min(cii_width), np.mean(cii_width)

velo = (cii_wl-cii_ref_wl)/cii_ref_wl*c / 1e5 # unit km/s

# projection and convert to arcsec
plot_ra = (ra-ra_cen)*np.cos(np.radians(dec_cen))*3600
plot_dec = (dec-dec_cen)*3600

# for i in range(len(velo)):
#   print plot_ra[i], plot_dec[i], velo[i]

# make the plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# ra_rebin = np.linspace(min(plot_ra), max(plot_ra), 100)
# dec_rebin = np.linspace(min(plot_dec), max(plot_dec), 100)
ra_rebin = np.linspace(20, -20, 100)
dec_rebin = np.linspace(-20, 20, 100)

z = griddata((plot_ra, plot_dec), velo, (ra_rebin[None,:], dec_rebin[:,None]), method='cubic')

# masked array for masking NaN values
masked_z = np.ma.array (z, mask=np.isnan(z))

# add the GREAT pixel configuration
great_config('LFA', ax, rot=-33)
# shift along the outflow direction
great_config('LFA', ax, rot=-33, shift=(6.9,0), color='m', linestyle='-')
great_config('LFA', ax, rot=-33, shift=(-6.9,0), color='cyan')

# plot the contour with color and lines
ax.contour(ra_rebin, dec_rebin, z, 10, linewidths=0.5,colors='k')
im = ax.imshow(z, cmap='jet', origin='lower', extent=[20,-20,-20,20])#,\
    # norm=LogNorm(vmin=velo.min(), vmax=velo.max()))
im.cmap.set_bad('w', 1.)
# Blues_r
ax.set_xticks(np.linspace(-20, 20, 5))
ax.set_xticklabels(np.linspace(-20, 20, 5))
ax.set_yticks(np.linspace(-20, 20, 5))
ax.set_yticklabels(np.linspace(-20, 20, 5))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{\Delta}v\,[km\,s^{-1}]$',fontsize=16)
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

fig.savefig('/Users/yaolun/test/L1551-IRS5_cii_velo_shift.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

# OI intensity map
ra = []
dec = []
oi_flux = []

# read in [OI] line from fitting
for i in range(1,26):
    oi_dum, oi_name_dum = read_fitting_oi(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt', 3)
    if len(oi_name_dum) > 0:
        oi_flux.append(oi_dum[4])
        ra.append(oi_dum[14])
        dec.append(oi_dum[15])
    else:
        data = ascii.read(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt')
        data = data[data['Line'] == 'OI3P1-3P2']
        oi_flux.append(0)
        ra.append(data['RA(deg)'].data[0])
        dec.append(data['Dec(deg)'].data[0])

oi_flux = np.array(oi_flux) * 1e7 / ((9.4/2)**2*np.pi)
ra = np.array(ra)
dec = np.array(dec)

# projection and convert to arcsec
plot_ra = (ra-ra_cen)*np.cos(np.radians(dec_cen))*3600
plot_dec = (dec-dec_cen)*3600

# make the plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# ra_rebin = np.linspace(min(plot_ra), max(plot_ra), 100)
# dec_rebin = np.linspace(min(plot_dec), max(plot_dec), 100)
ra_rebin = np.linspace(20, -20, 100)
dec_rebin = np.linspace(-20, 20, 100)

z = griddata((plot_ra, plot_dec), oi_flux, (ra_rebin[None,:], dec_rebin[:,None]), method='cubic') / 1e-14

# masked array for masking NaN values
masked_z = np.ma.array (z, mask=np.isnan(z))

# print out the averaged intensity per beam at specific position
pos = [[0,0], [5.78682692, 3.75800934], [-5.78682692, -3.75800934]]
beam_size = 6.3

for p in pos:
    total_I = 0
    pix_num = 0
    for rr in range(len(ra_rebin)):
        for dd in range(len(dec_rebin)):
            if (ra_rebin[rr]-p[0])**2+(dec_rebin[dd]-p[1])**2 <= (beam_size/2)**2:
                # the z array is inversed
                total_I = total_I+z[dd,rr]
                pix_num += 1
    print total_I/pix_num

# plot the contour with color and lines
levels = np.linspace(np.nanmin(z), np.nanmax(z), 10)[1:]
ax.contour(ra_rebin, dec_rebin, z, levels, linewidths=1,colors='k')
im = ax.imshow(z, cmap='Blues', origin='lower', extent=[20,-20,-20,20])#,\
    # norm=LogNorm(vmin=velo.min(), vmax=velo.max()))
im.cmap.set_bad('w', 1.)
# Blues_r
ax.set_xticks(np.linspace(-20, 20, 5))
ax.set_xticklabels(np.linspace(-20, 20, 5))
ax.set_yticks(np.linspace(-20 , 20 , 5))
ax.set_yticklabels(np.linspace(-20, 20, 5))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{F_{[OI]}\,[10^{-14}\,erg\,s^{-1}\,cm^{-2}\,arcsec^{-2}]}$',fontsize=16)
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

# add the GREAT pixel configuration
great_config('HFA', ax, rot=-33)
# shift along the outflow direction
great_config('HFA', ax, rot=-33, shift=(6.9,0), color='m', linestyle='-')
great_config('HFA', ax, rot=-33, shift=(-6.9,0), color='cyan')

fig.savefig('/Users/yaolun/test/L1551-IRS5_oi_intensity.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

# CII intensity map
ra = []
dec = []
cii_flux = []

# read in [OI] line from fitting
for i in range(1,26):
    cii_dum, cii_name_dum = read_fitting_cii(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt', 3)
    if len(cii_name_dum) > 0:
        cii_flux.append(cii_dum[4])
        ra.append(cii_dum[14])
        dec.append(cii_dum[15])
    else:
        data = ascii.read(indir+'L1551-IRS5_pacs_pixel'+str(i)+'_os8_sf7_lines.txt')
        data = data[data['Line'] == 'CII2P3_2-2P1_2']
        cii_flux.append(0)
        ra.append(data['RA(deg)'].data[0])
        dec.append(data['Dec(deg)'].data[0])

cii_flux = np.array(cii_flux) * 1e7 / ((9.4/2)**2*np.pi)
ra = np.array(ra)
dec = np.array(dec)

# projection and convert to arcsec
plot_ra = (ra-ra_cen)*np.cos(np.radians(dec_cen))*3600
plot_dec = (dec-dec_cen)*3600

# make the plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# ra_rebin = np.linspace(min(plot_ra), max(plot_ra), 100)
# dec_rebin = np.linspace(min(plot_dec), max(plot_dec), 100)
ra_rebin = np.linspace(20, -20, 100)
dec_rebin = np.linspace(-20, 20, 102)

z = griddata((plot_ra, plot_dec), cii_flux, (ra_rebin[None,:], dec_rebin[:,None]), method='cubic') / 1e-16

# masked array for masking NaN values
masked_z = np.ma.array (z, mask=np.isnan(z))

# print out the averaged intensity per beam at specific position
pos = [[0,0], [5.78682692, 3.75800934], [-5.78682692, -3.75800934]]
beam_size = 15.5


for p in pos:
    total_I = 0
    pix_num = 0
    for rr in range(len(ra_rebin)):
        for dd in range(len(dec_rebin)):
            if (ra_rebin[rr]-p[0])**2+(dec_rebin[dd]-p[1])**2 <= (beam_size/2)**2:
                # the z array is inversed
                total_I = total_I+z[dd,rr]
                pix_num += 1
    print total_I/pix_num


# plot the contour with color and lines
levels = np.linspace(np.nanmin(z), np.nanmax(z), 10)[1:]
ax.contour(ra_rebin, dec_rebin, z, levels, linewidths=1,colors='k')
im = ax.imshow(z, cmap='Blues', origin='lower', extent=[20,-20,-20,20])#,\
    # norm=LogNorm(vmin=velo.min(), vmax=velo.max()))
im.cmap.set_bad('w', 1.)
# Blues_r
ax.set_xticks(np.linspace(-20, 20, 5))
ax.set_xticklabels(np.linspace(-20, 20, 5))
ax.set_yticks(np.linspace(-20 , 20 , 5))
ax.set_yticklabels(np.linspace(-20, 20, 5))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.solids.set_edgecolor("face")
cb.ax.minorticks_on()
cb.ax.set_ylabel(r'$\rm{F_{[CII]}\,[10^{-16}\,erg\,s^{-1}\,cm^{-2}\,arcsec^{-2}]}$',fontsize=16)
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

# add the GREAT pixel configuration
great_config('LFA', ax, rot=-33)
great_config('LFA', ax, rot=-33, shift=(6.9,0), color='m', linestyle='-')
great_config('LFA', ax, rot=-33, shift=(-6.9,0), color='cyan')


fig.savefig('/Users/yaolun/test/L1551-IRS5_cii_intensity.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import os
from read_fitting import *
home = os.path.expanduser('~')

c = const.c.cgs.value
k = 1.380658e-16
obj = 'L1551-IRS5'
indir = '/data/digit_v65/'+obj+'/cube/'
lab = 63.1836709
obs = []
sig = []
coord = []
all_coord = []
velocity = []
all_pixel = []
Tr = []
omega_b = np.pi/4/np.log(2)*(9.4/3600./180*np.pi)**2
print omega_b
for pix in range(1,26):
	[data, header, name] = read_fitting(indir+obj+'_pacs_pixel'+str(pix)+'_os8_sf7_lines_localbaseline_fixwidth_global_noise.txt',0)
	#Plot the wavelength offest
	#wavelength version
	wl_lab = data[:,0].astype('float').T
	wl_obs = data[:,1].astype('float').T
	wl_obs_sig = data[:,2].astype('float').T
	flux = data[:,3].astype('float').T
	fwhm = data[:,5].astype('float').T
	snr_all = data[:,8].astype('float').T
	ra_dum = data[:,12].astype('float').T
	dec_dum = data[:,13].astype('float').T

	ind = np.nonzero((name=='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0))
	wl_lab = wl_lab[ind]
	wl_obs = wl_obs[ind]
	wl_obs_sig = wl_obs_sig[ind]
	flux = flux[ind]
	fwhm = fwhm[ind]
	snr_all = snr_all[ind]

	indd = np.nonzero((name=='OI3P1-3P2'))
	all_coord.append([float(ra_dum[indd])*np.cos(float(dec_dum[indd])*np.pi/180.),float(dec_dum[indd])])
	all_pixel.append(pix)
	if len(wl_obs) == 0:
		continue
	else:
		for iwl in range(0,len(wl_lab)):
			wl_lab[iwl] = 63.1836709
		obs.append(float(wl_obs))
		sig.append(float(wl_obs_sig))
		coord.append([float(ra_dum[ind])*np.cos(float(dec_dum[ind])*np.pi/180.),float(dec_dum[ind])])
		dwl = lab-wl_obs
		dv = -dwl/lab*c/1e5
		width = float(fwhm)/float(wl_obs)*c/1e5
		Tr.append((float(wl_obs)*1e-4)**3./2/k*float(flux)*1e7/omega_b/(width*1e5))
		velocity.append(dv)
	print 'Pixel: %2d, Velocity: %8f, Strength: %e' % (pix,dv,flux)
	print (float(wl_obs)*1e-4)**3./2/k*float(flux)*1e7/omega_b/(dv*1e5)
sig = np.array(sig)
obs = np.array(obs)
mean = np.sum((1/sig**2)*(lab-obs))/np.sum(1/sig**2)
# mean = np.sum(lab-obs)/len(obs)
sig = (1/np.sum(1/sig**2))**0.5

velo = mean/lab*c/1e5
sig_velo = sig/lab*c/1e5

print obj
print '%8f  %8f' % (velo, sig_velo)

# Plot the velocity shift spatiallly
# 
all_coord = np.array(all_coord)*3600
coord = np.array(coord)*3600
cen_ra = all_coord[12,0]
cen_dec = all_coord[12,1]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(cen_ra-all_coord[:,0],cen_dec-all_coord[:,1],'s',color='k')
ax.plot(cen_ra-coord[:,0],cen_dec-coord[:,1],'s',color='Red')

for i in range(0,len(all_coord[:,0])):
	ax.annotate('%d' % all_pixel[i],xy=(cen_ra-all_coord[i,0],cen_dec-all_coord[i,1]),xytext=(cen_ra-all_coord[i,0],cen_dec-all_coord[i,1]+1))
	# print 'Pixel %d, RA: %16f, Dec: %16f' % (i+1,all_coord[i,0]/3600,all_coord[i,1]/3600)

for i in range(0,len(velocity)):
	if velocity[i] < 0:
		color='Blue'
	elif velocity[i] > 0:
		color='Red'
	else:
		color='k'
	ax.annotate('%6.4f' % velocity[i],xy=(cen_ra-coord[i,0],cen_dec-coord[i,1]),xytext=(cen_ra-coord[i,0]+2,cen_dec-coord[i,1]-4),color=color)

# ax.set_xlim([max(coord[:,0])-cen_ra,min(coord[:,0])-cen_ra])
ax.set_xlim([30,-30])
ax.set_ylim([-30,30])
ax.set_xlabel(r'$\mathrm{RA~(arcsec)}$',fontsize=18)
ax.set_ylabel(r'$\mathrm{Dec~(arcsec)}$',fontsize=18)
ax.set_aspect('equal')
fig.savefig(home+indir+'OI_velocity_shift_spatial.eps',format='eps',dpi=300)
ax.cla()
fig.clf()

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(cen_ra-all_coord[:,0],cen_dec-all_coord[:,1],'s',color='k')
ax.plot(cen_ra-coord[:,0],cen_dec-coord[:,1],'s',color='Red')

for i in range(0,len(all_coord[:,0])):
	ax.annotate('%d' % all_pixel[i],xy=(cen_ra-all_coord[i,0],cen_dec-all_coord[i,1]),xytext=(cen_ra-all_coord[i,0],cen_dec-all_coord[i,1]+1))
	print 'Pixel %d, RA: %16f, Dec: %16f' % (i+1,all_coord[i,0]/3600,all_coord[i,1]/3600)

for i in range(0,len(Tr)):
	if velocity[i] < 0:
		color='Blue'
	elif velocity[i] > 0:
		color='Red'
	else:
		color='k'
	ax.annotate('%6.4f' % Tr[i],xy=(cen_ra-coord[i,0],cen_dec-coord[i,1]),xytext=(cen_ra-coord[i,0]+2,cen_dec-coord[i,1]-4),color=color)

# ax.set_xlim([max(coord[:,0])-cen_ra,min(coord[:,0])-cen_ra])
ax.set_xlim([30,-30])
ax.set_ylim([-30,30])
ax.set_xlabel(r'$\mathrm{RA~(arcsec)}$',fontsize=18)
ax.set_ylabel(r'$\mathrm{Dec~(arcsec)}$',fontsize=18)
ax.set_aspect('equal')
fig.savefig(home+indir+'OI_Tr_spatial.eps',format='eps',dpi=300)
ax.cla()
fig.clf()
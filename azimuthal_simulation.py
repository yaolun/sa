def azimuthal_simulation(rtout, beam_size, wave, dist=178., group=22):
	"""
	rtout: the filepath to the output file of Hyperion
	beam_size: the beam size used for the width of annulus
	dist: the physical distance to the source
	group: the group which contains image
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.constants as const
	from hyperion.model import ModelOutput

	# constant setup
	pc = const.pc.cgs.value
	au = const.au.cgs.value

	# Read in the Hyperion output file
	m = ModelOutput(rtout)
	# get image
	image = m.get_image(group=5, inclination=0, distance=dist*pc, units='Jy')

	# Calculate the image width in arcsec given the distance to the source
	rmax = max(m.get_quantities().r_wall)
	w = np.degrees(rmax / image.distance) * 3600
	# grid of radii of annulus
	annuli = np.linspace(beam_size/2., np.floor((w-beam_size/2.)/beam_size)*beam_size, np.floor((w-beam_size/2.)/beam_size))	# plot

	fig = plt.figure(figsize=(8,6))
	ax = fig.add_subplot(111)

	# iternate through wavelength
	if type(wave) == int or type(wave) == float:
		wave = [wave]
	color_list = plt.cm.viridis(np.linspace(0, 1, len(wave)+1))
	for i in range(len(wave)):
		wav = wave[i]
		# Find the closest wavelength
		iwav = np.argmin(np.abs(wav - image.wav))
		# avoid zero when log, and flip the image
		val = image.val[::-1, :, iwav]
		# determine the center of the image
		npix = len(val[:,0])
		center = np.array([npix/2. + 0.5, npix/2. + 0.5])
		scale = 2*rmax/npix
		# create index array of the image
		x = np.empty_like(val)
		for j in range(len(val[0,:])):
			x[:,j] = j

		flux_annuli = np.empty_like(annuli)
		for k in range(len(annuli)):
			flux_annuli[k] = np.sum(val[(((x-center[0])**2+(x.T-center[1])**2)**0.5*2*w/npix >= annuli[k]-beam_size/2.) & \
										(((x-center[0])**2+(x.T-center[1])**2)**0.5*2*w/npix < annuli[k]+beam_size/2.)])
		flux_annuli = flux_annuli/np.nanmax(flux_annuli)

		ax.plot(np.log10(annuli*dist), np.log10(flux_annuli), 'o-', color=color_list[i], \
				markersize=3, mec='None', label=r'$\rm{'+str(wav)+'\,\mu m}$')
	ax.axvline(np.log10((w-beam_size/2.)*dist), linestyle='--', color='k')
	ax.axvline(np.log10(w*dist), linestyle='-', color='k')

	ax.legend(loc='best', fontsize=12, numpoints=1, ncol=2)
	ax.set_xlabel(r'$\rm{log(Radius)\,[AU]}$', fontsize=18)
	ax.set_ylabel(r'${\rm log(}F/F_{\rm max})$', fontsize=18)
	fig.gca().set_ylim(top=0.1)
	[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
	ax.minorticks_on()
	ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=15,length=5)
	ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=15,length=2.5)

	fig.savefig('/Users/yaolun/test/annuli_profile.pdf', format='pdf', dpi=300, bbox_inches='tight')
	fig.clf()


rtout = '/Users/yaolun/test/model11.rtout'
beam_size = 7.9
wave = [300., 350., 400., 450., 500., 550., 600., 650.]
azimuthal_simulation(rtout, beam_size, wave)

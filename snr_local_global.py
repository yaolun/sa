def snr_local_global(pathlocal, pathglobal, cube=True, oned=False):
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.io.ascii as ascii
	import os
	home = os.path.expanduser('~')

	# Header of the all cube fitting results
	# ====================================================================================================
	# Object,   Line,			LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
	# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), Noise(W/cm2/um),SNR,		 	E_u(K),         A(s-1)        
	# g,        RA(deg),        Dec(deg),       Blend,          Validity
	# ====================================================================================================

	if cube == True:
		pacs_local = ascii.read(home+pathlocal+'CDF_archive_pacs_cube_lines.txt')
		pacs_local = pacs_local[(np.isnan(pacs_local['SNR'])!=True) & (pacs_local['Validity']==1)]
		
		spire_local = ascii.read(home+pathlocal+'CDF_archive_spire_cube_lines.txt')
		spire_local = spire_local[(np.isnan(spire_local['SNR'])!=True) & (spire_local['Validity']==1)]

		pacs_global = ascii.read(home+pathglobal+'CDF_archive_pacs_cube_lines.txt')
		pacs_global = pacs_global[(np.isnan(pacs_global['SNR'])!=True) & (pacs_global['Validity']==1)]
		
		spire_global = ascii.read(home+pathglobal+'CDF_archive_spire_cube_lines.txt')
		spire_global = spire_global[(np.isnan(spire_global['SNR'])!=True) & (spire_global['Validity']==1)]

	if oned == True:
		pacs_local = ascii.read(home+pathlocal+'CDF_archive_pacs_1d_lines.txt')
		pacs_local = pacs_local[(np.isnan(pacs_local['SNR'])!=True) & (pacs_local['Validity']==1)]
		
		spire_local = ascii.read(home+pathlocal+'CDF_archive_spire_1d_lines.txt')
		spire_local = spire_local[(np.isnan(spire_local['SNR'])!=True) & (spire_local['Validity']==1)]

		pacs_global = ascii.read(home+pathglobal+'CDF_archive_pacs_1d_lines.txt')
		pacs_global = pacs_global[(np.isnan(pacs_global['SNR'])!=True) & (pacs_global['Validity']==1)]
		
		spire_global = ascii.read(home+pathglobal+'CDF_archive_spire_1d_lines.txt')
		spire_global = spire_global[(np.isnan(spire_global['SNR'])!=True) & (spire_global['Validity']==1)]

	# PACS
	if len(pacs_local) > len(pacs_global):
		fewer_line_pacs = pacs_global
		more_line_pacs = pacs_local
	else:
		fewer_line_pacs = pacs_local
		more_line_pacs = pacs_global
	ind_pacs = np.empty(len(fewer_line_pacs))
	remove = []
	for i  in range(0, len(fewer_line_pacs)):
		if cube == False:
			if len(np.where((more_line_pacs['Object'] == fewer_line_pacs['Object'][i]) & (more_line_pacs['Line'] == fewer_line_pacs['Line'][i]))[0]) != 1:
				remove = [remove, i]
			else:
				ind_pacs[i] = int(np.where((more_line_pacs['Object'] == fewer_line_pacs['Object'][i]) & (more_line_pacs['Line'] == fewer_line_pacs['Line'][i]))[0])
		else:
			if len(np.where((more_line_pacs['Object'] == fewer_line_pacs['Object'][i]) & (more_line_pacs['Line'] == fewer_line_pacs['Line'][i]) & (more_line_pacs['Pixel_No.'] == fewer_line_pacs['Pixel_No.'][i]))[0]) != 1:
				remove = [remove, i]
			else:
				ind_pacs[i] = int(np.where((more_line_pacs['Object'] == fewer_line_pacs['Object'][i]) & (more_line_pacs['Line'] == fewer_line_pacs['Line'][i]) & (more_line_pacs['Pixel_No.'] == fewer_line_pacs['Pixel_No.'][i]))[0])
	for j in remove:
		fewer_line_pacs.remove_row(j)

	# SPIRE
	if len(spire_local) > len(spire_global):
		fewer_line_spire = spire_global
		more_line_spire = spire_local
	else:
		fewer_line_spire = spire_local
		more_line_spire = spire_global
	ind_spire = np.empty(len(fewer_line_spire))
	remove = []
	for i  in range(0, len(fewer_line_spire)):
		if cube == False:
			if len(np.where((more_line_spire['Object'] == fewer_line_spire['Object'][i]) & (more_line_spire['Line'] == fewer_line_spire['Line'][i]))[0]):
				remove = [remove, i]
			ind_spire[i] = int(np.where((more_line_spire['Object'] == fewer_line_spire['Object'][i]) & (more_line_spire['Line'] == fewer_line_spire['Line'][i]))[0])
		else:
			if len(np.where((more_line_spire['Object'] == fewer_line_spire['Object'][i]) & (more_line_spire['Line'] == fewer_line_spire['Line'][i]) & (more_line_spire['Pixel_No.'] == fewer_line_spire['Pixel_No.'][i]))[0]) != 1:
				remove = [remove, i]
			else:
				ind_spire[i] = int(np.where((more_line_spire['Object'] == fewer_line_spire['Object'][i]) & (more_line_spire['Line'] == fewer_line_spire['Line'][i]) & (more_line_spire['Pixel_No.'] == fewer_line_spire['Pixel_No.'][i]))[0])
	for j in remove:
		fewer_line_spire.remove_row(j)

	fig = plt.figure(figsize=(12,9))
	ax  = fig.add_subplot(111)

	# snr_ratio, = ax.plot(fewer_line_pacs['ObsWL(um)'], more_line_pacs['SNR'][ind_pacs]/fewer_line_pacs['SNR'], 'go', alpha=0.7)
	ax.plot(fewer_line_spire['ObsWL(um)'], more_line_spire['SNR'][ind_spire]/fewer_line_spire['SNR'], 'go', alpha=0.7)

	# lg = plt.legend([n_local, n_global],[r'$\mathrm{1^{st}~fitting~only}$',r'$\mathrm{Full~procedure}$'], fontsize=18, loc='best', numpoints=1)
	ax.set_xlabel(r'$\mathrm{Wavelength~[\mu m]}$', fontsize=20)
	ax.set_ylabel(r'$\mathrm{SNR_{full}/SNR_{1^{st}~only}}$', fontsize=20)
	ax.tick_params('both',labelsize=18,width=1.5,which='major')
	ax.tick_params('both',labelsize=18,width=1.5,which='minor')
	[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

	fig.savefig(home+'/Copy/snr_local_global.png', format='png', dpi=300, bbox_inches='tight')
	fig.clf()

snr_local_global('/test/local/', '/test/')
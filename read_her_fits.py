def unc_test(filepath,plotdir,noise_test=False,oned=True,zoom=False,spire=False,scatter_hist=False,png=False):
	from astropy.io import ascii
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	home = os.path.expanduser('~')

	filepath = home + filepath
	data = ascii.read(filepath)
	if oned == True:
		# Header of the 1-D fitting results
		# =========================================================================================
		# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), Noise(W/cm2/um),SNR,		 	E_u(K),         A(s-1)        
		# g,        RA(deg),        Dec(deg),       Blend,          Validity
		# =========================================================================================
		header = data.colnames
		data = data[np.isnan(data['SNR'])!=True]        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
		data = data[(data['Sig_Cen(um)']!=-999.0) & (data['Sig_FWHM(um)']!=-999.0)]
		snr = abs(data['SNR'][np.argsort(data['ObsWL(um)'])])
		snr_flux = (data['Str(W/cm2)']/data['Sig_str(W/cm2)'])[np.argsort(data['ObsWL(um)'])]
		wl = data['ObsWL(um)'][np.argsort(data['ObsWL(um)'])]
	else:
		# Header of the all cube fitting results
		# ====================================================================================================
		# Object,   Line,			LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), Noise(W/cm2/um),SNR,		 	E_u(K),         A(s-1)        
		# g,        RA(deg),        Dec(deg),       Blend,          Validity
		# ====================================================================================================
		header = data.colnames
		data = data[np.isnan(data['SNR'])!=True]        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
		data = data[(data['Sig_Cen(um)']!=-999.0) & (data['Sig_FWHM(um)']!=-999.0)]
		snr = abs(data['SNR'][np.argsort(data['ObsWL(um)'])])
		snr_flux = (data['Str(W/cm2)']/data['Sig_str(W/cm2)'])[np.argsort(data['ObsWL(um)'])]
		wl = data['ObsWL(um)'][np.argsort(data['ObsWL(um)'])]

	if noise_test == True:
		# Plot noise * FWHM verses Sig_str
		fig = plt.figure(figsize=(12,9))
		ax = fig.add_subplot(111)
		sig_str = data['Sig_str(W/cm2)'][np.argsort(data['ObsWL(um)'])]
		noise = data['Noise(W/cm2/um)'][np.argsort(data['ObsWL(um)'])]
		fwhm = data['FWHM(um)'][np.argsort(data['ObsWL(um)'])]
		if spire != True:
			low, = ax.plot(np.log10( noise[(snr >= 3.0) & (snr < 10.0)] * fwhm[(snr >= 3.0) & (snr < 10.0)] ), np.log10( sig_str[(snr >= 3.0) & (snr < 10.0)] ), 'go')
			high, = ax.plot(np.log10(noise[(snr >= 10.0)]*fwhm[(snr >= 10.0)]), np.log10(sig_str[(snr >= 10.0)]), 'ro')
		else:
			ssw, = ax.plot(np.log10( noise[(snr >= 3.0) & (data['ObsWL(um)'] <= 304)] * fwhm[(snr >= 3.0) & (data['ObsWL(um)'] <= 304)] ), np.log10( sig_str[(snr >= 3.0) & (data['ObsWL(um)'] <= 304)] ), 'go')
			slw, = ax.plot(np.log10( noise[(snr >= 3.0) & (data['ObsWL(um)'] > 304)] * fwhm[(snr >= 3.0) & (data['ObsWL(um)'] > 304)] ), np.log10( sig_str[(snr >= 3.0) & (data['ObsWL(um)'] > 304)] ), 'ro')
			# high, = ax.plot(np.log10(noise[(snr >= 10.0)]*fwhm[(snr >= 10.0)]), np.log10(sig_str[(snr >= 10.0)]), 'ro')
		# ax.plot(np.log10(sig_str[snr < 3.0]), np.log10(noise[snr < 3.0]*fwhm[snr < 3.0]), 'go', alpha = 0.5)
		# overplot the slope = 1 line
		slope, = ax.plot(np.hstack((np.log10(sig_str[snr >= 3.0]), np.log10(noise[snr >= 3.0]*fwhm[snr >= 3.0]))),np.hstack((np.log10(sig_str[snr >= 3.0]), np.log10(noise[snr >= 3.0]*fwhm[snr >= 3.0]))),'b-',linewidth=1.5)
		if spire != True:
			lg = plt.legend([low, high, slope],[r'$\mathrm{3<SNR<10}$',r'$\mathrm{10<SNR}$',r'$\mathrm{Equality}$'], loc='best', numpoints=1, fontsize=20)
		else:
			lg = plt.legend([ssw,slw,slope],[r'$\mathrm{SSW}$',r'$\mathrm{SLW}$',r'$\mathrm{Equality}$'], loc='best', numpoints=1, fontsize=20)
		ax.set_ylabel(r'$\mathrm{log(\sigma_{fit})}$', fontsize=22)
		ax.set_xlabel(r'$\mathrm{log(Noise \times FWHM)}$', fontsize=22)
		ax.tick_params('both',labelsize=18,width=1.5,which='major')
		ax.tick_params('both',labelsize=18,width=1.5,which='minor')
		# ax.set_ylim([-1e-19,5e-19])
		# ax.set_xlim([-1e-19,5e-19])
		addname = ''
		if zoom == True:
			ax.set_xlim([-25,-15])
			ax.set_ylim([-25,-15])
			addname = '_zoomin'
		[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
		if png == False:
			fig.savefig(home+plotdir+'unc_comparison'+addname+'.pdf',format='pdf',dpi=300,bbox_inches='tight')
		else:
			fig.savefig(home+plotdir+'unc_comparison'+addname+'.png',format='png',dpi=300,bbox_inches='tight')
		fig.clf()

		if scatter_hist == True:
			print 'Plotting scatter/histogram...'
			from matplotlib.ticker import NullFormatter

			nullfmt   = NullFormatter()         # no labels

			# data
			x = np.log10( noise[(snr >= 3.0) & (snr < 10.0)] * fwhm[(snr >= 3.0) & (snr < 10.0)] )
			y = np.log10( sig_str[(snr >= 3.0) & (snr < 10.0)] )

			xx = np.log10( noise[(snr >= 10.0)] * fwhm[(snr >= 10.0)] )
			yy = np.log10( sig_str[(snr >= 10.0)] )

			line = np.hstack((np.log10(sig_str[snr >= 3.0]), np.log10(noise[snr >= 3.0]*fwhm[snr >= 3.0])))

			# definitions for the axes
			left, width = 0.1, 0.65
			bottom, height = 0.1, 0.65
			bottom_h = left_h = left+width+0.02

			rect_scatter = [left, bottom, width, height]
			rect_histx = [left, bottom_h, width, 0.2]
			rect_histy = [left_h, bottom, 0.2, height]

			# start with a rectangular Figure
			fig = plt.figure(1,figsize=(12,12))

			axScatter = fig.add_axes(rect_scatter)
			axHistx = fig.add_axes(rect_histx)
			axHisty = fig.add_axes(rect_histy)

			# set labels for scatter plot
			axScatter.set_xlabel(r'$\mathrm{log(Noise\times FWHM)~(W~cm^{-2})}$', fontsize=22)
			axScatter.set_ylabel(r'$\mathrm{log(\sigma_{fit})~(W~cm^{-2})}$', fontsize=22)

			# no labels
			axHistx.xaxis.set_major_formatter(nullfmt)
			axHisty.yaxis.set_major_formatter(nullfmt)

			# no tick labels
			axHistx.yaxis.set_ticks([])
			axHisty.xaxis.set_ticks([])

			# the scatter plot:
			low = axScatter.scatter(x, y, c='g', s=27, lw=0.5)
			high = axScatter.scatter(xx, yy, c='r', s=27, lw=0.5)
			slope, = axScatter.plot(line, line, 'b-', linewidth=1.5)
			lg = plt.legend([low, high, slope],[r'$\mathrm{3<SNR<10}$',r'$\mathrm{10<SNR}$',r'$\mathrm{Equality}$'],\
							loc='best', bbox_to_anchor=[0.35,1],bbox_transform=axScatter.transAxes, numpoints=1, scatterpoints=1, fontsize=18)
			
			# now determine nice limits by hand:
			binwidth = 0.05
			xymax = np.max( [np.max(x), np.max(y)] )
			xymin = np.min( [np.min(x), np.min(y)] )
			ulim = ( int(xymax/binwidth) + 1) * binwidth + 0.5
			llim = ( int(xymin/binwidth) + 1) * binwidth - 0.5

			axScatter.set_xlim( (llim, ulim) )
			axScatter.set_ylim( (llim, ulim) )

			bins = np.arange(llim, ulim + binwidth , binwidth)
			axHistx.hist(x, bins=bins, histtype='step', edgecolor='Green',hatch='////')
			axHisty.hist(y, bins=bins, histtype='step', edgecolor='Green',hatch='////', orientation='horizontal')

			axHistx.hist(xx, bins=bins, histtype='step', edgecolor='Red', hatch='\\\\\\\\',)
			axHisty.hist(yy, bins=bins, histtype='step', edgecolor='Red', hatch='\\\\\\\\', orientation='horizontal')

			axHistx.set_xlim( axScatter.get_xlim() )
			axHisty.set_ylim( axScatter.get_ylim() )

			# Tweak the axes thickness
			axes = [axHistx,axHisty]
			for ax in axes:
				ax.tick_params('both',labelsize=14,width=1.5,which='major')
				ax.tick_params('both',labelsize=14,width=1.5,which='minor')
				[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
			ax = axScatter
			ax.tick_params('both',labelsize=14,width=1.5,which='major')
			ax.tick_params('both',labelsize=14,width=1.5,which='minor')
			[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

			if png == False:
				plt.savefig(home+plotdir+'scatter_hist.pdf',format='pdf',dpi=300,bbox_inches='tight')
			else:
				plt.savefig(home+plotdir+'scatter_hist.png',format='png',dpi=300,bbox_inches='tight')
			plt.cla()
			plt.clf()


	# fig = plt.figure(figsize=(12,9))
	# ax = fig.add_subplot(111)

	# SNR, = ax.plot(wl[snr >= 3.0], snr[snr >= 3.0]/snr_flux[snr >= 3.0], 'go', linewidth=1.5)
	# ax.plot(wl[snr < 3.0], snr[snr < 3.0]/snr_flux[snr < 3.0], 'go', alpha = 0.5, linewidth=1.5)
	# # SNR_flux, = ax.plot(wl[snr_flux >= 3.0], snr_flux[snr_flux >= 3.0], 'go', linewidth=1.5)
	# # ax.plot(wl[snr_flux < 3.0], snr_flux[snr_flux < 3.0], 'go', alpha = 0.5, linewidth=1.5)
	# ax.set_xlabel(r'$\mathrm{Wavelength~(\mu m)}$', fontsize=20)
	# ax.set_ylabel(r'$\mathrm{SNR}$', fontsize=20)
	# # lg = ax.legend([SNR, SNR_flux], [r'$\mathrm{SNR}$',r'$\mathrm{F_{line}/\sigma_{line}}$'], fontsize=18)
	# lg = ax.legend([SNR], [r'$\mathrm{\frac{SNR}{F_{line}/\sigma_{line}}}$'], numpoints=1, fontsize=18)
	# ax.tick_params('both',labelsize=18,width=1.5,which='major')
	# ax.tick_params('both',labelsize=18,width=1.5,which='minor')
	# [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
	# ax.set_ylim([min(snr/snr_flux),max(snr/snr_flux)])
	# fig.savefig(home+plotdir+'snr_comparison.pdf',format='pdf',dpi=300,bbox_inches='tight')
	# fig.clf()
zoom = False
scatter_hist = True
png = True
# filepath = '/FWD_archive_pacs_1d_lines.txt'
# plotdir = '/pacs_1d_'
# unc_test(filepath,plotdir,noise_test=True,zoom=zoom)
# filepath = '/FWD_archive_spire_1d_lines.txt'
# plotdir = '/spire_1d_'
# unc_test(filepath,plotdir,noise_test=True)
filepath = '/FWD_archive_pacs_cube_lines.txt'
plotdir = '/pacs_cube_'
unc_test(filepath,plotdir,noise_test=True,zoom=zoom,scatter_hist=scatter_hist, png=png)
filepath = '/FWD_archive_spire_cube_lines.txt'
plotdir = '/spire_cube_'
unc_test(filepath,plotdir,noise_test=True,zoom=zoom,scatter_hist=scatter_hist, png=png)
# filepath = '/test/fitting_test/WL12_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
# plotdir = '/WL12_pacs_'
# unc_test(filepath,plotdir,noise_test=True)
# filepath = '/test/fitting_test/BHR71_spire_corrected_lines.txt'
# plotdir = '/BHR71_spire_'
# unc_test(filepath,plotdir,noise_test=True)

# # Using the uncertainty from the fits file
# filepath = '/bhr71/fitting/unc_test/pacs/advanced_products/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
# plotname = '/bhr71/plots/snr_comparison_unc'
# unc_test(filepath,plotname)
# # Doesn't input the uncertainty
# filepath = '/bhr71/fitting/1107/pacs/advanced_products/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
# plotname = '/bhr71/plots/snr_comparison_zero'
# unc_test(filepath,plotname)

def fitting_check(indir,obj=None,outdir=False):
	from astropy.io import ascii
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	home = os.path.expanduser('~')

	# Child function to find the files that match the given pattern
	def find(pattern, path):
		import os, fnmatch
		result = []
		for root, dirs, files in os.walk(path):
			for name in files:
				if fnmatch.fnmatch(name, pattern):
					result.append(os.path.join(root, name))
		return result

	# Search for the fitting results tables under the primary directory
	# If the object name is specified, then only search in the the object directory
	if obj == None:
		# objdir = [x[0] for x in os.walk(home+indir)][1:]
		objdir = os.walk(home+indir).next()[1]
		pacspath = []
		spirepath = []
		num_pacs = 0.0
		num_spire = 0.0
		for o in objdir:
			# PACS search
			path_dum = find('*_lines.txt', home+indir+'/'+o+'/pacs/')
			if len(path_dum) > 0:
				pacspath.extend(path_dum)
				num_pacs  += 1
			# SPIRE search
			path_dum = find('*_lines.txt', home+indir+'/'+o+'/spire/')
			if len(path_dum) > 0:
				spirepath.extend(path_dum)
				num_spire += 1
	else:
		num_pacs = 1.0
		num_spire = 1.0
		pacspath = find('*_lines.txt', home+indir+'/'+obj+'/pacs/')
		spirepath = find('*_lines.txt', home+indir+'/'+obj+'/spire/')

	num_fit = 0.0
	num_line = 0.0
	num1 = 0.0
	num2 = 0.0
	num3 = 0.0
	num4 = 0.0
	num5 = 0.0
	num6 = 0.0
	num7 = 0.0
	num8 = 0.0
	num9 = 0.0
	num_line2 = 0.0
	num_dg = 0.0
	num_dg_line = 0.0

	# PACS statistic
	for path in pacspath:
		data = ascii.read(path)
		# Header of the 1-D fitting results
		# =========================================================================================
		# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,		 	E_u(K),         A(s-1)        
		# g,        RA(deg),        Dec(deg),       Blend,          Validity
		# =========================================================================================
		header = data.colnames
		# If there is a missing segment in the spectra, the SNR will return NaN.
		num1 = len(data[np.isnan(data['SNR'])==True]) + num1
		# Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
		data = data[np.isnan(data['SNR'])!=True]
		num_fit = len(data['SNR']) + num_fit
		num_line = len(data[data['SNR'] >= 3]) + num_line
		num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
		num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
		num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
		num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
		num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
		num7 = len(data[data['Sig_str(W/cm2)'] == 0.0]) + num7
		num8 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) + num8
		num9 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
		num_line2 = len(data[(data['SNR'] >=3) & (data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & \
							 (data['Sig_str(W/cm2)'] != 0.0) & (data['Validity'] ==1)]) + num_line2
		num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
		num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

		# Print out the detail information of the line that has zero in line strength uncertainty.
		if len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) != 0:
			print path
			print data['Line'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['ObsWL(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['Sig_Cen(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['Str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['Sig_str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
			print data['Sig_FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]


	# Print out the statistic of the pacs fitting results
	print '<PACS>'
	print '\t Number of object: %d ' % num_pacs
	print '\t %d lines fitted, %.2f lines fitted per object' % (num_fit,num_fit/num_pacs)
	print '\t %d detections, %.2f detections per object.' % (num_line,num_line/num_pacs)
	print '\t %d lines fitted with blend Gaussian, %d lines detections among them.' % (num_dg,num_dg_line)
	print '\t <<Anomaly>>'
	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
	print '\t \t Zeros in line centroid uncertainty: %d and %d with detections.' % (num2,num3)
	print '\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.' % (num4,num5,num6)
	print '\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian' % (num7,num8,num9)
	print '\t %d detections without anomalous, and %.2f lines per object.' % (num_line2,num_line2/num_pacs)

	# ============================================================================================================================
	# ============================================================================================================================
	print '=============================================================================================================='
	num_fit = 0.0
	num_line = 0.0
	num1 = 0.0
	num2 = 0.0
	num3 = 0.0
	num4 = 0.0
	num5 = 0.0
	num6 = 0.0
	num7 = 0.0
	num8 = 0.0
	num9 = 0.0
	num_line2 = 0.0
	num_dg = 0.0
	num_dg_line = 0.0

	# SPIRE statistic
	for path in spirepath:
		data = ascii.read(path)
		# Header of the 1-D fitting results
		# =========================================================================================
		# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,		 	E_u(K),         A(s-1)        
		# g,        RA(deg),        Dec(deg),       Blend,          Validity
		# =========================================================================================
		header = data.colnames
		# If there is a missing segment in the spectra, the SNR will return NaN.
		num1 = len(data[np.isnan(data['SNR'])==True]) + num1
		# Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
		data = data[np.isnan(data['SNR'])!=True]
		num_fit = len(data['SNR']) + num_fit
		num_line = len(data[data['SNR'] >= 3]) + num_line
		num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
		num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
		num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
		num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
		num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
		num7 = len(data[data['Sig_str(W/cm2)'] == 0.0]) + num7
		num8 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) + num8
		num9 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
		num_line2 = len(data[(data['SNR'] >=3) & (data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & \
							 (data['Sig_str(W/cm2)'] != 0.0) & (data['Validity'] == 1)]) + num_line2
		num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
		num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

	# Print out the statistic of the pacs fitting results
	print '<SPIRE>'
	print '\t Number of object: %d ' % num_spire
	print '\t %d lines fitted, %.2f lines fitted per object' % (num_fit,num_fit/num_spire)
	print '\t %d detections, %.2f detections per object.' % (num_line,num_line/num_spire)
	print '\t %d lines fitted with blend Gaussian, %d lines detections among them.' % (num_dg,num_dg_line)
	print '\t <<Anomaly>>'
	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
	print '\t \t Zeros in line centroid uncertainty: %d and %d with detections.' % (num2,num3)
	print '\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.' % (num4,num5,num6)
	print '\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian' % (num7,num8,num9)
	print '\t %d detections without anomalous, and %.2f lines per object.' % (num_line2,num_line2/num_spire)


# indir = '/bhr71/fitting/'
# fitting_check(indir)

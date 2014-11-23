def unc_test(filepath,plotname):
	from astropy.io import ascii
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	home = os.path.expanduser('~')

	filepath = home + filepath
	data = ascii.read(filepath)
	# Header of the 1-D fitting results
	# =========================================================================================
	# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
	# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,		 	E_u(K),         A(s-1)        
	# g,        RA(deg),        Dec(deg),       Blend,          Validity
	# =========================================================================================
	header = data.colnames
	data = data[np.isnan(data['SNR'])!=True]        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
	data = data[(data['Sig_Cen(um)']!=-999.0) & (data['Sig_FWHM(um)']!=-999.0)]
	snr = abs(data['SNR'][np.argsort(data['ObsWL(um)'])])
	snr_flux = (data['Str(W/cm2)']/data['Sig_str(W/cm2)'])[np.argsort(data['ObsWL(um)'])]
	wl = data['ObsWL(um)'][np.argsort(data['ObsWL(um)'])]

	fig = plt.figure(figsize=(12,9))
	ax = fig.add_subplot(111)

	SNR, = ax.plot(wl[snr >= 3.0], snr[snr >= 3.0]/snr_flux[snr >= 3.0], 'go', linewidth=1.5)
	ax.plot(wl[snr < 3.0], snr[snr < 3.0]/snr_flux[snr < 3.0], 'go', alpha = 0.5, linewidth=1.5)
	# SNR_flux, = ax.plot(wl[snr_flux >= 3.0], snr_flux[snr_flux >= 3.0], 'go', linewidth=1.5)
	# ax.plot(wl[snr_flux < 3.0], snr_flux[snr_flux < 3.0], 'go', alpha = 0.5, linewidth=1.5)
	ax.set_xlabel(r'$\mathrm{Wavelength~(\mu m)}$', fontsize=20)
	ax.set_ylabel(r'$\mathrm{SNR}$', fontsize=20)
	# lg = ax.legend([SNR, SNR_flux], [r'$\mathrm{SNR}$',r'$\mathrm{F_{line}/\sigma_{line}}$'], fontsize=18)
	lg = ax.legend([SNR], [r'$\mathrm{\frac{SNR}{F_{line}/\sigma_{line}}}$'], numpoints=1, fontsize=18)
	ax.tick_params('both',labelsize=18,width=1.5,which='major')
	ax.tick_params('both',labelsize=18,width=1.5,which='minor')
	[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
	ax.set_ylim([min(snr/snr_flux),max(snr/snr_flux)])
	fig.savefig(home+plotname+'.pdf',format='pdf',dpi=300,bbox_inches='tight')
	fig.clf()

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


	# Print out the statistic of the pacs fitting results
	print '<PACS>'
	print '\t Number of object: %d ' % num_pacs
	print '\t %d lines are fitted, %.2f lines are fitted in each object on average.' % (num_fit,num_fit/num_pacs)
	print '\t %d lines are detected, %.2f lines are fitted in each objects on average.' % (num_line,num_line/num_pacs)
	print '\t %d lines are fitted with double Gaussian, %d lines are detected among them.' % (num_dg,num_dg_line)
	print '\t <<Anomaly>>'
	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
	print '\t \t Zeros in line centroid uncertainty: %d and %d with detection.' % (num2,num3)
	print '\t \t Zeros in FWHM uncertainty: %d, %d with detection, and %d with detection and double Gaussian.' % (num4,num5,num6)
	print '\t \t Zeros in line strength uncertainty: %d, %d with detection, and %d with detection and double Gaussian' % (num7,num8,num9)
	print '\t %d lines are detected without anomalous, and %.2f lines on average.' % (num_line2,num_line2/num_pacs)

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
							 (data['Sig_str(W/cm2)'] != 0.0) & (data['Validity'] ==1)]) + num_line2
		num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
		num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

	# Print out the statistic of the pacs fitting results
	print '<SPIRE>'
	print '\t Number of object: %d ' % num_spire
	print '\t %d lines are fitted, %.2f lines are fitted in each object on average.' % (num_fit,num_fit/num_spire)
	print '\t %d lines are detected, %.2f lines are fitted in each objects on average.' % (num_line,num_line/num_spire)
	print '\t %d lines are fitted with double Gaussian, %d lines are detected among them.' % (num_dg,num_dg_line)
	print '\t <<Anomaly>>'
	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
	print '\t \t Zeros in line centroid uncertainty: %d and %d with detection.' % (num2,num3)
	print '\t \t Zeros in FWHM uncertainty: %d, %d with detection, and %d with detection and double Gaussian.' % (num4,num5,num6)
	print '\t \t Zeros in line strength uncertainty: %d, %d with detection, and %d with detection and double Gaussian' % (num7,num8,num9)
	print '\t %d lines are detected without anomalous, and %.2f lines on average.' % (num_line2,num_line2/num_spire)



# indir = '/bhr71/fitting/'
# fitting_check(indir)

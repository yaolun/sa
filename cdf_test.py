def obj_com(indir):
	import numpy as np
	import os
	home = os.path.expanduser('~')

	objdir = os.walk(home+indir).next()[1]
	print objdir
obj_com('/')

def fits_com(indir, objname):

def strong_line(indir,obj=None):

def unc_test(filepath, outdir):

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

	num_test = 0.0

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

		num_test = len(data[data['Sig_FWHM(um)'] == 0.0]) + num_test

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

	print num_test

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

	num_test = 0.0

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

		num_test = len(data[data['Sig_FWHM(um)'] == 0.0]) + num_test


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

	print num_test

def cdf_test(indir,outdir):

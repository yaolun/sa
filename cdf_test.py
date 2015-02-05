def obj_com(indir):
	import os
	home = os.path.expanduser('~')

	# pre-define the full object list

	pacsobj = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
			   'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
			   'HD97048','HD98922','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
			   'L1489','L1527','L1551-IRS5','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RULup','RYLup','SCra','SR21',\
			   'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']


	spireobj = ['B1-a','B1-c','B335','BHR71','Ced110-IRS4','FUOri','GSS30-IRS1','HH46','HH100','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','L1014',\
			   'L1157','L1455-IRS3','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO91','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg',\
			   'V1515Cyg','V1735Cyg','VLA1623','WL12']

	objlist = [pacsobj, spireobj]

	objlist = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','Ced110-IRS4','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
			   'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
			   'HD97048','HD98922','HH46','HH100','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
			   'L1489','L1527','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RNO91','RULup','RYLup','SCra','SR21',\
			   'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']

	objdir = os.walk(home+indir).next()[1]
	objdir = objdir[objdir != 'contour']
	err = 0
	for o in objlist:
		if len(objdir[objdir == o]) != 1:
			print 'Cannot find ', o
			err += 1
		if o in pacsobj:
			# Check 1d and cube fitting results
			if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt') == False:
				err += 1
				print 'Missing PACS 1d fitting on ', o
			if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/cube/'+o+'_pacs_pixel13_os8_sf7_lines.txt') == False:
				err += 1
				print 'Missing PACS cube fitting on ', o
		if o in spireobj:
			# Check 1d and cube fitting results
			if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/'+o+'_spire_corrected_lines.txt') == False:
				err += 1
				print 'Missing SPIRE 1d fitting on ', o
			if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_SLWC3_lines.txt') == False:
				err += 1
				print 'Missing SPIRE-SLW cube fitting on ', o
			if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_SSWD4_lines.txt') == False:
				err += 1
				print 'Missing SPIRE-SSW cube fitting on ', o
	if err == 0:
		print 'Passed the object test!'

# obj_com('/data/FWD_bettyjo/FWD_archive_slim/')

def fits_com(indir):
	import os
	home = os.path.expanduser('~')

	# pre-define OBSID info

	objlist = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','Ced110-IRS4','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
			   'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
			   'HD97048','HD98922','HH46','HH100','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
			   'L1489','L1527','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RNO91','RULup','RYLup','SCra','SR21',\
			   'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']

	obsid = [['ABAur','1342217842','1342217843','0'],\
			 ['AS205','1342215737','1342215738','0'],\
			 ['B1-a','1342216182','1342216183','1342249475'],\
			 ['B1-c','1342216213','1342216214','1342249476'],\
			 ['B335','1342208889','1342208888','1342253652'],\
			 ['BHR71','1342212230','1342212231','1342248249'],\
			 ['Ced110-IRS4','0','0','1342248246'],\
			 ['DGTau','1342225730','1342225731','0'],\
			 ['EC82','1342192975','1342219435','0'],\
			 ['Elias29','1342228519','1342228520','0'],\
			 ['FUOri','1342250907','1342250908','1342230412'],\
			 ['GSS30-IRS1','1342215678','1342215679','1342251286'],\
			 ['HD100453','1342211695','1342211696','0'],\
			 ['HD100546','1342188037','1342188038','0'],\
			 ['HD104237','1342207819','1342207820','0'],\
			 ['HD135344B','1342213921','1342213922','0'],\
			 ['HD139614','1342215683','1342215684','0'],\
			 ['HD141569','1342213913','0','0'],\
			 ['HD142527','1342216174','1342216175','0'],\
			 ['HD142666','1342213916','0','0'],\
			 ['HD144432','1342213919','0','0'],\
			 ['HD144668','1342215641','1342215642','0'],\
			 ['HD150193','1342227068','0','0'],\
			 ['HD163296','1342217819','1342217820','0'],\
			 ['HD169142','1342206987','1342206988','0'],\
			 ['HD179218','1342208884','1342208885','0'],\
			 ['HD203024','1342206975','0','0'],\
			 ['HD245906','1342228582','0','0'],\
			 ['HD35187','1342217846','0','0'],\
			 ['HD36112','1342228247','1342228248','0'],\
			 ['HD38120','1342226212','1342226213','0'],\
			 ['HD50138','1342206991','1342206992','0'],\
			 ['HD97048','1342199412','1342199413','0'],\
			 ['HD98922','1342210385','0','0'],\
			 ['HH46','0','0','1342245084'],\
			 ['HH100','0','0','1342252897'],\
			 ['HTLup','1342213920','0','0'],\
			 ['IRAM04191','1342216654','1342216655','0'],\
			 ['IRAS03245','1342214677','1342214676','1342249053'],\
			 ['IRAS03301','1342215668','1342216181','1342249477'],\
			 ['IRAS12496','1342188039','1342188040','1342254037'],\
			 ['IRAS15398','0','0','1342250515'],\
			 ['IRS46','1342228474','1342228475','1342251289'],\
			 ['IRS48','1342227069','1342227070','0'],\
			 ['L1014','1342208911','1342208912','1342245857'],\
			 ['L1157','1342208909','1342208908','1342247625'],\
			 ['L1448-MM','1342213683','1342214675','0'],\
			 ['L1455-IRS3','1342204122','1342204123','1342249474'],\
			 ['L1489','1342216216','1342216215','0'],\
			 ['L1527','1342192981','1342192982','0'],\
			 ['L1551-IRS5','1342192805','1342229711','1342249470'],\
			 ['L483','0','0','1342253649'],\
			 ['L723-MM','0','0','1342245094'],\
			 ['RCrA-IRS5A','1342207806','1342207805','1342253646'],\
			 ['RCrA-IRS7B','1342207807','1342207808','1342242620'],\
			 ['RCrA-IRS7C','1342206990','1342206989','1342242621'],\
			 ['RNO90','1342228206','0','0'],\
			 ['RNO91','0','0','1342251285'],\
			 ['RULup','1342215682','0','0'],\
			 ['RYLup','1342216171','0','0'],\
			 ['SCra','1342207809','1342207810','0'],\
			 ['SR21','1342227209','1342227210','0'],\
			 ['Serpens-SMM3','1342193216','1342193214','0'],\
			 ['Serpens-SMM4','1342193217','1342193215','0'],\
			 ['TMC1','1342225803','1342225804','1342250512'],\
			 ['TMC1A','1342192987','1342192988','1342250510'],\
			 ['TMR1','1342192985','1342192986','1342250509'],\
			 ['V1057Cyg','1342235853','1342235852','1342221695'],\
			 ['V1331Cyg','1342233446','1342233445','1342221694'],\
			 ['V1515Cyg','1342235691','1342235690','1342221685'],\
			 ['V1735Cyg','1342235849','1342235848','1342219560'],\
			 ['VLA1623','1342213918','1342213917','1342251287'],\
			 ['WL12','1342228187','1342228188','1342251290']]
	err = 0
	for i in range(0, len(obsid)):
		for j in range(1,3):
			if obsid[i][j] != '0':
				# check FITS files for PACS 1d and cube
				# 1d - 9Spx
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_blue_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_red_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				# 1d - 3x3NO
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				# 1d - 3x3YES
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits') == False:
					print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits' % (obsid[i][0], obsid[i][j], obsid[i][0])
					err += 1
				# cube - finalcube
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_finalcubes_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_blue_finalcubes_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_finalcubes_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_red_finalcubes_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
				# cube - noda
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnoda_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_blue_rebinnedcubesnoda_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnoda_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_red_rebinnedcubesnoda_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
				# cube - nodb
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnodb_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_blue_rebinnedcubesnodb_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
				if os.path.exists(home+indir+'/'+obsid[i][0]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnodb_os8_sf7.fits') == False:
					print '%s missing OBSID_%s_red_rebinnedcubesnodb_os8_sf7.fits' % (obsid[i][0], obsid[i][j])
					err += 1
		if obsid[i][3] != '0':
			# check FITS files for SPIRE 1d and cube
			if (os.path.exists(home+indir+'/'+obsid[i][0]+'/spire/data/fits/OBSID_'+obsid[i][j]+'_spire_corrected.fits') or \
				os.path.exists(home+indir+'/'+obsid[i][0]+'/spire/data/fits/OBSID_'+obsid[i][j].lower()+'_spire_corrected.fits')) == False:
				print '%s missing OBSID_%s_spire_corrected.fits' % (obsid[i][0], obsid[i][j])
				err += 1
			if os.path.exists(home+indir+'/'+obsid[i][0]+'/spire/data/fits/'+obsid[i][j]+'_spectrum_extended_HR_aNB_15.fits') == False:
				print '%s missing %s_spectrum_extended_HR_aNB_15.fits' % (obsid[i][0], obsid[i][j])
				err += 1
	if err == 0:
		print 'passed the FITS test!'

# fits_com('/data/FWD_bettyjo/FWD_archive_slim/')

# def strong_line(indir,obj=None):

# def unc_test(filepath, outdir):

# def fitting_check(indir,obj=None,outdir=False):
# 	from astropy.io import ascii
# 	import numpy as np
# 	import matplotlib.pyplot as plt
# 	import os
# 	home = os.path.expanduser('~')

# 	# Child function to find the files that match the given pattern
# 	def find(pattern, path):
# 		import os, fnmatch
# 		result = []
# 		for root, dirs, files in os.walk(path):
# 			for name in files:
# 				if fnmatch.fnmatch(name, pattern):
# 					result.append(os.path.join(root, name))
# 		return result

# 	# Search for the fitting results tables under the primary directory
# 	# If the object name is specified, then only search in the the object directory
# 	if obj == None:
# 		# objdir = [x[0] for x in os.walk(home+indir)][1:]
# 		objdir = os.walk(home+indir).next()[1]
# 		pacspath = []
# 		spirepath = []
# 		num_pacs = 0.0
# 		num_spire = 0.0
# 		for o in objdir:
# 			# PACS search
# 			path_dum = find('*_lines.txt', home+indir+'/'+o+'/pacs/')
# 			if len(path_dum) > 0:
# 				pacspath.extend(path_dum)
# 				num_pacs  += 1
# 			# SPIRE search
# 			path_dum = find('*_lines.txt', home+indir+'/'+o+'/spire/')
# 			if len(path_dum) > 0:
# 				spirepath.extend(path_dum)
# 				num_spire += 1
# 	else:
# 		num_pacs = 1.0
# 		num_spire = 1.0
# 		pacspath = find('*_lines.txt', home+indir+'/'+obj+'/pacs/')
# 		spirepath = find('*_lines.txt', home+indir+'/'+obj+'/spire/')

# 	num_fit = 0.0
# 	num_line = 0.0
# 	num1 = 0.0
# 	num2 = 0.0
# 	num3 = 0.0
# 	num4 = 0.0
# 	num5 = 0.0
# 	num6 = 0.0
# 	num7 = 0.0
# 	num8 = 0.0
# 	num9 = 0.0
# 	num_line2 = 0.0
# 	num_dg = 0.0
# 	num_dg_line = 0.0

# 	num_test = 0.0

# 	# PACS statistic
# 	for path in pacspath:
# 		data = ascii.read(path)
# 		# Header of the 1-D fitting results
# 		# =========================================================================================
# 		# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
# 		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,		 	E_u(K),         A(s-1)        
# 		# g,        RA(deg),        Dec(deg),       Blend,          Validity
# 		# =========================================================================================
# 		header = data.colnames
# 		# If there is a missing segment in the spectra, the SNR will return NaN.
# 		num1 = len(data[np.isnan(data['SNR'])==True]) + num1
# 		# Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
# 		data = data[np.isnan(data['SNR'])!=True]
# 		num_fit = len(data['SNR']) + num_fit
# 		num_line = len(data[data['SNR'] >= 3]) + num_line
# 		num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
# 		num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
# 		num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
# 		num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
# 		num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
# 		num7 = len(data[data['Sig_str(W/cm2)'] == 0.0]) + num7
# 		num8 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) + num8
# 		num9 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
# 		num_line2 = len(data[(data['SNR'] >=3) & (data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & \
# 							 (data['Sig_str(W/cm2)'] != 0.0) & (data['Validity'] ==1)]) + num_line2
# 		num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
# 		num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

# 		num_test = len(data[data['Sig_FWHM(um)'] == 0.0]) + num_test

# 		# Print out the detail information of the line that has zero in line strength uncertainty.
# 		if len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) != 0:
# 			print path
# 			print data['Line'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['ObsWL(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['Sig_Cen(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['Str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['Sig_str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
# 			print data['Sig_FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]


# 	# Print out the statistic of the pacs fitting results
# 	print '<PACS>'
# 	print '\t Number of object: %d ' % num_pacs
# 	print '\t %d lines fitted, %.2f lines fitted per object' % (num_fit,num_fit/num_pacs)
# 	print '\t %d detections, %.2f detections per object.' % (num_line,num_line/num_pacs)
# 	print '\t %d lines fitted with blend Gaussian, %d lines detections among them.' % (num_dg,num_dg_line)
# 	print '\t <<Anomaly>>'
# 	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
# 	print '\t \t Zeros in line centroid uncertainty: %d and %d with detections.' % (num2,num3)
# 	print '\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.' % (num4,num5,num6)
# 	print '\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian' % (num7,num8,num9)
# 	print '\t %d detections without anomalous, and %.2f lines per object.' % (num_line2,num_line2/num_pacs)

# 	print num_test

# 	# ============================================================================================================================
# 	# ============================================================================================================================
# 	print '=============================================================================================================='
# 	num_fit = 0.0
# 	num_line = 0.0
# 	num1 = 0.0
# 	num2 = 0.0
# 	num3 = 0.0
# 	num4 = 0.0
# 	num5 = 0.0
# 	num6 = 0.0
# 	num7 = 0.0
# 	num8 = 0.0
# 	num9 = 0.0
# 	num_line2 = 0.0
# 	num_dg = 0.0
# 	num_dg_line = 0.0

# 	num_test = 0.0

# 	# SPIRE statistic
# 	for path in spirepath:
# 		data = ascii.read(path)
# 		# Header of the 1-D fitting results
# 		# =========================================================================================
# 		# Line,		LabWL(um),		ObsWL(um),		Sig_Cen(um),	Str(W/cm2),		Sig_str(W/cm2)
# 		# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,		 	E_u(K),         A(s-1)        
# 		# g,        RA(deg),        Dec(deg),       Blend,          Validity
# 		# =========================================================================================
# 		header = data.colnames
# 		# If there is a missing segment in the spectra, the SNR will return NaN.
# 		num1 = len(data[np.isnan(data['SNR'])==True]) + num1
# 		# Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
# 		data = data[np.isnan(data['SNR'])!=True]
# 		num_fit = len(data['SNR']) + num_fit
# 		num_line = len(data[data['SNR'] >= 3]) + num_line
# 		num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
# 		num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
# 		num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
# 		num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
# 		num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
# 		num7 = len(data[data['Sig_str(W/cm2)'] == 0.0]) + num7
# 		num8 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) + num8
# 		num9 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
# 		num_line2 = len(data[(data['SNR'] >=3) & (data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & \
# 							 (data['Sig_str(W/cm2)'] != 0.0) & (data['Validity'] == 1)]) + num_line2
# 		num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
# 		num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

# 		num_test = len(data[data['Sig_FWHM(um)'] == 0.0]) + num_test


# 	# Print out the statistic of the pacs fitting results
# 	print '<SPIRE>'
# 	print '\t Number of object: %d ' % num_spire
# 	print '\t %d lines fitted, %.2f lines fitted per object' % (num_fit,num_fit/num_spire)
# 	print '\t %d detections, %.2f detections per object.' % (num_line,num_line/num_spire)
# 	print '\t %d lines fitted with blend Gaussian, %d lines detections among them.' % (num_dg,num_dg_line)
# 	print '\t <<Anomaly>>'
# 	print '\t \t SNR anomalies due to the missing spectra: %d' % num1
# 	print '\t \t Zeros in line centroid uncertainty: %d and %d with detections.' % (num2,num3)
# 	print '\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.' % (num4,num5,num6)
# 	print '\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian' % (num7,num8,num9)
# 	print '\t %d detections without anomalous, and %.2f lines per object.' % (num_line2,num_line2/num_spire)

# 	print num_test

# def cdf_test(indir,outdir):

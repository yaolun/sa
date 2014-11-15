def get_name(objname):
	if objname == 'b1-a_spire': name = 'B1-a'
	if objname == 'b1-c_spire': name = 'B1-c'
	if objname == 'b335_spire': name = 'B335'
	if objname == 'bhr71': name = 'BHR71'
	if objname == 'ced110_spire': name = 'Ced110-IRS4'
	if objname == 'dkcha_spire': name = 'DKCha'
	if objname == 'fuori_spire': name = 'FU Ori'
	if objname == 'gss30_spire': name = 'GSS30-IRS1'
	if objname == 'hbc722r1': name = 'HBC722'
	if objname == 'hd100546_spire': name = 'HD100546'
	if objname == 'hh100': name = 'HH100'
	if objname == 'hh46_spire': name = 'HH46'
	if objname == 'iras03245_spire': name = 'IRAS03245'
	if objname == 'iras03301_spire': name = 'IRAS03301'
	if objname == 'iras15398_spire': name = 'IRAS15398'
	if objname == 'irs44': name = 'IRS44'
	if objname == 'irs46': name = 'IRS46'
	if objname == 'l1157_spire': name = 'L1157'
	if objname == 'l1455-irs3': name = 'L1455-IRS3'
	if objname == 'l1489_spire': name = 'L1489'
	if objname == 'l1551-irs5_spire': name = 'L1551-IRS5'
	if objname == 'l483_spire': name = 'L483'
	if objname == 'l723-mm_spire': name = 'L723-MM'
	if objname == 'rcra-irs5a_spire': name = 'RCrA-IRS5A'
	if objname == 'rcra-irs7b_spire': name = 'RCrA-IRS7B'
	if objname == 'rcra-irs7c_spire': name = 'RCrA-IRS7C'
	if objname == 'rno91_spire': name = 'RNO91'
	if objname == 'tmc1a_spire': name = 'TMC1A'
	if objname == 'tmc1_spire': name = 'TMC1'
	if objname == 'tmr1_spire': name = 'TMR1'
	if objname == 'v1057cyg_spire': name = 'V1057 Cyg'
	if objname == 'v1331cyg_spire': name = 'V1331 Cyg'
	if objname == 'v1515cyg_spire': name = 'V1515 Cyg'
	if objname == 'v1735cyg': name = 'V1735 Cyg'
	if objname == 'vla1623_spire': name = 'VLA1623-243'
	if objname == 'wl12_spire': name = 'WL12'
	return name

def comparison(objname):
	import numpy as np
	import scipy as sp
	from scipy.interpolate import interp1d
	import matplotlib
	import matplotlib.pyplot as plt
	import os
	home = os.path.expanduser('~')
	pixel = ['SLWA1','SLWA2','SLWB1','SLWB2','SLWB3','SLWC2','SLWC3','SLWC4','SLWC5','SLWD2','SLWD3','SLWD4','SLWE2','SLWE3']
	pixel = ['center_ssw', 'center_slw']
	pixel = [objname+'_slw_lines', objname+'_ssw_lines', 'pacs_pixel13_lines']
	pixel = ['SLWC3lines','SSWD4lines','pacs_pixel13_lines']
	name_tot = []
	wl_tot = []
	wl_sig_tot = []
	flux_tot = []
	flux_sig_tot = []
	for i in pixel:
		#filename = '/Users/yaolun/Rebecca/data/'+i+'lines.txt'
		#filename = '/Users/yaolun/bhr71/data/'+i+'.txt'
		filename = home+'/data/spire_corrected_joel/line_fitting/'+i+'.txt'
		data = open(filename, "r").readlines()
		name = []
		wl = []
		wl_sig = []
		flux = []
		flux_sig = []
		snr = []
		data.pop(0)
		wl_below = []
		wl_sig_below = []
		flux_below = []
		flux_sig_below = []
		snr_below = []
		for line in data:
			t = line.split()
			if t[1] <> 'ERROR:':
				if float(t[8]) >= 3:
					#print t[7]
					name.append(t[0])
					wl.append(float(t[1]))
					wl_sig.append(float(t[2]))
					flux.append(float(t[3])/1e-22)
					flux_sig.append(float(t[4])/1e-22)
					snr.append(float(t[8]))
				else:
					wl_below.append(float(t[1]))
					wl_sig_below.append(float(t[2]))
					flux_below.append(float(t[3])/1e-22)
					flux_sig_below.append(float(t[4])/1e-22)
					snr_below.append(float(t[8]))
	#	filename_becca = '/Users/yaolun/Rebecca/L1455/L1455-IRS3_'+i+'_'+j+'.txt'
	#		if os.path.isfile(filename_becca) == True:
	#			data_b = open(filename_becca, "r").readlines()
	#			data_b.pop(0)
	#			for line in data_b:
	#				t = line.split()
	#				name_b.append(t[0]+t[1])
	#				wl_b.append(float(t[3]))
	#				flux_b.append(float(t[5]))
	#	print i
	#	print name
	#	print wl_sig
	#	print flux_sig
		#print wl
		#print flux
	#	if wl <> []:
	#		plt.plot(wl, flux, 'bo')
	#		plt.errorbar(wl, flux, xerr=wl_sig, yerr=flux_sig, fmt='bo')
	#	if wl_b <> []:
			#plt.plot(wl_b, flux_b, 'ro')
		#plt.xlim([200, 670])
		#flux.extend(flux_b)
		#plt.ylim([min(flux)-1, max(flux)+1])
		#plt.xlabel('Wavelength ($\mu$m)')
		#plt.ylabel('Flux Intensity (10$^{-22}$ W/cm$^{2}$/$\mu$m)')
		#plt.title(i)
		#plt.savefig(filename=home+'/data/spire_corrected_joel/plots/comparison_'+i+'.eps', format='eps', dpi=300)
		#plt.clf()
		#plt.cla()
		name_tot.extend(name)
		wl_tot.extend(wl)
		wl_sig_tot.extend(wl_sig)
		flux_tot.extend(flux)
		flux_sig_tot.extend(flux_sig)
	print wl_tot
	name_b = []
	wl_b = []
	wl_b_sig = []
	flux_b = []
	flux_b_sig = []
	filename_becca = home+'/data/cops_linelist/'+objname+'_line_joel.txt'
	data_b = open(filename_becca, "r").readlines()
	#data_b.pop(0)
	for line in data_b:
		if len(line) > 1:
			t = line.split()
			#print len(line)
			name_b.append(t[0])
			wl_b.append(float(t[1]))
			flux_b.append(float(t[3]))
			flux_b_sig.append(float(t[4]))
	#plot the difference of two reductions
	wl_frac = []
	flux_frac = []
	for line_name in name_b:
		if (line_name in name_tot) == True:
			ind_tot = name_tot.index(line_name)
			ind_b = name_b.index(line_name)
			wl_frac.append(wl_tot[ind_tot])
			flux_frac.append(flux_tot[ind_tot]-flux_b[ind_b])
	#plot the comparison of 12CO and 13CO
	wl_tot_co = []
	flux_tot_co = []
	wl_b_co = []
	flux_b_co = []
	wl_tot_13co = []
	flux_tot_13co = []
	wl_b_13co = []
	flux_b_13co = []
	for line_name in name_b:
		if '13CO' in line_name:
			ind_b = name_b.index(line_name)
			wl_b_13co.append(wl_b[ind_b])
			flux_b_13co.append(flux_b[ind_b])
		elif line_name in ['CO13-12','CO12-11','CO11-10','CO10-9','CO9-8','CO8-7','CO7-6','CO6-5','CO5-4','CO4-3']:
			ind_b = name_b.index(line_name)
			wl_b_co.append(wl_b[ind_b])
			flux_b_co.append(flux_b[ind_b])
	for line_name in name_tot:
		if '13CO' in line_name:
			ind_tot = name_tot.index(line_name)
			wl_tot_13co.append(wl_tot[ind_tot])
			flux_tot_13co.append(flux_tot[ind_tot])
		elif line_name in ['CO13-12','CO12-11','CO11-10','CO10-9','CO9-8','CO8-7','CO7-6','CO6-5','CO5-4','CO4-3']:
			ind_tot = name_tot.index(line_name)
			wl_tot_co.append(wl_tot[ind_tot])
			flux_tot_co.append(flux_tot[ind_tot])
	plt.plot(wl_tot, flux_tot, 'bo')
	plt.errorbar(wl_tot, flux_tot, xerr=wl_sig_tot, yerr=flux_sig_tot, fmt='bo')
	#plt.plot(wl_b, flux_b, 'r+')
	#plt.errorbar(wl_b, flux_b, yerr=flux_b_sig, fmt='r+')
	plt.xlim([55, 670])
	plt.ylim([0, max(flux_tot)*1.1])
	plt.xlabel('Wavelength ($\mu$m)')
	plt.ylabel('Flux Intensity (10$^{-22}$ W/cm$^{2}$/$\mu$m)')
	plt.text(550, max(flux_tot), get_name(objname), fontsize=16)
	#plt.title('Extended Correction')
	plt.savefig(filename=home+'/data/spire_corrected_joel/plots/comparison/'+objname+'_comparison_Ext_corrected.eps', format='eps', dpi=300)
	plt.clf()
	plt.cla()
	if flux_frac != []:
		plt.plot(wl_frac, flux_frac, 'bo')
		plt.xlim([55,670])
		plt.xlabel('Wavelength ($\mu$m)')
		plt.ylabel('Difference (new vs old)')
		plt.text(550, max(flux_frac)*0.96, get_name(objname), fontsize=16)
		plt.savefig(filename=home+'/data/spire_corrected_joel/plots/comparison/'+objname+'_fraction_Ext_corrected.eps', format='eps', dpi=300)
		plt.clf()
		plt.cla()
	plt.plot(wl_tot_co, flux_tot_co, 'bo')
	plt.plot(wl_tot_13co, flux_tot_13co, 'b^')
	plt.plot(wl_b_co, flux_b_co, 'ro')
	plt.plot(wl_b_13co, flux_b_13co, 'r^')
	plt.xlim([55,670])
	plt.xlabel('Wavelength ($\mu$m)')
	plt.ylabel('Flux Intensity (10$^{-22}$ W/cm$^{2}$/$\mu$m)')
	plt.text(550, max(flux_tot_co)*0.96, get_name(objname), fontsize=16)
	plt.savefig(filename=home+'/data/spire_corrected_joel/plots/comparison/'+objname+'CO_diff_Ext_corrected.eps',format='eps',dpi=300)
	plt.cla()
	plt.clf()
	#interpol_fit = interp1d(wl_tot, flux_tot)
	#interpol_old = interp1d(wl_b, flux_b)
	#wl_tot_new = np.linspace(min(wl_tot), max(wl_tot), 4700)
	#wl_b_new = np.linspace(min(wl_b), max(wl_b), 4700)
	#plt.plot(wl_tot_new, interpol_fit(wl_tot_new)/interpol_old(wl_b_new))
	#plt.xlim([180, 670])
	#plt.ylim([min(interpol_fit(wl_tot_new)/interpol_old(wl_b_new)), max(intepol_fit(wl_tot_new)/interpol_old(wl_b_new))])
	#plt.xlabel('Wavelength ($\mu$m)')
	#plt.ylabel('Fitting/Old')
	#plt.text(550, max(interpol_fit(wl_tot_new)/interpol_old(wl_b_new)), get_name(objname), fontsize=16)
	#plt.savefig(filename=home+'/data/spire_corrected_joel/plots/comparison/'+objname+'_fraction_Ext_corrected.eps', format='eps', dpi=300)
	#plt.clf()
	#plt.cla()
	print objname

obj = ['b1-a_spire','hh46_spire','l1551-irs5_spire','tmc1a_spire','b1-c_spire','iras03245_spire','l483_spire','tmr1_spire','b335_spire','iras03301_spire','l723-mm_spire','v1057cyg_spire','ced110_spire','iras15398_spire','rcra-irs5a_spire','v1331cyg_spire','dkcha_spire','rcra-irs7b_spire','v1515cyg_spire','fuori_spire','l1157_spire','rcra-irs7c_spire','vla1623_spire','gss30_spire','l1455-irs3','rno91_spire','wl12_spire','hd100546_spire','l1489_spire','tmc1_spire','bhr71']   #No correspond line list from Joel for l1014 right now.
obj = ['bhr71']

for objname in obj:
	comparison(objname)

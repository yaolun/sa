def wl_offset(pacs_path=None,spire_path=None,outdir=None,mag=1.5,obj=False,cube=False):

	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from read_fitting import read_fitting
	import scipy.optimize as opt

	# Define the file I/O
	#
	if outdir == None:
		plot_dir = os.path.expanduser('~')+os.path.dirname(pacs_path)+'/'

	# Setting for Gaussian fitting
	#
	# Equation for Gaussian
	def f(x, a, b, c):
		return a * np.exp(-(x - b)**2.0 / (2 * c**2))

	# Read the resolution table
	#
	c = 2.998e10
	# PACS
	home = os.path.expanduser('~')
	[wl1,res1] = np.genfromtxt(home+'/bhr71/data/spectralresolution_order1.txt').T
	[wl2,res2] = np.genfromtxt(home+'/bhr71/data/spectralresolution_order2.txt').T
	# [wl3,res3] = np.genfromtxt(home+'/bhr71/data/spectralresolution_order3.txt').T
	wl_pacs_res = np.concatenate((wl1[wl1 > max(wl2)],wl2))
	res_pacs = wl_pacs_res/np.concatenate((res1[wl1 > max(wl2)],res2))/2.354

	res_pacs = c/(wl_pacs_res*1e-4)**2*res_pacs*1e-4
	freq_pacs_res = c*1e4/wl_pacs_res

	res_pacs = res_pacs[np.argsort(freq_pacs_res)]/1e9
	freq_pacs_res = freq_pacs_res[np.argsort(freq_pacs_res)]/1e9
	# SPIRE
	wl_spire_res = np.arange(200,670,0.01)
	res_spire = 1.4472/2.998e8*1e3*(wl_spire_res**2)/2.354

	res_spire = c/(wl_spire_res*1e-4)**2*res_spire*1e-4
	freq_spire_res = c*1e4/wl_spire_res

	res_spire = res_spire[np.argsort(freq_spire_res)]/1e9
	freq_spire_res = freq_spire_res[np.argsort(freq_spire_res)]/1e9

	# Read the full fitting text file
	#
	# header of the fitting results
	# Line				LabWL(um)	ObsWL(um)	Sig_Cen(um)	Str(W/cm2)	Sig_str(W/cm2)	FWHM(um)	Sig_FWHM(um)
	# Base(W/cm2/um)	SNR			E_u(K)		A(s-1)		g           RA(deg)          Dec(deg)	Blend
	# Validity
	#
	fig = plt.figure(figsize=(8*mag,6*mag))
	# PACS
	if pacs_path != None:
		# fig = plt.figure(figsize=(8*mag,6*mag))
		ax_pacs = fig.add_subplot(111)
		[data,header] = read_fitting(pacs_path,3)
		lab_wl = data['LabWL(um)']		# np.array(data[:,1]).astype('float')
		obs_wl = data['ObsWL(um)']		# np.array(data[:,2]).astype('float')
		sig_wl = data['Sig_Cen(um)'] 	# np.array(data[:,3]).astype('float')
		snr    = data['SNR']			# np.array(data[:,9]).astype('float')

		ind = np.where(sig_wl != -999)
		indd = np.where(sig_wl == -999)

		lab_freq = c*1e4/lab_wl[ind]/1e9
		obs_freq = c*1e4/obs_wl[ind]/1e9
		sig_freq = c*1e4/(obs_wl[ind])*(sig_wl[ind])/(obs_wl[ind])/1e9
		snr      = snr[ind]
		for dum in snr:
			if np.isnan(dum) == True:
				print dum
		delta_freq = lab_freq-obs_freq

		lab_freqq = c*1e4/lab_wl[indd]/1e9
		obs_freqq = c*1e4/obs_wl[indd]/1e9
		sig_freqq = c*1e4/(obs_wl[indd])*(sig_wl[indd])/(obs_wl[indd])/1e9
		delta_freqq = lab_freqq-obs_freqq

		# obs_wl = obs_wl[(sig_wl != 0) & (sig_wl < 5)]
		# delta_wl = delta_wl[(sig_wl != 0) & (sig_wl < 5)]
		# sig_wl = sig_wl[(sig_wl != 0) & (sig_wl < 5)]

		delta = np.sum(delta_freq/sig_freq**2)/np.sum(1/sig_freq**2)
		sig_delta = (1/np.sum(1/sig_freq**2))**0.5
		print delta,sig_delta

		delta_strong = np.sum(delta_freq[snr>5]/sig_freq[snr>5]**2)/np.sum(1/sig_freq[snr>5]**2)
		sig_delta_strong = (1/np.sum(1/sig_freq[snr>5]**2))**0.5
		print delta_strong,sig_delta_strong

		# Plot the range of line centroids that put into the fitting constraint
		#
		fill = ax_pacs.fill_between(freq_pacs_res[freq_pacs_res>min(obs_freq)],res_pacs[freq_pacs_res>min(obs_freq)]*(-4),res_pacs[freq_pacs_res>min(obs_freq)]*4,facecolor='Orange',alpha=0.3,interpolate=True)
		# Plot the detection from the fitting results
		#
		ax_pacs.plot(obs_freq,delta_freq,'o',mec='None',mfc='Blue',alpha=0.2)
		ax_pacs.errorbar(obs_freq,delta_freq,linestyle='None',yerr=sig_freq,color='b',alpha=0.2)
		ax_pacs.plot(obs_freqq,delta_freqq,'o',mec='None',mfc='Gray',alpha=0.2)
		#
		ax_pacs.plot(freq_pacs_res[freq_pacs_res>min(obs_freq)],res_pacs[freq_pacs_res>min(obs_freq)],color='k',linestyle='--',linewidth=2)
		ax_pacs.plot(freq_pacs_res[freq_pacs_res>min(obs_freq)],-res_pacs[freq_pacs_res>min(obs_freq)],color='k',linestyle='--',linewidth=2)
		#
		ax_pacs.axhline(linewidth=2,color='k',linestyle='--')
		ax_pacs.set_xlim([min(obs_freq),max(obs_freq)])
		ax_pacs.set_xlabel(r'$\mathrm{Line~Centroid~(GHz)}$',fontsize=18*mag)
		ax_pacs.set_ylabel(r'$\mathrm{\Delta \nu~(GHz)}$',fontsize=18*mag)
		ax_pacs.set_ylim([min(delta_freq)*1.5,-min(delta_freq)*1.5])
		[ax_pacs.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
		ax_pacs.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
		ax_pacs.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)
		pacs_text = ax_pacs.text(0.35,0.9,r'$\mathrm{<\Delta\nu>= %5.2f \times 10^{%d} \pm %5.2f \times 10^{%d}~GHz}$' % (delta/10**np.floor(np.log10(abs(delta))),np.floor(np.log10(abs(delta))),sig_delta/10**np.floor(np.log10(sig_delta)),np.floor(np.log10(sig_delta))),fontsize=mag*16,transform=ax_pacs.transAxes)
		fig.savefig(plot_dir+'wl_offset_pacs.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)
		pacs_text.remove()
		# Plot the frequency offset distribution
		#
		fig_hist = plt.figure(figsize=(8*mag,6*mag))
		ax_hist = fig_hist.add_subplot(111)
		num_bins = 50.0
		# All detections
		n, bins, patches = ax_hist.hist(delta_freq/obs_freq*c/1e5, bins=np.linspace(min(delta_freq/obs_freq*c/1e5),max(delta_freq/obs_freq*c/1e5),num_bins),\
										facecolor='green', alpha=0.5)
		x = (bins[1:]+bins[:-1])/2.0
		y = n
		popt, pcov = opt.curve_fit(f, x, y,p0=[max(n),0,max(res_pacs)/6000*c/1e5])
		fwhm = abs(popt[2])*2.354
		x_fit = np.linspace(x[0], x[-1], 500)
		y_fit = f(x_fit, *popt)
		ax_hist.plot(x_fit,y_fit,linewidth=2,color='Blue')
		# Detections with S/N > 5
		n, bins, patches = ax_hist.hist(delta_freq[snr>5]/obs_freq[snr>5]*c/1e5, bins=np.linspace(min(delta_freq/obs_freq*c/1e5),max(delta_freq/obs_freq*c/1e5),num_bins),\
										facecolor='green', alpha=1)
		ax_hist.text(0.6,0.9,r'$\mathrm{Center= %6.3f~km/s}$' % popt[1],fontsize=mag*16,transform=ax_hist.transAxes)
		ax_hist.text(0.6,0.8,r'$\mathrm{FWHM= %6.3f~km/s}$' % fwhm,fontsize=mag*16,transform=ax_hist.transAxes)
		ax_hist.set_xlabel(r'$\mathrm{Velocity~(km/s)}$',fontsize=18*mag)
		ax_hist.set_ylabel(r'$\mathrm{Number~of~detections}$',fontsize=18*mag)
		[ax_hist.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
		ax_hist.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
		ax_hist.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)
		ax_hist.set_xlim([-600,600])
		ax_hist.plot([-550,-550+max(res_pacs[freq_pacs_res>min(obs_freq)])/max(freq_pacs_res)*c/1e5],[max(n)*0.2,max(n)*0.2],linewidth=4,color='Green')
		ax_hist.text(-550,max(n)*0.12,r'$\mathrm{%6.3f~km/s}$' % (max(res_pacs[freq_pacs_res>min(obs_freq)])/max(freq_pacs_res)*c/1e5),fontsize=mag*16)
		fig_hist.savefig(plot_dir+'wl_offset_distribution_pacs.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)
		plt.close(fig_hist)

	if spire_path != None:
		# fig = plt.figure(figsize=(8*mag,6*mag))
		ax_spire = fig.add_subplot(111)
		[data,header] = read_fitting(spire_path,3)
		lab_wl = data['LabWL(um)']		# np.array(data[:,1]).astype('float')
		obs_wl = data['ObsWL(um)']		# np.array(data[:,2]).astype('float')
		sig_wl = data['Sig_Cen(um)'] 	# np.array(data[:,3]).astype('float')
		snr    = data['SNR']			# np.array(data[:,9]).astype('float')
		line   = data['Line']			# np.array(data[:,0]).astype('str')

		ind = np.where(sig_wl != -999)
		indd = np.where(sig_wl == -999)

		lab_freq = c*1e4/lab_wl[ind]/1e9
		obs_freq = c*1e4/obs_wl[ind]/1e9
		sig_freq = c*1e4/(obs_wl[ind])*(sig_wl[ind])/(obs_wl[ind])/1e9
		line     = line[ind]
		snr      = snr[ind]
		delta_freq = lab_freq-obs_freq

		lab_freqq = c*1e4/lab_wl[indd]/1e9
		obs_freqq = c*1e4/obs_wl[indd]/1e9
		sig_freqq = c*1e4/(obs_wl[indd])*(sig_wl[indd])/(obs_wl[indd])/1e9
		delta_freqq = lab_freqq-obs_freqq

		delta = np.sum(delta_freq/sig_freq**2)/np.sum(1/sig_freq**2)
		sig_delta = (1/np.sum(1/sig_freq**2))**0.5
		print delta,sig_delta

		delta_strong = np.sum(delta_freq[snr>5]/sig_freq[snr>5]**2)/np.sum(1/sig_freq[snr>5]**2)
		sig_delta_strong = (1/np.sum(1/sig_freq[snr>5]**2))**0.5
		print delta_strong,sig_delta_strong


		velo = delta_freq/obs_freq*c/1e5
		# if ('Object' in data.colnames) == True:
		# 	name = (data['Object'])[ind]
		# 	print velo[velo>1000],name[velo>1000],line[velo>1000]

		# Plot the range of line centroids that put into the fitting constraint
		#
		fill = ax_spire.fill_between(freq_spire_res[freq_spire_res>min(obs_freq)],res_spire[freq_spire_res>min(obs_freq)]*(-4),res_spire[freq_spire_res>min(obs_freq)]*4,facecolor='Green',alpha=0.3,interpolate=True)
		# Plot the detection from the fitting results
		#
		ax_spire.plot(obs_freq,delta_freq,'o',mec='None',mfc='Red',alpha=0.2)
		ax_spire.errorbar(obs_freq,delta_freq,linestyle='None',yerr=sig_freq,color='r',alpha=0.2)
		ax_spire.plot(obs_freqq,delta_freqq,'o',mec='None',mfc='Gray',alpha=0.2)
		#
		ax_spire.plot(freq_spire_res[freq_spire_res>min(obs_freq)],res_spire[freq_spire_res>min(obs_freq)],color='k',linestyle='--',linewidth=2)
		ax_spire.plot(freq_spire_res[freq_spire_res>min(obs_freq)],-res_spire[freq_spire_res>min(obs_freq)],color='k',linestyle='--',linewidth=2)
		#
		ax_spire.axhline(linewidth=2,color='k',linestyle='--')
		ax_spire.set_xlim([c*1e4/670/1e9,c*1e4/200/1e9])
		# ax_spire.set_xlim([165,175])
		ax_spire.set_xlabel(r'$\mathrm{Line~Centroid~(GHz)}$',fontsize=18*mag)
		ax_spire.set_ylabel(r'$\mathrm{\Delta \nu~(GHz)}$',fontsize=18*mag)
		ax_spire.set_ylim([min(delta_freq)*2,-min(delta_freq)*2])
		ax_spire.set_xlabel(r'$\mathrm{Line~Centroid~(GHz)}$',fontsize=18*mag)
	 	ax_spire.set_ylabel(r'$\mathrm{\Delta \nu~(GHz)}$',fontsize=18*mag)
	 	[ax_spire.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
		ax_spire.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
		ax_spire.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)
		spire_text = ax_spire.text(0.35,0.9,r'$\mathrm{<\Delta\nu>= %5.2f \times 10^{%d} \pm %5.2f \times 10^{%d}~GHz}$' % (delta/10**np.floor(np.log10(abs(delta))),np.floor(np.log10(abs(delta))),sig_delta/10**np.floor(np.log10(sig_delta)),np.floor(np.log10(sig_delta))),fontsize=mag*16,transform=ax_spire.transAxes)
		ax_spire.set_ylim([min(delta_freq)*1.5,-min(delta_freq)*1.5])
	 	fig.savefig(plot_dir+'wl_offset_spire.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)
	 	spire_text.remove()
		# Plot the frequency offset distribution
		#
		fig_hist = plt.figure(figsize=(8*mag,6*mag))
		ax_hist = fig_hist.add_subplot(111)
		num_bins = 50
		# All detections
		n, bins, patches = ax_hist.hist(delta_freq/obs_freq*c/1e5, bins=np.linspace(min(delta_freq/obs_freq*c/1e5),max(delta_freq/obs_freq*c/1e5),num_bins),\
										facecolor='Orange', alpha=0.5)
		x = (bins[1:]+bins[:-1])/2.0
		y = n
		popt, pcov = opt.curve_fit(f, x, y,p0=[max(n),0,max(res_spire)/1000*c/1e5])
		fwhm = abs(popt[2])*2.354
		x_fit = np.linspace(x[0], x[-1], 500)
		y_fit = f(x_fit, *popt)
		ax_hist.plot(x_fit,y_fit,linewidth=2,color='Blue')
		# Detections with S/N > 5
		n, bins, patches = ax_hist.hist(delta_freq[snr>5]/obs_freq[snr>5]*c/1e5, bins=np.linspace(min(delta_freq/obs_freq*c/1e5),max(delta_freq/obs_freq*c/1e5),num_bins),\
										facecolor='Orange', alpha=1)
		ax_hist.plot([-1300,-1300+max(res_spire[freq_spire_res>min(obs_freq)])/max(freq_spire_res)*c/1e5],[max(n)*0.2,max(n)*0.2],linewidth=4,color='Orange')
		ax_hist.text(-1300,max(n)*0.12,r'$\mathrm{%6.3f~km/s}$' % (max(res_spire[freq_spire_res>min(obs_freq)])/max(freq_spire_res)*c/1e5),fontsize=16*mag)
		ax_hist.text(0.6,0.9,r'$\mathrm{Center= %6.3f~km/s}$' % popt[1],fontsize=mag*16,transform=ax_hist.transAxes)
		ax_hist.text(0.6,0.8,r'$\mathrm{FWHM= %6.3f~km/s}$' % fwhm,fontsize=mag*16,transform=ax_hist.transAxes)
		ax_hist.set_xlabel(r'$\mathrm{Velocity~(km/s)}$',fontsize=18*mag)
		ax_hist.set_ylabel(r'$\mathrm{Number~of~detections}$',fontsize=18*mag)
		[ax_hist.spines[axis].set_linewidth(1.5*mag) for axis in ['top','bottom','left','right']]
		ax_hist.tick_params('both',labelsize=mag*18,width=1.5*mag,which='major',pad=15,length=5*mag)
		ax_hist.tick_params('both',labelsize=mag*18,width=1.5*mag,which='minor',pad=15,length=2.5*mag)
		ax_hist.set_xlim([-1500,1500])
		fig_hist.savefig(plot_dir+'wl_offset_distribution_spire.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)
		plt.close(fig_hist)

	if (pacs_path != None) and (spire_path != None):
		plt.xlim([c*1e4/670/1e9,c*1e4/50/1e9])
		plt.ylim([-6,6])
		fig.savefig(plot_dir+'wl_offset.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)

pacs_path = '/data/FWD_bettyjo/FWD_archive_slim/FWD_archive_pacs_cube_lines.txt'
spire_path = '/data/FWD_bettyjo/FWD_archive_slim/FWD_archive_spire_cube_lines.txt'
# pacs_path = '/data/FWD_archive/FWD_archive_pacs_1d_lines.txt'
# spire_path = None
# pacs_path = '/data/FWD_archive/TMC1/pacs/advanced_products/TMC1_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
wl_offset(pacs_path,spire_path,cube=True)

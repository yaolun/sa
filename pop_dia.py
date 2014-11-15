def pop_dia_co(indir,filename,plotdir,spire=None):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from read_fitting import read_fitting_co
	from leastsqfit import lin_leastsqfit
	home = os.path.expanduser('~')

	filepath = indir+filename+'.txt'
	[co_data,co_data_name] = read_fitting_co(filepath,3)
	#Check the fitting see if the spectrum has enough co detection for rotational temperature fitting
	if len(co_data_name)<=2:
		return None
	#co_data structure: wl, wl_sig, flux(W/cm2), flux_sig(W/cm2), E_u, A, g
	c = 2.99792458e10
	h = 6.6260755e-27
	k = 1.380658e-16
	B = 1.9225#192.25 m-1
	pc = 3.086e18
	v = c/(co_data[:,0]*1e-4)
	N = 4*np.pi*co_data[:,2]*1e7*(178*pc)**2/(co_data[:,5]*h*v)
	N_sigma = 4*np.pi*co_data[:,3]*1e7*(178*pc)**2/(co_data[:,5]*h*v)
	x = co_data[:,4]
	y = np.log10(N/co_data[:,6])
	yerr_hi = np.log10((N+N_sigma)/co_data[:,6])-np.log10(N/co_data[:,6])
	yerr_low = np.log10(N/co_data[:,6])-np.log10((N-N_sigma)/co_data[:,6])
	y_sig = y*0
	for i in range(0,len(y)):
	 		y_sig[i] = max(yerr_hi[i], yerr_low[i])
	ind = np.argsort(x)
	x = x[ind]
	y = y[ind]
	y_sig = y_sig[ind]

	#Single temperature fitting
	if len(x)>2:
		[yfit,yerr,t_rot,sig_t_rot,s_min,yoff] = lin_leastsqfit(x, y, y_sig)
		Q = float(k*t_rot/h/c/B)
		N_fit = Q*10**(float(yoff))
		x = x.reshape(len(x))
		y = y.reshape(len(y))
		y_sig = y_sig.reshape(len(y_sig))
		fig_rot_dia = plt.figure()
		ax_rot_dia = fig_rot_dia.add_subplot(111)
		data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
		ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

		fit, = ax_rot_dia.plot(x,yfit,color='DarkMagenta')
		ax_rot_dia.plot(x,yfit+yerr,'--',color='Magenta')
		ax_rot_dia.plot(x,yfit-yerr,'--',color='Magenta')
		ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=14)
		ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=14)
		ax_rot_dia.set_xlim([500,4500])
		ax_rot_dia.legend([fit],[r'$\mathrm{T_{rot}= %8.4f \pm %8.6f~K,~\mathcal{N}= %.3e~g/cm^{2}}$' % (t_rot,sig_t_rot,N_fit)],numpoints=1,loc='upper right',fontsize=12)
		fig_rot_dia.savefig(home+plotdir+filename+'_co_rot_single.eps',format='eps',dpi=300)
		ax_rot_dia.cla()
		fig_rot_dia.clf()
		print filename
		print 'T_rot: %8.6f K' % t_rot
		s_min_single = s_min
	#Two temperature fitting
		if len(x)>=6:
			best_fit = []
			s_min = []
			for i in range(int(x[2]+1),int(x[-3]-1),10):
				turning_pt = i
				[yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
				[yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
				best_fit.append(i)
				s_min.append(s_min_warm+s_min_cool)
			best_fit = np.array(best_fit)
			s_min = np.array(s_min)
			# fig_s_min = plt.figure()
			# ax_s_min = fig_s_min.add_subplot(111)
			# ax_s_min.plot(best_fit,s_min)
			# ax_s_min.set_xlabel(r'$\mathrm{Turning~points}$')
			# ax_s_min.set_ylabel(r'$\mathrm{S_{min}}$')
			# fig_s_min.savefig(home+'/bhr71/plots/rotational_diagram/s_min.eps',format='eps',dpi=300)
			# ax_s_min.cla()
			# fig_s_min.clf()

			turning_pt = np.mean(best_fit[s_min == min(s_min)])
			[yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
			[yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
			if (s_min_cool+s_min_warm)<s_min_single:
				Q_warm = float(k*t_rot_warm/h/c/B)
				Q_cool = float(k*t_rot_cool/h/c/B)
				N_warm_fit = Q_warm*10**(float(yoff_warm))
				N_cool_fit = Q_cool*10**(float(yoff_cool))
				fig_rot_dia = plt.figure()
				ax_rot_dia = fig_rot_dia.add_subplot(111)
				data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
				ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

				fit_warm, = ax_rot_dia.plot(x[x>=turning_pt],yfit_warm,color='DarkMagenta')
				ax_rot_dia.plot(x[x>=turning_pt],yfit_warm+yerr_warm,'--',color='Magenta')
				ax_rot_dia.plot(x[x>=turning_pt],yfit_warm-yerr_warm,'--',color='Magenta')

				fit_cool, = ax_rot_dia.plot(x[x<turning_pt],yfit_cool,color='Blue')
				ax_rot_dia.plot(x[x<turning_pt],yfit_cool+yerr_cool,'--',color='MediumBlue')
				ax_rot_dia.plot(x[x<turning_pt],yfit_cool-yerr_cool,'--',color='MediumBlue')

				ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=14)
				ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=14)
				ax_rot_dia.set_xlim([500,4500])

				ax_rot_dia.legend([fit_warm,fit_cool],[r'$\mathrm{T_{rot,warm}= %8.4f \pm %8.6f~K,~\mathcal{N}= %.3e~g/cm^{2}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit),r'$\mathrm{T_{rot,cool}= %8.4f \pm %8.6f~K,~\mathcal{N}= %.3e~g/cm^{2}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit)],numpoints=1,loc='upper right',fontsize=12)
				fig_rot_dia.savefig(home+plotdir+filename+'_co_rot_two.eps',format='eps',dpi=300)
				ax_rot_dia.cla()
				fig_rot_dia.clf()
				print filename
				print 'T_rot(warm): %8.6f K, T_rot(cool): %8.6f K' % (t_rot_warm,t_rot_cool)
	# Add SPIRE data
	if spire != None:
		[co_spire,co_spire_name] = read_fitting_co(spire,3)




filename = 'centralSpaxels_correctedYES_lines_localbaseline_fixwidth_global_noise'
pop_dia_co('/bhr71/data/latest/',filename,'/bhr71/plots/rotational_diagram/')
#PACS
for i in range(1,26):
	filename = 'pacs_pixel'+str(i)+'_lines_localbaseline_fixwidth_global_noise'
	pop_dia_co('/bhr71/data/latest/',filename,'/bhr71/plots/rotational_diagram/')
#SPIRE

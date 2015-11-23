import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys
import os
from phot_reader import phot_reader
from phot_filter import phot_filter

objlist = ['B1-a','B1-c','B335','BHR71','FUOri','GSS30-IRS1','HH46','IRAS03245',\
		   'IRAS03301','IRAS12496','IRAS15398','L1014','L1157','L1455-IRS3','L1551-IRS5',\
		   'RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','TMC1','TMC1A','TMR1','V1331Cyg','V1515Cyg',\
		   'VLA1623','WL12'] #HD135344
# objlist = ['B1-a','B335','BHR71']

# print objlist
# objlist = ['AS205','B1-a','B1-c','B335','BHR71','Elias29','FUOri','GSS30-IRS1','HD100546']

photdir = '/Users/yaolun/data/herschel_phot/'
CDFdir = '/Users/yaolun/data/CDF_archive/'
archival_dir = '/Users/yaolun/data/herschel_archival/'

SLW_beam = np.pi/4.*34**2
SSW_beam = np.pi/4.*19**2

# option to plot different bands in different symbol
diff = False

data_dict = {'object': [], 'phot': [], 'spec_phot': [], 'archival_spec_phot': [], 'wave': []}

def herschel_spec_phot(wl, flux, pacs=False, spire=True, filter_func=False):
	import numpy as np
	from scipy.interpolate import interp1d
	from phot_filter import phot_filter
	phot_wl = []
	if pacs == True:
		phot_wl.extend([70.,100.]) # no enough point in longer wavelength to make an accurate 160um photometry
	if spire == True:
		phot_wl.extend([250.,350.,500.])
	phot_wl = np.array(phot_wl)
	phot_flux = np.empty_like(phot_wl)

	# clean up NaN value in observation
	wl = wl[np.isnan(flux) == False]
	flux = flux[np.isnan(flux) == False]

	for i in range(len(phot_wl)):
		if filter_func == False:
			res = 3     # temp. resolution to mimic the photometry filters
			ind = np.where((wl < phot_wl[i]*(1+1./res)) & (wl > phot_wl[i]*(1-1./res)))
			if len(ind[0]) != 0:
				phot_flux[i] = np.nanmean(flux[ind])
		else:
			# apply the filter function
			# decide the filter name
			if phot_wl[i] == 70:
				fil_name = 'Herschel PACS 70um'
			elif phot_wl[i] == 100:
				fil_name = 'Herschel PACS 100um'
			elif phot_wl[i] == 160:
				fil_name = 'Herschel PACS 160um'
			elif phot_wl[i] == 250:
				fil_name = 'Herschel SPIRE 250um'
			elif phot_wl[i] == 350:
				fil_name = 'Herschel SPIRE 350um'
			elif phot_wl[i] == 500:
				fil_name = 'Herschel SPIRE 500um'

			filter_func = phot_filter(fil_name)

			# trim the filter function
			if phot_wl[i] in [70.,100.,160.]:
				filter_func = filter_func[(filter_func['wave']/1e4 >= max(54.8,min(wl)))*((filter_func['wave']/1e4 <= 95.05)+(filter_func['wave']/1e4 >=103))*(filter_func['wave']/1e4 <= min(190.31,max(wl)))]
			elif phot_wl[i] in [250.,350.,500.]:
				filter_func = filter_func[(filter_func['wave']/1e4 >= 195) & (filter_func['wave']/1e4 <= 670)]

			f = interp1d(wl, flux)
			# print fil_name
			# print filter_func['wave']/1e4
			# print min(wl), max(wl)
			phot_flux[i] = np.trapz(filter_func['wave']/1e4, f(filter_func['wave']/1e4)*filter_func['transmission'])/np.trapz(filter_func['wave']/1e4, filter_func['transmission'])

	return phot_wl, phot_flux

no_match_obj = []
# for determining which one has smaller std
total_phot = []
delta_spec_phot = []
delta_archival_spec_phot = []
# pdb.set_trace()
for i in range(len(objlist)):
	print 'Processing %s' % objlist[i]
	# get the photometry data
	phot_dum = phot_reader(photdir, objlist[i], spire=True)

	# get the archival data
	# read in SLWC3 and SSWD4 spectra and combine them with the same trimming scheme in the CDF archive
	if not (os.path.exists(CDFdir+objlist[i]+'/spire/advanced_products/cube/'+objlist[i]+'_SLWC3_continuum.txt') and \
			os.path.exists(CDFdir+objlist[i]+'/spire/advanced_products/cube/'+objlist[i]+'_SSWD4_continuum.txt')):
		CDF_slw = ascii.read(CDFdir+objlist[i]+'/spire/data/cube/'+objlist[i]+'_SLWC3.txt',\
					header_start=None, data_start=1, names=['wave','intensity'], fill_values=('NaN',np.nan))
		CDF_ssw = ascii.read(CDFdir+objlist[i]+'/spire/data/cube/'+objlist[i]+'_SSWD4.txt',\
					header_start=None, data_start=1, names=['wave','intensity'], fill_values=('NaN',np.nan))
		print '--- either SLWC3 or SSWD4 continuum not found, original data is read instead.'
	else:
		CDF_slw = ascii.read(CDFdir+objlist[i]+'/spire/advanced_products/cube/'+objlist[i]+'_SLWC3_continuum.txt',\
					header_start=None, data_start=1, names=['wave','intensity'], fill_values=('NaN',np.nan))
		CDF_ssw = ascii.read(CDFdir+objlist[i]+'/spire/advanced_products/cube/'+objlist[i]+'_SSWD4_continuum.txt',\
					header_start=None, data_start=1, names=['wave','intensity'], fill_values=('NaN',np.nan))
	# trim and combine two spectra
	spec_archival = {'wave': np.hstack((CDF_ssw['wave'][(CDF_ssw['wave'] >= 195.) & (CDF_ssw['wave'] <= 310. )], \
								   CDF_slw['wave'][CDF_slw['wave'] > 310.])),
				'flux': np.hstack((CDF_ssw['intensity'][(CDF_ssw['wave'] >= 195.) & (CDF_ssw['intensity'] <= 310. )]*SSW_beam, \
								   		CDF_slw['intensity'][CDF_slw['wave'] > 310.]*SLW_beam))}

	archival_phot_wl, archival_phot_flux = herschel_spec_phot(spec_archival['wave'], spec_archival['flux'], filter_func=True)

	# get the spectra photometry from the CDF archive 
	if not os.path.exists(CDFdir+objlist[i]+'/spire/advanced_products/'+objlist[i]+'_spire_corrected_continuum.txt'):
		spec_CDF = ascii.read(CDFdir+objlist[i]+'/spire/data/'+objlist[i]+'_spire_corrected.txt',\
							  header_start=None, data_start=1, names=['wave','flux'])
		print '--- the continuum of extended-corrected spectrum is not found, original data is read instead'
	else:
		spec_CDF = ascii.read(CDFdir+objlist[i]+'/spire/advanced_products/'+objlist[i]+'_spire_corrected_continuum.txt',\
							  header_start=None, data_start=1, names=['wave','flux'])
	# get the filter-convolved spectro-photometry
	spec_phot_wl, spec_phot_flux = herschel_spec_phot(spec_CDF['wave'].data, spec_CDF['flux'].data, filter_func=True)

	# find the matched photometry and store them into data_dict
	mutual_wl = []
	phot = []
	spec_phot = []
	archival_spec_phot = []

	for wl in phot_dum['wave']:
		if (wl in spec_phot_wl) and (wl in archival_phot_wl):
			mutual_wl.append(wl)
			phot.append(float(phot_dum['flux'][phot_dum['wave'] == wl]))
			spec_phot.append(float(spec_phot_flux[spec_phot_wl == wl]))
			archival_spec_phot.append(float(archival_phot_flux[archival_phot_wl == wl]))

			total_phot.append(float(phot_dum['flux'][phot_dum['wave'] == wl]))
			delta_spec_phot.append(float(spec_phot_flux[spec_phot_wl == wl])-float(phot_dum['flux'][phot_dum['wave'] == wl]))
			delta_archival_spec_phot.append(float(archival_phot_flux[archival_phot_wl == wl])-float(phot_dum['flux'][phot_dum['wave'] == wl]))
	if len(mutual_wl) != 0:
		data_dict['object'].append(np.array(objlist[i]))
		data_dict['phot'].append(np.array(phot))
		data_dict['spec_phot'].append(np.array(spec_phot))
		data_dict['archival_spec_phot'].append(np.array(archival_spec_phot))
		data_dict['wave'].append(np.array(mutual_wl))
	else:
		print '--- No matched photometry data.  Object not included!'
		no_match_obj.append(objlist[i])

from pprint import pprint
print no_match_obj
print 'number of objects: ', len(data_dict['object'])-len(no_match_obj)
# calculate the standard deviation of two products
total_phot = np.array(total_phot)
delta_spec_phot = np.array(delta_spec_phot)
delta_archival_spec_phot = np.array(delta_archival_spec_phot)
mean_phot = np.mean(total_phot)
std_spec_phot = np.std(delta_spec_phot)/mean_phot
std_archival_spec_phot = np.std(delta_archival_spec_phot)/mean_phot

# print data_dict['phot'][:]

# plot!

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# phot = np.array([])
# spec_phot = np.array([])
phot = []
spec_phot = []
archival_spec_phot = []

for i in range(len(data_dict['object'])):
	if diff:
		# 250 um
		cdf, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 250.], \
					   data_dict['spec_phot'][i][data_dict['wave'][i] == 250.], 'o', color='Blue', mec='None', alpha=0.7)
		archiv, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 250.], \
						  data_dict['archival_spec_phot'][i][data_dict['wave'][i] == 250.], 'o', color='Red', mec='None', alpha=0.7)
		# 350 um
		cdf_350, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 350.], \
					   data_dict['spec_phot'][i][data_dict['wave'][i] == 350.], 's', color='Blue', mec='None', alpha=0.7)
		archiv_350, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 350.], \
						  data_dict['archival_spec_phot'][i][data_dict['wave'][i] == 350.], 's', color='Red', mec='None', alpha=0.7)
		# 500 um
		cdf_500, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 500.], \
					   data_dict['spec_phot'][i][data_dict['wave'][i] == 500.], '^', color='Blue', mec='None', alpha=0.7)
		archiv_500, = ax.plot(data_dict['phot'][i][data_dict['wave'][i] == 500.], \
						  data_dict['archival_spec_phot'][i][data_dict['wave'][i] == 500.], '^', color='Red', mec='None', alpha=0.7)
	else:
		cdf, = ax.plot(data_dict['phot'][i], \
					   data_dict['spec_phot'][i], 'o', color='Blue', mec='None', alpha=0.7)
		archiv, = ax.plot(data_dict['phot'][i], \
						  data_dict['archival_spec_phot'][i], 'o', color='Red', mec='None', alpha=0.7)
	# print data_dict['phot'][i], data_dict['spec_phot'][i], data_dict['archival_spec_phot'][i]
	if (np.mean(np.array(data_dict['phot'][i])/np.array(data_dict['spec_phot'][i])) > 10.) or (np.mean(np.array(data_dict['phot'][i])/np.array(data_dict['spec_phot'][i])) < 0.1):
		continue
	phot.extend(data_dict['phot'][i])
	spec_phot.extend(data_dict['spec_phot'][i])
	archival_spec_phot.extend(data_dict['archival_spec_phot'][i])
	# np.hstack((phot, np.array(data_dict['phot'][i])))
	# np.hstack((spec_phot, np.array(data_dict['spec_phot'][i])))

# fit the cdf-only spectrophotometric data
phot = np.array(phot)
spec_phot = np.array(spec_phot)
archival_spec_phot = np.array(archival_spec_phot)

fit_para = np.polyfit(np.log10(phot), np.log10(spec_phot), 1)
fit_para_arc = np.polyfit(np.log10(phot), np.log10(archival_spec_phot), 1)
cdf_fit = fit_para[0]*np.log10(phot) + fit_para[1]
arc_fit = fit_para_arc[0]*np.log10(phot) + fit_para_arc[1]
print fit_para
print fit_para_arc

# cdf, = ax.plot(phot, spec_phot, 'o', color='Blue', mec='None', alpha=0.7)
fit, = ax.plot(phot, 10**cdf_fit, color='Blue', alpha=0.7, linewidth=1.5)
fit_arc, = ax.plot(phot, 10**arc_fit, color='Red', alpha=0.7, linewidth=1.5)
equal, = ax.plot([min(phot), max(phot)], [min(phot), max(phot)], '-', color='k', linewidth=1.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.3,1000])
ax.set_ylim([0.3,1000])

ax.legend([cdf, archiv, fit, fit_arc], [r'$\rm{CDF\,(\sigma/<F_{phot.}>=%2.2f)}$' % std_spec_phot, \
	r'$\rm{HSA\,(HIPE\,11)\,(\sigma/<F_{phot.}>=%2.2f)}$' % std_archival_spec_phot, r'$\rm{CDF\,fit}$', r'$\rm{HSA\,fit}$'],\
	numpoints=1, fontsize=14, loc='upper left', framealpha=0.5)
ax.set_xlabel(r'$\rm{log(F_{photometry})\,[Jy]}$', fontsize=18)
ax.set_ylabel(r'$\rm{log(F_{spec.\,phot})\,[Jy]}$', fontsize=18)
[ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
ax.minorticks_on()
ax.tick_params('both',labelsize=18,width=1.5,which='major',pad=10,length=5)
ax.tick_params('both',labelsize=18,width=1.5,which='minor',pad=10,length=2.5)

fig.savefig('/Users/yaolun/test/phot_com_spire.pdf', format='pdf', dpi=300, bbox_inches='tight')
fig.clf()

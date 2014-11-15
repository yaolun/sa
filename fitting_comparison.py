def read_fitting(filepath,noiselevel):
    import numpy as np
    import os
    home = os.path.expanduser('~')
    data = np.genfromtxt(home+filepath, dtype=None)
    header = data[0]
    data = data[1:,:]
    #data = data.astype(float)
    detection = np.empty(len(data[0,:]))
    i=0
    while i != len(data[:,0]):
        if float(data[i,9])< noiselevel:
            data = np.delete(data, i, 0)
        else:
            i += 1
    name = data[:,0]
    data = data[:,1:].astype(float)
    return data, header, name


#different data stored in column as the following order:
# wl, sig_wl, str, sig_str, fwhm, sig_fwhm, base_str, snr, e_u, A, g, ra, dec


#This part is for comparing the flux difference among different version of reduction
import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')
noiselevel=3
# for i in range(1,26):
#     filepath = '/bhr71/data/pacs_pixel'+str(i)+'_lines.txt'
#     [data, header] = read_fitting(filepath, noiselevel)
#     file_com = '/bhr71/data/pacs_pixel'+str(i)+'_lines copy.txt'
#     [data_com, header_com] = read_fitting(file_com, noiselevel)
#     if len(data[:,0]) > 0:
#         result, = plt.plot(data[:,0], data[:,2]/1e-22, 'bo')
#         plt.errorbar(data[:,0], data[:,2]/1e-22, yerr=data[:,3]/1e-22, linestyle='None', color = 'blue')
#         #snr
#         #result, = plt.plot(data[:,0], data[:,7],'bo')
#     if len(data_com[:,0]) > 0:
#         ref, = plt.plot(data_com[:,0], data_com[:,2]/1e-22, 'ro')
#         plt.errorbar(data_com[:,0], data_com[:,2]/1e-22, yerr=data_com[:,3]/1e-22, linestyle='None', color='red')
#         #snr
#         #ref, = plt.plot(data_com[:,0], data_com[:,7], 'ro')
#     lg = plt.legend([result, ref],['V65(os8 sf7) '+str(len(data[:,0]))+' lines','os2 '+str(len(data_com[:,0]))+' lines'], loc='upper right')
#     plt.xlabel('Wavelength ($\mu$m)')
#     plt.ylabel('Flux (10$^{-22}$ W/cm$^{2}$)')
#     #plt.ylabel('S/N')
#     plt.xlim([50,210])
#     plt.savefig(home+'/bhr71/plots/flux_comparison/pix'+str(i)+'.eps',format='eps',dpi=300)
#     #plt.savefig(home+'/bhr71/plots/flux_comparison/snr_pix'+str(i)+'.eps',format='eps',dpi=300)
#     plt.clf()
#     plt.cla()
#     print 'Finish pixel '+str(i)


###########
#  9Spax  #
###########
#Fixed width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/central9Spaxels_lines_localbaseline_fixwidth_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([9.96,0.0,0.0,0.67,0.75,0.42,0.6,0.67,0.19,0.88,0.77,0.3])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
# for line in wl_ref:
#     ind = (np.abs(data[:,0].T-line)).argmin()
#     print line, data[ind,0]
#     wl.append(data[ind,0])
#     flux.append(data[ind,2]/1e-20)
#     flux_sig.append(data[ind,3]/1e-20)
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T
print flux[snr>=3]
result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/central9Spaxels_localbaseline_fixwidth_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 9Spax fixed width'

#flexible width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/central9Spaxels_lines_localbaseline_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([9.96,0.0,0.0,0.67,0.75,0.42,0.6,0.67,0.19,0.88,0.77,0.3])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
# for line in wl_ref:
#     ind = (np.abs(data[:,0].T-line)).argmin()
#     print line, data[ind,0]
#     wl.append(data[ind,0])
#     flux.append(data[ind,2]/1e-20)
#     flux_sig.append(data[ind,3]/1e-20)
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T
result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/central9Spaxels_localbaseline_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 9Spax'

##########
# 3x3yes #
##########
#Fixed width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedYES_lines_localbaseline_fixwidth_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([7.88,0.30,0.43,0.84,0.88,0.35,0.58,0.88,0.25,0.73,0.46,0.44])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T

result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedYES_localbaseline_fixwidth_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 3x3yes fixed width'

#Flexible width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedYES_lines_localbaseline_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([7.88,0.30,0.43,0.84,0.88,0.35,0.58,0.88,0.25,0.73,0.46,0.44])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T

result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedYES_localbaseline_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 3x3yes'

######################
#1x1(central no corr)#
######################
#Fixed width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedNO_lines_localbaseline_fixwidth_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([7.28,0.21,0.34,0.71,0.74,0.27,0.45,0.63,0.18,0.49,0.30,0.29])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T
result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedNO_localbaseline_fixwidth_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 1x1 (central Spx no corr) fixed width'

#Flexible width
[data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedNO_lines_localbaseline_global_noise.txt',0)
wl_ref = np.array([63.1840,81.8060,82.0315,84.4110,84.6,108.0732,108.7630,124.1930,125.3537,145.5250,157.74,179.5267])
flux_ref = np.array([7.28,0.21,0.34,0.71,0.74,0.27,0.45,0.63,0.18,0.49,0.30,0.29])
name_ref = np.array(['OI3P1-3P2','CO32-31','o-H2O6_16-5_05','CO31-30','OH9-3','o-H2O2_21-1_10','CO24-23','CO21-20','p-H2O4_04-3_13','OI3P0-3P1','CII2P3_2-2P1_2','o-H2O2_12-1_01'])
wl = []
flux = []
flux_sig = []
snr = []
for ref in name_ref:
    for i in range(0,len(data[:,0])):
        if name[i] == ref:
            wl.append(data[i,1])
            flux.append(data[i,3]/1e-20)
            flux_sig.append(data[i,4]/1e-20)
            snr.append(data[i,8])
            break
wl = np.array(wl).T
flux = np.array(flux).T
flux_sig = np.array(flux_sig).T
snr = np.array(snr).T
result, = plt.plot(wl[snr>=3], flux_ref[snr>=3]/flux[snr>=3], 'bo')
result_non, = plt.plot(wl[snr<3], flux_ref[snr<3]/flux[snr<3], 'bo',mfc='None',mec='blue')
sigma = flux_ref/flux**2*flux_sig
plt.errorbar(wl, flux_ref/flux, yerr=sigma, linestyle='None', color='blue')
plt.axhline(y=1,color='k',ls='dashed')
plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
plt.ylabel('F$_{agata}$/F$_{yaolun}$',fontsize=14)
plt.title('Mean: '+str(np.mean(flux_ref[snr>=3]/flux[snr>=3]))+' $\sigma$: '+str(np.sum(sigma[snr>=3]**2)**0.5/len(sigma[snr>=3])))
lg = plt.legend([result, result_non],['S/N $\geq$ 3$\sigma$','S/N $<$ 3$\sigma$'],numpoints=1)
plt.xlim([50,200])
plt.ylim([-0.1,2])
plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/centralSpaxels_correctedNO_localbaseline_global_noise.eps', format='eps',dpi=300)
plt.cla()
plt.clf()
print 'Finish 1x1 (central Spx no corr)'
print 'Done'
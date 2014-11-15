import numpy as np
import matplotlib.pyplot as plt
import os
from read_fitting import *
home = os.path.expanduser('~')
noiselevel=3

#indir = '/bhr71/data/BHR71_V65/'
indir = '/tmc1/data/'
#plotdir = '/bhr71/plots/'
plotdir = '/tmc1/plots/'

obj = 'TMC1'

fig_freq = plt.figure()
ax_freq = fig_freq.add_subplot(111)
fig_wl = plt.figure()
ax_wl = fig_wl.add_subplot(111)
fig_wl_diff = plt.figure()
ax_wl_diff = fig_wl_diff.add_subplot(111)
########################
#5x5 summed with jitter#
########################
[data, header, name] = read_fitting(indir+obj+'_lines_localbaseline_fixwidth_global_noise.txt',0)
#Plot the wavelength offest
#For plot in frequency
freq_lab = 2.997e5/data[:,0].T
freq_obs = 2.997e5/data[:,1].T
freq_obs_sig = (freq_obs**2/2.997e5)*data[:,2].T
snr_all = data[:,8].T
cube_summed, = ax_freq.plot(freq_lab[snr_all>=3], freq_lab[snr_all>=3]-freq_obs[snr_all>=3],'bo')
ax_freq.errorbar(freq_lab[snr_all>=3], freq_lab[snr_all>=3]-freq_obs[snr_all>=3], yerr=freq_obs_sig[snr_all>=3],linestyle='None',color='b')
#Calculate the mean and uncertainty of the wavelength offset. And exclude the OI3P1-3P2 line since it is possibly a high velocity gas
ind = np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(freq_obs_sig!=0))
freq_lab = freq_lab[ind]
freq_obs = freq_obs[ind]
freq_all = snr_all[ind]
freq_obs_sig = freq_obs_sig[ind]

mean_low = np.sum((1/freq_obs_sig[freq_lab<3000]**2)*(freq_lab[freq_lab<3000]-freq_obs[freq_lab<3000]))/np.sum(1/freq_obs_sig[freq_lab<3000]**2)
mean_hi = np.sum((1/freq_obs_sig[freq_lab>=3000]**2)*(freq_lab[freq_lab>=3000]-freq_obs[freq_lab>=3000]))/np.sum(1/freq_obs_sig[freq_lab>=3000]**2)
sig_low = (1/np.sum(1/freq_obs_sig[freq_lab<3000]**2))**0.5
sig_hi = (1/np.sum(1/freq_obs_sig[freq_lab>=3000]**2))**0.5

label_cube_summed_freq = '5x5 cube summed (R/B:~%6.4f~+/-~%6.4f,~%6.4f~+/-~%6.4f)' % (mean_low,sig_low,mean_hi,sig_hi)

#wavelength version
wl_lab = data[:,0].T
wl_obs = data[:,1].T
wl_obs_sig = data[:,2].T
snr_all = data[:,8].T
cube_summed_wl, = ax_wl.plot(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3],'bo')
ax_wl.errorbar(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3], yerr=wl_obs_sig[snr_all>=3],linestyle='None',color='b')

ind = np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0))
wl_lab = wl_lab[ind]
wl_obs = wl_obs[ind]
snr_all = snr_all[ind]
wl_obs_sig = wl_obs_sig[ind]

mean_blue = np.sum((1/wl_obs_sig[wl_lab<101]**2)*(wl_lab[wl_lab<101]-wl_obs[wl_lab<101]))/np.sum(1/wl_obs_sig[wl_lab<101]**2)
mean_red = np.sum((1/wl_obs_sig[wl_lab>=101]**2)*(wl_lab[wl_lab>=101]-wl_obs[wl_lab>=101]))/np.sum(1/wl_obs_sig[wl_lab>=101]**2)
sig_blue = (1/np.sum(1/wl_obs_sig[wl_lab<101]**2))**0.5
sig_red = (1/np.sum(1/wl_obs_sig[wl_lab>=101]**2))**0.5

label_cube_summed_wl = 'B/R:~%6.4f~+/-~%6.4f,~%6.4f~+/-~%6.4f' % (mean_blue,sig_blue,mean_red,sig_red)

#Do the calculation of difference of wavelength fitted in different reduction
name_5x5 = name[ind]
wl_5x5 = wl_obs
wl_sig_5x5 = wl_obs_sig
###########
#  9Spax  #
###########
[data, header, name] = read_fitting(indir+'central9Spaxels_lines_localbaseline_fixwidth_global_noise.txt',0)
#Plot the wavelength offest
#For plot in frequency
freq_lab = 2.997e5/data[:,0].T
freq_obs = 2.997e5/data[:,1].T
freq_obs_sig = (freq_obs**2/2.997e5)*data[:,2].T
snr_all = data[:,8].T

central9, = ax_freq.plot(freq_lab[snr_all>=3], freq_lab[snr_all>=3]-freq_obs[snr_all>=3],'ro')
ax_freq.errorbar(freq_lab[snr_all>=3], freq_lab[snr_all>=3]-freq_obs[snr_all>=3], yerr=freq_obs_sig[snr_all>=3],linestyle='None',color='r')
#Calculate the mean and uncertainty of the wavelength offset. And exclude the OI3P1-3P2 line since it is possibly a high velocity gas
ind = np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(freq_obs_sig!=0))
freq_lab = freq_lab[ind]
freq_obs = freq_obs[ind]
freq_all = snr_all[ind]
freq_obs_sig = freq_obs_sig[ind]

mean_low = np.sum((1/freq_obs_sig[freq_lab<3000]**2)*(freq_lab[freq_lab<3000]-freq_obs[freq_lab<3000]))/np.sum(1/freq_obs_sig[freq_lab<3000]**2)
mean_hi = np.sum((1/freq_obs_sig[freq_lab>=3000]**2)*(freq_lab[freq_lab>=3000]-freq_obs[freq_lab>=3000]))/np.sum(1/freq_obs_sig[freq_lab>=3000]**2)
sig_low = (1/np.sum(1/freq_obs_sig[freq_lab<3000]**2))**0.5
sig_hi = (1/np.sum(1/freq_obs_sig[freq_lab>=3000]**2))**0.5

label_central9_freq = 'Central9Spaxels~~~~(R/B:~%6.4f~+/-~%6.4f,~%6.4f~+/-~%6.4f)' % (mean_low,sig_low,mean_hi,sig_hi)

#wavelength version
wl_lab = data[:,0].T
wl_obs = data[:,1].T
wl_obs_sig = data[:,2].T
snr_all = data[:,8].T
central9_wl, = ax_wl.plot(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3],'ro')
ax_wl.errorbar(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3], yerr=wl_obs_sig[snr_all>=3],linestyle='None',color='r')

ind = np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0))
wl_lab = wl_lab[ind]
wl_obs = wl_obs[ind]
snr_all = snr_all[ind]
wl_obs_sig = wl_obs_sig[ind]

mean_blue = np.sum((1/wl_obs_sig[wl_lab<101]**2)*(wl_lab[wl_lab<101]-wl_obs[wl_lab<101]))/np.sum(1/wl_obs_sig[wl_lab<101]**2)
mean_red = np.sum((1/wl_obs_sig[wl_lab>=101]**2)*(wl_lab[wl_lab>=101]-wl_obs[wl_lab>=101]))/np.sum(1/wl_obs_sig[wl_lab>=101]**2)
sig_blue = (1/np.sum(1/wl_obs_sig[wl_lab<101]**2))**0.5
sig_red = (1/np.sum(1/wl_obs_sig[wl_lab>=101]**2))**0.5

label_central9_wl = 'B/R:~%6.4f~+/-~%6.4f,~%6.4f~+/-~%6.4f' % (mean_blue,sig_blue,mean_red,sig_red)

#Do the calculation of difference of wavelength fitted in different reduction
name_9spx = name[ind]
wl_9spx = wl_obs
wl_sig_9spx = wl_obs_sig
#Find the same lines in two reductions
wl = []
wl_diff = []
wl_sig_diff = []
for ref in name_5x5:
    for i in range(0,len(name_9spx)):
        if name_9spx[i] == ref:
            wl.append(float(wl_5x5[name_5x5==ref]))
            wl_diff.append(float(wl_5x5[name_5x5==ref]-wl_9spx[i]))
            wl_sig_diff.append(float((wl_sig_5x5[name_5x5==ref]**2+wl_sig_9spx[i]**2)**0.5))
            break
#Plot the wavelength difference between two reductions and the wavelength
wl = np.array(wl)
wl_diff = np.array(wl_diff)
wl_sig_diff = np.array(wl_sig_diff)

ax_wl_diff.plot(wl,wl_diff,'o',color='DarkGreen')
ax_wl_diff.errorbar(wl,wl_diff,yerr=wl_sig_diff,linestyle='None',color='DarkGreen')
ax_wl_diff.axhline(y=0,color='k',ls='dashed')
ax_wl_diff.set_ylim([-0.06,0.01])
ax_wl_diff.set_xlabel(r'$\mathrm{Wavelength~(\mu m)}$',fontsize=14)
ax_wl_diff.set_ylabel(r'$\mathrm{\lambda_{5x5,summed}-\lambda_{9Spx}~(\mu m)}$',fontsize=14)
fig_wl_diff.savefig(home+plotdir+'wl_diff_5x5-9spx.eps',format='eps',dpi=300)

ax_freq.legend([cube_summed,central9],[label_cube_summed_freq,label_central9_freq],numpoints=1,loc='lower center')
ax_wl.legend([cube_summed_wl,central9_wl],[label_cube_summed_wl,label_central9_wl],numpoints=1,loc='upper center')

ax_freq.axhline(y=0,color='k',ls='dashed')
ax_freq.set_xlabel('Frequency (GHz)')
ax_freq.set_ylabel('Frequency Offset (GHz)')
ax_freq.set_ylim([-3,5])
fig_freq.savefig(home+plotdir+'freq_offset_comparison.eps',format='eps',dpi=300)

ax_wl.axhline(y=0,color='k',ls='dashed')
ax_wl.set_xlabel(r'$\mathrm{Wavelength~(\mu m)}$',fontsize=14)
ax_wl.set_ylabel(r'$\mathrm{Wavelength~Offset~(\mu m)}$',fontsize=14)
ax_wl.set_ylim([-0.25,0.15])
fig_wl.savefig(home+plotdir+'wl_offset_comparison.eps',format='eps',dpi=300)

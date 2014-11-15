import numpy as np
import matplotlib.pyplot as plt
import os

def read_fitting(filepath,noiselevel):
    home = os.path.expanduser('~')
    data = np.genfromtxt(home+filepath, dtype=None)
    header = data[0]
    data = data[1:,:]
    #data = data.astype(float)
    detection = np.empty(len(data[0,:]))
    i=0
    while i != len(data[:,0]):
        if float(data[i,8])< noiselevel:
            data = np.delete(data, i, 0)
        else:
            i += 1
    name = data[:,0]
    data = data[:,1:].astype(float)
    return data, header, name

home = os.path.expanduser('~')
noiselevel=3

objname = ['tmr1','tmc1a','iras03301','iras03245','gss30','dkcha','WL12','VLA1623','TMC1','RCrA-IRS7B','L1551-IRS5','L1489','L1448-MM','L1157','Elias29','B335','B1-c','B1-a','l1527','HD100546']#'Serpens-SMM3',
indir = '/data/digit_range_scan_v65/'
for obj in objname:
    [data, header, name] = read_fitting(indir+obj+'/data/'+obj+'_pacs_v65_trim_lines_localbaseline_fixwidth_global_noise.txt',0)
    wl_lab = data[:,0].T
    wl_obs = data[:,1].T
    wl_obs_sig = data[:,2].T
    snr_all = data[:,8].T
    wloffset = plt.figure()
    plt.plot(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3],'bo')
    plt.errorbar(wl_lab[snr_all>=3], wl_lab[snr_all>=3]-wl_obs[snr_all>=3], yerr=wl_obs_sig[snr_all>=3],linestyle='None')
    #Calculate the mean and uncertainty of the wavelength offset. And exclude the OI3P1-3P2 line since it is possibly a high velocity gas
    wl_lab = data[np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0)),0].T
    wl_obs = data[np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0)),1].T
    wl_obs_sig = data[np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0)),2].T
    snr_all = data[np.nonzero((name!='OI3P1-3P2')&(snr_all >=3)&(wl_obs_sig!=0)),8].T
    mean_blue = np.sum((1/wl_obs_sig[wl_lab<101]**2)*(wl_lab[wl_lab<101]-wl_obs[wl_lab<101]))/np.sum(1/wl_obs_sig[wl_lab<101]**2)
    mean_red = np.sum((1/wl_obs_sig[wl_lab>=101]**2)*(wl_lab[wl_lab>=101]-wl_obs[wl_lab>=101]))/np.sum(1/wl_obs_sig[wl_lab>=101]**2)
    sig_blue = (1/np.sum(1/wl_obs_sig[wl_lab<101]**2))**0.5
    sig_red = (1/np.sum(1/wl_obs_sig[wl_lab>=101]**2))**0.5
    plt.title('Blue/Red: %8.6f +/- %8.6f, %8.6f +/- %8.6f' % (mean_blue,sig_blue,mean_red,sig_red))
    plt.axhline(y=0,color='k',ls='dashed')
    #print 'Mean Blue: %8.4f' % mean_blue
    #print 'Mean Red : %8.4f' % mean_red
    plt.xlabel('Lab Wavelength ($\mu$m)')
    plt.ylabel('Wavelength Offset ($\mu$m)')
    plt.xlim([50,200])
    plt.savefig(home+indir+obj+'/plots/'+obj+'_wl_offset_localbaseline_fixwidth_global_noise.eps',format='eps',dpi=300)
    plt.cla()
    plt.clf()
    print 'Finish'+ obj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')

objname = 'bhr71'
#filename = ['pacs_pixel13_lines.txt','SSWD4lines.txt','SLWC3lines.txt']
filename = []
for j in range(1, 26):
    filename.append('pacs_pixel'+str(j)+'_lines.txt')
name_tot = []
wl_tot = []
wl_sig_tot = []
flux_tot = []
flux_sig_tot = []
fwhm_tot = []
fwhm_sig_tot = []
pix = 0
for i in filename:
    pix = 1+pix
    #[name, wl, wl_sig, flux, flux_sig, fwhm, fwh m_sig, basestr, snr, E_u, A, g, ra, dec] = np.loadtxt(home+'/'+objname+'/data/'+i).T
    data = open(home+'/'+objname+'/data/'+i).readlines()
    name = []
    wl = []
    wl_sig = []
    flux = []
    flux_sig = []
    fwhm = []
    fwhm_sig = []
    snr = []
    data.pop(0)
    wl_below = []
    wl_sig_below = []
    flux_below = []
    flux_sig_below = []
    fwhm_below = []
    fwhm_sig_below = []
    snr_below = []
    for line in data:
        t = line.split()
        #print t
        if len(t) > 1:
            if float(t[8]) >= 3:
                name.append(t[0])                    
                wl.append(float(t[1]))
                wl_sig.append(float(t[2]))
                flux.append(float(t[3])/1e-22)
                flux_sig.append(float(t[4])/1e-22)
                fwhm.append(float(t[5]))
                fwhm_sig.append(float(t[6]))
                snr.append(float(t[8]))
            else:
                wl_below.append(float(t[1]))
                wl_sig_below.append(float(t[2]))
                flux_below.append(float(t[3])/1e-22)
                flux_sig_below.append(float(t[4])/1e-22)
                fwhm_below.append(float(t[5]))
                fwhm_sig_below.append(float(t[6]))
                snr_below.append(float(t[8]))
    name_tot.extend(name)
    wl_tot.extend(wl) 
    wl_sig_tot.extend(wl_sig)
    flux_tot.extend(flux)
    flux_sig_tot.extend(flux_sig)
    fwhm_tot.extend(fwhm)
    fwhm_sig_tot.extend(fwhm_sig)
    #Read the theoretical fwhm value
    #order1 = open(home+'/'+objname+'/data/spectralresolution_order1.txt','r').readlines()
    #order2 = open(home+'/'+objname+'/data/spectralresolution_order2.txt','r').readlines()
    #order3 = open(home+'/'+objname+'/data/spectralresolution_order3.txt','r').readlines()
    [wl_order1, res1] = np.loadtxt(home+'/'+objname+'/data/spectralresolution_order1.txt').T
    [wl_order2, res2] = np.loadtxt(home+'/'+objname+'/data/spectralresolution_order2.txt').T
    [wl_order3, res3] = np.loadtxt(home+'/'+objname+'/data/spectralresolution_order3.txt').T
    #make a plot
    if len(wl) > 0:
        result, = plt.plot(wl, fwhm, 'rx')
        plt.errorbar(wl, fwhm, yerr=fwhm_sig, linestyle='None',color='red')
        band1, = plt.plot(wl_order1, wl_order1/res1, 'b-')
        band2, = plt.plot(wl_order2, wl_order2/res2, 'g-')
        band3, = plt.plot(wl_order3, wl_order3/res3, 'y-')
        lg = plt.legend([result, band1, band2, band3], ['Data','order1','order2','order3'], loc='upper right')
    if len(wl) == 0:
        band1, = plt.plot(wl_order1, wl_order1/res1, 'b-')
        band2, = plt.plot(wl_order2, wl_order2/res2, 'g-')
        band3, = plt.plot(wl_order3, wl_order3/res3, 'y-')
        lg = plt.legend([band1, band2, band3], ['order1','order2','order3'], loc='upper right')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('FWHM ($\mu$m)')
    plt.xlim([50, 210])
    plt.ylim([0,0.5])
    plt.text(60, 0.4, str(len(wl))+' lines detected')
    plt.savefig(home+'/'+objname+'/plots/fwhm_resolution_comparison_pacs_pixel'+str(pix)+'_latest.eps', format='eps', dpi=300)
    print 'Finish pixel '+str(pix)+' comparison. There are '+str(len(wl))+' lines detected.'
    plt.cla()
    plt.clf()

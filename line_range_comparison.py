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
import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')
#range scan and line scan comparison
#TMC-1
reduction = ['_central9Spaxels_PointSourceCorrected_slice_0','_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_0','_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_0']
reduction_r = ['central9Spaxels','centralSpaxels_correctedNO','centralSpaxels_correctedYES']
aor_name = ['1342225830_TMC_1_', '1342225831_TMC_1_', '1342225830_TMC_1_', '1342225831_TMC_1_']
band = ['blue','blue','red','red']
suffix = '_os8sf7_lines_LS_localbaseline_global_noise'
line_slice = [[0,2],[0,5],[1,3],[0,7]]
for i_rec in range(0, len(reduction)):
    name_l = np.array([])
    wl_l = np.array([])
    flux_l = np.array([])
    flux_sig_l = np.array([])
    snr_l = np.array([])
    for iband in range(0, len(band)):
        for i in range(line_slice[iband][0], line_slice[iband][1]+1):
            filename = 'linescan_OBSID_' + aor_name[iband] + band[iband] + reduction[i_rec] + str(i) + suffix
            [data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+filename+'.txt', 1)
            name_l = np.append(name_l, name[:])
            wl_l = np.append(wl_l, data[:,0])
            flux_l = np.append(flux_l, data[:,2])
            flux_sig_l = np.append(flux_sig_l, data[:,3])
            snr_l = np.append(snr_l, data[:,7])
    [data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+reduction_r[i_rec]+'_lines_localbaseline_global_noise.txt', 1)
    wl = np.array(data[:,0]).T
    flux = np.array(data[:,2]).T
    flux_sig = np.array(data[:,3]).T
    snr = np.array(data[:,7]).T
    #Find the corresponding lines
    ind = np.arange(0, len(name)-1)
    ind_linr = []
    ind_rinl = []
    for i in range(0, len(name_l)):
        for j in ind:
            if name_l[i] == name[j]:
                #ind = np.delete(ind, j)
                ind = ind[ind != j]
                ind_linr.append(j)
                ind_rinl.append(i)
                break
    wl_l = wl_l[ind_rinl]
    flux_l = flux_l[ind_rinl]
    flux_sig_l = flux_sig_l[ind_rinl]
    snr_l = snr_l[ind_rinl]
    wl = wl[ind_linr]
    flux = flux[ind_linr]
    flux_sig = flux_sig[ind_linr]
    snr = snr[ind_linr]
    #Plot 
    #range_scan, = plt.semilogy(wl[ind_linr], flux[ind_linr], 'ro', markersize=4)
    #line_scan, = plt.semilogy(wl_l[ind_rinl], flux_l[ind_rinl], 'bo', markersize=4)
    sigma = ((flux_l/flux)**2*((flux_sig_l/flux_l)**2+(flux_sig/flux)**2))**0.5
    ratio, = plt.plot(wl_l[snr_l >= 3], flux_l[snr_l >= 3]/flux[snr_l >= 3], 'bo', markersize=5,mec='blue')
    plt.errorbar(wl_l[snr_l >= 3], flux_l[snr_l >= 3]/flux[snr_l >= 3], yerr=sigma[snr_l >= 3], linestyle='None', color='blue')
    ratio_non, = plt.plot(wl_l[snr_l < 3], flux_l[snr_l < 3]/flux[snr_l < 3], 'bo', markersize=5,mfc='None',mec='blue')
    plt.errorbar(wl_l[snr_l < 3], flux_l[snr_l < 3]/flux[snr_l < 3], yerr=sigma[snr_l < 3], linestyle='None',color='blue')
    plt.axhline(y=1,color='k',ls='dashed')
    ax = plt.gca()
    #ax.set_yscale('log')
    #plt.errorbar(wl_l, flux_l, yerr=flux_sig_l, linestyle='None', color='blue')
    #plt.errorbar(wl, flux, yerr=flux_sig, linestyle='None', color='red')
    #lg = plt.legend([line_scan, range_scan], ['line scan','range scan'], loc='upper right', numpoints=1)
    plt.title('Mean: '+str(np.mean(flux_l/flux))+' $\sigma$: '+str(np.sum(sigma**2)**0.5/len(sigma)))
    plt.xlabel('Wavelength ($\mu$m)',fontsize=16)
    plt.ylabel('F$_{line}$/F$_{range}$',fontsize=16)
    plt.ylim([-1,5])
    #plt.ylabel('Flux')
    plt.savefig(home+'/line_fitting_comparison/fitting_result/tmc1/'+reduction_r[i_rec]+'.eps',format='eps',dpi=300)
    plt.cla()
    plt.clf()
    print reduction_r[i_rec]
    print 'There are '+str(len(wl_l))+' lines,',str(len(wl_l[snr_l < 3]))+' of them are below 3 sigma'
print 'Finish range scan vs line scan comparison'

for rec in reduction_r:
    [data, header, name] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+rec+'_lines.txt',1)
    [data_lb, header, name_lb] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+rec+'_lines_localbaseline.txt', 1)
    [data_lb_fw, header, name_lb_fw] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+rec+'_lines_localbaseline_fixwidth.txt',1)
    [data_lb_gn, header, name_lb_gn] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+rec+'_lines_localbaseline_global_noise.txt',1)
    [data_lb_fw_gn, header, name_lb_fw_gn] = read_fitting('/line_fitting_comparison/fitting_result/tmc1/'+rec+'_lines_localbaseline_fixwidth_global_noise.txt',1)
    
    [wl, flux, flux_sig, snr] = [data[:,0].T, data[:,2].T, data[:,3].T, data[:,7].T]
    [wl_lb, flux_lb, flux_sig_lb, snr_lb] = [data_lb[:,0].T, data_lb[:,2].T, data_lb[:,3].T, data_lb[:,7].T]
    [wl_lb_fw, flux_lb_fw, flux_sig_lb_fw, snr_lb_fw] = [data_lb_fw[:,0].T, data_lb_fw[:,2].T, data_lb_fw[:,3].T, data_lb_fw[:,7].T]
    [wl_lb_gn, flux_lb_gn, flux_sig_lb_gn, snr_lb_gn] = [data_lb_gn[:,0].T, data_lb_gn[:,2].T, data_lb_gn[:,3].T, data_lb_gn[:,7].T]
    [wl_lb_fw_gn, flux_lb_fw_gn, flux_sig_lb_fw_gn, snr_lb_fw_gn] = [data_lb_fw_gn[:,0].T, data_lb_fw_gn[:,2].T, data_lb_fw_gn[:,3].T, data_lb_fw_gn[:,7].T]

    data = np.array(zip(wl[wl.argsort()], flux[wl.argsort()], flux_sig[wl.argsort()], snr[wl.argsort()]))
    data_lb = np.array(zip(wl_lb[wl_lb.argsort()], flux_lb[wl_lb.argsort()], flux_sig_lb[wl_lb.argsort()], snr_lb[wl_lb.argsort()]))
    data_lb_fw = np.array(zip(wl_lb_fw[wl_lb_fw.argsort()], flux_lb_fw[wl_lb_fw.argsort()], flux_sig_lb_fw[wl_lb_fw.argsort()], snr_lb_fw[wl_lb_fw.argsort()]))
    data_lb_gn = np.array(zip(wl_lb_gn[wl_lb_gn.argsort()], flux_lb_gn[wl_lb_gn.argsort()], flux_sig_lb_gn[wl_lb_gn.argsort()], snr_lb_gn[wl_lb_gn.argsort()]))
    data_lb_fw_gn = np.array(zip(wl_lb_fw_gn[wl_lb_fw_gn.argsort()], flux_lb_fw_gn[wl_lb_fw_gn.argsort()], flux_sig_lb_fw_gn[wl_lb_fw_gn.argsort()], snr_lb_fw_gn[wl_lb_fw_gn.argsort()]))

    #Plot
    #Compare the effect of localbaseline
    plt.figure()
    ax = plt.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')
    normal, = plt.plot(data[snr >= 3,0],data[snr >=3,1],'bo',markersize=4, mec='None')
    plt.errorbar(data[snr >= 3,0],data[snr >=3,1], yerr=data[snr>=3,2],color='blue',linestyle='None')
    normal_below, = plt.plot(data[snr<3,0],data[snr<3,1],'o',markersize=4,mfc='None',mec='blue')
    plt.errorbar(data[snr<3,0],data[snr<3,1], yerr=data[snr<3,2], color='blue',linestyle='None')

    lb, = plt.plot(data_lb[snr_lb>=3,0], data_lb[snr_lb>=3,1],'ro',markersize=4,mec='None')
    plt.errorbar(data_lb[snr_lb>=3,0], data_lb[snr_lb>=3,1], yerr=data_lb[snr_lb>=3,2],color='red',linestyle='None')
    lb_below, = plt.plot(data_lb[snr_lb<3,0], data_lb[snr_lb<3,1],'o',markersize=4,mfc='None',mec='red')
    plt.errorbar(data_lb[snr_lb<3,0], data_lb[snr_lb<3,1], yerr=data_lb[snr_lb<3,2],color='red',linestyle='None')

    lg = plt.legend([normal, normal_below, lb, lb_below],['w/o local baseline ('+str(len(data[snr>=3,0]))+')','w/o local baseline $<$ 3 sigma ('+str(len(data[snr<3,0]))+')','local baseline ('+str(len(data_lb[snr_lb>=3,0]))+')','local baseline $<$ 3 sigma ('+str(len(data_lb[snr_lb<3,0]))+')'],numpoints=1,loc='upper left')

    plt.xlabel('Wavelength ($\mu$m)',fontsize=16)
    plt.ylabel('Flux',fontsize=16)
    plt.savefig(home + '/line_fitting_comparison/fitting_result/tmc1/'+rec+'_normal_localbaseline.eps',format='eps',dpi=300)
    plt.cla()
    plt.clf()
    print 'Finish '+rec+' normal vs baseline'

    #Compare the effect of fixed width
    plt.figure()
    ax = plt.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')

    lb, = plt.plot(data_lb[snr_lb>=3,0], data_lb[snr_lb>=3,1],'bo',markersize=4,mec='None')
    plt.errorbar(data_lb[snr_lb>=3,0], data_lb[snr_lb>=3,1], yerr=data_lb[snr_lb>=3,2],color='blue',linestyle='None')
    lb_below, = plt.plot(data_lb[snr_lb<3,0], data_lb[snr_lb<3,1],'o',markersize=4,mfc='None',mec='blue')
    plt.errorbar(data_lb[snr_lb<3,0], data_lb[snr_lb<3,1], yerr=data_lb[snr_lb<3,2],color='blue',linestyle='None')

    lb_fw, = plt.plot(data_lb_fw[snr_lb_fw>=3,0], data_lb_fw[snr_lb_fw>=3,1],'ro',markersize=4,mec='None')
    plt.errorbar(data_lb_fw[snr_lb_fw>=3,0], data_lb_fw[snr_lb_fw>=3,1], yerr=data_lb_fw[snr_lb_fw>=3,2],color='red',linestyle='None')
    lb_fw_below, = plt.plot(data_lb_fw[snr_lb_fw<3,0], data_lb_fw[snr_lb_fw<3,1],'o',markersize=4,mfc='None',mec='red')
    plt.errorbar(data_lb_fw[snr_lb_fw<3,0], data_lb_fw[snr_lb_fw<3,1], yerr=data_lb_fw[snr_lb_fw<3,2],color='red',linestyle='None')

    lg = plt.legend([lb, lb_below, lb_fw, lb_fw_below],['w/o fixed width ('+str(len(data_lb[snr_lb>=3,0]))+')','w/o fixed width $<$ 3 sigma ('+str(len(data_lb[snr_lb<3,0]))+')','fixed width ('+str(len(data_lb_fw[snr_lb_fw>=3,0]))+')','fixed width $<$ 3 sigma ('+str(len(data_lb_fw[snr_lb_fw<3,0]))+')'],numpoints=1,loc='upper left')

    plt.xlabel('Wavelength ($\mu$m)',fontsize=16)
    plt.ylabel('Flux',fontsize=16)
    plt.savefig(home + '/line_fitting_comparison/fitting_result/tmc1/'+rec+'_fixwidth.eps',format='eps',dpi=300)
    plt.cla()
    plt.clf()
    print 'Finish '+rec+' fixed width'

    #Compare the effect of global noise
    plt.figure()
    ax = plt.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')

    lb_fw, = plt.plot(data_lb_fw[snr_lb_fw>=3,0], data_lb_fw[snr_lb_fw>=3,1],'bs',markersize=8,mec='None')
    #plt.errorbar(data_lb_fw[snr_lb_fw>=3,0], data_lb_fw[snr_lb_fw>=3,1], yerr=data_lb_fw[snr_lb_fw>=3,2],color='blue',linestyle='None')
    lb_fw_below, = plt.plot(data_lb_fw[snr_lb_fw<3,0], data_lb_fw[snr_lb_fw<3,1],'s',markersize=8,mfc='None',mec='blue')
    #plt.errorbar(data_lb_fw[snr_lb_fw<3,0], data_lb_fw[snr_lb_fw<3,1], yerr=data_lb_fw[snr_lb_fw<3,2],color='blue',linestyle='None')

    lb_fw_gn, = plt.plot(data_lb_fw_gn[snr_lb_fw_gn>=3,0], data_lb_fw_gn[snr_lb_fw_gn>=3,1],'rD',markersize=6,mec='None')
    #plt.errorbar(data_lb_fw_gn[snr_lb_fw_gn>=3,0], data_lb_fw_gn[snr_lb_fw_gn>=3,1], yerr=data_lb_fw_gn[snr_lb_fw_gn>=3,2],color='red',linestyle='None')
    lb_fw_gn_below, = plt.plot(data_lb_fw_gn[snr_lb_fw_gn<3,0], data_lb_fw_gn[snr_lb_fw_gn<3,1],'D',markersize=6,mfc='None',mec='red')
    #plt.errorbar(data_lb_fw_gn[snr_lb_fw_gn<3,0], data_lb_fw_gn[snr_lb_fw_gn<3,1], yerr=data_lb_fw_gn[snr_lb_fw_gn<3,2],color='red',linestyle='None')

    lg = plt.legend([lb_fw, lb_fw_below, lb_fw_gn, lb_fw_gn_below],['w/o global noise ('+str(len(data_lb_fw[snr_lb_fw>=3,0]))+')','w/o global noise $<$ 3 sigma ('+str(len(data_lb_fw[snr_lb_fw<3,0]))+')','global noise ('+str(len(data_lb_fw_gn[snr_lb_fw_gn>=3,0]))+')','global noise $<$ 3 sigma ('+str(len(data_lb_fw_gn[snr_lb_fw_gn<3,0]))+')'],numpoints=1,loc='upper left')

    plt.xlabel('Wavelength ($\mu$m)',fontsize=16)
    plt.ylabel('Flux',fontsize=16)
    plt.savefig(home + '/line_fitting_comparison/fitting_result/tmc1/'+rec+'_global_noise.eps',format='eps',dpi=300)
    plt.cla()
    plt.clf()
    print 'Finish '+rec+' global noise'

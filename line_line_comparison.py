import numpy as np
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')
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

def read_linescan_greg(indir,reduction,aor_name,line_slice):
    #range scan and line scan comparison
    #TMC-1
    #Greg's reduction
    #reduction = ['_central9Spaxels_PointSourceCorrected_slice_0','_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_0','_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_0']
    #reduction_r = ['central9Spaxels','centralSpaxels_correctedNO','centralSpaxels_correctedYES']
    #aor_name = ['1342225830_TMC_1_', '1342225831_TMC_1_', '1342225830_TMC_1_', '1342225831_TMC_1_']
    band = ['blue','blue','red','red']
    suffix = '_lines_LS_localbaseline_fixwidth'
    #line_slice = [[0,2],[0,5],[1,3],[0,7]]
    name_l = np.array([])
    wl_l = np.array([])
    flux_l = np.array([])
    flux_sig_l = np.array([])
    snr_l = np.array([])
    for iband in range(0, len(band)):
        for i in range(line_slice[iband][0], line_slice[iband][1]+1):
            filename = 'linescan_' + aor_name[iband] +'v13os4_'+ band[iband] +'norm'+ reduction +'_sl'+ str(i) + suffix
            [data, header, name] = read_fitting(indir+filename+'.txt', 0)
            name_l = np.append(name_l, name[:])
            wl_l = np.append(wl_l, data[:,1])
            flux_l = np.append(flux_l, data[:,3])
            flux_sig_l = np.append(flux_sig_l, data[:,4])
            snr_l = np.append(snr_l, data[:,8])
    return name_l, wl_l, flux_l, flux_sig_l, snr_l
wl_ref = [63.1840,72.8430,78.9282,79.12,79.18,79.36,81.8060,84.4110,84.6,87.19,89.9883,90.163,108.0732,108.763,118.581,119.23,119.44,138.5278,144.784,145.525,157.74,162.812,163.12,163.4,174.6259,179.5267,185.999]
wl_ref5 = [63.1840,84.4110,84.6,108.763,118.581,119.23,119.44,138.5278,144.784,145.525,157.74,162.812,174.6259,179.5267,185.999]
wl_ref3 = [63.1840,79.12,79.18,84.4110,84.6,108.0732,108.763,118.581,119.23,119.44,138.5278,144.784,145.525,157.74,162.812,174.6259,179.5267,185.999]
flux_ref5 = [12.63,0.85,0.57,0.64,0.38,0.27,0.50,0.24,0.65,0.93,1.05,0.84,0.26,0.25,0.61]
flux_ref3 = [10.98,0.35,0.37,0.87,0.74,0.26,0.49,0.52,0.42,0.55,0.23,0.77,0.72,0.66,0.73,0.30,0.27,0.59]
#9Spax
#[name_l, wl_l, flux_l, flux_sig_l, snr_l] = read_linescan_greg('/tmc1/v13wishauto/data/','central9Spaxels_PointSourceCorrected',['1342225830_TMC_1_', '1342225831_TMC_1_', '1342225830_TMC_1_', '1342225831_TMC_1_'],[[0,2],[0,5],[1,3],[0,7]])
#[name_l, wl_l, flux_l, flux_sig_l, snr_l] = read_linescan_greg('/tmc1/v13wishauto/data/','centralSpaxel_PointSourceCorrected_Corrected3x3YES',['1342225830_TMC_1_', '1342225831_TMC_1_', '1342225830_TMC_1_', '1342225831_TMC_1_'],[[0,2],[0,5],[1,3],[0,7]])
#[name_l, wl_l, flux_l, flux_sig_l, snr_l] = read_linescan_greg('/tmc1/v13wishauto/data/','centralSpaxel_PointSourceCorrected_Corrected3x3NO',['1342225830_TMC_1_', '1342225831_TMC_1_', '1342225830_TMC_1_', '1342225831_TMC_1_'],[[0,2],[0,5],[1,3],[0,7]])
[name_l, wl_l, flux_l, flux_sig_l, snr_l] = read_linescan_greg('/tmc1/data/v13wish/','',['1342225830', '1342225831', '1342225830', '1342225831'],[[0,2],[0,7],[1,3],[0,7]])

for i in range(0, len(name_l)):
    print '%16s %16.6f %16.6f %6.3f' % (name_l[i],wl_l[i],flux_l[i]/1e-20,snr_l[i])
import numpy as np
from astropy.io import ascii
import os
home = os.path.expanduser('~')
# Header of the 1-D fitting results
# =========================================================================================
# Line,     LabWL(um),      ObsWL(um),      Sig_Cen(um),    Str(W/cm2),     Sig_str(W/cm2)
# FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,            E_u(K),         A(s-1)
# g,        RA(deg),        Dec(deg),       Blend,          Validity
# =========================================================================================

def read_fitting(filepath,noiselevel,obj=False,cube=False):
    data = ascii.read(home+filepath)
    # data = np.genfromtxt(home+filepath, skip_header=1, dtype=str)
    header = data.colnames
    data = data[np.isnan(data['SNR'])!=True]        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
    data = data[(data['SNR']>noiselevel) & (data['Validity']!=0) & (data['Str(W/cm2)']!=0)]
    return data, header
def read_fitting_co(filepath,noiselevel,save_txt=True):
    from astropy.io import ascii
    co_label = []
    upper_level = 48
    lower_level = 4
    for level in range(0,upper_level-lower_level):
        co_label.append('CO'+str(lower_level+level)+'-'+str(lower_level+level-1))
    data = ascii.read(home+filepath)
    # data = np.genfromtxt(home+filepath, skip_header=1, dtype=str)
    # header = data[0]
    # data = data[1:,:]
    header = data.colnames
    data = data[(data['SNR']>noiselevel) & (data['Validity']!=0) & (data['Str(W/cm2)']>0)]
    #data = data.astype(float)

    ind_co = []
    for i in range(0, len(data['Line'])):
        if data['Line'][i] in co_label:
            #wl, wl_sig, flux, flux_sig, E_u, A, g
            ind_co.append(i)
    co_data = data[ind_co]
    return co_data, co_data['Line']
def read_fitting_h2o(filepath,noiselevel,positive=False):
    label = ['p-H2O','o-H2O']
    data = np.genfromtxt(home+filepath, skip_header=1, dtype=str)
    # header = data[0]
    # data = data[1:,:]
    header = []
    data = data
    #data = data.astype(float)
    detection = np.empty(len(data[0,:]))
    i=0
    while i != len(data[:,0]):
        if abs(float(data[i,9])) < noiselevel:
            data = np.delete(data, i, 0)
        else:
            i += 1
    if positive == True:
        while i != len(data[:,0]):
            if abs(float(data[i,4])) <= 0:
                data = np.delete(data, i, 0)
            else:
                i += 1
    name = data[:,0]
    data = data[:,1:]#.astype(float)
    ph2o_data = []
    ph2o_data_name = []
    oh2o_data = []
    oh2o_data_name = []
    for line in name:
        line = str(line)
        if line.find(label[0]) != -1:
            #wl, wl_sig, flux, flux_sig, E_u, A, g
            if np.reshape(data[np.where(name == line),3],1).astype('float') > 0:
                ph2o_data_name.append(line)
                ph2o_data.append([data[np.where(name == line),1],data[np.where(name == line),2],data[np.where(name == line),3],data[np.where(name == line),4],data[np.where(name == line),9],data[np.where(name == line),10],data[np.where(name == line),11]])
        elif line.find(label[1]) != -1:
            #wl, wl_sig, flux, flux_sig, E_u, A, g
            if np.reshape(data[np.where(name == line),3],1).astype('float') > 0:
                oh2o_data_name.append(line)
                oh2o_data.append([data[np.where(name == line),1],data[np.where(name == line),2],data[np.where(name == line),3],data[np.where(name == line),4],data[np.where(name == line),9],data[np.where(name == line),10],data[np.where(name == line),11]])
    ph2o_data = np.array(ph2o_data)
    oh2o_data = np.array(oh2o_data)
    ph2o_data = ph2o_data.reshape(len(ph2o_data_name),7)
    oh2o_data = oh2o_data.reshape(len(oh2o_data_name),7)
    return ph2o_data, ph2o_data_name, oh2o_data, oh2o_data_name
#Uncomment below to test the performance
# filepath = '/bhr71/data/latest/pacs/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines_localbaseline_fixwidth_global_noise.txt'
# [p,p_name, o, o_name] = read_fitting_h2o(filepath,3)
# [co_data, co_data_name] = read_fitting_co(filepath,3)
# print co_data

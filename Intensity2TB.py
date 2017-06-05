def FluxDensity2Intensity(nu, F_v, beamsize, nunit='Hz', funit='Jy'):
    """
    nu: frequency in Hz
    Convert Intensity (Jy) to brightness temperature (K)
    beamsize: diameter of the beam in arcsec
    """
    import numpy as np
    import astropy.constants as const
    c = const.c.cgs.value
    h = const.h.cgs.value
    k = const.k_B.cgs.value

    # other unit of frequency
    if nunit == 'um':
        nu = c/nu*1e4

    # other unit of flux density -> all covert to Jy
    if funit == 'W/cm2/um':
        F_v = F_v*1e4*c/nu**2*1e7*1e23

    # Flux density (Jy) to cgs (erg/s/cm2/Hz/sr)
    I_v = F_v*1e-23/(np.pi*beamsize**2/4)*4.25e10

    return I_v

def Planck(nu, T):
    import astropy.constants as const
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    c = const.c.cgs.value

    # for given temperature, calculate the corresponding B_v
    B_v = 2*h*nu**3/c**2*(np.exp(h*nu/k/T)-1)**-1

    return B_v

def InversePlanck(nu, Bv):
    import astropy.constants as const
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    c = const.c.cgs.value

    # For given value of B_v, calculate the corresponding temperature
    T = (h*nu/k)*(np.log(2*h*nu**3/c**2/Bv+1))**-1

    return T

from astropy.io import ascii
from scipy.interpolate import interp1d
import numpy as np
import astropy.constants as const
h = const.h.cgs.value
k = const.k_B.cgs.value
c = const.c.cgs.value

fitting = ascii.read('/Volumes/SD-Mac/CDF_archive_v2/CDF_archive_v2_lines.txt')
tau4TB = ascii.read('/Volumes/SD-Mac/Dropbox/real_cops-spire/tau_for_TB.txt')

T_rot_list = {'Ced110': [58.4],
              'GSS30-IRS1': [449.4, 121.6, 96.0, 56.3],
              'VLA1623': [348.1, 125.0, 79.3, 43.8],
              'WL12': [589.5, 162.1, 101.7, 35.9],
              'L1551-IRS5': [484.8, 99.7, 39.1],
              'RCrA-IRS5A': [318.2, 121.1, 77.7, 43.0],
              'RCrA-IRS7B': [367.5, 122.8, 92.9, 61.6],
              'RCrA-IRS7C': [394.9, 117.0, 98.1, 51.5],
              'BHR71': [451.3, 145.1, 87.7, 37.1]}

# get the object list for tau4TB, where the 13CO is detected
obj_list = list(set(tau4TB['Object']))

# set up the beam profiles for SLW and SSW
# read in beam profile of spire
ssw_beam = ascii.read('/Users/yaolun/data/SSW_beam_profile.txt', names=['wavelength(um)', 'size(arcsec)'])
slw_beam = ascii.read('/Users/yaolun/data/SLW_beam_profile.txt', names=['wavelength(um)', 'size(arcsec)'])
# trimming
trimmer_ssw = (ssw_beam['wavelength(um)'] < 310) & (ssw_beam['wavelength(um)'] >= 195)
trimmer_slw = (slw_beam['wavelength(um)'] > 310)
# function to interpolate the beam size at given wavelength
f_ssw_beam = interp1d(ssw_beam['wavelength(um)'][trimmer_ssw], ssw_beam['size(arcsec)'][trimmer_ssw])
f_slw_beam = interp1d(slw_beam['wavelength(um)'][trimmer_slw], slw_beam['size(arcsec)'][trimmer_slw])

# use the aperture for extracting 1D SPIRE spectrum
fitted_size = {'Ced110': 30.5, 'GSS30-IRS1': 38.5, 'VLA1623': 25.75, 'WL12': 41.0,
               'L1551-IRS5': 14.25, 'RCrA-IRS5A': 40.0, 'RCrA-IRS7B': 37.0, 'RCrA-IRS7C': 38.0,
               'BHR71': 15.5}

# loop through each pair of object and J_up found in tau4TB file
for i, obj in enumerate(obj_list):

    if obj not in T_rot_list.keys():
        continue
    print(obj)

    # get the J_up list
    Jup_list = tau4TB['J_up'][tau4TB['Object'] == obj]

    for j in Jup_list:
        # get the data
        data = fitting[(fitting['Object'] == obj) & (fitting['SNR'] >= 4) & (fitting['Validity'] == 1) &\
                       (fitting['Pixel_No.'] == 'c') & (fitting['Line'] == 'CO'+str(j)+'-'+str(j-1))]

        # fail safe
        if len(data) == 0:
            print('No data found for J_up='+str(j)+' in'+obj)
            continue

        # use the averaged line width of 5 km/s
        F_v = data['Str(W/cm2)'].data/(13.8e5/c*data['ObsWL(um)'].data)/1.064
        nu = c/data['ObsWL(um)'].data*1e4

        if data['ObsWL(um)'] >= 310:
            # SLW
            T_B = InversePlanck(nu, FluxDensity2Intensity(nu, F_v,
                                                          (f_slw_beam(data['ObsWL(um)'])**2+fitted_size[obj]**2)**0.5, funit='W/cm2/um'))
        else:
            # SSW
            T_B = InversePlanck(nu, FluxDensity2Intensity(nu, F_v,
                                                          (f_ssw_beam(data['ObsWL(um)'])**2+fitted_size[obj]**2)**0.5, funit='W/cm2/um'))

        print('The brightness temperature is ', T_B[0])
        print('The optical depth is ', tau4TB['tau'][(tau4TB['Object'] == obj) & (tau4TB['J_up'] == j)].data[0])
        # calculate the filling factor for given rotational temperature fitted from rotational diagrams

        # The rotational temperatures fitted from rotational diagram
        try:
            T_rot = T_rot_list[obj]
        except KeyError:
            continue
        for t in T_rot:
            print(t, j)

            # print('f - opt. thick', Planck(nu, T_B[0])/Planck(nu, t))
            print('f - full', Planck(nu, T_B[0])/\
                              Planck(nu, t)/(1-np.exp(-tau4TB['tau'][(tau4TB['Object'] == obj) & (tau4TB['J_up'] == j)].data[0])))


# # sort by wavelength
# data = data[np.argsort(data['ObsWL(um)'])]
#
# data_trim_ssw = (data['ObsWL(um)'] < 310) & (data['ObsWL(um)'] >= 195)
# data_trim_slw = (data['ObsWL(um)'] > 310)
#
#
# # use the averaged line width of 5 km/s
# F_v = data['Str(W/cm2)'].data/(50e5/c*data['ObsWL(um)'].data)/1.064
#
# # SSW
# nu_ssw = c/data[data_trim_ssw]['ObsWL(um)'].data*1e4
# T_B_ssw = InversePlanck(nu_ssw, FluxDensity2Intensity(nu_ssw, F_v[data_trim_ssw], f_ssw_beam(data['ObsWL(um)'][data_trim_ssw]), funit='W/cm2/um'))
#
# # SLW
# nu_slw = c/data[data_trim_slw]['ObsWL(um)'].data*1e4
# T_B_slw = InversePlanck(nu_slw, FluxDensity2Intensity(nu_slw, F_v[data_trim_slw], f_slw_beam(data['ObsWL(um)'][data_trim_slw]), funit='W/cm2/um'))
#
# # print the brightness temperatures of all lines
# # [print(data['Line'][data_trim_ssw].data[i], T_B_ssw[i]) for i in range(len(T_B_ssw))]
# # [print(data['Line'][data_trim_slw].data[i], T_B_slw[i]) for i in range(len(T_B_slw))]
#
# # print(T_B_slw[data['Line'][data_trim_slw] == 'o-H2O1_10-1_01'])
# # print(data['SNR'][data['Line'][data_trim_slw] == 'o-H2O1_10-1_01'])
#
# # get the indice for selecting CO data
# ind_co_ssw = []
# ind_co_slw = []
#
# for i in range(len(data[data_trim_ssw])):
#     if len(data['Line'][data_trim_ssw][i].split('CO')[0]) == 0:
#         ind_co_ssw.append(i)
# for i in range(len(data[data_trim_slw])):
#     if len(data['Line'][data_trim_slw][i].split('CO')[0]) == 0:
#         ind_co_slw.append(i)
#
# # print the beam sizes for CO lines
# print('Beam size')
# print(f_ssw_beam(data['ObsWL(um)'][data_trim_ssw][ind_co_ssw]))
# print(f_slw_beam(data['ObsWL(um)'][data_trim_slw][ind_co_slw]))
#
# # print the calculated Intensity for CO lines
# print('I_v')
# print(FluxDensity2Intensity(nu_ssw[ind_co_ssw], F_v[data_trim_ssw][ind_co_ssw], f_ssw_beam(data['ObsWL(um)'][data_trim_ssw][ind_co_ssw]), funit='W/cm2/um'))
# print(FluxDensity2Intensity(nu_slw[ind_co_slw], F_v[data_trim_slw][ind_co_slw], f_slw_beam(data['ObsWL(um)'][data_trim_slw][ind_co_slw]), funit='W/cm2/um'))
#
# # optical depth as a function of upper energy
# def tau12(E_u):
#     p = [1.43984169, -0.00235609]
#     return 10**(p[1]*E_u+p[0])
#
# # calculate the filling factor for given rotational temperature fitted from rotational diagrams
#
# # print the brightness temperatures of all CO lines
# print('Line', 'Str(W/cm2)', 'T_B(K)')
# [print(data['Line'][data_trim_ssw][ind_co_ssw].data[i], data['Str(W/cm2)'][data_trim_ssw][ind_co_ssw].data[i], T_B_ssw[ind_co_ssw][i]) for i in range(len(T_B_ssw[ind_co_ssw]))]
# [print(data['Line'][data_trim_slw][ind_co_slw].data[i], data['Str(W/cm2)'][data_trim_slw][ind_co_slw].data[i], T_B_slw[ind_co_slw][i]) for i in range(len(T_B_slw[ind_co_slw]))]
#
# print('Line', 'ObsWL(um)', 'FWHM(um)')
# [print(data['Line'][data_trim_ssw][ind_co_ssw].data[i], data['ObsWL(um)'][data_trim_ssw][ind_co_ssw].data[i], data['FWHM(um)'][data_trim_ssw][ind_co_ssw].data[i]) for i in range(len(T_B_ssw[ind_co_ssw]))]
# [print(data['Line'][data_trim_slw][ind_co_slw].data[i], data['ObsWL(um)'][data_trim_slw][ind_co_slw].data[i], data['FWHM(um)'][data_trim_slw][ind_co_slw].data[i]) for i in range(len(T_B_slw[ind_co_slw]))]
# # print(f_slw_beam(data['ObsWL(um)'][data_trim_slw][ind_co_slw][-1]))
#
# # print the optical depth of all CO lines
# print('Tau12')
# print(tau12(data[data_trim_ssw][ind_co_ssw]['E_u(K)'].data))
# print(tau12(data[data_trim_slw][ind_co_slw]['E_u(K)'].data))
#
# # The rotational temperatures fitted from rotational diagram
# T_rot = np.array([35.0, 77.4, 102.6])
# for i in range(len(T_rot)):
#     print(T_rot[i])
#
#     print('f - opt. thick', Planck(nu_ssw[ind_co_ssw], T_B_ssw[ind_co_ssw])/Planck(nu_ssw[ind_co_ssw], T_rot[i]))
#     print('f - opt. thick', Planck(nu_slw[ind_co_slw], T_B_slw[ind_co_slw])/Planck(nu_slw[ind_co_slw], T_rot[i]))
#     print('f - full', Planck(nu_ssw[ind_co_ssw], T_B_ssw[ind_co_ssw])/\
#                       Planck(nu_ssw[ind_co_ssw], T_rot[i])/(1-np.exp(-tau12(data[data_trim_ssw][ind_co_ssw]['E_u(K)'].data))))
#     print('f - full', Planck(nu_slw[ind_co_slw], T_B_slw[ind_co_slw])/\
#                       Planck(nu_slw[ind_co_slw], T_rot[i])/(1-np.exp(-tau12(data[data_trim_slw][ind_co_slw]['E_u(K)'].data))))

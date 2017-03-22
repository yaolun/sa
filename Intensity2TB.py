def Intensity2TB(nu, F_v, beamsize, nunit='Hz', funit='Jy'):
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

    T_B = h*nu/k*(np.log(2*h*nu**3/c**2/I_v+1))**-1

    return T_B

from astropy.io import ascii
from scipy.interpolate import interp1d
import numpy as np

fitting = ascii.read('/Volumes/SD-Mac/CDF_archive_v2/CDF_archive_v2_lines.txt')

data = fitting[(fitting['Object'] == 'BHR71') & (fitting['SNR'] >= 5) & (fitting['Validity'] == 1) &\
               (fitting['Pixel_No.'] == 'c')]

# sort by wavelength
data = data[np.argsort(data['ObsWL(um)'])]

# read in beam profile of spire
ssw_beam = ascii.read('/Users/yaolun/data/SSW_beam_profile.txt', names=['wavelength(um)', 'size(arcsec)'])
slw_beam = ascii.read('/Users/yaolun/data/SLW_beam_profile.txt', names=['wavelength(um)', 'size(arcsec)'])
# trimming
trimmer_ssw = (ssw_beam['wavelength(um)'] < 310) & (ssw_beam['wavelength(um)'] >= 195)
trimmer_slw = (slw_beam['wavelength(um)'] > 310)

data_trim_ssw = (data['ObsWL(um)'] < 310) & (data['ObsWL(um)'] >= 195)
data_trim_slw = (data['ObsWL(um)'] > 310)

# function to interpolate the beam size at given wavelength
f_ssw_beam = interp1d(ssw_beam['wavelength(um)'][trimmer_ssw], ssw_beam['size(arcsec)'][trimmer_ssw])
f_slw_beam = interp1d(slw_beam['wavelength(um)'][trimmer_slw], slw_beam['size(arcsec)'][trimmer_slw])

nu = data['ObsWL(um)']
F_v = data['Str(W/cm2)']/data['FWHM(um)']/1.064

# SSW
T_B_ssw = Intensity2TB(nu[data_trim_ssw].data, F_v[data_trim_ssw].data, (15.5**2+(f_ssw_beam(nu[data_trim_ssw].data)**2))**0.5,
                       nunit='um', funit='W/cm2/um')
# SLW
T_B_slw = Intensity2TB(nu[data_trim_slw].data, F_v[data_trim_slw].data, (15.5**2+(f_slw_beam(nu[data_trim_slw].data)**2))**0.5,
                       nunit='um', funit='W/cm2/um')

# [print(data['Line'][data_trim_ssw].data[i], T_B_ssw[i]) for i in range(len(T_B_ssw))]
# [print(data['Line'][data_trim_slw].data[i], T_B_slw[i]) for i in range(len(T_B_slw))]

# get the CO data
ind_co_ssw = []
ind_co_slw = []

for i in range(len(data[data_trim_ssw])):
    if len(data['Line'][data_trim_ssw][i].split('CO')[0]) == 0:
        ind_co_ssw.append(i)
for i in range(len(data[data_trim_slw])):
    if len(data['Line'][data_trim_slw][i].split('CO')[0]) == 0:
        ind_co_slw.append(i)


# optical depth as a function of upper energy
def tau12(E_u):
    p = [1.43984169, -0.00235609]
    return 10**(p[1]*E_u+p[0])

# function to derive T_ex
# def T_ex(T_B, tau, nu, unit='Hz'):
#     import numpy as np
#     import astropy.constants as const
#     h = const.h.cgs.value
#     k = const.k_B.cgs.value
#     c = const.c.cgs.value
#
#     if unit == 'um':
#         nu = c/nu*1e4
#
#     return (h*nu/k)*(np.log(1+h*nu/k/T_B*(1-np.exp(-tau))))**-1
#
#
#
# T_ex_slw = T_ex(T_B_slw[ind_co_slw],
#                 tau12(data[data_trim_slw][ind_co_slw]['E_u(K)'].data),
#                 data[data_trim_slw][ind_co_slw]['ObsWL(um)'].data)
# T_ex_ssw = T_ex(T_B_ssw[ind_co_ssw],
#                 tau12(data[data_trim_ssw][ind_co_ssw]['E_u(K)'].data),
#                 data[data_trim_ssw][ind_co_ssw]['ObsWL(um)'].data)

# print(T_ex_slw)
# print(T_ex_ssw)

# calculate the filling factor for given rotational temperature fitted from rotational diagrams
import astropy.constants as const
h = const.h.cgs.value
k = const.k_B.cgs.value
c = const.c.cgs.value

def J_v(nu, T, unit='Hz'):
    import astropy.constants as const
    import numpy as np
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    c = const.c.cgs.value

    if unit == 'um':
        nu = c/nu*1e4

    return h*nu/k/(np.exp(h*nu/k/T)-1)

[print(data['Line'][data_trim_ssw][ind_co_ssw].data[i], T_B_ssw[ind_co_ssw][i]) for i in range(len(T_B_ssw[ind_co_ssw]))]
[print(data['Line'][data_trim_slw][ind_co_slw].data[i], T_B_slw[ind_co_slw][i]) for i in range(len(T_B_slw[ind_co_slw]))]

T_rot = np.array([35.0, 77.4, 102.6])
f_slw = np.empty_like(T_rot)
f_ssw = np.empty_like(T_rot)
for i in range(len(T_rot)):
    print(T_rot[i])
    print('RJ-limit',h*(c/data[data_trim_slw][ind_co_slw]['ObsWL(um)'].data*1e4)/k/T_rot[i])
    print('RJ-limit',h*(c/data[data_trim_ssw][ind_co_ssw]['ObsWL(um)'].data*1e4)/k/T_rot[i])
    print('f (full)',T_B_slw[ind_co_slw]/J_v(data[data_trim_slw][ind_co_slw]['ObsWL(um)'].data, T_rot[i], unit='um')/(1-np.exp(-tau12(data[data_trim_slw][ind_co_slw]['E_u(K)'].data))))
    print('f (full)',T_B_slw[ind_co_ssw]/J_v(data[data_trim_ssw][ind_co_ssw]['ObsWL(um)'].data, T_rot[i], unit='um')/(1-np.exp(-tau12(data[data_trim_ssw][ind_co_ssw]['E_u(K)'].data))))

# if T_ex satisfies RJ-limit
for i in range(len(T_rot)):
    print(T_rot[i])
    print('f (RJ)',T_B_slw[ind_co_slw]/(T_rot[i])/(1-np.exp(-tau12(data[data_trim_slw][ind_co_slw]['E_u(K)'].data))))
    print('f (RJ)',T_B_ssw[ind_co_ssw]/(T_rot[i])/(1-np.exp(-tau12(data[data_trim_ssw][ind_co_ssw]['E_u(K)'].data))))

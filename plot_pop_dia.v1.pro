pro plot_pop_dia
;label = ['SLWA1','SLWA2','SLWB1','SLWB2','SLWB3','SLWC2','SLWC3','SLWC4','SLWC5','SLWD2','SLWD3','SLWD4','SLWE2','SLWE3','SSWA1','SSWA2','SSWA3',$
;         'SSWA4','SSWB1','SSWB2','SSWB3','SSWB4','SSWB5','SSWC1','SSWC2','SSWC3','SSWC4','SSWC5','SSWC6','SSWD1','SSWD2','SSWD3','SSWD4','SSWD6','SSWD7','SSWE1','SSWE2','SSWE3',$
;         'SSWE4','SSWE5','SSWE6','SSWF1','SSWF2','SSWF3','SSWF5','SSWG1','SSWG2','SSWG3','SSWG4']
label = [['SSWD4', 'SLWC3'], ['SSWE2', 'SLWD2'], ['SSWF3', 'SLWC2'], ['SSWE5', 'SLWB2'], ['SSWC5', 'SLWB3'], ['SSWB3', 'SLWC4'], ['SSWC2', 'SLWD3']]
peak = [[654.2, 645.8], [523.15, 516.85], [436.125, 430.875], [373.075, 369.925], [326.5626, 323.9375], [290.05, 287.95], [260.9, 259.5], [237.23, 235.97], [217.635, 216.165]]
A = [6.13E-06, 1.22E-05, 2.14E-05, 3.42E-05, 5.13E-05, 7.33E-05, 1.01E-04, 1.34E-04, 1.75E-04]
utran = [4,5,6,7,8,9,10,11,12]
E_u = [55.32, 82.97, 116.16, 154.87, 199.11, 248.88, 304.16, 364.97, 431.29]
for i = 0, 6 do begin
    get_spire, label[0,i], label[1,i], wl, flux
    wl_fit = dblarr(9)
    flux_fit = dblarr(9) & sig_flux = dblarr(9)
    for j = 0, n_elements(peak[0,*])-1 do begin
        print, label[*,i], utran[j]
        ind = where(wl gt peak[1,j] and wl le peak[0,j])
        wll = wl[ind] & fluxx = flux[ind]
        fit_line, label[0,i]+label[1,i]+'_'+strtrim(utran[j],1), wll, fluxx, status, errmsg, cen_wl, sig_cen_wl, str, sig_str, fwhm, sig_fwhm, base_str, rms
        wl_fit[j] = cen_wl
        flux_fit[j] = str & sig_flux[j] = sig_str
    endfor
    pop_dia_co, label[0,i]+label[1,i], wl_fit, flux_fit, sig_flux, A, E_u, utran
endfor
end

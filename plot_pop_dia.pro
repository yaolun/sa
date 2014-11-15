pro plot_pop_diav2, indir, outdir, objname
legned = objname
;---------------------------------
angle = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, 1*COS(angle), 1*SIN(angle), /FILL, color = 80
;---------------------------------
;SPIRE
pixelname = ['SLWA1','SLWA2','SLWB1','SLWB2','SLWB3','SLWC2','SLWC3','SLWC4','SLWC5','SLWD2','SLWD3','SLWD4','SLWE2','SLWE3','SSWA1','SSWA2','SSWA3','SSWA4','SSWB1','SSWB2','SSWB3','SSWB4','SSWB5','SSWC1','SSWC2','SSWC3','SSWC4','SSWC5','SSWC6','SSWD1','SSWD2','SSWD3','SSWD4','SSWD6','SSWD7','SSWE1','SSWE2','SSWE3','SSWE4','SSWE5','SSWE6','SSWF1','SSWF2','SSWF3','SSWF5','SSWG1','SSWG2','SSWG3','SSWG4']

label = [['SSWD4', 'SLWC3'], ['SSWE2', 'SLWD2'], ['SSWF3', 'SLWC2'], ['SSWE5', 'SLWB2'], ['SSWC5', 'SLWB3'], ['SSWB3', 'SLWC4'], ['SSWC2', 'SLWD3']]
;peak = [[654.2, 645.8], [523.15, 516.85], [436.125, 430.875], [373.075, 369.925], [326.5626, 323.9375], [290.05, 287.95], [260.9, 259.5], [237.23, 235.97], [217.635, 216.165]]
A = [6.13E-06, 1.22E-05, 2.14E-05, 3.42E-05, 5.13E-05, 7.33E-05, 1.01E-04, 1.34E-04, 1.75E-04, 2.20E-04]
utran = [4,5,6,7,8,9,10,11,12,13]
g = [9,11,13,15,17,19,21,23,25,27]
E_u = [55.32, 82.97, 116.16, 154.87, 199.11, 248.88, 304.16, 364.97, 431.29, 503.13]

plot_ra = fltarr(7)
plot_dec = fltarr(7)
plot_t_rot = fltarr(7)
for i = 0, 6 do begin
	filename_s = indir+'/data/'+label[0,i]+'lines.txt' & filename_l = indir+'/data/'+label[1,i]+'lines.txt'
	readcol, filename_s, format='A,D,D,D,D,D,D,D,D,D,D,D,D,D',line_s, wl_s, sig_wl_s, str_s, sig_str_s, fwhm_s, sig_fwhm_s, base_str_s,snr_s, E_u_s, A_s, g_s, ra_s, dec_s, skipline=1
	readcol, filename_l, format='A,D,D,D,D,D,D,D,D,D,D,D,D,D',line_l, wl_l, sig_wl_l, str_l, sig_str_l, fwhm_l, sig_fwhm_l, base_str_l,snr_l, E_u_l, A_l, g_l, ra_l, dec_l, skipline=1
	line = [line_s, line_l]
	wl = [wl_s, wl_l] & sig_wl = [sig_wl_s, sig_wl_l]
	str = [str_s, str_l] & sig_str = [sig_str_s, sig_str_l]
	fwhm = [fwhm_s, fwhm_l] & sig_fwhm = [sig_fwhm_s, sig_fwhm_l]
	snr = [snr_s, snr_l]
	ra = [ra_s, ra_l]
	dec = [dec_s, dec_l]
	plot_ra[i] = ra_l[0]
	plot_dec[i] = dec_l[0]
	i4 = where(line eq 'CO4-3')
	i5 = where(line eq 'CO5-4')
	i6 = where(line eq 'CO6-5')
	i7 = where(line eq 'CO7-6')
	i8 = where(line eq 'CO8-7')
	i9 = where(line eq 'CO9-8')
	i10 = where(line eq 'CO10-9')
	i11 = where(line eq 'CO11-10')
	i12 = where(line eq 'CO12-11')
	i13 = where(line eq 'CO13-12')
	ind = [i4,i5,i6,i7,i8,i9,i10,i11,i12,i13]
	wl_co = fltarr(10)
	sig_wl_co = fltarr(10)
	str_co = fltarr(10)
	sig_str_co = fltarr(10)
	for j = 0, 9 do begin
		if ind[j] ge 0 then begin
			wl_co[j] = wl[ind[j]]
			sig_wl_co[j] = sig_wl[ind[j]]
			str_co[j] = str[ind[j]]
			sig_str_co[j] = sig_str[ind[j]]
		endif
	endfor
	pop_dia_co, outdir, label[0,i]+'_'+label[1,i], wl_co, str_co, sig_str_co, A, E_u, g, t_rot,legend='BHR71'
	plot_t_rot[i] = t_rot
endfor
ra_cen = 180.3998
dec_cen = -65.14834
;plot_ra = plot_ra - ra_cen
;plot_dec = plot_dec - dec_cen
print, plot_ra
print, plot_dec
print, plot_t_rot
set_plot, 'ps'
!p.font = 0
device, filename = outdir+'/plots/CO_rotational_temperature_spire.eps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 16, decomposed = 0, /color, xoffset = 200
loadct, 13
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
  plot_ra = (plot_ra-180.3998)*3600*cos(plot_dec*!pi/180.)
  plot_dec = (plot_dec+65.14834)*3600
  plot, plot_ra, plot_dec, psym = 8, xtitle='RA offset (arcsec)', ytitle='Dec offset (arcsec)', xstyle = 2, ystyle = 2, position=aspect(1.);, xrange=[-0.05,0.05]
  for i = 0, 6 do begin
     xyouts, plot_ra[i]+2, plot_dec[i], strtrim(string(floor(plot_t_rot[i])),1)+'K', charsize=0.7
  endfor
device, /close_file, decomposed = 1
!p.multi = 0

;PACS
plot_ra = fltarr(25)
plot_dec = fltarr(25)
plot_t_rot = fltarr(25)
for i = 0, 24 do begin
    filename = indir+'/data/pacs_pixel'+strtrim(string(i+1),1)+'_lines.txt'
    readcol, filename, format='A,D,D,D,D,D,D,D,D,D,D,D',line, wl, sig_wl, str, sig_str, fwhm, sig_fwhm, base_str,snr, E_u, A, g, ra, dec, skipline=1
    plot_ra[i] = ra[0]
    plot_dec[i] = dec[0]
    i14 = where(line eq 'CO14-13')
    i15 = where(line eq 'CO15-14')
    i16 = where(line eq 'CO16-15')
    i17 = where(line eq 'CO17-16')
    i18 = where(line eq 'CO18-17')
    i19 = where(line eq 'CO19-18')
    i20 = where(line eq 'CO20-19')
    i21 = where(line eq 'CO21-20')
    i22 = where(line eq 'CO22-21')
    i23 = where(line eq 'CO23-22')
    i24 = where(line eq 'CO24-23')
    i25 = where(line eq 'CO25-24')
    i26 = where(line eq 'CO26-25')
    i27 = where(line eq 'CO27-26')
    i28 = where(line eq 'CO28-27')
    i29 = where(line eq 'CO29-28')
    i30 = where(line eq 'CO30-29')
    i31 = where(line eq 'CO31-30')
    i32 = where(line eq 'CO32-31')
    i33 = where(line eq 'CO33-32')
    i34 = where(line eq 'CO34-33')
    i35 = where(line eq 'CO35-34')
    i36 = where(line eq 'CO36-35')
    i37 = where(line eq 'CO37-36')
    i38 = where(line eq 'CO38-37')
    i39 = where(line eq 'CO39-38')
    i40 = where(line eq 'CO40-39')
    ind = [i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35,i36,i37,i38,i39,i40]
    wl_co = []
    sig_wl_co = []
    str_co = []
    sig_str_co = []
    use_ind = []
    A_j = []
    utran_j = []
    g_j = []
    E_u_j = []
    A = [2.739e-4,3.354e-4,4.05e-4,4.829e-4,5.695e-4, 6.65e-4,7.695e-4,8.833e-4,1.006e-3,1.139e-3,1.281e-3,1.432e-3,1.592e-3,1.761e-3,1.94e-3,2.126e-3,2.321e-3,2.524e-3,2.735e-3,2.952e-3,3.175e-3,$
    3.404e-3,3.638e-3,3.878e-3,4.12e-3,4.365e-3,4.613e-3]
    utran = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
    g = [29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81]
    E_u = [580.49,663.35,751.72,845.59,944.97,1049.84,1160.2,1276.05,1397.38,1524.19,1656.47,1794.23,1937.44,2086.12,2240.24,2399.82,2564.83,2735.28,2911.15,3092.45,3279.15,3471.27,3668.78,3871.69,$
      4079.98,4293.64,4512.67]
    for j = 0, 26 do begin
        if ind[j] ge 0 then begin
            if snr[ind[j]] ge 3 then begin
                wl_co = [wl_co,wl[ind[j]]]
                sig_wl_co = [sig_wl_co,sig_wl[ind[j]]]
                str_co = [str_co,str[ind[j]]]
                sig_str_co = [sig_str_co,sig_str[ind[j]]]
                use_ind = [use_ind,ind[j]]
                A_j = [A_j,A[j]]
                utran_j = [utran_j,utran[j]]
                g_j = [g_j,g[j]]
                E_u_j = [E_u_j,E_u[j]]
            endif
        endif
    endfor
    if n_elements(str_co) ge 3 then begin
	openw, lun, '~/bhr71/data/pacs_co_rot_pixel'+strtrim(string(i+1),1)+'.txt',/get_lun
	for line = 0, n_elements(wl_co)-1 do begin
	printf, lun, format='(6(g10.4,2x))',wl_co[line], str_co[line], sig_str_co[line], A_j[line], E_u_j[line], g_j[line]
	endfor
	close, lun
	free_lun, lun
        pop_dia_co, outdir, 'pacs_pixel'+strtrim(string(i+1),1), wl_co, str_co, sig_str_co, A_j, E_u_j, g_j, t_rot,legend=objname,two_comp=1800
        plot_t_rot[i] = t_rot
    endif
endfor
ra_cen = 180.3998
dec_cen = -65.14834
;plot_ra = plot_ra - ra_cen
;plot_dec = plot_dec - dec_cen
;print, plot_ra
;print, plot_dec
;print, plot_t_rot
set_plot, 'ps'
!p.font = 0
device, filename = outdir+'/plots/CO_rotational_temperature_pacs.eps', /portrait, /helvetica, /encapsulated, isolatin = 1, font_size = 14, decomposed = 0, /color, xoffset = 200
loadct, 6
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
  plot_ra = (plot_ra-180.3998)*3600*cos(plot_dec*!pi/180.)
  plot_dec = (plot_dec+65.14834)*3600
  plot, plot_ra, plot_dec, psym = 8, xtitle='RA Offset (arcsec)', ytitle='Dec Offset (arcsec)', xstyle = 2, ystyle = 2, position=aspect(1.), charthick=3, xrange=[40,-40],/nodata
  for pix = 0, 24 do begin
    if plot_t_rot[pix] ge 1000 then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=65
    if (plot_t_rot[pix] ge 500) and (plot_t_rot[pix] lt 1000) then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=127
    if (plot_t_rot[pix] lt 500) and (plot_t_rot[pix] gt 0) then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=200
    if plot_t_rot[pix] eq 0 then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=0
  endfor
  for i = 0, 24 do begin
     xyouts, plot_ra[i]-0.7, plot_dec[i]+0.9, strtrim(string(floor(plot_t_rot[i])),1), charsize=0.5
  endfor
device, /close_file, decomposed = 1
!p.multi = 0
set_plot,'x'
;Make the rotational temperature compare to the emissions lines contour.
noise=[1,1,1]
plot_contour, data_slw, data_ssw, data_pacs, noise, /no_plot
;SPIRE
;SLW
line_name = ['p-H2O2_02-1_11','CO8-7','13CO8-7','CO7-6','CI370','13CO7-6','p-H2O2_11-2_02','CO6-5','13CO6-5','HCO+P7-6','CO5-4','o-H2O1_10-1_01','13CO5-4','CI610','CO4-3']
line_name = ['CO4-3']
for i = 0, n_elements(line_name)-1 do begin
    wl = []
    flux = []
    flux_sig = []
    ra = []
    ra_tot = []
    dec = []
    dec_tot = []
    for pix = 0, n_elements(data_slw[*].ra)-1 do begin
        data_ind = where(line_name[i] eq data_slw[pix].line)

        wl = [wl, data_slw[pix].wl[data_ind]]
        flux = [flux, data_slw[pix].flux[data_ind]]
        flux_sig = [flux_sig, data_slw[pix].flux_sig[data_ind]]
        ra = [ra, data_slw[pix].ra[data_ind]]
        dec = [dec, data_slw[pix].dec[data_ind]]
    endfor
    if n_elements(flux) ge 3 then begin
        set_plot, 'ps'
        !p.font = 0
        device, filename = '/Users/yaolun/bhr71/plots/CO_rotational_temperature_'+line_name[i]+'_contour.eps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 14, decomposed = 0, /color
        loadct,6
        !p.thick = 4 & !x.thick = 5 & !y.thick = 5
        flux_sig[where(flux_sig gt 1)] = 0
        noisefloor = median(flux_sig[where(flux_sig gt 0)])*noise[2]
        level = indgen(11)/10.*(max(flux)-max([min(flux),noisefloor]))+max([min(flux),noisefloor]); > noisefloor
        level = level[sort(level)]
        ra = (ra-180.3998)*3600*cos(dec*!pi/180.)
        dec = (dec+65.14834)*3600
        flux_smooth = min_curve_surf(flux,ra,dec)
        ra_smooth = min_curve_surf(ra,ra,dec)
        dec_smooth = min_curve_surf(dec,ra,dec)
        plotposition = aspect(1.)
        plot, plot_ra, plot_dec, psym=8,xtitle='RA offset (arcsec)', ytitle='Dec offset (arcsec)', xrange=[120,-120], yrange=[-120,120],position=plotposition, /nodata
        contour, flux_smooth, ra_smooth, dec_smooth, levels=level, /irregular, /overplot, position=plotposition, color=250
        for pix = 0, 24 do begin
            if plot_t_rot[pix] ge 1000 then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=65
            if (plot_t_rot[pix] ge 500) and (plot_t_rot[pix] lt 1000) then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=127
            if (plot_t_rot[pix] lt 500) and (plot_t_rot[pix] gt 0) then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=200
            ;if plot_t_rot[pix] eq 0 then oplot, [plot_ra[pix]-0.006], [plot_dec[pix]+0.8], psym=8, color=0
        endfor
        xyouts, -10, 125, title_name(line_name[i])
        device, /close_file, decomposed = 1
        !p.multi = 0
        ;stop
        set_plot,'x'
    endif
endfor
end

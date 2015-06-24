pro get_oi_temp,filename,outdir,cube=cube,oi=oi,cii=cii
	c = 2.99792458e10
	k = 1.380658e-16
	;all_digit_jitter_centralyes_lines_localbaseline_fixwidth_global_noise
	if not keyword_set(cube) then begin
		readcol,filename,format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,A,D',$
			object, line_name, lab_wl, cen_wl, sig_cen_wl, str, sig_str, fwhm, sig_fwhm, base_str, snr, noise, E_u, A, g, ra, dec, blend, valid, skipline=1,/silent
		msg = ''
	endif else begin
		readcol,filename,format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,A,A,D',$
			object, line_name, lab_wl, cen_wl, sig_cen_wl, str, sig_str, fwhm, sig_fwhm, base_str, snr, noise, E_u, A, g, ra, dec, pix, blend, valid, skipline=1,/silent
		msg = '_cube_'
	endelse
	if keyword_set(oi) then begin
		ind = where(line_name eq 'OI3P1-3P2')
		name = 'oi'
	endif
	if keyword_set(cii) then begin
		ind = where(line_name eq 'CII2P3_2-2P1_2')
		name = 'cii'
	endif
	;pix_size = 15.3^2  ; pixel size of the SOFIA-GREAT: 15.3 arcsec
	pix_size = 9.4^2
	openw, lun, outdir+'brightness_temperature_'+name+'_for_GREAT'+msg+'.txt',/get_lun
	if not keyword_set(cube) then begin
		printf, lun, format='(a10,2x,a14,2x,a10,2x,a14,2x,a8,2x,a10)','Object', 'Line','Wavelength','Flux(W/cm2)','Velocity','T_rad'
	endif else begin
		printf, lun, format='(a10,2x,a8,2x,a14,2x,a10,2x,a14,2x,a8,2x,a10)','Object', 'Pixel', 'Line','Wavelength','Flux(W/cm2)','Velocity','T_rad'
	endelse
	for i = 0, n_elements(ind)-1 do begin
		j = ind[i]
		obj = object[j]
		if obj ne 'L1551-IRS5' then continue
		width = fwhm[j];/2.354						; um
		velo = width/cen_wl[j]*c/1e5				; km/s
		; assume all velocity should not be greater than 100 km/s
		velo = 20.0
		omega_b = !pi/4/alog(2)*pix_size*(1/3600.0/180*!pi)^2		; sr
		Tr = (cen_wl[j]*1e-4)^3/2/k*(str[j]*1e7)/omega_b/(velo*1e5)
		;if object[j] eq 'L1551-IRS5' then print, width,omega_b
		if not keyword_set(cube) then begin
			print, format='(a10,2x,a14,2x,g10.8,2x,g14.8,2x,g8.5,2x,g10.4)',object[j], line_name[j],cen_wl[j],str[j],velo,Tr
			printf, lun, format='(a10,2x,a14,2x,g10.8,2x,g14.8,2x,g8.5,2x,g10.4)',object[j], line_name[j],cen_wl[j],str[j],velo,Tr
		endif else begin
			print, format='(a10,2x,a8,2x,a14,2x,g10.8,2x,g14.8,2x,g8.5,2x,g10.4)',object[j], pix[j], line_name[j],cen_wl[j],str[j],velo,Tr
			printf, lun, format='(a10,2x,a8,2x,a14,2x,g10.8,2x,g14.8,2x,g8.5,2x,g10.4)',object[j], pix[j], line_name[j],cen_wl[j],str[j],velo,Tr
		endelse
	endfor
	free_lun,lun
	close,lun
end

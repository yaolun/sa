pro pop_dia_co, outdir,name, wl, flux, sigma, A, E_u, g, t_rot,no_fit=no_fit,legend=legend,two_comp=two_comp
;
;
;---------------------------------
angle = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, 0.5*COS(angle), 0.5*SIN(angle), /FILL
;---------------------------------
;
;plotsym, 0,/fill, color = 100
;units are in CGS unit.
k=1.38d-16
h=6.626d-27
c=3d10
v=c/(wl*1d-4)
p = 3.086
;name = ssw_name+'_'+slw_name

flux = flux*1e7 & sigma = sigma*1e7

N = 4*!PI*flux*(178*p)^2/(A*h*v)
N_sigma = 4*!PI*sigma*(178*p)^2/(A*h*v)
N_sig_hi = alog10((N+N_sigma)/g)-alog10(N/g)
n_sig_low = alog10(N/g)-alog10((N+N_sigma)/g)
;N1=8*!PI*k*(v^2)*flux/(h*(c^3)*A)
;N1_sigma=8*!PI*k*(v^2)*sigma/(h*(c^3)*A)
;g = (ind_utran)*2+1
if not keyword_set(legend) then legend=''
if not keyword_set(no_fit) then begin
;print, E_u[*], alog10(N[*]/g[*])+36
;Single temperature fitting
    start = fltarr(2)
    start[0] = (max(alog10(N/g))-min(alog10(N/g)))/(max(E_u)-min(E_u))
    start[1] = max(alog10(N/g))
    result = mpfitfun('base1d', E_u, alog10(N/g), N_sigma/(g*N), start, /quiet, perror=sigma, status=status, errmsg=errmsg)
    fit = result[0]*E_u+result[1]
    t_rot = -1/result[1]*alog10(exp(1)) & sig_t_rot = -(t_rot)*sigma[0]/result[0]*alog10(exp(1))
    two=0
    if keyword_set(two_comp) then begin
            ind1 = where(E_u gt two_comp) & ind2 = where(E_u le two_comp)
	if (n_elements(N[ind1]) ge 2) and (n_elements(N[ind2]) ge 2) then begin
	    ;result1 = linfit(E_u[ind1], alog10(N[ind1]/g[ind1]),measure_errors=N_sigma[ind1]/g[ind1]*N[ind1],/double, sigma=sigma1,chisq=chisq,yfit=fit1)
	    ;result2 = linfit(E_u[ind2], alog10(N[ind2]/g[ind2]),measure_errors=N_sigma[ind2]/g[ind2]*N[ind2],/double, sigma=sigma2,chisq=chisq,yfit=fit2)
	    start1 = fltarr(2) & start2 = fltarr(2)
	    start1[0] = (max(alog10(N[ind1]/g[ind1]))-min(alog10(N[ind1]/g[ind1])))/(max(E_u[ind1])-min(E_u[ind1]))
	    start1[1] = max(alog10(N[ind1]/g[ind1]))
	    start2[0] = (max(alog10(N[ind2]/g[ind2]))-min(alog10(N[ind2]/g[ind2])))/(max(E_u[ind2])-min(E_u[ind2]))
	    start2[1] = max(alog10(N[ind2]/g[ind2]))
	    y1 = alog10(N[ind1]/g[ind1]) & sig1 = N_sigma[ind1]/(g[ind1]*N[ind1])
	    y2 = alog10(N[ind2]/g[ind2]) & sig2 = N_sigma[ind2]/(g[ind2]*N[ind2])
	    result1 = mpfitfun('base1d', E_u[ind1], y1, sig1, start1, /quiet, perror=sigma1, status = status1, errmsg = errmsg1)
	    result2 = mpfitfun('base1d', E_u[ind2], y2, sig2, start2, /quiet, perror=sigma2, status = status2, errmsg = errmsg2)
	    fit1 = result1[0]*E_u[ind1]+result1[1]
	    fit2 = result2[0]*E_u[ind2]+result2[1]
	    t_rot1 = -1/result1[0]*alog10(exp(1)) & sig_t_rot1 = -(t_rot1)*sigma1[0]/result1[0]*alog10(exp(1))
	    t_rot2 = -1/result2[0]*alog10(exp(1)) & sig_t_rot2 = -(t_rot2)*sigma2[0]/result2[0]*alog10(exp(1))
	    ;t_rot1 = -1/result1[1] & t_rot2 = -1/result2[1]
	    ;stop
	    two=1
	endif
    endif
endif
set_plot, 'ps'
!p.font = -1
device, filename = outdir+'/plots/'+name+'_co.eps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 16, decomposed = 0, /color
loadct, 13
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
if not keyword_set(no_fit) then begin
  plot, E_u,alog10(N/g)+36, xtitle = 'E!du!n/k (K)', ytitle = 'log[!13N!3/g]', psym = 8, yrange = [43,50],xrange=[0,4500], charthick=4
  oploterror,E_u,alog10(N/g)+36,N_sig_hi, /hibar, psym = 8
  oploterror,E_u,alog10(N/g)+36,N_sig_low, /lowbar, psym=8
  oplot, E_u, fit+36, color = 80
  xyouts, 0.8, 0.85, legend, /normal, charthick=4
  xyouts, 0.77, 0.8, name, /normal,charsize=0.7, charthick=3
  xyouts, 0.7, 0.7, 'T!drot!n = '+strtrim(string(t_rot),1)+'K', /normal, charthick=3,charsize=0.7
endif else begin
  plot, E_u,alog10(N/g)+36, xtitle = 'E!du!n/k (K)', ytitle = 'log[!13N!3/g]', psym = 8, yrange = [43,50],xrange=[0,4500], charthick=4
  oploterror,E_u,alog10(N/g)+36,N_sig_hi, /hibar, psym = 8
  oploterror,E_u,alog10(N/g)+36,N_sig_low, /lowbar, psym=8
  xyouts, 0.8, 0.85, legend, /normal, charthick=3
  xyouts, 0.77, 0.8, name, /normal,charsize=0.7, charthick=3
endelse
device, /close_file, decomposed = 1
!p.multi = 0
set_plot, 'x'
if keyword_set(two_comp) and two eq 1 then begin
    set_plot,'ps'
    !p.font = -1
    device,filename=outdir+'/plots/'+name+'_two_co.eps',/helvetica,/portrait,/encapsulated,isolatin=1,font_size=16,decomposed=0,/color
    loadct,13
    !p.thick=4 & !x.thick=3 & !y.thick=3
    plot, E_u,alog10(N/g)+36,xtitle='E!du!n/k (K)',ytitle='log[!13N!3/g]',psym=8,yrange=[43,50],xrange=[0,4500],charthick=4
    oploterror,E_u,alog10(N/g)+36,N_sig_hi, /hibar, psym = 8
    oploterror,E_u,alog10(N/g)+36,N_sig_low, /lowbar, psym=8
    oplot, E_u[ind1], fit1+36, color=80
    oplot, E_u[ind2], fit2+36, color=80
    xyouts, 0.8, 0.85, legend, /normal, charthick=4
    xyouts, 0.77, 0.8, name, /normal, charsize=0.7,charthick=3
    t_hot = max([t_rot1,t_rot2]) & t_warm = min([t_rot1,t_rot2])
    xyouts, 0.7, 0.7, 'T(hot) ='+strtrim(string(floor(t_hot)),1)+'K+/-'+strtrim(string(floor(sig_t_rot1)),1)+'K', charsize=0.7,charthick=3,/normal
    xyouts, 0.7, 0.65, 'T(warm) ='+strtrim(string(floor(t_warm)),1)+'K+/-'+strtrim(string(floor(sig_t_rot2)),1)+'K', charsize=0.7, charthick=3,/normal
    device, /close_file, decomposed=1
    !p.multi=0
    set_plot, 'x'
endif
end

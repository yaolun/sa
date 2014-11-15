def match_spec(wl,flux,target):
	import numpy as np
	wl = np.array(wl).astype('float')
	flux = np.array(flux).astype('float')
	target = np.array(target).astype('float')
	tar_flux = np.empty([len(target)])
	for iwl in target:
		ind = np.where(target == iwl)
		wl_dum = wl - iwl
		if min(wl_dum) == 0:
			tar_flux[ind] = flux[wl_dum == 0]
		else:
			ind_dum = np.argsort(abs(wl_dum))[0:2]
			print ind_dum
			tar_flux[ind] = flux[wl == wl[ind_dum[0]]] + \
			(flux[wl == wl[ind_dum[1]]]-flux[wl == wl[ind_dum[0]]])/(wl[ind_dum[1]]-wl[ind_dum[0]])*(iwl-wl[ind_dum[0]])
	return tar_flux
wl = [1,2,3,4,5,6]
flux = [23,54,1,25,54,54]
target = [1,2.5,3.5,4,5.6]
print match_spec(wl,flux,target)

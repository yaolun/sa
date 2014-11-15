def read_sedfit(filename,outdir,num):
	import numpy as np
	import os
	data = np.genfromtxt(filename,skip_header=1,dtype=None)
	for i in range(0,num):
		print 'Model '+str(data[i][1])
		foo = open(outdir+'model_'+str(data[i][1])+'.dat','w')
		foo.write('# Parameters Setup \n# \n')
		foo.write('tstar \t \t %e \t # Stellar temperature         [K]\n' % data[i][8])
		foo.write('mstar \t \t %e \t # Stellar mass                [Solar Mass]\n' % data[i][6])
		foo.write('rstar \t \t %e \t # Stellar radius              [Solar radius]\n' % data[i][7])
		foo.write('M_env_dot \t %e \t # Envelope accretion rate     [Solar mass/yr]\n' % data[i][9])
		foo.write('M_disk_dot \t %e \t # Disk accretion rate         [Solar mass/yr]\n' % data[i][23])
		foo.write('R_env_max \t %e \t # Envelope outer radius       [AU]\n' % data[i][10])
		foo.write('R_env_min \t %e \t # Envelope inner radius       [AU]\n' % data[i][12])
		foo.write('theta_cav \t %f \t # Outflow cavity opening angle[deg]\n' % data[i][11])
		foo.write('R_disk_max \t %e \t # Disk outer radius           [AU]\n' % data[i][14])
		foo.write('R_disk_min \t %e \t # Disk inner radius           [AU]\n' % data[i][15])
		foo.write('M_disk \t \t %e \t # Disk mass                   [Solar mass]\n' % data[i][13])
		foo.write('beta \t \t %f \t # Disk flare factor           []\n' % data[i][19])
		foo.write('h100 \t \t %f \t # Disk scale height at 100 AU [AU]\n' % data[i][27])
		foo.write('rho_cav \t %e \t # Outflow cavity density      [g/cm3]\n' % data[i][21])
		foo.write('\n')
		foo.close()
# filename = '/home/ylyang/bhr71/data/results_full.txt'
filename = '/home/ylyang/bhr71/set_of_env_outer_radius.txt'
outdir = '/home/ylyang/bhr71/radmc3d_params/r_env_outer_v2/'
num = 4
read_sedfit(filename,outdir,num)

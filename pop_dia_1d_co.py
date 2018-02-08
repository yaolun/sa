def pop_dia_1d(objname,plotdir,dstar,fitting_table,
               single=False, opt_correction=None, fitting=True, full=False,
               spire_snr_threshold=4, fixed=False, write_out=False):
    """
    opt_correction is the largest J-level that need to be corrected from optical depth
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from read_fitting import read_fitting_co
    from leastsqfit import lin_leastsqfit
    import astropy.table as astal
    import astropy.constants as const
    import pprint
    from astropy.io import ascii

    # functions for further extract useful information for rotational diagram analysis
    def rot_lin_leastsqfit(x, y, y_sig):
        [yfit ,yerr ,a_hat, cov_hat, s_min] = lin_leastsqfit(x, y, y_sig)

        t_rot = -1/a_hat[1]*np.log10(np.e)
        sig_t_rot = -t_rot*cov_hat[1,1]**0.5/a_hat[1]*np.log10(np.e)
        yoff = a_hat[0]
        sig_yoff = cov_hat[0,0]**0.5

        return yfit ,yerr ,t_rot ,sig_t_rot ,s_min ,yoff, sig_yoff

    # Constants Setup
    c = const.c.cgs.value
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    B = 1.9225  # 192.25 m-1
    pc = const.pc.cgs.value

    # find CO fitting result
    data = ascii.read(fitting_table)
    if objname != 'RCrA-IRS7C':
        data = data[(data['Object'] == objname) & (data['Pixel_No.'] == 'c') & \
                    (data['Validity'] == 1)]
        ind_co = []

        for i in range(len(data)):
            if len(data['Line'][i].split('CO')[0]) == 0:
                # Adopt different SNR threshold for PACS and SPIRE data
                if data['ObsWL(um)'][i] > 195.05:
                    if data['SNR'][i] >= spire_snr_threshold:
                        ind_co.append(i)
                else:
                    if data['SNR'][i] >= 3:
                        ind_co.append(i)

        co_data = data[ind_co]
        co_data_name = data['Line'][ind_co]

    # average the CO line fluxes of both RCrA-IRS7B and 7C
    else:
        selector1 = (data['Object'] == objname) & (data['Pixel_No.'] == 'c') & \
                    (data['Validity'] == 1)
        selector2 = (data['Object'] == 'RCrA-IRS7B') & (data['Pixel_No.'] == 'c') & \
                    (data['Validity'] == 1)
        data = data[selector1+selector2]
        ind_co = []

        for i in range(len(data)):
            if len(data['Line'][i].split('CO')[0]) == 0:
                # Adopt different SNR threshold for PACS and SPIRE data
                if data['ObsWL(um)'][i] > 195.05:
                    if data['SNR'][i] >= spire_snr_threshold:
                        ind_co.append(i)
                else:
                    if data['SNR'][i] >= 3:
                        ind_co.append(i)

        # co_data = data[ind_co]
        data_dum = data[ind_co]
        co_data_name = np.array(list(set(data['Line'][ind_co].data)))

        # average the CO data
        co_data = {'Line': co_data_name,
                   'ObsWL(um)': np.empty(len(co_data_name)),
                   'Str(W/cm2)': np.empty(len(co_data_name)),
                   'FWHM(um)': np.empty(len(co_data_name)),
                   'Noise(W/cm2/um)': np.empty(len(co_data_name)),
                   'A(s-1)': np.empty(len(co_data_name)),
                   'E_u(K)': np.empty(len(co_data_name)),
                   'g': np.empty(len(co_data_name))}
        for i, co_line in enumerate(co_data['Line']):
            co_data['ObsWL(um)'][i] = np.mean(data_dum['ObsWL(um)'][data_dum['Line'] == co_line])
            co_data['Str(W/cm2)'][i] = np.mean(data_dum['Str(W/cm2)'][data_dum['Line'] == co_line])
            co_data['FWHM(um)'][i] = np.mean(data_dum['FWHM(um)'][data_dum['Line'] == co_line])
            co_data['Noise(W/cm2/um)'][i] = np.mean(data_dum['Noise(W/cm2/um)'][data_dum['Line'] == co_line])
            co_data['A(s-1)'][i] = np.mean(data_dum['A(s-1)'][data_dum['Line'] == co_line])
            co_data['E_u(K)'][i] = np.mean(data_dum['E_u(K)'][data_dum['Line'] == co_line])
            co_data['g'][i] = np.mean(data_dum['g'][data_dum['Line'] == co_line])


    if 'co_data_name' not in locals():
        return None
    if len(co_data_name) <= 2:
        return None

    # Calculate the N/g and Eu from the data
    v = c/(co_data['ObsWL(um)']*1e-4)
    N = 4*np.pi*co_data['Str(W/cm2)']*1e7*(dstar*pc)**2/(co_data['A(s-1)']*h*v)
    N_sigma = 4*np.pi*(1.064*co_data['FWHM(um)']*co_data['Noise(W/cm2/um)'])*1e7*(dstar*pc)**2/(co_data['A(s-1)']*h*v)


    if opt_correction is not None:
        # option for correcting optical depth if opt_correction is given with a J_up level
        # levels below (including the given one) will be performed correction

        # correct the 12CO line strength by applying optical depth derived from the 12CO/13CO analysis
        def tau_12co(E_u):
            p = [1.36547367, -0.00261421]
            return 10**(p[1]*E_u+p[0])*0.5

        cor_select = (co_data['g'] <= 2*opt_correction+1)
        str_corrected = co_data['Str(W/cm2)'][cor_select]*tau_12co(co_data['E_u(K)'][cor_select])/(1-np.exp(-tau_12co(co_data['E_u(K)'][cor_select])))
        str_unc_corrected = co_data['Noise(W/cm2/um)'][cor_select]*tau_12co(co_data['E_u(K)'][cor_select])/(1-np.exp(-tau_12co(co_data['E_u(K)'][cor_select])))

        N[cor_select] = 4*np.pi*str_corrected*1e7*(dstar*pc)**2/\
                        (co_data['A(s-1)'][cor_select]*h*v[cor_select])
        N_sigma[cor_select] = 4*np.pi*(1.064*co_data['FWHM(um)'][cor_select]*str_unc_corrected)*1e7*\
                              (dstar*pc)**2/(co_data['A(s-1)'][cor_select]*h*v[cor_select])

    # add systematic uncertainties from the difference between spectroscopy and photometry on top of that
    N_sigma[co_data['ObsWL(um)'] <= 85] = (N_sigma[co_data['ObsWL(um)'] <= 85]**2 + (0.16*N[co_data['ObsWL(um)'] <= 85])**2)**0.5
    N_sigma[(co_data['ObsWL(um)'] > 85) & (co_data['ObsWL(um)'] <= 125)] = \
            (N_sigma[(co_data['ObsWL(um)'] > 85) & (co_data['ObsWL(um)'] <= 125)]**2 + \
            (0.08*N[(co_data['ObsWL(um)'] > 85) & (co_data['ObsWL(um)'] <= 125)])**2)**0.5
    N_sigma[(co_data['ObsWL(um)'] > 125) & (co_data['ObsWL(um)'] <= 200)] = \
            (N_sigma[(co_data['ObsWL(um)'] > 125) & (co_data['ObsWL(um)'] <= 200)]**2 + \
            (0.10*N[(co_data['ObsWL(um)'] > 125) & (co_data['ObsWL(um)'] <= 200)])**2)**0.5
    N_sigma[(co_data['ObsWL(um)'] > 200) & (co_data['ObsWL(um)'] <= 300)] = \
            (N_sigma[(co_data['ObsWL(um)'] > 200) & (co_data['ObsWL(um)'] <= 300)]**2 + \
            (0.08*N[(co_data['ObsWL(um)'] > 200) & (co_data['ObsWL(um)'] <= 300)])**2)**0.5
    N_sigma[(co_data['ObsWL(um)'] > 300) & (co_data['ObsWL(um)'] <= 428)] = \
            (N_sigma[(co_data['ObsWL(um)'] > 300) & (co_data['ObsWL(um)'] <= 428)]**2 + \
            (0.03*N[(co_data['ObsWL(um)'] > 300) & (co_data['ObsWL(um)'] <= 428)])**2)**0.5
    N_sigma[co_data['ObsWL(um)'] > 428] = (N_sigma[co_data['ObsWL(um)'] > 428]**2 + (0.24*N[co_data['ObsWL(um)'] > 428])**2)**0.5


    x = co_data['E_u(K)']
    y = np.log10(N/co_data['g'])
    yerr_hi = np.log10((N+N_sigma)/co_data['g'])-np.log10(N/co_data['g'])
    yerr_low = np.log10(N/co_data['g'])-np.log10((N-N_sigma)/co_data['g'])
    y_sig = y*0
    for i in range(0,len(y)):
            y_sig[i] = max(yerr_hi[i], yerr_low[i])
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    y_sig = y_sig[ind]

    if not fitting:
        return x, y, y_sig

    # option to use a set of fixed break points
    if fixed:
        turning_pt = [199.11+1, 503.13+1, 1937.44+1]  # J_up = 26, 13, and 8.  Define the highest J element in a segament.

        # plot the data
        fig_rot_dia = plt.figure()
        ax_rot_dia = fig_rot_dia.add_subplot(111)
        data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
        ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

        # fit for four temperatures
        if len(x[x >= turning_pt[2]]) > 2:
            [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
            Q_hot  = float(k*t_rot_hot/h/c/B)
            N_hot_fit  = Q_hot *10**(float(yoff_hot))
            N_hot_fit_sig = [N_hot_fit-Q_hot*10**(float(yoff_hot-sig_yoff_hot)), Q_hot*10**(float(yoff_hot+sig_yoff_hot))-N_hot_fit]
            fit_hot,  = ax_rot_dia.plot(x[x>turning_pt[2]],yfit_hot,linewidth=1.5, color='Red')
            ax_rot_dia.fill_between(x[x>turning_pt[2]], yfit_hot-yerr_hot, yfit_hot+yerr_hot, color='Red', edgecolor='None', alpha=0.5)
        else:
            t_rot_hot = 0
            sig_t_rot_hot = 0
            N_hot_fit = 1
            N_hot_fit_sig = [1,1]

        if len(x[(x<turning_pt[2]) & (x>=turning_pt[1])]) > 2:
            [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
            Q_warm = float(k*t_rot_warm/h/c/B)
            N_warm_fit = Q_warm*10**(float(yoff_warm))
            N_warm_fit_sig = [N_warm_fit-Q_warm*10**(float(yoff_warm-sig_yoff_warm)), Q_warm*10**(float(yoff_warm+sig_yoff_warm))-N_warm_fit]
            fit_warm, = ax_rot_dia.plot(x[(x>=turning_pt[1]) & (x<turning_pt[2])],yfit_warm,linewidth=1.5, color='Orange')
            ax_rot_dia.fill_between(x[(x>=turning_pt[1]) & (x<turning_pt[2])], yfit_warm-yerr_warm, yfit_warm+yerr_warm, color='Orange', edgecolor='None', alpha=0.5)
        else:
            t_rot_warm = 0
            sig_t_rot_warm = 0
            N_warm_fit = 1
            N_warm_fit_sig = [1,1]

        if len(x[(x<turning_pt[1]) & (x>=turning_pt[0])]) > 2:
            [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
            Q_cool = float(k*t_rot_cool/h/c/B)
            N_cool_fit = Q_cool*10**(float(yoff_cool))
            N_cool_fit_sig = [N_cool_fit-Q_cool*10**(float(yoff_cool-sig_yoff_cool)), Q_cool*10**(float(yoff_cool+sig_yoff_cool))-N_cool_fit]
            fit_cool, = ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_cool,linewidth=1.5, color='DarkMagenta')
            ax_rot_dia.fill_between(x[(x>=turning_pt[0]) & (x<turning_pt[1])], yfit_cool-yerr_cool, yfit_cool+yerr_cool, color='DarkMagenta', edgecolor='None', alpha=0.5)
        else:
            t_rot_cool = 0
            sig_t_rot_cool = 0
            N_cool_fit = 1
            N_cool_fit_sig = [1,1]

        if len(x[x<turning_pt[0]]) > 2:
            [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
            Q_cold = float(k*t_rot_cold/h/c/B)
            N_cold_fit = Q_cold*10**(float(yoff_cold))
            N_cold_fit_sig = [N_cold_fit-Q_cold*10**(float(yoff_cold-sig_yoff_cold)), Q_cold*10**(float(yoff_cold+sig_yoff_cold))-N_cool_fit]
            fit_cold, = ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cold,linewidth=1.5, color='Blue')
            ax_rot_dia.fill_between(x[x<turning_pt[0]], yfit_cold-yerr_cold, yfit_cold+yerr_cold, color='Blue', edgecolor='None', alpha=0.5)
        else:
            t_rot_cold = 0
            sig_t_rot_cold = 0
            N_cold_fit = 1
            N_cold_fit_sig = [1,1]

        ax_rot_dia.set_xlabel(r'$\rm{E_{u}\,(K)}$',fontsize=18)
        ax_rot_dia.set_ylabel(r'$\rm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
        ax_rot_dia.set_xlim([0,x.max()+200])
        ax_rot_dia.set_ylim([42,52])

        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
        [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

        # ax_rot_dia.legend([fit_hot,fit_warm,fit_cool,fit_cold],\
        #     [r'$\rm{T_{rot,hot}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
        #      r'$\rm{T_{rot,warm}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
        #      r'$\rm{T_{rot,cool}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit))),\
        #      r'$\rm{T_{rot,cold}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cold,sig_t_rot_cold,N_cold_fit/10**np.floor(np.log10(N_cold_fit)),np.floor(np.log10(N_cold_fit)))],\
        #      numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)

        if full:
            # uncertainty on N
            ax_rot_dia.text(0.1,0.4, r'$\Delta \mathcal{N}_{\rm hot}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                            (N_hot_fit_sig[0]/10**np.floor(np.log10(N_hot_fit_sig[0])), np.floor(np.log10(N_hot_fit_sig[0])),\
                            N_hot_fit_sig[1]/10**np.floor(np.log10(N_hot_fit_sig[1])), np.floor(np.log10(N_hot_fit_sig[1]))),
                            transform=ax_rot_dia.transAxes, fontsize=10)
            ax_rot_dia.text(0.1,0.3, r'$\Delta \mathcal{N}_{\rm warm}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                            (N_warm_fit_sig[0]/10**np.floor(np.log10(N_warm_fit_sig[0])), np.floor(np.log10(N_warm_fit_sig[0])),\
                            N_warm_fit_sig[1]/10**np.floor(np.log10(N_warm_fit_sig[1])), np.floor(np.log10(N_warm_fit_sig[1]))),
                            transform=ax_rot_dia.transAxes, fontsize=10)
            ax_rot_dia.text(0.1,0.2, r'$\Delta \mathcal{N}_{\rm cool}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                            (N_cool_fit_sig[0]/10**np.floor(np.log10(N_cool_fit_sig[0])), np.floor(np.log10(N_cool_fit_sig[0])),\
                            N_cool_fit_sig[1]/10**np.floor(np.log10(N_cool_fit_sig[1])), np.floor(np.log10(N_cool_fit_sig[1]))),
                            transform=ax_rot_dia.transAxes, fontsize=10)
            ax_rot_dia.text(0.1,0.1, r'$\Delta \mathcal{N}_{\rm cold}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                            (N_cold_fit_sig[0]/10**np.floor(np.log10(N_cold_fit_sig[0])), np.floor(np.log10(N_cold_fit_sig[0])),\
                            N_cold_fit_sig[1]/10**np.floor(np.log10(N_cold_fit_sig[1])), np.floor(np.log10(N_cold_fit_sig[1]))),
                            transform=ax_rot_dia.transAxes, fontsize=10)

        fig_rot_dia.savefig(plotdir+objname+'_co_rot_fixfour.pdf',format='pdf',dpi=300, bbox_inches='tight')
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        print 'T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K, T_rot(cold): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool,t_rot_cold)
        output = (float(t_rot_hot), float(t_rot_warm), float(t_rot_cool), float(t_rot_cold),
                  float(sig_t_rot_hot), float(sig_t_rot_warm), float(sig_t_rot_cool), float(sig_t_rot_cold),
                  N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1], N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1],
                  N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1], N_cold_fit, N_cold_fit_sig[0], N_cold_fit_sig[1])


        # print out the fitted temperatures
        if opt_correction != None:
            suffix = '_optcor'
        else:
            suffix = ''
        if os.path.exists(plotdir+'co_rotational_fixtemp'+suffix+'.txt'):
            foo = open(plotdir+'co_rotational_fixtemp'+suffix+'.txt', 'a')
        else:
            foo = open(plotdir+'co_rotational_fixtemp'+suffix+'.txt', 'w')
            colhead = ['Object', 't_hot','t_warm','t_cool','t_cold','sig_t_hot','sig_t_warm','sig_t_cool','sig_t_cold',
                       'N_hot', 'N_hot_lowerr', 'N_hot_hierr', 'N_warm', 'N_warm_lowerr', 'N_warm_hierr',
                       'N_cool','N_cool_lowerr','N_cool_hierr','N_cold','N_cold_lowerr','N_cold_hierr']
            strformat = len(colhead)*'{:<13s}  '
            foo.write(strformat.format(*colhead))
            foo.write('\n')

        foo.write('{:<13s}  '.format(objname))
        foo.write('{:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}\n'.format(*output))
        foo.close()

        return None


    # collect the min-chisq for different kinds of temperature fitting
    s_min_total = []

    # Single temperature fitting
    #
    if len(x)>2:
        # [yfit, yerr, a_hat, cov_hat, s_min] = lin_leastsqfit(x, y, y_sig)
        # t_rot, sig_t_rot = getT_rot(a_hat, cov_hat)
        # yoff = a_hat[0]
        [yfit,yerr,t_rot,sig_t_rot,s_min,yoff,sig_yoff] = rot_lin_leastsqfit(x, y, y_sig)
        Q = float(k*t_rot/h/c/B)
        N_fit = Q*10**(float(yoff))
        N_fit_sig = [N_fit-Q*10**(float(yoff-sig_yoff)), Q*10**(float(yoff+sig_yoff))-N_fit]
        x = x.reshape(len(x))
        y = y.reshape(len(y))
        y_sig = y_sig.reshape(len(y_sig))
        s_min_total.append(s_min)

        fig_rot_dia = plt.figure()
        ax_rot_dia = fig_rot_dia.add_subplot(111)
        data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
        ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

        fit, = ax_rot_dia.plot(x,yfit,color='DarkMagenta', linewidth=1.5)
        ax_rot_dia.fill_between(x, yfit-yerr, yfit+yerr, facecolor='DarkMagenta', edgecolor='None', alpha=0.5)
        ax_rot_dia.set_xlabel(r'$\rm{E_{u}\,(K)}$',fontsize=18)
        ax_rot_dia.set_ylabel(r'$\rm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
        [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax_rot_dia.set_xlim([0,x.max()+200])
        ax_rot_dia.set_ylim([42,52])
        ax_rot_dia.legend([fit],[r'$\rm{T_{rot}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot,sig_t_rot,N_fit/10**np.floor(np.log10(N_fit)),np.floor(np.log10(N_fit)))],\
            numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
        if full:
            ax_rot_dia.text(0.1,0.1, r'$\Delta \mathcal{N}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                            (N_fit_sig[0]/10**np.floor(np.log10(N_fit_sig[0])), np.floor(np.log10(N_fit_sig[0])),\
                            N_fit_sig[1]/10**np.floor(np.log10(N_fit_sig[1])), np.floor(np.log10(N_fit_sig[1]))),
                            transform=ax_rot_dia.transAxes, fontsize=10)
        fig_rot_dia.savefig(plotdir+objname+'_co_rot_single.pdf',format='pdf',dpi=300, bbox_inches='tight')
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        print 'T_rot: %8.6f K' % t_rot
        # print t_rot, sig_t_rot, N_fit/10**np.floor(np.log10(N_fit)), np.floor(np.log10(N_fit))
        output = [0.0,0.0,0.0,float(t_rot.A1),
                  0.0,0.0,0.0,float(sig_t_rot.A1),
                  0.0,0.0,0.0,0.0,0.0,0.0,
                  0.0,0.0,0.0, N_fit, N_fit_sig[0], N_fit_sig[1]]
        s_min_single = s_min

    # Two temperature fitting
    #
        if len(x)>=8:
            best_fit = []
            s_min = []
            for i in range(2, len(x)-4):
                turning_pt = x[i]+1
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
                best_fit.append(turning_pt)
                s_min_dum = (s_min_warm*(len(x[x >= turning_pt])+2)+s_min_cool*(len(x[x < turning_pt])+2)) / (len(x)+5)
                s_min.append(s_min_dum)
            best_fit = np.array(best_fit)
            s_min = np.array(s_min)

            turning_pt = np.mean(best_fit[s_min == min(s_min)])
            [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
            [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
            s_min_double = (s_min_warm*(len(x[x >= turning_pt])+2)+s_min_cool*(len(x[x < turning_pt])+2)) / (len(x)+5)
            s_min_total.append(s_min_double)

            if s_min_double < s_min_single:
                Q_warm = float(k*t_rot_warm/h/c/B)
                Q_cool = float(k*t_rot_cool/h/c/B)
                N_warm_fit = Q_warm*10**(float(yoff_warm))
                N_cool_fit = Q_cool*10**(float(yoff_cool))
                N_warm_fit_sig = [N_warm_fit-Q_warm*10**(float(yoff_warm-sig_yoff_warm)), Q_warm*10**(float(yoff_warm+sig_yoff_warm))-N_warm_fit]
                N_cool_fit_sig = [N_cool_fit-Q_cool*10**(float(yoff_cool-sig_yoff_cool)), Q_cool*10**(float(yoff_cool+sig_yoff_cool))-N_cool_fit]

                fig_rot_dia = plt.figure()
                ax_rot_dia = fig_rot_dia.add_subplot(111)
                data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
                ax_rot_dia.errorbar(x,y,yerr=(yerr_low,yerr_hi),linestyle='None',color='DarkGreen')

                fit_warm, = ax_rot_dia.plot(x[x>=turning_pt],yfit_warm,color='DarkMagenta', linewidth=1.5)
                ax_rot_dia.fill_between(x[x>=turning_pt], yfit_warm-yerr_warm, yfit_warm+yerr_warm, facecolor='DarkMagenta', edgecolor='None', alpha=0.5)

                fit_cool, = ax_rot_dia.plot(x[x<turning_pt],yfit_cool,color='Blue', linewidth=1.5)
                ax_rot_dia.fill_between(x[x<turning_pt], yfit_cool-yerr_cool, yfit_cool+yerr_cool, facecolor='Blue', edgecolor='None', alpha=0.5)

                ax_rot_dia.set_xlabel(r'$\rm{E_{u}\,(K)}$',fontsize=18)
                ax_rot_dia.set_ylabel(r'$\rm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
                ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
                ax_rot_dia.set_xlim([0,x.max()+200])
                ax_rot_dia.set_ylim([42,52])

                ax_rot_dia.legend([fit_warm,fit_cool],\
                    [r'$\rm{T_{rot,warm}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\rm{T_{rot,cool}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
                if full:
                    # uncertainty on N
                    ax_rot_dia.text(0.1,0.2, r'$\Delta \mathcal{N}_{\rm warm}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                    (N_warm_fit_sig[0]/10**np.floor(np.log10(N_warm_fit_sig[0])), np.floor(np.log10(N_warm_fit_sig[0])),\
                                    N_warm_fit_sig[1]/10**np.floor(np.log10(N_warm_fit_sig[1])), np.floor(np.log10(N_warm_fit_sig[1]))),
                                    transform=ax_rot_dia.transAxes, fontsize=10)
                    ax_rot_dia.text(0.1,0.1, r'$\Delta \mathcal{N}_{\rm cool}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                    (N_cool_fit_sig[0]/10**np.floor(np.log10(N_cool_fit_sig[0])), np.floor(np.log10(N_cool_fit_sig[0])),\
                                    N_cool_fit_sig[1]/10**np.floor(np.log10(N_cool_fit_sig[1])), np.floor(np.log10(N_cool_fit_sig[1]))),
                                    transform=ax_rot_dia.transAxes, fontsize=10)

                fig_rot_dia.savefig(plotdir+objname+'_co_rot_two.pdf',format='pdf',dpi=300, bbox_inches='tight')
                ax_rot_dia.cla()
                fig_rot_dia.clf()
                print 'T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_warm,t_rot_cool)
                output = [0.0,0.0,float(t_rot_warm.A1), float(t_rot_cool.A1),
                          0.0,0.0,float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1),
                          0.0,0.0,0.0,0.0,0.0,0.0,
                          N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]
            else:
                print 'No smaller chi2 found for two temperature fitting.'
            # Three temperature fitting
            #
            if x.min() < 500 and len(x) > 12:
                best_fit = []
                s_min = []
                for i in range(2, len(x)-4):
                    for j in range(2, len(x)-4):
                        if (j-i) < 3:
                            continue
                        turning_pt = [x[i]+1,x[j]+1]

                        [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                        [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                        [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])

                        # make sure the temperatures are monotonically increasing
                        temp_list = [t_rot_hot, t_rot_warm, t_rot_cool]
                        if not all(xx > yy for xx, yy in zip(temp_list, temp_list[1:])):
                            continue

                        best_fit.append(turning_pt)
                        s_min_dum = (s_min_hot*(len(x[x >= turning_pt[1]])+2) +\
                                     s_min_warm*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) +\
                                     s_min_cool*(len(x[x < turning_pt[0]])+2))/(len(x)+8)
                        s_min.append(s_min_dum)

                best_fit = np.array(best_fit)
                s_min = np.array(s_min)
                turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1])]
                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                s_min_triple = (s_min_hot*(len(x[x >= turning_pt[1]])+2) +\
                                s_min_warm*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) +\
                                s_min_cool*(len(x[x < turning_pt[0]])+2))/(len(x)+8)
                s_min_total.append(s_min_triple)

                if s_min_triple < s_min_double:
                    Q_hot  = float(k*t_rot_hot/h/c/B)
                    Q_warm = float(k*t_rot_warm/h/c/B)
                    Q_cool = float(k*t_rot_cool/h/c/B)
                    N_hot_fit  = Q_hot *10**(float(yoff_hot))
                    N_warm_fit = Q_warm*10**(float(yoff_warm))
                    N_cool_fit = Q_cool*10**(float(yoff_cool))

                    N_hot_fit_sig = [N_hot_fit-Q_hot*10**(float(yoff_hot-sig_yoff_hot)), Q_hot*10**(float(yoff_hot+sig_yoff_hot))-N_hot_fit]
                    N_warm_fit_sig = [N_warm_fit-Q_warm*10**(float(yoff_warm-sig_yoff_warm)), Q_warm*10**(float(yoff_warm+sig_yoff_warm))-N_warm_fit]
                    N_cool_fit_sig = [N_cool_fit-Q_cool*10**(float(yoff_cool-sig_yoff_cool)), Q_cool*10**(float(yoff_cool+sig_yoff_cool))-N_cool_fit]

                    fig_rot_dia = plt.figure()
                    ax_rot_dia = fig_rot_dia.add_subplot(111)
                    data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
                    ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

                    fit_hot,  = ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot,color='Red', linewidth=1.5)
                    ax_rot_dia.fill_between(x[x>turning_pt[1]], yfit_hot-yerr_hot, yfit_hot+yerr_hot, facecolor='Red', edgecolor='None', alpha=0.5)

                    fit_warm, = ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm,color='DarkMagenta', linewidth=1.5)
                    ax_rot_dia.fill_between(x[(x>=turning_pt[0]) & (x<turning_pt[1])], yfit_warm-yerr_warm, yfit_warm+yerr_warm, facecolor='DarkMagenta', edgecolor='None', alpha=0.5)

                    fit_cool, = ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool,color='Blue', linewidth=1.5)
                    ax_rot_dia.fill_between(x[x<turning_pt[0]], yfit_cool-yerr_cool, yfit_cool+yerr_cool, facecolor='Blue', edgecolor='None', alpha=0.5)

                    ax_rot_dia.set_xlabel(r'$\rm{E_{u}\,(K)}$',fontsize=18)
                    ax_rot_dia.set_ylabel(r'$\rm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
                    ax_rot_dia.set_xlim([0,x.max()+200])
                    ax_rot_dia.set_ylim([42,52])

                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                    [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

                    ax_rot_dia.legend([fit_hot,fit_warm,fit_cool],\
                        [r'$\rm{T_{rot,warm}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
                         r'$\rm{T_{rot,cool}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\rm{T_{rot,cold}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)

                    if full:
                        # uncertainty on N
                        ax_rot_dia.text(0.1,0.3, r'$\Delta \mathcal{N}_{\rm hot}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                        (N_hot_fit_sig[0]/10**np.floor(np.log10(N_hot_fit_sig[0])), np.floor(np.log10(N_hot_fit_sig[0])),\
                                        N_hot_fit_sig[1]/10**np.floor(np.log10(N_hot_fit_sig[1])), np.floor(np.log10(N_hot_fit_sig[1]))),
                                        transform=ax_rot_dia.transAxes, fontsize=10)
                        ax_rot_dia.text(0.1,0.2, r'$\Delta \mathcal{N}_{\rm warm}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                        (N_warm_fit_sig[0]/10**np.floor(np.log10(N_warm_fit_sig[0])), np.floor(np.log10(N_warm_fit_sig[0])),\
                                        N_warm_fit_sig[1]/10**np.floor(np.log10(N_warm_fit_sig[1])), np.floor(np.log10(N_warm_fit_sig[1]))),
                                        transform=ax_rot_dia.transAxes, fontsize=10)
                        ax_rot_dia.text(0.1,0.1, r'$\Delta \mathcal{N}_{\rm cool}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                        (N_cool_fit_sig[0]/10**np.floor(np.log10(N_cool_fit_sig[0])), np.floor(np.log10(N_cool_fit_sig[0])),\
                                        N_cool_fit_sig[1]/10**np.floor(np.log10(N_cool_fit_sig[1])), np.floor(np.log10(N_cool_fit_sig[1]))),
                                        transform=ax_rot_dia.transAxes, fontsize=10)

                    fig_rot_dia.savefig(plotdir+objname+'_co_rot_three.pdf',format='pdf',dpi=300, bbox_inches='tight')
                    ax_rot_dia.cla()
                    fig_rot_dia.clf()
                    print 'T_rot(hot): %5.1f K, T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_hot,t_rot_warm,t_rot_cool)
                    output = [0.0, float(t_rot_hot.A1), float(t_rot_warm.A1), float(t_rot_cool.A1),
                              0.0, float(sig_t_rot_hot.A1), float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1),
                              0.0,0.0,0.0, N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1],
                              N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]

                else:
                    print 'No smaller chi2 found for three temperature fitting.'

                # four temperature fitting
                if (x.min() < 500) and (len(x) > 16):
                    best_fit = []
                    s_min = []
                    for i in range(2, len(x)-4):
                        for j in range(2, len(x)-4):
                            for jj in range(2, len(x)-4):
                                if (j-i < 3) or (jj-j < 3) or (jj-i < 6):
                                    continue
                                turning_pt = [x[i]+1, x[j]+1, x[jj]+1]
                                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
                                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
                                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                                [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])

                                # make sure the temperatures are monotonically increasing
                                temp_list = [t_rot_hot, t_rot_warm, t_rot_cool, t_rot_cold]
                                if not all(xx > yy for xx, yy in zip(temp_list, temp_list[1:])):
                                    continue
                                best_fit.append(turning_pt)
                                s_min_dum = (s_min_hot*(len(x[x >= turning_pt[2]])+2) + \
                                             s_min_warm*(len(x[(x < turning_pt[2]) & (x >= turning_pt[1])])+2) + \
                                             s_min_cool*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) + \
                                             s_min_cold*(len(x[x < turning_pt[0]])+2)) / (len(x)+11)
                                s_min.append(s_min_dum)

                    best_fit = np.array(best_fit)
                    s_min = np.array(s_min)
                    # if objname == 'BHR71':
                    #     print s_min.min()
                    #     dum = np.argsort(s_min)
                    #     print s_min[dum[0:6]]
                    #     print best_fit[s_min == min(s_min),0],best_fit[s_min == min(s_min),1],best_fit[s_min == min(s_min),2]
                    #     print best_fit[s_min == s_min[dum[1]],0],best_fit[s_min == s_min[dum[1]],1],best_fit[s_min == s_min[dum[1]],2]
                    #     print best_fit[s_min == s_min[dum[2]],0],best_fit[s_min == s_min[dum[2]],1],best_fit[s_min == s_min[dum[2]],2]
                    #
                    #     turning_pt = [ 155.87, 846.59, 1795.23]
                    #     [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
                    #     [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
                    #     [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                    #     [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                    #     print s_min_hot+s_min_warm+s_min_cool+s_min_cold
                    #     print t_rot_hot, t_rot_warm, t_rot_cool, t_rot_cold

                    turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1]),np.mean(best_fit[s_min == min(s_min),2])]
                    [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
                    [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
                    [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                    [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                    s_min_quad = (s_min_hot*(len(x[x >= turning_pt[2]])+2) + \
                                  s_min_warm*(len(x[(x < turning_pt[2]) & (x >= turning_pt[1])])+2) + \
                                  s_min_cool*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) + \
                                  s_min_cold*(len(x[x < turning_pt[0]])+2)) / (len(x)+11)
                    s_min_total.append(s_min_quad)

                    if (t_rot_hot < 0) or (t_rot_warm < 0) or (t_rot_cool < 0) or (t_rot_cold < 0):
                        print 'negative temperature found'

                        if write_out:
                            # print out the fitted temperatures
                            if opt_correction != None:
                                suffix = '_optcor'
                            else:
                                suffix = ''
                            if os.path.exists(plotdir+'co_rotational_temp'+suffix+'.txt'):
                                foo = open(plotdir+'co_rotational_temp'+suffix+'.txt', 'a')
                            else:
                                foo = open(plotdir+'co_rotational_temp'+suffix+'.txt', 'w')
                                colhead = ['Object', 't_hot','t_warm','t_cool','t_cold','sig_t_hot','sig_t_warm','sig_t_cool','sig_t_cold',
                                           'N_hot', 'N_hot_lowerr', 'N_hot_hierr', 'N_warm', 'N_warm_lowerr', 'N_warm_hierr',
                                           'N_cool','N_cool_lowerr','N_cool_hierr','N_cold','N_cold_lowerr','N_cold_hierr']
                                strformat = len(colhead)*'{:<13s}  '
                                foo.write(strformat.format(*colhead))
                                foo.write('\n')

                            foo.write('{:<13s}  '.format(objname))
                            foo.write('{:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}\n'.format(*output))
                            foo.close()

                            # print out the break points
                            if 'turning_pt' in locals():
                                if os.path.exists(plotdir+'co_breakpoints'+suffix+'.txt'):
                                    foo = open(plotdir+'co_breakpoints'+suffix+'.txt', 'a')
                                else:
                                    foo = open(plotdir+'co_breakpoints'+suffix+'.txt', 'w')

                                if type(turning_pt) != list:
                                    foo.write(str(turning_pt)+', ')
                                else:
                                    for pt in list(turning_pt):
                                        foo.write(str(pt)+', ')
                                foo.write('\n')
                                foo.close()

                        return None

                    if s_min_quad < s_min_triple:
                        Q_hot  = float(k*t_rot_hot/h/c/B)
                        Q_warm = float(k*t_rot_warm/h/c/B)
                        Q_cool = float(k*t_rot_cool/h/c/B)
                        Q_cold = float(k*t_rot_cold/h/c/B)
                        N_hot_fit  = Q_hot *10**(float(yoff_hot))
                        N_warm_fit = Q_warm*10**(float(yoff_warm))
                        N_cool_fit = Q_cool*10**(float(yoff_cool))
                        N_cold_fit = Q_cold*10**(float(yoff_cold))

                        N_hot_fit_sig = [N_hot_fit-Q_hot*10**(float(yoff_hot-sig_yoff_hot)), Q_hot*10**(float(yoff_hot+sig_yoff_hot))-N_hot_fit]
                        N_warm_fit_sig = [N_warm_fit-Q_warm*10**(float(yoff_warm-sig_yoff_warm)), Q_warm*10**(float(yoff_warm+sig_yoff_warm))-N_warm_fit]
                        N_cool_fit_sig = [N_cool_fit-Q_cool*10**(float(yoff_cool-sig_yoff_cool)), Q_cool*10**(float(yoff_cool+sig_yoff_cool))-N_cool_fit]
                        N_cold_fit_sig = [N_cold_fit-Q_cold*10**(float(yoff_cold-sig_yoff_cold)), Q_cold*10**(float(yoff_cold+sig_yoff_cold))-N_cool_fit]

                        fig_rot_dia = plt.figure()
                        ax_rot_dia = fig_rot_dia.add_subplot(111)
                        data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
                        ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')

                        fit_hot,  = ax_rot_dia.plot(x[x>turning_pt[2]],yfit_hot,linewidth=1.5, color='Red')
                        ax_rot_dia.fill_between(x[x>turning_pt[2]], yfit_hot-yerr_hot, yfit_hot+yerr_hot, color='Red', edgecolor='None', alpha=0.5)

                        fit_warm, = ax_rot_dia.plot(x[(x>=turning_pt[1]) & (x<turning_pt[2])],yfit_warm,linewidth=1.5, color='Orange')
                        ax_rot_dia.fill_between(x[(x>=turning_pt[1]) & (x<turning_pt[2])], yfit_warm-yerr_warm, yfit_warm+yerr_warm, color='Orange', edgecolor='None', alpha=0.5)

                        fit_cool, = ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_cool,linewidth=1.5, color='DarkMagenta')
                        ax_rot_dia.fill_between(x[(x>=turning_pt[0]) & (x<turning_pt[1])], yfit_cool-yerr_cool, yfit_cool+yerr_cool, color='DarkMagenta', edgecolor='None', alpha=0.5)

                        fit_cold, = ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cold,linewidth=1.5, color='Blue')
                        ax_rot_dia.fill_between(x[x<turning_pt[0]], yfit_cold-yerr_cold, yfit_cold+yerr_cold, color='Blue', edgecolor='None', alpha=0.5)

                        ax_rot_dia.set_xlabel(r'$\rm{E_{u}\,(K)}$',fontsize=18)
                        ax_rot_dia.set_ylabel(r'$\rm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
                        ax_rot_dia.set_xlim([0,x.max()+200])
                        ax_rot_dia.set_ylim([42,52])

                        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                        [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

                        ax_rot_dia.legend([fit_hot,fit_warm,fit_cool,fit_cold],\
                            [r'$\rm{T_{rot,hot}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
                             r'$\rm{T_{rot,warm}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                             r'$\rm{T_{rot,cool}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit))),\
                             r'$\rm{T_{rot,cold}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cold,sig_t_rot_cold,N_cold_fit/10**np.floor(np.log10(N_cold_fit)),np.floor(np.log10(N_cold_fit)))],\
                             numpoints=1,loc='upper right',fontsize=12,framealpha=0.3)

                        if full:
                            # uncertainty on N
                            ax_rot_dia.text(0.1,0.4, r'$\Delta \mathcal{N}_{\rm hot}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                            (N_hot_fit_sig[0]/10**np.floor(np.log10(N_hot_fit_sig[0])), np.floor(np.log10(N_hot_fit_sig[0])),\
                                            N_hot_fit_sig[1]/10**np.floor(np.log10(N_hot_fit_sig[1])), np.floor(np.log10(N_hot_fit_sig[1]))),
                                            transform=ax_rot_dia.transAxes, fontsize=10)
                            ax_rot_dia.text(0.1,0.3, r'$\Delta \mathcal{N}_{\rm warm}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                            (N_warm_fit_sig[0]/10**np.floor(np.log10(N_warm_fit_sig[0])), np.floor(np.log10(N_warm_fit_sig[0])),\
                                            N_warm_fit_sig[1]/10**np.floor(np.log10(N_warm_fit_sig[1])), np.floor(np.log10(N_warm_fit_sig[1]))),
                                            transform=ax_rot_dia.transAxes, fontsize=10)
                            ax_rot_dia.text(0.1,0.2, r'$\Delta \mathcal{N}_{\rm cool}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                            (N_cool_fit_sig[0]/10**np.floor(np.log10(N_cool_fit_sig[0])), np.floor(np.log10(N_cool_fit_sig[0])),\
                                            N_cool_fit_sig[1]/10**np.floor(np.log10(N_cool_fit_sig[1])), np.floor(np.log10(N_cool_fit_sig[1]))),
                                            transform=ax_rot_dia.transAxes, fontsize=10)
                            ax_rot_dia.text(0.1,0.1, r'$\Delta \mathcal{N}_{\rm cold}=- %3.2f \times 10^{%d}/+ %3.2f \times 10^{%d}$' % \
                                            (N_cold_fit_sig[0]/10**np.floor(np.log10(N_cold_fit_sig[0])), np.floor(np.log10(N_cold_fit_sig[0])),\
                                            N_cold_fit_sig[1]/10**np.floor(np.log10(N_cold_fit_sig[1])), np.floor(np.log10(N_cold_fit_sig[1]))),
                                            transform=ax_rot_dia.transAxes, fontsize=10)

                        fig_rot_dia.savefig(plotdir+objname+'_co_rot_four.pdf',format='pdf',dpi=300, bbox_inches='tight')
                        ax_rot_dia.cla()
                        fig_rot_dia.clf()
                        print 'T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K, T_rot(cold): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool,t_rot_cold)
                        output = (float(t_rot_hot.A1), float(t_rot_warm.A1), float(t_rot_cool.A1), float(t_rot_cold.A1),
                                  float(sig_t_rot_hot.A1), float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1), float(sig_t_rot_cold.A1),
                                  N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1], N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1],
                                  N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1], N_cold_fit, N_cold_fit_sig[0], N_cold_fit_sig[1])
                    else:
                        print 'No smaller chi2 found for four temperature fitting.'

        print s_min_total

        if write_out:
            # print out the fitted temperatures
            if opt_correction != None:
                suffix = '_optcor'
            else:
                suffix = ''
            if os.path.exists(plotdir+'co_rotational_temp'+suffix+'.txt'):
                foo = open(plotdir+'co_rotational_temp'+suffix+'.txt', 'a')
            else:
                foo = open(plotdir+'co_rotational_temp'+suffix+'.txt', 'w')
                colhead = ['Object', 't_hot','t_warm','t_cool','t_cold','sig_t_hot','sig_t_warm','sig_t_cool','sig_t_cold',
                           'N_hot', 'N_hot_lowerr', 'N_hot_hierr', 'N_warm', 'N_warm_lowerr', 'N_warm_hierr',
                           'N_cool','N_cool_lowerr','N_cool_hierr','N_cold','N_cold_lowerr','N_cold_hierr']
                strformat = len(colhead)*'{:<13s}  '
                foo.write(strformat.format(*colhead))
                foo.write('\n')

            foo.write('{:<13s}  '.format(objname))
            foo.write('{:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13f}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}  {:<13.4e}\n'.format(*output))
            foo.close()

            # print out the break points
            if 'turning_pt' in locals():
                if os.path.exists(plotdir+'co_breakpoints'+suffix+'.txt'):
                    foo = open(plotdir+'co_breakpoints'+suffix+'.txt', 'a')
                else:
                    foo = open(plotdir+'co_breakpoints'+suffix+'.txt', 'w')

                if type(turning_pt) != list:
                    foo.write(str(turning_pt)+', ')
                else:
                    for pt in list(turning_pt):
                        foo.write(str(pt)+', ')
                foo.write('\n')
                foo.close()


        return None

import numpy as np
import matplotlib.pyplot as plt
import os
from read_fitting import read_fitting_co
from leastsqfit import lin_leastsqfit
home = os.path.expanduser('~')

pacs = '/bhr71/best_calibrated/fitting/pacs/BHR71_pacs_weighted_lines.txt'
spire = '/bhr71/best_calibrated/fitting/spire/BHR71_spire_corrected_lines.txt'
# fitting_table = '/Users/yaolun/data/CDF_archive_v2/CDF_archive_v2_lines.txt'
fitting_table = '/Volumes/SD-Mac/CDF_archive_v2/CDF_archive_v2_lines.txt'
spire_snr_threshold = 4

# get the entire object list
from astropy.io import ascii
data = ascii.read(fitting_table)
obj_list = list(set(data['Object'].data))

dist = ascii.read('/Users/yaolun/data/cops-spire_distance.txt')

# loop through all objects
# with optical depth correction
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_rotational_temp_optcor.txt'):
    # os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_rotational_temp_optcor.txt')
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_breakpoints.txt'):
#     os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_breakpoints.txt')
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_breakpoints_optcor.txt'):
    # os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_breakpoints_optcor.txt')
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_rotational_temp.txt'):
#     os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_rotational_temp.txt')
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_rotational_fixtemp.txt'):
#     os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/co_rotational_fixtemp.txt')
# if os.path.exists('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_rotational_fixtemp_optcor.txt'):
#     os.remove('/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/co_rotational_fixtemp_optcor.txt')


for o in obj_list:
    print o
    pop_dia_1d(o, '/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/80percent/',
               dist['distance'][dist['object'] == o],
               fitting_table, opt_correction=13, spire_snr_threshold=spire_snr_threshold, write_out=True)
# # without optical depth correction
# for o in obj_list:
#     print o
#     pop_dia_1d(o, '/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/',
#                dist['distance'][dist['object'] == o],
#                fitting_table, opt_correction=None, spire_snr_threshold=spire_snr_threshold, write_out=True)
# # fix the breakpoints and with the optical correction
# for o in obj_list:
#     print o
#     pop_dia_1d(o, '/Volumes/SD-Mac/research/cops-spire/rotational_diagrams/corrected/',
#                dist['distance'][dist['object'] == o],
#                fitting_table, opt_correction=13, spire_snr_threshold=spire_snr_threshold, fixed=True, write_out=True)



# pop_dia_1d('BHR71','/test/',200.,pacs=pacs,spire=spire)
# pop_dia_1d('BHR71', '/test/', 200., fitting_table)
# pacs_cube = '/bhr71/data/HSA/cube/BHR71_pacs_pixel'
# for i in range(1,26):
#     print i
#     pop_dia_1d('BHR71_pacs_pixel'+str(i),'/bhr71/plots/',pacs=pacs_cube+str(i)+'_os8_sf7_lines.txt')

# pop_dia_h2o_1d('BHR71','/bhr71/plots/',pacs=pacs,spire=spire)



# pacs = '/bhr71/fitting/latest/pacs/advanced_products/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
# [co_pacs,co_name_pacs] = read_fitting_co(pacs,3)
# spire = '/bhr71/fitting/latest/spire/advanced_products/BHR71_spire_corrected_lines.txt'
# [co_spire,co_name_spire] = read_fitting_co(spire,3)

def pop_dia_1d(objname,plotdir,dstar,line_data, species_info,
               opt_correction=None, fitting=True, polyfit=False):
    """
    line_data: selected line fitting results with typical header as the fitting table
    species_info: A dictionary contains the neccessary info for constructin rotational diagram
                  This script currently only deals with linear molecules such as CO.
                  Thus, only the rotational constant, B, is accepted.  However, more information
                  provided won't have any negative effect.
    opt_correction should be a function that will give the proper correction factor as a function
    of wavelength.  This has not been implemented yet.
    polyfit: use the numpy.polyfit for least-sqaure fitting.  There is a case when I get too small
             uncertainty from my own version of fitting.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from leastsqfit import lin_leastsqfit
    import astropy.table as astal
    import astropy.constants as const
    import pprint

    from astropy.io import ascii

    # functions for further extract useful information for rotational diagram analysis
    def rot_lin_leastsqfit(x, y, y_sig, polyfit=False):
        if not polyfit:
            [yfit ,yerr ,a_hat, cov_hat, s_min] = lin_leastsqfit(x, y, y_sig)
            t_rot = -1/a_hat[1]*np.log10(np.e)
            sig_t_rot = -t_rot*cov_hat[1,1]**0.5/a_hat[1]*np.log10(np.e)
            yoff = a_hat[0]
            sig_yoff = cov_hat[0,0]**0.5
        else:
            p, cov = np.polyfit(x.data.flatten(), y.data, 1, w=1/y_sig.data, cov=True)
            yfit = p[1]+p[0]*x.data.flatten()
            yerr = (abs(cov[1,1])+2*cov[0,1]*x.data.flatten()+abs(cov[0,0])*x.data.flatten()**2)**0.5
            s_min = np.sum((y-yfit)**2)
            t_rot = -1/p[0]*np.log10(np.e)
            sig_t_rot = -t_rot*abs(cov[0,0])**0.5/p[0]*np.log10(np.e)
            yoff = p[1]
            sig_yoff = abs(cov[1,1])**0.5

        return yfit ,yerr ,t_rot ,sig_t_rot ,s_min ,yoff, sig_yoff

    home = os.path.expanduser('~')

    # Constants Setup
    c = const.c.cgs.value
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    pc = const.pc.cgs.value

    # species information
    B = species_info['B']
    species_name = species_info['name']

    # Calculate the N/g and Eu from the data
    v = c/(line_data['ObsWL(um)']*1e-4)
    N = 4*np.pi*line_data['Str(W/cm2)']*1e7*(dstar*pc)**2/(line_data['A(s-1)']*h*v)
    N_sigma = 4*np.pi*(1.064*line_data['FWHM(um)']*line_data['Noise(W/cm2/um)'])*1e7*(dstar*pc)**2\
              /(line_data['A(s-1)']*h*v)

    if opt_correction is not None:

        str_corrected = line_data['Str(W/cm2)']*opt_correction(line_data['Wavelength(um)'])/(1-np.exp(-opt_correction(line_data['Wavelength(um)'])))
        str_unc_corrected = line_data['Noise(W/cm2/um)']*opt_correction(line_data['Wavelength(um)'])/(1-np.exp(-opt_correction(line_data['Wavelength(um)'])))

        N = 4*np.pi*str_corrected*1e7*(dstar*pc)**2/\
                        (line_data['A(s-1)']*h*v)
        N_sigma = 4*np.pi*(1.064*line_data['FWHM(um)']*str_unc_corrected)*1e7*\
                              (dstar*pc)**2/(line_data['A(s-1)']*h*v)

    # add systematic uncertainties from the difference between spectroscopy and photometry on top of that
    N_sigma[line_data['ObsWL(um)'] <= 85] = (N_sigma[line_data['ObsWL(um)'] <= 85]**2 + (0.16*N[line_data['ObsWL(um)'] <= 85])**2)**0.5
    N_sigma[(line_data['ObsWL(um)'] > 85) & (line_data['ObsWL(um)'] <= 125)] = \
            (N_sigma[(line_data['ObsWL(um)'] > 85) & (line_data['ObsWL(um)'] <= 125)]**2 + \
            (0.08*N[(line_data['ObsWL(um)'] > 85) & (line_data['ObsWL(um)'] <= 125)])**2)**0.5
    N_sigma[(line_data['ObsWL(um)'] > 125) & (line_data['ObsWL(um)'] <= 200)] = \
            (N_sigma[(line_data['ObsWL(um)'] > 125) & (line_data['ObsWL(um)'] <= 200)]**2 + \
            (0.10*N[(line_data['ObsWL(um)'] > 125) & (line_data['ObsWL(um)'] <= 200)])**2)**0.5
    N_sigma[(line_data['ObsWL(um)'] > 200) & (line_data['ObsWL(um)'] <= 300)] = \
            (N_sigma[(line_data['ObsWL(um)'] > 200) & (line_data['ObsWL(um)'] <= 300)]**2 + \
            (0.08*N[(line_data['ObsWL(um)'] > 200) & (line_data['ObsWL(um)'] <= 300)])**2)**0.5
    N_sigma[(line_data['ObsWL(um)'] > 300) & (line_data['ObsWL(um)'] <= 428)] = \
            (N_sigma[(line_data['ObsWL(um)'] > 300) & (line_data['ObsWL(um)'] <= 428)]**2 + \
            (0.03*N[(line_data['ObsWL(um)'] > 300) & (line_data['ObsWL(um)'] <= 428)])**2)**0.5
    N_sigma[line_data['ObsWL(um)'] > 428] = (N_sigma[line_data['ObsWL(um)'] > 428]**2 + (0.24*N[line_data['ObsWL(um)'] > 428])**2)**0.5

    x = line_data['E_u(K)']
    y = np.log10(N/line_data['g'])
    yerr_hi = np.log10((N+N_sigma)/line_data['g'])-np.log10(N/line_data['g'])
    yerr_low = np.log10(N/line_data['g'])-np.log10((N-N_sigma)/line_data['g'])
    y_sig = y*0
    for i in range(0,len(y)):
        y_sig[i] = max(yerr_hi[i], yerr_low[i])
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    y_sig = y_sig[ind]

    if not fitting:
        return x, y, y_sig

    # collect the min-chisq for different kinds of temperature fitting
    s_min_total = []

    # Single temperature fitting
    #
    if len(x) > 2:
        [yfit, yerr, t_rot, sig_t_rot, s_min, yoff, sig_yoff] = rot_lin_leastsqfit(x, y, y_sig, polyfit=polyfit)

        Q = float(k*t_rot/h/B)
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
        # ax_rot_dia.set_xlim([0,x.max()+200])
        # ax_rot_dia.set_ylim([42,50])
        ax_rot_dia.legend([fit],
        [r'$\rm{T_{rot}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot,sig_t_rot,N_fit/10**np.floor(np.log10(N_fit)),np.floor(np.log10(N_fit)))],\
            numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
        fig_rot_dia.savefig(home+plotdir+objname+'_'+species_name+'_rot_single.pdf',format='pdf',dpi=300, bbox_inches='tight')
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        print('T_rot: %8.6f K' % t_rot)
        s_min_single = s_min

        if polyfit:
            output = [0.0,0.0,0.0,float(t_rot),
                      0.0,0.0,0.0,float(sig_t_rot),
                      0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0, N_fit, N_fit_sig[0], N_fit_sig[1]]
        else:
            output = [0.0,0.0,0.0,float(t_rot.A1),
                      0.0,0.0,0.0,float(sig_t_rot.A1),
                      0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0, N_fit, N_fit_sig[0], N_fit_sig[1]]

    # Two temperature fitting
    #
        if len(x) >= 8:
            best_fit = []
            s_min = []
            for i in range(2, len(x)-4):
                turning_pt = x[i]+1
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt], polyfit=polyfit)
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt], polyfit=polyfit)
                best_fit.append(turning_pt)
                s_min_dum = (s_min_warm*(len(x[x >= turning_pt])+2)+s_min_cool*(len(x[x < turning_pt])+2)) / (len(x)+5)
                s_min.append(s_min_dum)
            best_fit = np.array(best_fit)
            s_min = np.array(s_min)

            turning_pt = np.mean(best_fit[s_min == min(s_min)])

            [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt], polyfit=polyfit)
            [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt], polyfit=polyfit)
            s_min_double = (s_min_warm*(len(x[x >= turning_pt])+2)+s_min_cool*(len(x[x < turning_pt])+2)) / (len(x)+5)
            s_min_total.append(s_min_double)

            if (s_min_double)<s_min_single:
                Q_warm = float(k*t_rot_warm/h/B)
                Q_cool = float(k*t_rot_cool/h/B)
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

                ax_rot_dia.legend([fit_warm,fit_cool],\
                    [r'$\rm{T_{rot,warm}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\rm{T_{rot,cool}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
                fig_rot_dia.savefig(home+plotdir+objname+'_'+species_name+'_rot_two.pdf',format='pdf',dpi=300, bbox_inches='tight')
                ax_rot_dia.cla()
                fig_rot_dia.clf()
                print('T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_warm,t_rot_cool))

                if polyfit:
                    output = [0.0,0.0,float(t_rot_warm), float(t_rot_cool),
                              0.0,0.0,float(sig_t_rot_warm), float(sig_t_rot_cool),
                              0.0,0.0,0.0,0.0,0.0,0.0,
                              N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]
                else:
                    output = [0.0,0.0,float(t_rot_warm.A1), float(t_rot_cool.A1),
                              0.0,0.0,float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1),
                              0.0,0.0,0.0,0.0,0.0,0.0,
                              N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]

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

                        [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]], polyfit=polyfit)
                        [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])], polyfit=polyfit)
                        [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]], polyfit=polyfit)

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
                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]], polyfit=polyfit)
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])], polyfit=polyfit)
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]], polyfit=polyfit)
                s_min_triple = (s_min_hot*(len(x[x >= turning_pt[1]])+2) +\
                                s_min_warm*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) +\
                                s_min_cool*(len(x[x < turning_pt[0]])+2))/(len(x)+8)
                s_min_total.append(s_min_triple)

                if s_min_triple < s_min_double:
                    Q_hot  = float(k*t_rot_hot/h/B)
                    Q_warm = float(k*t_rot_warm/h/B)
                    Q_cool = float(k*t_rot_cool/h/B)
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

                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                    [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

                    ax_rot_dia.legend([fit_hot,fit_warm,fit_cool],\
                        [r'$\rm{T_{rot,warm}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
                         r'$\rm{T_{rot,cool}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\rm{T_{rot,cold}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
                    fig_rot_dia.savefig(home+plotdir+objname+'_'+species_name+'_rot_three.pdf',format='pdf',dpi=300, bbox_inches='tight')
                    ax_rot_dia.cla()
                    fig_rot_dia.clf()
                    print('T_rot(hot): %5.1f K, T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_hot,t_rot_warm,t_rot_cool))

                    if polyfit:
                        output = [0.0, float(t_rot_hot), float(t_rot_warm), float(t_rot_cool),
                                  0.0, float(sig_t_rot_hot), float(sig_t_rot_warm), float(sig_t_rot_cool),
                                  0.0,0.0,0.0, N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1],
                                  N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]
                    else:
                        output = [0.0, float(t_rot_hot.A1), float(t_rot_warm.A1), float(t_rot_cool.A1),
                                  0.0, float(sig_t_rot_hot.A1), float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1),
                                  0.0,0.0,0.0, N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1],
                                  N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1], N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1]]

                else:
                    print('No smaller chi2 found for three temperature fitting.')


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
                                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]], polyfit=polyfit)
                                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])], polyfit=polyfit)
                                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])], polyfit=polyfit)
                                [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]], polyfit=polyfit)

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
                    turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1]),np.mean(best_fit[s_min == min(s_min),2])]
                    [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot,sig_yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]], polyfit=polyfit)
                    [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm,sig_yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])], polyfit=polyfit)
                    [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool,sig_yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])], polyfit=polyfit)
                    [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold,sig_yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]], polyfit=polyfit)
                    s_min_quad = (s_min_hot*(len(x[x >= turning_pt[2]])+2) + \
                                  s_min_warm*(len(x[(x < turning_pt[2]) & (x >= turning_pt[1])])+2) + \
                                  s_min_cool*(len(x[(x < turning_pt[1]) & (x >= turning_pt[0])])+2) + \
                                  s_min_cold*(len(x[x < turning_pt[0]])+2)) / (len(x)+11)
                    s_min_total.append(s_min_quad)

                    if (t_rot_hot < 0) or (t_rot_warm < 0) or (t_rot_cool < 0) or (t_rot_cold < 0):
                        print('negative temperature found')

                        # print out the fitted temperatures
                        if os.path.exists(home+plotdir+species_info['name']+'_rotational_temp.txt'):
                            foo = open(home+plotdir+species_info['name']+'_rotational_temp.txt', 'a')
                        else:
                            foo = open(home+plotdir+species_info['name']+'_rotational_temp.txt', 'w')
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
                            if os.path.exists(home+plotdir+species_info['name']+'_breakpoints.txt'):
                                foo = open(home+plotdir+species_info['name']+'_breakpoints.txt', 'a')
                            else:
                                foo = open(home+plotdir+species_info['name']+'_breakpoints.txt', 'w')

                            if type(turning_pt) != list:
                                foo.write(str(turning_pt)+', ')
                            else:
                                for pt in list(turning_pt):
                                    foo.write(str(pt)+', ')
                            foo.write('\n')
                            foo.close()

                        return None


                    if s_min_quad < s_min_triple:
                        Q_hot  = float(k*t_rot_hot/h/B)
                        Q_warm = float(k*t_rot_warm/h/B)
                        Q_cool = float(k*t_rot_cool/h/B)
                        Q_cold = float(k*t_rot_cold/h/B)
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
                        # ax_rot_dia.set_xlim([0,x.max()+200])
                        # ax_rot_dia.set_ylim([42,50])

                        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                        [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

                        ax_rot_dia.legend([fit_hot,fit_warm,fit_cool,fit_cold],\
                            [r'$\rm{T_{rot,hot}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
                             r'$\rm{T_{rot,warm}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                             r'$\rm{T_{rot,cool}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit))),\
                             r'$\rm{T_{rot,cold}\,\,=\,\, %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cold,sig_t_rot_cold,N_cold_fit/10**np.floor(np.log10(N_cold_fit)),np.floor(np.log10(N_cold_fit)))],\
                             numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
                        fig_rot_dia.savefig(home+plotdir+objname+'_'+species_name+'_rot_four.pdf',format='pdf',dpi=300, bbox_inches='tight')
                        ax_rot_dia.cla()
                        fig_rot_dia.clf()
                        print('T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K, T_rot(cold): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool,t_rot_cold))

                        if polyfit:
                            output = (float(t_rot_hot), float(t_rot_warm), float(t_rot_cool), float(t_rot_cold),
                                      float(sig_t_rot_hot), float(sig_t_rot_warm), float(sig_t_rot_cool), float(sig_t_rot_cold),
                                      N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1], N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1],
                                      N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1], N_cold_fit, N_cold_fit_sig[0], N_cold_fit_sig[1])
                        else:
                            output = (float(t_rot_hot.A1), float(t_rot_warm.A1), float(t_rot_cool.A1), float(t_rot_cold.A1),
                                      float(sig_t_rot_hot.A1), float(sig_t_rot_warm.A1), float(sig_t_rot_cool.A1), float(sig_t_rot_cold.A1),
                                      N_hot_fit, N_hot_fit_sig[0], N_hot_fit_sig[1], N_warm_fit, N_warm_fit_sig[0], N_warm_fit_sig[1],
                                      N_cool_fit, N_cool_fit_sig[0], N_cool_fit_sig[1], N_cold_fit, N_cold_fit_sig[0], N_cold_fit_sig[1])

                    else:
                        print('No smaller chi2 found for four temperature fitting.')

        print(s_min_total)

        # print out the fitted temperatures
        if os.path.exists(home+plotdir+species_info['name']+'_rotational_temp.txt'):
            foo = open(home+plotdir+species_info['name']+'_rotational_temp.txt', 'a')
        else:
            foo = open(home+plotdir+species_info['name']+'_rotational_temp.txt', 'w')
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
            if os.path.exists(home+plotdir+species_info['name']+'_breakpoints.txt'):
                foo = open(home+plotdir+species_info['name']+'_breakpoints.txt', 'a')
            else:
                foo = open(home+plotdir+species_info['name']+'_breakpoints.txt', 'w')

            if type(turning_pt) != list:
                foo.write(str(turning_pt)+', ')
            else:
                for pt in list(turning_pt):
                    foo.write(str(pt)+', ')
            foo.write('\n')
            foo.close()

    return x, y, y_sig

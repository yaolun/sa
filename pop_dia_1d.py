def pop_dia_1d(objname,plotdir,dstar,line_data, species_info,
               opt_correction=None, fitting=True):
    """
    line_data: selected line fitting results with typical header as the fitting table
    species_info: A dictionary contains the neccessary info for constructin rotational diagram
                  This script currently only deals with linear molecules such as CO.
                  Thus, only the rotational constant, B, is accepted.  However, more information
                  provided won't have any negative effect.
    opt_correction is the largest J-level that need to be corrected from optical depth
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
    def rot_lin_leastsqfit(x, y, y_sig):
        [yfit ,yerr ,a_hat, cov_hat, s_min] = lin_leastsqfit(x, y, y_sig)

        t_rot = -1/a_hat[1]*np.log10(np.e)
        sig_t_rot = -t_rot*cov_hat[1,1]**0.5/a_hat[1]*np.log10(np.e)
        yoff = a_hat[0]

        return yfit ,yerr ,t_rot ,sig_t_rot ,s_min ,yoff

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
        [yfit,yerr,t_rot,sig_t_rot,s_min,yoff] = rot_lin_leastsqfit(x, y, y_sig)
        Q = float(k*t_rot/h/B)
        N_fit = Q*10**(float(yoff))
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
        print 'T_rot: %8.6f K' % t_rot
        s_min_single = s_min

    # Two temperature fitting
    #
        if len(x) >= 8:
            best_fit = []
            s_min = []
            for i in range(3, len(x)-4):
                turning_pt = x[i]+1

                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
                best_fit.append(turning_pt)
                s_min.append(s_min_warm+s_min_cool)
            best_fit = np.array(best_fit)
            s_min = np.array(s_min)

            turning_pt = np.mean(best_fit[s_min == min(s_min)])

            [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
            [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
            s_min_double = s_min_cool+s_min_warm
            s_min_total.append(s_min_double)
            if (s_min_cool+s_min_warm)<s_min_single:
                Q_warm = float(k*t_rot_warm/h/B)
                Q_cool = float(k*t_rot_cool/h/B)
                N_warm_fit = Q_warm*10**(float(yoff_warm))
                N_cool_fit = Q_cool*10**(float(yoff_cool))
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
                # ax_rot_dia.set_xlim([0,x.max()+200])
                # ax_rot_dia.set_ylim([42,50])

                ax_rot_dia.legend([fit_warm,fit_cool],\
                    [r'$\rm{T_{rot,warm}= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\rm{T_{rot,cool}\,\,= %5.1f \pm %5.1f\,K,\,\mathcal{N}= %3.2f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=10,framealpha=0.3)
                fig_rot_dia.savefig(home+plotdir+objname+'_'+species_name+'_rot_two.pdf',format='pdf',dpi=300, bbox_inches='tight')
                ax_rot_dia.cla()
                fig_rot_dia.clf()
                print 'T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_warm,t_rot_cool)

            # Three temperature fitting
            #
            if x.min() < 500 and len(x) > 12:
                best_fit = []
                s_min = []
                for i in range(3, len(x)-4):
                    for j in range(3, len(x)-4):
                        if (j-i) <4:
                            continue
                        turning_pt = [x[i]+1,x[j]+1]

                        [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                        [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                        [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                        best_fit.append(turning_pt)
                        s_min.append(s_min_hot+s_min_warm+s_min_cool)

                best_fit = np.array(best_fit)
                s_min = np.array(s_min)
                turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1])]
                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                s_min_triple = s_min_cool+s_min_warm+s_min_hot
                s_min_total.append(s_min_triple)

                if (s_min_cool+s_min_warm+s_min_hot) < s_min_double:
                    Q_hot  = float(k*t_rot_hot/h/B)
                    Q_warm = float(k*t_rot_warm/h/B)
                    Q_cool = float(k*t_rot_cool/h/B)
                    N_hot_fit  = Q_hot *10**(float(yoff_hot))
                    N_warm_fit = Q_warm*10**(float(yoff_warm))
                    N_cool_fit = Q_cool*10**(float(yoff_cool))
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
                    # ax_rot_dia.set_xlim([0,x.max()+200])
                    # ax_rot_dia.set_ylim([42,50])

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
                    print 'T_rot(hot): %5.1f K, T_rot(warm): %5.1f K, T_rot(cool): %5.1f K' % (t_rot_hot,t_rot_warm,t_rot_cool)

                # four temperature fitting
                if (x.min() < 500) and (len(x) > 16):
                    best_fit = []
                    s_min = []
                    for i in range(3, len(x)-4):
                        for j in range(3, len(x)-4):
                            for jj in range(3, len(x)-4):
                                if (j-i<4) or (jj-j<4) or (jj-i<8):
                                    continue
                                turning_pt = [x[i]+1, x[j]+1, x[jj]+1]
                                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
                                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
                                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                                [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                                best_fit.append(turning_pt)
                                s_min.append(s_min_hot+s_min_warm+s_min_cool+s_min_cold)

                    best_fit = np.array(best_fit)
                    s_min = np.array(s_min)
                    turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1]),np.mean(best_fit[s_min == min(s_min),2])]
                    [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = rot_lin_leastsqfit(x[x>=turning_pt[2]], y[x>=turning_pt[2]], y_sig[x>=turning_pt[2]])
                    [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = rot_lin_leastsqfit(x[(x<turning_pt[2]) & (x>=turning_pt[1])], y[(x<turning_pt[2]) & (x>=turning_pt[1])], y_sig[(x<turning_pt[2]) & (x>=turning_pt[1])])
                    [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = rot_lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                    [yfit_cold,yerr_cold,t_rot_cold,sig_t_rot_cold,s_min_cold,yoff_cold] = rot_lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                    s_min_total.append(s_min_hot+s_min_warm+s_min_cool+s_min_cold)

                    if (t_rot_hot < 0) or (t_rot_warm < 0) or (t_rot_cool < 0) or (t_rot_cold < 0):
                        print 'negative temperature found'
                        return None

                    if (s_min_cold+s_min_cool+s_min_warm+s_min_hot) < s_min_triple:
                        Q_hot  = float(k*t_rot_hot/h/B)
                        Q_warm = float(k*t_rot_warm/h/B)
                        Q_cool = float(k*t_rot_cool/h/B)
                        Q_cold = float(k*t_rot_cold/h/B)
                        N_hot_fit  = Q_hot *10**(float(yoff_hot))
                        N_warm_fit = Q_warm*10**(float(yoff_warm))
                        N_cool_fit = Q_cool*10**(float(yoff_cool))
                        N_cold_fit = Q_cold*10**(float(yoff_cold))
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
                        print 'T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K, T_rot(cold): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool,t_rot_cold)

    return x, y, y_sig



import numpy as np
import os
home = os.path.expanduser('~')

fitting_table = '/Volumes/SD-Mac/CDF_archive_v2/CDF_archive_v2_lines.txt'

# get the entire object list
from astropy.io import ascii
data = ascii.read(fitting_table)
obj_list = list(set(data['Object'].data))
# get distance info
dist = ascii.read('/Users/yaolun/data/cops-spire_distance.txt')

# species info
# # HCO+
B_hco = 44.594262e9 # Hz
species_info = {'name': 'HCO+', 'B': B_hco}
#
# # 13CO
# B_13co = 55.1010138e9 # Hz
# species_info = {'name': '13CO', 'B': B_13co}

# 12CO
# B_12co = 57.63596828e9 # Hz
# species_info = {'name': 'CO', 'B': B_12co}

# line data
fitting = ascii.read(fitting_table)

# # focus on GSS30-IRS1
# obj_list = ['GSS30-IRS1']
# # additional slicer for choosing the same J-lines as 13CO
# slicer = (fitting['ObsWL(um)'] >= 350) & (fitting['ObsWL(um)'] <= 600)

# loop through all objects
for o in obj_list:
    print o
    # select line data
    data_dum = fitting[(fitting['Object'] == o) & (fitting['Pixel_No.'] == 'c') & \
                       (fitting['Validity'] == 1) & (fitting['SNR'] >= 4)]
    ind = []
    for i in range(len(data_dum)):
        if len(data_dum['Line'][i].split(species_info['name'])[0]) == 0:
            ind.append(i)
    if len(ind) == 0: continue
    line_data = data_dum[ind]

    # rotational diagram with fitting
    pop_dia_1d(o, '/test/', dist['distance'][dist['object'] == o],
               line_data, species_info)

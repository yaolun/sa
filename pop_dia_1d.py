# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def pop_dia_1d(objname,plotdir,pacs=None,spire=None,single=False):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from read_fitting import read_fitting_co
    from leastsqfit import lin_leastsqfit
    home = os.path.expanduser('~')
    
    # Constants Setup
    c = 2.99792458e10
    h = 6.6260755e-27
    k = 1.380658e-16
    B = 1.9225#192.25 m-1
    pc = 3.086e18
    
    if (pacs != None) & (spire == None):
        [co_pacs,co_name_pacs] = read_fitting_co(pacs,3)
        co_data = co_pacs.astype('float')
        co_data_name = co_name_pacs
    if (spire != None) & (pacs == None):
        [co_spire,co_name_spire] = read_fitting_co(spire,3)
        co_data = co_spire.astype('float')
        co_data_name = co_name_spire
    if (pacs != None) & (spire != None):
        [co_pacs,co_name_pacs] = read_fitting_co(pacs,3)
        [co_spire,co_name_spire] = read_fitting_co(spire,3)
        co_data = np.concatenate((co_pacs,co_spire)).astype('float')
        co_data_name = np.concatenate((co_name_pacs,co_name_spire))
    if 'co_data_name' not in locals():
        return None
    if len(co_data_name) <= 2:
        return None
    
    # Calculate the N/g and Eu from the data
    v = c/(co_data[:,0]*1e-4)
    N = 4*np.pi*co_data[:,2]*1e7*(178*pc)**2/(co_data[:,5]*h*v)
    N_sigma = 4*np.pi*co_data[:,3]*1e7*(178*pc)**2/(co_data[:,5]*h*v)
    x = co_data[:,4]
    y = np.log10(N/co_data[:,6])
    yerr_hi = np.log10((N+N_sigma)/co_data[:,6])-np.log10(N/co_data[:,6])
    yerr_low = np.log10(N/co_data[:,6])-np.log10((N-N_sigma)/co_data[:,6])
    y_sig = y*0
    for i in range(0,len(y)):
            y_sig[i] = max(yerr_hi[i], yerr_low[i])
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    y_sig = y_sig[ind]
    
    
    # Single temperature fitting
    #
    if len(x)>2:
        [yfit,yerr,t_rot,sig_t_rot,s_min,yoff] = lin_leastsqfit(x, y, y_sig)
        Q = float(k*t_rot/h/c/B)
        N_fit = Q*10**(float(yoff))
        x = x.reshape(len(x))
        y = y.reshape(len(y))
        y_sig = y_sig.reshape(len(y_sig))
        
        fig_rot_dia = plt.figure()
        ax_rot_dia = fig_rot_dia.add_subplot(111)
        data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
        ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')
        
        fit, = ax_rot_dia.plot(x,yfit,color='DarkMagenta')
        ax_rot_dia.plot(x,yfit+yerr,'--',color='Magenta')
        ax_rot_dia.plot(x,yfit-yerr,'--',color='Magenta')
        ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=18)
        ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
        ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
        [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax_rot_dia.set_xlim([50,4500])
        ax_rot_dia.legend([fit],[r'$\mathrm{T_{rot}= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot,sig_t_rot,N_fit/10**np.floor(np.log10(N_fit)),np.floor(np.log10(N_fit)))],numpoints=1,loc='upper right',fontsize=14)
        fig_rot_dia.savefig(home+plotdir+objname+'_co_rot_single.pdf',format='pdf',dpi=300)
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        print 'T_rot: %8.6f K' % t_rot
        s_min_single = s_min
        
    # Two temperature fitting
    #
        if len(x)>=6:
            best_fit = []
            s_min = []
            for i in range(int(x[2]+1),int(x[-3]-1),10):
                turning_pt = i
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
                best_fit.append(i)
                s_min.append(s_min_warm+s_min_cool)
            best_fit = np.array(best_fit)
            s_min = np.array(s_min)
            # fig_s_min = plt.figure()
            # ax_s_min = fig_s_min.add_subplot(111)
            # ax_s_min.plot(best_fit,s_min)
            # ax_s_min.set_xlabel(r'$\mathrm{Turning~points}$')
            # ax_s_min.set_ylabel(r'$\mathrm{S_{min}}$')
            # fig_s_min.savefig(home+'/bhr71/plots/rotational_diagram/s_min.eps',format='eps',dpi=300)
            # ax_s_min.cla()
            # fig_s_min.clf()
    
            turning_pt = np.mean(best_fit[s_min == min(s_min)])
            [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[x>=turning_pt], y[x>=turning_pt], y_sig[x>=turning_pt])
            [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt], y[x<turning_pt], y_sig[x<turning_pt])
            if (s_min_cool+s_min_warm)<s_min_single:
                Q_warm = float(k*t_rot_warm/h/c/B)
                Q_cool = float(k*t_rot_cool/h/c/B)
                N_warm_fit = Q_warm*10**(float(yoff_warm))
                N_cool_fit = Q_cool*10**(float(yoff_cool))
                s_min_double = s_min_cool+s_min_warm
                fig_rot_dia = plt.figure()
                ax_rot_dia = fig_rot_dia.add_subplot(111)
                data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
                ax_rot_dia.errorbar(x,y,yerr=(yerr_low,yerr_hi),linestyle='None',color='DarkGreen')
    
                fit_warm, = ax_rot_dia.plot(x[x>=turning_pt],yfit_warm,color='DarkMagenta')
                ax_rot_dia.plot(x[x>=turning_pt],yfit_warm+yerr_warm,'--',color='Magenta')
                ax_rot_dia.plot(x[x>=turning_pt],yfit_warm-yerr_warm,'--',color='Magenta')
    
                fit_cool, = ax_rot_dia.plot(x[x<turning_pt],yfit_cool,color='Blue')
                ax_rot_dia.plot(x[x<turning_pt],yfit_cool+yerr_cool,'--',color='MediumBlue')
                ax_rot_dia.plot(x[x<turning_pt],yfit_cool-yerr_cool,'--',color='MediumBlue')
    
                ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=18)
                ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
                ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
                ax_rot_dia.set_xlim([50,4500])
    
                ax_rot_dia.legend([fit_warm,fit_cool],\
                    [r'$\mathrm{T_{rot,warm}= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\mathrm{T_{rot,cool}~~= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],\
                         numpoints=1,loc='upper right',fontsize=14)
                fig_rot_dia.savefig(home+plotdir+objname+'_co_rot_two.pdf',format='pdf',dpi=300)
                ax_rot_dia.cla()
                fig_rot_dia.clf()
                print 'T_rot(warm): %8.6f K, T_rot(cool): %8.6f K' % (t_rot_warm,t_rot_cool)
            # Three temperature fitting
            #
            if spire != None and len(x) > 9:
                best_fit = []
                s_min = []
                for i in range(int(x[2]+1),int(x[-6]-1),10):
                    for j in range(int(x[5]+1),int(x[-3]-1),10):
                        if len(x[(x>=i) & (x<=j)])< 3:
                            break
                        turning_pt = [i,j]
                        [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                        [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                        [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                        best_fit.append(turning_pt)
                        s_min.append(s_min_hot+s_min_warm+s_min_cool)
                        
                best_fit = np.array(best_fit)
                s_min = np.array(s_min)
                turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1])]
                [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
                [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
                [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                
                if (s_min_cool+s_min_warm+s_min_hot) < s_min_double:
                    Q_hot  = float(k*t_rot_hot/h/c/B)
                    Q_warm = float(k*t_rot_warm/h/c/B)
                    Q_cool = float(k*t_rot_cool/h/c/B)
                    N_hot_fit  = Q_hot *10**(float(yoff_hot))
                    N_warm_fit = Q_warm*10**(float(yoff_warm))
                    N_cool_fit = Q_cool*10**(float(yoff_cool))
                    fig_rot_dia = plt.figure()
                    ax_rot_dia = fig_rot_dia.add_subplot(111)
                    data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
                    ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')
                    
                    fit_hot,  = ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot,color='Red')
                    ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot+yerr_hot,'--',color='Tomato')
                    ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot-yerr_hot,'--',color='Tomato')
        
                    fit_warm, = ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm,color='DarkMagenta')
                    ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm+yerr_warm,'--',color='Magenta')
                    ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm-yerr_warm,'--',color='Magenta')
        
                    fit_cool, = ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool,color='Blue')
                    ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool+yerr_cool,'--',color='MediumBlue')
                    ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool-yerr_cool,'--',color='MediumBlue')
        
                    ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=18)
                    ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
                    ax_rot_dia.set_xlim([50,4500])
                    
                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
                    ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
                    [ax_rot_dia.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
                    # ax_rot_dia.minorticks_on()

                    ax_rot_dia.legend([fit_hot,fit_warm,fit_cool],\
                        [r'$\mathrm{T_{rot,warm}= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
                         r'$\mathrm{T_{rot,cool}~~= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
                         r'$\mathrm{T_{rot,cold}~~=~~ %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],numpoints=1,loc='upper right',fontsize=14)
                    fig_rot_dia.savefig(home+plotdir+objname+'_co_rot_three.pdf',format='pdf',dpi=300)
                    ax_rot_dia.cla()
                    fig_rot_dia.clf()
                    print 'T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool)
            # four temperature fitting
            # if (spire != None) and (len(x) > 12):
            #     best_fit = []
            #     s_min = []
            #     for i in range(int(x[2]+1),int(x[-6]-1),10):
            #         for j in range(int(x[5]+1),int(x[-3]-1),10):
            #             if len(x[(x>=i) & (x<=j)])< 3:
            #                 break
            #             turning_pt = [i,j]
            #             [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
            #             [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
            #             [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
            #             best_fit.append(turning_pt)
            #             s_min.append(s_min_hot+s_min_warm+s_min_cool)
                        
            #     best_fit = np.array(best_fit)
            #     s_min = np.array(s_min)
            #     turning_pt = [np.mean(best_fit[s_min == min(s_min),0]),np.mean(best_fit[s_min == min(s_min),1])]
            #     [yfit_hot,yerr_hot,t_rot_hot,sig_t_rot_hot,s_min_hot,yoff_hot] = lin_leastsqfit(x[x>=turning_pt[1]], y[x>=turning_pt[1]], y_sig[x>=turning_pt[1]])
            #     [yfit_warm,yerr_warm,t_rot_warm,sig_t_rot_warm,s_min_warm,yoff_warm] = lin_leastsqfit(x[(x<turning_pt[1]) & (x>=turning_pt[0])], y[(x<turning_pt[1]) & (x>=turning_pt[0])], y_sig[(x<turning_pt[1]) & (x>=turning_pt[0])])
            #     [yfit_cool,yerr_cool,t_rot_cool,sig_t_rot_cool,s_min_cool,yoff_cool] = lin_leastsqfit(x[x<turning_pt[0]], y[x<turning_pt[0]],y_sig[x<turning_pt[0]])
                
            #     if (s_min_cool+s_min_warm+s_min_hot) < s_min_double:
            #         Q_hot  = float(k*t_rot_hot/h/c/B)
            #         Q_warm = float(k*t_rot_warm/h/c/B)
            #         Q_cool = float(k*t_rot_cool/h/c/B)
            #         N_hot_fit  = Q_hot *10**(float(yoff_hot))
            #         N_warm_fit = Q_warm*10**(float(yoff_warm))
            #         N_cool_fit = Q_cool*10**(float(yoff_cool))
            #         fig_rot_dia = plt.figure()
            #         ax_rot_dia = fig_rot_dia.add_subplot(111)
            #         data, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
            #         ax_rot_dia.errorbar(x,y,yerr=y_sig,linestyle='None',color='DarkGreen')
                    
            #         fit_hot,  = ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot,color='Red')
            #         ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot+yerr_hot,'--',color='Tomato')
            #         ax_rot_dia.plot(x[x>turning_pt[1]],yfit_hot-yerr_hot,'--',color='Tomato')
        
            #         fit_warm, = ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm,color='DarkMagenta')
            #         ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm+yerr_warm,'--',color='Magenta')
            #         ax_rot_dia.plot(x[(x>=turning_pt[0]) & (x<turning_pt[1])],yfit_warm-yerr_warm,'--',color='Magenta')
        
            #         fit_cool, = ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool,color='Blue')
            #         ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool+yerr_cool,'--',color='MediumBlue')
            #         ax_rot_dia.plot(x[x<turning_pt[0]],yfit_cool-yerr_cool,'--',color='MediumBlue')
        
            #         ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=18)
            #         ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=18)
            #         ax_rot_dia.set_xlim([50,4500])
                    
            #         ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='major')
            #         ax_rot_dia.tick_params('both',labelsize=16,width=1.5,which='minor')
            #         # ax_rot_dia.minorticks_on()

            #         ax_rot_dia.legend([fit_hot,fit_warm,fit_cool],\
            #             [r'$\mathrm{T_{rot,warm}= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_hot,sig_t_rot_hot,N_hot_fit/10**np.floor(np.log10(N_hot_fit)),np.floor(np.log10(N_hot_fit))),\
            #              r'$\mathrm{T_{rot,cool}~~= %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_warm,sig_t_rot_warm,N_warm_fit/10**np.floor(np.log10(N_warm_fit)),np.floor(np.log10(N_warm_fit))),\
            #              r'$\mathrm{T_{rot,cold}~~=~~ %8.4f \pm %8.4f~K,~\mathcal{N}= %4.3f \times 10^{%d}}$' % (t_rot_cool,sig_t_rot_cool,N_cool_fit/10**np.floor(np.log10(N_cool_fit)),np.floor(np.log10(N_cool_fit)))],numpoints=1,loc='upper right',fontsize=14)
            #         fig_rot_dia.savefig(home+plotdir+objname+'_co_rot_three.eps',format='eps',dpi=300)
            #         ax_rot_dia.cla()
            #         fig_rot_dia.clf()
            #         print 'T_rot(hot): %8.6f K, T_rot(warm): %8.6f K, T_rot(cool): %8.6f K' % (t_rot_hot,t_rot_warm,t_rot_cool)


def pop_dia_h2o_1d(objname,plotdir,pacs=None,spire=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from read_fitting import read_fitting_h2o
    from leastsqfit import lin_leastsqfit
    home = os.path.expanduser('~')
    
    # Constants Setup
    c = 2.99792458e10
    h = 6.6260755e-27
    k = 1.380658e-16
    pc = 3.086e18

    # Partition function for H2O adopt from Green 2013
    # Q = -4.9589 + 0.28147*T_rot + 0.0012848*T_rot^2 - 5.8343e-7*T_rot^3
    def Q_h2o(T_rot):
        Q = -4.9589 + 0.28147*T_rot + 0.0012848*T_rot**2 - 5.8343e-7*T_rot**3
        return Q
    
    if (pacs != None) & (spire == False):
        [ph2o_pacs,ph2o_name_pacs,oh2o_pacs,oh2o_name_pacs] = read_fitting_h2o(pacs,3,positive=True)
        ph2o_data = ph2o_pacs.astype('float')
        ph2o_data_name = ph2o_name_pacs
        oh2o_data = oh2o_pacs.astype('float')
        oh2o_data_name = oh2o_name_pacs
    if (spire != None) & (pacs == False):
        [ph2o_spire,ph2o_name_spire,oh2o_spire,oh2o_name_spire] = read_fitting_h2o(spire,3,positive=True)
        ph2o_data = ph2o_spire.astype('float')
        ph2o_data_name = ph2o_name_spire
        oh2o_data = oh2o_spire.astype('float')
        oh2o_data_name = oh2o_name_spire
    if (pacs != None) & (spire != None):
        [ph2o_pacs,ph2o_name_pacs,oh2o_pacs,oh2o_name_pacs] = read_fitting_h2o(pacs,3,positive=True)
        [ph2o_spire,ph2o_name_spire,oh2o_spire,oh2o_name_spire] = read_fitting_h2o(spire,3,positive=True)
        ph2o_data = np.concatenate((ph2o_pacs,ph2o_spire)).astype('float')
        ph2o_data_name = np.concatenate((ph2o_name_pacs,ph2o_name_spire))
        oh2o_data = np.concatenate((oh2o_pacs,oh2o_spire)).astype('float')
        oh2o_data_name = np.concatenate((oh2o_name_pacs,oh2o_name_spire))

    # Iterate through ph2o and oh2o
    water_data = [ph2o_data,oh2o_data]
    water_data_name = [ph2o_data_name,oh2o_data_name]


    fig_rot_dia = plt.figure()
    ax_rot_dia = fig_rot_dia.add_subplot(111)
    for j in range(0,2):
        water = water_data[j]
        water_name = water_data_name[j]
        
        # Calculate the N/g and Eu from the data
        v = c/(water[:,0]*1e-4)
        N = 4*np.pi*water[:,2]*1e7*(178*pc)**2/(water[:,5]*h*v)
        N_sigma = 4*np.pi*water[:,3]*1e7*(178*pc)**2/(water[:,5]*h*v)
        x = water[:,4]
        y = np.log10(N/water[:,6])
        yerr_hi = np.log10((N+N_sigma)/water[:,6])-np.log10(N/water[:,6])
        yerr_low = abs(np.log10(N/water[:,6])-np.log10((abs(N-N_sigma))/water[:,6]))
        y_sig = y*0
        for i in range(0,len(y)):
                y_sig[i] = max(yerr_hi[i], yerr_low[i])
        ind = np.argsort(x)
        x = x[ind]
        y = y[ind]
        y_sig = y_sig[ind]

        # Plot the data
        if j == 0:
            ph2o, = ax_rot_dia.plot(x,y,'o',color='DarkGreen',markersize=6)
            ax_rot_dia.errorbar(x,y,yerr=(yerr_low,yerr_hi),linestyle='None',color='DarkGreen')
        else:
            oh2o, = ax_rot_dia.plot(x,y,'o',color='Crimson',markersize=6)
            ax_rot_dia.errorbar(x,y,yerr=(yerr_low,yerr_hi),linestyle='None',color='Crimson')

    # If the spectrum doesn't have enough water lines detection for the fitting, then just plot the result.
    water_data = np.concatenate((ph2o_data,oh2o_data))
    water_name = np.concatenate((ph2o_data_name,oh2o_data_name))

    # Calculate the N/g and Eu from the data
    v = c/(water_data[:,0]*1e-4)
    N = 4*np.pi*water_data[:,2]*1e7*(178*pc)**2/(water_data[:,5]*h*v)
    N_sigma = 4*np.pi*water_data[:,3]*1e7*(178*pc)**2/(water_data[:,5]*h*v) 

    x = water_data[:,4]
    y = np.log10(N/water_data[:,6])
    yerr_hi = np.log10((N+N_sigma)/water_data[:,6])-np.log10(N/water_data[:,6])
    yerr_low = abs(np.log10(N/water_data[:,6])-np.log10((abs(N-N_sigma))/water_data[:,6]))
    y_sig = y*0
    for i in range(0,len(y)):
            y_sig[i] = max(yerr_hi[i], yerr_low[i])
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    y_sig = y_sig[ind]
    for t in range(0,len(x)):
        print '%s %8.6f %6.2f %8.6e %8.6e' % (water_name[t],water_data[t,0],x[t],y[t],y_sig[t])

    if len(water_name) <= 2:
        fig_rot_dia.savefig(home+plotdir+objname+'_h2o_rot_dia.eps',format='eps',dpi=300)
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        return None
    # Single temperature fitting
    #
    if len(x)>2:
        # Fitting
        #
        [yfit,yerr,t_rot,sig_t_rot,s_min,yoff] = lin_leastsqfit(x, y, y_sig)
        Q = Q_h2o(t_rot)
        N_fit = Q*10**(float(yoff))
        x = x.reshape(len(x))
        y = y.reshape(len(y))
        y_sig = y_sig.reshape(len(y_sig))
        
        fit, = ax_rot_dia.plot(x,yfit,color='DarkMagenta')
        ax_rot_dia.plot(x,yfit+yerr,'--',color='Magenta')
        ax_rot_dia.plot(x,yfit-yerr,'--',color='Magenta')
        ax_rot_dia.set_xlabel(r'$\mathrm{E_{u}~(K)}$',fontsize=14)
        ax_rot_dia.set_ylabel(r'$\mathrm{log(\mathcal{N}_{J}/g_{J})}$',fontsize=14)
        ax_rot_dia.set_xlim([0,1600])
        lg_fit = ax_rot_dia.legend([fit],[r'$\mathrm{T_{rot}= %8.4f \pm %8.6f~K,~\mathcal{N}= %.3e~g/cm^{2}}$' % (t_rot,sig_t_rot,N_fit)],numpoints=1,loc='lower left',fontsize=12)
        lg_data = ax_rot_dia.legend([oh2o,ph2o],[r'$\mathrm{o-H_{2}O}$',r'$\mathrm{p-H_{2}O}$'],loc='lower right',fontsize=14,numpoints=1)
        plt.gca().add_artist(lg_fit)
        fig_rot_dia.savefig(home+plotdir+objname+'_h2o_rot_single.eps',format='eps',dpi=300)
        ax_rot_dia.cla()
        fig_rot_dia.clf()
        print 'T_rot: %8.6f K' % t_rot
        s_min_single = s_min
# <codecell>

import numpy as np
import matplotlib.pyplot as plt
import os
from read_fitting import read_fitting_co
from leastsqfit import lin_leastsqfit
home = os.path.expanduser('~')

pacs = '/bhr71/fitting/unc_test/pacs/advanced_products/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
spire = '/bhr71/fitting/1107/spire/advanced_products/BHR71_spire_corrected_lines.txt'
pop_dia_1d('BHR71','/bhr71/plots/',pacs=pacs,spire=spire)
pacs_cube = '/bhr71/fitting/unc_test/pacs/advanced_products/cube/BHR71_pacs_pixel'
# pacs_cube = [pacs_cube+str(i)+'_os8_sf7_lines.txt' for i in range(1,26)]
for i in range(1,26):
    print i
    pop_dia_1d('BHR71_pacs_pixel'+str(i),'/bhr71/plots/',pacs=pacs_cube+str(i)+'_os8_sf7_lines.txt')

# pop_dia_h2o_1d('BHR71','/bhr71/plots/',pacs=pacs,spire=spire)
# <codecell>



pacs = '/bhr71/fitting/unc_test/pacs/advanced_products/BHR71_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt'
[co_pacs,co_name_pacs] = read_fitting_co(pacs,3)
spire = '/bhr71/fitting/1107/spire/advanced_products/BHR71_spire_corrected_lines.txt'
[co_spire,co_name_spire] = read_fitting_co(spire,3)

# # <codecell>

# yfit = np.array(yfit)

# # <codecell>

# plt.plot(x,np.squeeze(yfit))

# # <codecell>



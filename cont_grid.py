def cont_grid(indir, objlist, outdir):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import astropy.constants as const
    from matplotlib.ticker import MaxNLocator
    home = os.path.expanduser('~')

    # constant setip
    c = const.c.cgs.value

    row = int(np.ceil(len(objlist)/2.))
    col = 2

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row')

    for i in range(0, row):
        for j in range(0, col):
            num = i + j*row
            ax = axarr[i,j]
            if (num+1) <= len(objlist):
                [wl_cont_pacs,flux_cont_pacs,unc_cont_pacs] = np.genfromtxt(indir+objlist[num]+'/pacs/advanced_products/'+objlist[num]+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',skip_header=1).T
                [wl_cont_spire,flux_cont_spire] = np.genfromtxt(indir+objlist[num]+'/spire/advanced_products/'+objlist[num]+'_spire_corrected_continuum.txt',skip_header=1).T

                cont, = ax.plot(np.log10(wl_cont_pacs[wl_cont_pacs<100]),np.log10(flux_cont_pacs[wl_cont_pacs<100]*(c/wl_cont_pacs[wl_cont_pacs<100]*1e4)),\
                                '-',color='Blue',linewidth=1)
                ax.plot(np.log10(wl_cont_pacs[wl_cont_pacs>100]),np.log10(flux_cont_pacs[wl_cont_pacs>100]*(c/wl_cont_pacs[wl_cont_pacs>100]*1e4)),\
                        '-',color='Blue',linewidth=1)
                ax.plot(np.log10(wl_cont_spire),np.log10(flux_cont_spire*(c/wl_cont_spire*1e4)),'r-',linewidth=1)

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on()
            ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            if (i != 0) & (j != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='upper'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))
            # ax.set_xlim([50,670])
            # ax.set_xlabel(r'$\mathrm{log\lambda~({\mu}m)}$',fontsize=20)
            # ax.set_ylabel(r'$\mathrm{log\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$',fontsize=20)

    fig.text(0.5, 0, r'$\mathrm{log~\lambda~({\mu}m)}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=20, va='center', rotation='vertical')
    fig.subplots_adjust(hspace=0,wspace=0)

    fig.savefig(outdir+'cont_grid.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)


# indir = '/home/bettyjo/yaolun/FWD_archive/CDF_archive_feb15/'    
# objlist = ['B1-a','B1-c','B335','BHR71','IRAS12496','FUOri','GSS30-IRS1','IRAS03245','IRAS03301','L1455-IRS3','L1157','L1551-IRS5',\
           # 'RCrA-IRS7B','TMC1A','TMC1','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']
           # jitter correction failed
           # IRS46, L1014, RCrA-IRS5A, RCrA-IRS7C
# sorted by Tbol

# testing
# indir = '/Users/yaolun/test/'
# objlist = ['BHR71','BHR71','BHR71','BHR71','BHR71','BHR71','BHR71']
# cont_grid(indir, objlist)
# Tbol = []
    
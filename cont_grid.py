def cont_grid(indir, objlist, outdir, Tbol=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import astropy.constants as const
    from matplotlib.ticker import MaxNLocator
    home = os.path.expanduser('~')

    objlist = np.array(objlist)
    Tbol = np.array(Tbol)

    # constant setip
    c = const.c.cgs.value

    # sorted by Tbol if provided
    if Tbol != None:
        objlist = objlist[np.argsort(Tbol)]

    row = int(np.ceil(len(objlist)/2.))
    col = 2

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,18))

    for i in range(0, row):
        for j in range(0, col):
            num = j + i*col
            ax = axarr[i,j]
            if (num+1) <= len(objlist):
                [wl_cont_pacs,flux_cont_pacs,unc_cont_pacs] = np.genfromtxt(indir+objlist[num]+'/pacs/advanced_products/'+objlist[num]+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',skip_header=1).T
                [wl_cont_spire,flux_cont_spire] = np.genfromtxt(indir+objlist[num]+'/spire/advanced_products/'+objlist[num]+'_spire_corrected_continuum.txt',skip_header=1).T

                cont, = ax.plot(np.log10(wl_cont_pacs[wl_cont_pacs<100]),np.log10(flux_cont_pacs[wl_cont_pacs<100]*(c/wl_cont_pacs[wl_cont_pacs<100]*1e4)),\
                                '-',color='Blue',linewidth=1)
                ax.plot(np.log10(wl_cont_pacs[wl_cont_pacs>100]),np.log10(flux_cont_pacs[wl_cont_pacs>100]*(c/wl_cont_pacs[wl_cont_pacs>100]*1e4)),\
                        '-',color='Blue',linewidth=1)
                ax.plot(np.log10(wl_cont_spire),np.log10(flux_cont_spire*(c/wl_cont_spire*1e4)),'r-',linewidth=1)

            if i == 1:
                ax.yaxis.tick_right()

            [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
            ax.minorticks_on()
            ax.tick_params('both',labelsize=12,width=1.5,which='major',pad=15,length=5)
            ax.tick_params('both',labelsize=12,width=1.5,which='minor',pad=15,length=2.5)

            # fix the overlap tick labels
            x_nbins = len(ax.get_xticklabels())
            y_nbins = len(ax.get_yticklabels())
            ax.text(0.7,0.8,r'$\mathrm{'+objlist[num]+'}$', fontsize=14, transform=ax.transAxes)
            if (i != 0) & (j != 0):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune='lower'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))
            # ax.set_xlim([50,670])
            # ax.set_xlabel(r'$\mathrm{log\lambda~({\mu}m)}$',fontsize=20)
            # ax.set_ylabel(r'$\mathrm{log\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$',fontsize=20)

    fig.text(0.5, 0.07, r'$\mathrm{log~\lambda~({\mu}m)}$', fontsize=20, ha='center')
    fig.text(0.05, 0.5, r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=20, va='center', rotation='vertical')
    fig.subplots_adjust(hspace=0,wspace=0)

    fig.savefig(outdir+'cont_grid.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)


# indir = '/home/bettyjo/yaolun/FWD_archive/CDF_archive_feb15/' 
indir = '/Users/yaolun/data/CDF_archive/'   
objlist = ['B1-a','B1-c','B335','BHR71','IRAS12496','FUOri','GSS30-IRS1','IRAS03245','IRAS03301','L1455-IRS3','L1157','L1551-IRS5',\
           'RCrA-IRS7B','TMC1A','TMC1','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','VLA1623','WL12'] # ,'V1735Cyg'
           # jitter correction failed
           # IRS46, L1014, RCrA-IRS5A, RCrA-IRS7C
# sorted by Tbol
Tbol = [86,51,41,47,579,2440,146,45,350,171,37,106,53,167,159,136,1070,1130,1390,30,235]

objlist = ['B1-a','B1-c','B335','BHR71','IRAS12496','FUOri','GSS30-IRS1','IRAS03245','IRAS03301','L1455-IRS3','L1157',\
           'RCrA-IRS7B','TMC1A','TMC1','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','VLA1623','WL12'] # ,'V1735Cyg'
           # jitter correction failed
           # IRS46, L1014, RCrA-IRS5A, RCrA-IRS7C
# sorted by Tbol
Tbol = [86,51,41,47,579,2440,146,45,350,171,37,53,167,159,136,1070,1130,1390,30,235]
cont_grid(indir, objlist, '/Users/yaolun/test/', Tbol=Tbol)

# testing
# indir = '/Users/yaolun/test/'
# objlist = ['BHR71','BHR71','BHR71','BHR71','BHR71','BHR71','BHR71']
# cont_grid(indir, objlist)
# Tbol = []
    
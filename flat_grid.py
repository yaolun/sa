def flat_grid(indir, objlist, outdir, Tbol=None):
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

    row = len(objlist)
    col = 1

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row', figsize=(12,24))

    for i in range(0, row):
        ax = axarr[i]

        [wl_flat_pacs,flux_flat_pacs,unc_flat_pacs] = np.genfromtxt(indir+objlist[i]+'/pacs/advanced_products/'+objlist[i]+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_flat_spectrum.txt',skip_header=1).T
        [wl_flat_spire,flux_flat_spire] = np.genfromtxt(indir+objlist[i]+'/spire/advanced_products/'+objlist[i]+'_spire_corrected_flat_spectrum.txt',skip_header=1).T

        pacs_spec, = ax.plot(c*1e-5/wl_flat_pacs[wl_flat_pacs<100],flux_flat_pacs[wl_flat_pacs<100],'-',color='Green',linewidth=1)
        ax.plot(c*1e-5/wl_flat_pacs[wl_flat_pacs>100],flux_flat_pacs[wl_flat_pacs>100],'-',color='Green',linewidth=1)
        spire_spec, = ax.plot(c*1e-5/wl_flat_spire,flux_flat_spire,'r-',linewidth=1)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5)
        ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

        ax.text(0.7,0.8,r'$\mathrm{'+objlist[i]+'}$', fontsize=14, transform=ax.transAxes)
        # fix the overlap tick labels
        y_nbins = len(ax.get_yticklabels())
        if (i != 0):
            ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    
    fig.text(0.5, 0, r'$\mathrm{Frequency (Hz)}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{S_{\nu}~(erg~s^{-1}~cm^{-2}~Hz^{-1})}$', fontsize=20, va='center', rotation='vertical')
    fig.subplots_adjust(hspace=0,wspace=0)

    fig.savefig(outdir+'flat_grid.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)


indir = '/Users/yaolun/data/CDF_archive/'   
objlist = ['B1-a','B1-c','B335','BHR71','IRAS12496','FUOri','GSS30-IRS1','IRAS03245','IRAS03301','L1455-IRS3','L1157','L1551-IRS5',\
           'RCrA-IRS7B','TMC1A','TMC1','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','VLA1623','WL12'] # ,'V1735Cyg'
           # jitter correction failed
           # IRS46, L1014, RCrA-IRS5A, RCrA-IRS7C
# sorted by Tbol
Tbol = [86,51,41,47,579,2440,146,45,350,171,37,106,53,167,159,136,1070,1130,1390,30,235]
# L1551 :106
flat_grid(indir, objlist, '/Users/yaolun/test/', Tbol=Tbol)
    
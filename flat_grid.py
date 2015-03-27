def flat_grid(indir, objlist):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import astropy.constants as const
    from matplotlib.ticker import MaxNLocator
    home = os.path.expanduser('~')

    # constant setip
    c = const.c.cgs.value

    row = len(objlist)
    col = 1

    fig, axarr = plt.subplots(row, col, sharex='col', sharey='row')

    for i in range(0, row):
        ax = axarr[i,j]

        [wl_cont_pacs,flux_cont_pacs,unc_cont_pacs] = np.genfromtxt(indir+objlist[i]+'/pacs/advanced_products/'+objlist[i]+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_continuum.txt',skip_header=1).T
        [wl_cont_spire,flux_cont_spire] = np.genfromtxt(indir+objlist[i]+'/spire/advanced_products/'+objlist[i]+'_spire_corrected_continuum.txt',skip_header=1).T
        
        pacs_spec, = ax.plot(c*1e-5/wl_flat_pacs[wl_pacs<100],flux_flat_pacs[wl_pacs<100],'-',color='Green',linewidth=1)
        ax.plot(c*1e-5/wl_flat_pacs[wl_pacs>100],flux_flat_pacs[wl_pacs>100],'-',color='Green',linewidth=1)
        spire_spec, = ax.plot(c*1e-5/wl_flat_spire,flux_flat_spire,'r-',linewidth=1)

        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
        ax.minorticks_on() 
        ax.tick_params('both',labelsize=14,width=1.5,which='major',pad=15,length=5);;llk,,n nczxvcz2        /;z
        ax.tick_params('both',labelsize=14,width=1.5,which='minor',pad=15,length=2.5)

        # fix the overlap tick labels
        y_nbins = len(ax.get_yticklabels())
        if (i != 0) & (j != 0):
            ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune='upper'))

    
    fig.text(0.5, 0, r'$\mathrm{log~\lambda~({\mu}m)}$', fontsize=20, ha='center')
    fig.text(0, 0.5, r'$\mathrm{log~\nu S_{\nu}~(erg~s^{-1}~cm^{-2})}$', fontsize=20, va='center', rotation='vertical')
    fig.subplots_adjust(hspace=0,wspace=0)

    fig.savefig('/Users/yaolun/test/flat_grid.pdf',format='pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)
    
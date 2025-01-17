def cdf_pacs_1d(osbid, objname, outdir, fits_for_header, line_fitting_dir):
    """
    cubedir: the directory contains FITS rebinned cube
    outcubedir: the directory contains ASCII style spectra for each spaxel
    out1ddir: the directory contains ASCII style 1-D spectrum
    """

    from pacs_weight import pacs_weight
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii, fits
    import pidly
    idl = pidly.IDL('/Applications/exelis/idl83/bin/idl')

    # outdir = '/Users/yaolun/bhr71/best_calibrated/'
    # cubedir = '/Users/yaolun/bhr71/data/HSA/'
    # photpath = '/Users/yaolun/bhr71/best_calibrated/bhr71.txt'

    cubedir = outdir+'data/fits/'

    cubefile = [cubedir+'OBSID_'+obsid[0]+'_blue_finalcubes_slice00_os8_sf7.fits',
                cubedir+'OBSID_'+obsid[0]+'_red_finalcubes_slice00_os8_sf7.fits',
                cubedit+'OBSID_'+obsid[1]+'_blue_finalcubes_slice00_os8_sf7.fits',
                cubedir+'OBSID_'+obsid[1]+'_red_finalcubes_slice00_os8_sf7.fits',]

    idl('.r '+line_fitting_dir+'get_pacs.pro')
    idl.pro('get_pacs', outdir=outdir+'data/cube/', objname=objname, filename=cubefile, suffix='os8_sf7')

    wl, flux = pacs_weight(outdir+'data/cube/', objname, 31.8, outdir+'data/', fits_for_header, suffix='os8_sf7')

    idl('.r '+line_fitting_dir+'gauss.pro')
    idl('.r '+line_fitting_dir+'extract_pacs.pro')

    idl.pro('extract_pacs', indir=outdir+'data/', filename=objname+'_pacs_weighted',
            outdir=outdir+'advanced_products/',
            plotdir=outdir+'advanced_products/plots/', noiselevel=3, ra=0, dec=0, global_noise=20,
            localbaseline=10, opt_width=1, continuum=1, flat=1, object=objname, double_gauss=1, fixed_width=1)

obsid_list = [['AB_Aur','1342217842','1342217843','0'],\
              ['AS205','1342215737','1342215738','0'],\
              ['B1-a','1342216182','1342216183','1342249475'],\
              ['B1-c','1342216213','1342216214','1342249476'],\
              ['B335','1342208889','1342208888','1342253652'],\
              ['BHR71','1342212230','1342212231','1342248249'],\
              ['Ced110','0','0','1342248246'],\
              ['DG_Tau','1342225730','1342225731','0'],\
              ['EC82','1342192975','1342219435','0'],\
              ['Elias29','1342228519','1342228520','0'],\
              ['FU_Orionis','1342250907','1342250908','0'],\
              ['FUOri','0','0','1342230412'],\
              ['GSS30-IRS1','1342215678','1342215679','0'],\
              ['GSS30','0','0','1342251286'],\
              ['HD100453','1342211695','1342211696','0'],\
              ['HD100546','1342188037','1342188038','0'],\
              ['HD104237','1342207819','1342207820','0'],\
              ['HD135344B-1','1342213921','1342213922','0'],\
              ['HD139614','1342215683','1342215684','0'],\
              ['HD141569','1342213913','0','0'],\
              ['HD142527','1342216174','1342216175','0'],\
              ['HD142666','1342213916','0','0'],\
              ['HD144432','1342213919','0','0'],\
              ['HD144668','1342215641','1342215642','0'],\
              ['HD150193','1342227068','0','0'],\
              ['HD163296','1342217819','1342217820','0'],\
              ['HD169142','1342206987','1342206988','0'],\
              ['HD179218','1342208884','1342208885','0'],\
              ['HD203024','1342206975','0','0'],\
              ['HD245906','1342228528','0','0'],\
              ['HD35187','1342217846','0','0'],\
              ['HD36112','1342228247','1342228248','0'],\
              ['HD38120','1342226212','1342226213','0'],\
              ['HD50138','1342206991','1342206992','0'],\
              ['HD97048','1342199412','1342199413','0'],\
              ['HD98922','1342210385','0','0'],\
              ['HH46','0','0','1342245084'],\
              ['HH100','0','0','1342252897'],\
              ['HT_Lup','1342213920','0','0'],\
              ['IRAM04191','1342216654','1342216655','0'],\
              ['IRAS03245+3002','1342214677','1342214676','0'],\
              ['IRAS03245','0','0','1342249053'],\
              ['IRAS03301+3111','1342215668','1342216181','0'],\
              ['IRAS03301','0','0','1342249477'],\
              ['IRAS12496','1342188039','1342188040','0'],\
              ['DKCha','0','0','1342254037'],\
              ['IRAS15398','0','0','1342250515'],\
              ['IRS46','1342228474','1342228475','1342251289'],\
              ['IRS48','1342227069','1342227070','0'],\
              ['IRS63','1342228473','1342228472','0'],\
              ['L1014','1342208911','1342208912','1342245857'],\
              ['L1157','1342208909','1342208908','1342247625'],\
              ['L1448-MM','1342213683','1342214675','0'],\
              ['L1455-IRS3','1342204122','1342204123','0'],\
              ['L1455-IRS','0','0','1342249474'],\
              ['L1489','1342216216','1342216215','0'],\
              ['L1527','1342192981','1342192982','0'],\
              ['L1551-IRS5','1342192805','1342229711','1342249470'],\
              ['L483','0','0','1342253649'],\
              ['L723-MM','0','0','1342245094'],\
              ['RCrA-IRS5A','1342207806','1342207805','1342253646'],\
              ['RCrA-IRS7B','1342207807','1342207808','1342242620'],\
              ['RCrA-IRS7C','1342206990','1342206989','1342242621'],\
              ['RNO90','1342228206','0','0'],\
              ['RNO91','0','0','1342251285'],\
              ['RU_Lup','1342215682','0','0'],\
              ['RY_Lup','1342216171','0','0'],\
              ['S_Cra','1342207809','1342207810','0'],\
              ['SR21','1342227209','1342227210','0'],\
              ['Serpens-SMM3','1342193216','1342193214','0'],\
              ['Serpens-SMM4','1342193217','1342193215','0'],\
              ['TMC1','1342225803','1342225804','1342250512'],\
              ['TMC1A','1342192987','1342192988','1342250510'],\
              ['TMR1','1342192985','1342192986','1342250509'],\
              ['V1057_Cyg-1','1342235853','1342235852','0'],\
              ['V1057Cyg','0','0','1342221695'],\
              ['V1331_Cyg-1','1342233446','1342233445','0'],\
              ['V1331Cyg','0','0','1342221694'],\
              ['V1515_Cyg-1','1342235691','1342235690','0'],\
              ['V1515Cyg','0','0','1342221685'],\
              ['V1735_Cyg-1','1342235849','1342235848','0'],\
              ['V1735Cyg','0','0','1342219560'],\
              ['VLA1623-243','1342213918','1342213917','0'],\
              ['VLA1623','0','0','1342251287'],\
              ['WL12','1342228187','1342228188','1342251290']]

line_fitting_dir = '/home/bettyjo/yaolun/programs/line_fitting/'
outdir =

for obsid in obsid_list:
    if (obsid[1] == '0') or (obsid[2] == '0'):
        continue
    fits_for_header = '/scratch/CDF_PACS_HSA/'+obsid[1]+'/herschel.pacs.signal.PacsRebinnedCube/hpacs'+obsid[1]+'_20hps3drbs_00fits'
    cdf_pacs_1d(osbid[1:3], obsid[0], outdir, fits_for_header, line_fitting_dir)

def obj_com(indir, noise_fix=False):
    """
    Test the completeness of the objects directories and the fitting results of both 1d and cube
    Usage:
        obj_com('/data/FWD_bettyjo/FWD_archive_slim/')
    """
    import os
    home = os.path.expanduser('~')
    # temp.
    import numpy as np
    import sys
    sys.path.append('/home/bettyjo/yaolun/programs/line_fitting')
    from extract_noise import extract_noise

    # pre-define the full object list

    pacsobj = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
               'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
               'HD97048','HD98922','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
               'L1489','L1527','L1551-IRS5','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RULup','RYLup','SCra','SR21',\
               'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']


    spireobj = ['B1-a','B1-c','B335','BHR71','Ced110-IRS4','FUOri','GSS30-IRS1','HH46','HH100','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','L1014',\
               'L1157','L1455-IRS3','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO91','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg',\
               'V1515Cyg','V1735Cyg','VLA1623','WL12']

    objlist = [pacsobj, spireobj]

    objlist = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','Ced110-IRS4','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
               'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
               'HD97048','HD98922','HH46','HH100','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
               'L1489','L1527','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RNO91','RULup','RYLup','SCra','SR21',\
               'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']

    objdir = os.walk(home+indir).next()[1]
    objdir = objdir[objdir != 'contour']
    err = 0
    # temp
    ra_std = np.empty(len(pacsobj))
    dec_std = np.empty(len(pacsobj))
    diff = []
    for o in objlist:
        if len(objdir[objdir == o]) != 1:
            print 'Cannot find ', o
            err += 1
        if o in pacsobj:
            # Check 1d and cube fitting results
            if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_lines.txt') == False:
                err += 1
                print 'Missing PACS 1d fitting on ', o
            # # temp. test for oversampling rate
            # else:
            #   print o, len(open(home+indir+'/'+o+'/pacs/data/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt','r').readlines())
            if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/cube/'+o+'_pacs_pixel13_os8_sf7_lines.txt') == False:
                err += 1
                print 'Missing PACS cube fitting on ', o
            else:
                wl, ra, dec = np.genfromtxt(home+indir+'/'+o+'/pacs/data/cube/'+o+'_pacs_pixel13_os8_sf7_coord.txt',skip_header=1).T
                # print o, np.std(ra*np.cos(np.mean(dec)*np.pi/180.)*3600), np.std(dec*3600)
                ra_std[pacsobj.index(o)] = np.std(ra*np.cos(np.mean(dec)*np.pi/180.)*3600)
                dec_std[pacsobj.index(o)] = np.std(dec*3600)
            #   print o, len(open(home+indir+'/'+o+'/pacs/data/cube/'+o+'_pacs_pixel13_os8_sf7.txt','r').readlines())
            if noise_fix == True:
                extract_noise(home+indir+'/'+o+'/pacs/advanced_products/', o, pacs=True)
                for i in range(1, 26):
                    extract_noise(home+indir+'/'+o+'/pacs/advanced_products/cube/', o, cube=str(i), pacs=True)
            else:
                if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_noise.txt'):
                    os.remove(home+indir+'/'+o+'/pacs/advanced_products/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim_noise.txt')
                for i in range(1, 26):
                    if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products/cube/'+o+'_pacs_pixel'+str(i)+'_os8_sf7_noise.txt'):
                    os.remove(home+indir+'/'+o+'/pacs/advanced_products/cube/'+o+'_pacs_pixel'+str(i)+'_os8_sf7_noise.txt')

        if o in spireobj:
            # Check 1d and cube fitting results
            if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/'+o+'_spire_corrected_lines.txt') == False:
                if (o != 'IRS46') and (o != 'HH100') and (o != 'V1735Cyg'):
                    err += 1
                    print 'Missing SPIRE 1d fitting on ', o
            if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_SLWC3_lines.txt') == False:
                err += 1
                print 'Missing SPIRE-SLW cube fitting on ', o
            if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_SSWD4_lines.txt') == False:
                err += 1
                print 'Missing SPIRE-SSW cube fitting on ', o
            # temp.
            # check the completeness of SPIRE spectra
            if os.path.exists(home+indir+'/'+o+'/spire/data/'+o+'_spire_corrected.txt') == False:
                if (o != 'IRS46') and (o != 'HH100') and (o != 'V1735Cyg'):
                    err += 1
                    print 'Cannot find SPIRE 1-D ASCII spectra.  Bug in the fitting routine.'
            else:
                (wl, flux) = np.genfromtxt(home+indir+'/'+o+'/spire/data/'+o+'_spire_corrected.txt', skip_header=1).T
                if min(wl) > 220:
                    err += 1
                    print 'SSW spectra is not included.'
            # print the offset between PACS and SPIRE
            if (o in pacsobj) and (o != 'IRS46') and (o != 'HH100') and (o != 'V1735Cyg'):
                (wl_pacs, flux_pacs, unc_pacs) = np.genfromtxt(home+indir+'/'+o+'/pacs/data/'+o+'_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt', skip_header=1).T
                (wl_spire, flux_spire) = np.genfromtxt(home+indir+'/'+o+'/spire/data/'+o+'_spire_corrected.txt', skip_header=1).T
                print flux_pacs[wl_pacs == max(wl_pacs)], flux_spire[wl_spire == min(wl_spire)], flux_spire[wl_spire == min(wl_spire)]-flux_pacs[wl_pacs == max(wl_pacs)]
                diff.append(float(flux_spire[wl_spire == min(wl_spire)]-flux_pacs[wl_pacs == max(wl_pacs)]))


            if noise_fix == True:
                spaxel = ['SLWA1','SLWA2','SLWA3','SLWB1','SLWB2','SLWB3','SLWB4','SLWC1','SLWC2','SLWC3','SLWC4','SLWC5','SLWD1','SLWD2','SLWD3','SLWD4','SLWE1','SLWE2','SLWE3',\
                          'SSWA1','SSWA2','SSWA3','SSWA4','SSWB1','SSWB2','SSWB3','SSWB4','SSWB5','SSWC1','SSWC2','SSWC3','SSWC4','SSWC5','SSWC6','SSWD1','SSWD2','SSWD3','SSWD4',\
                          'SSWD6','SSWD7','SSWE1','SSWE2','SSWE3','SSWE4','SSWE5','SSWE6','SSWF1','SSWF2','SSWF3','SSWF5','SSWG1','SSWG2','SSWG3','SSWG4']
                extract_noise(home+indir+'/'+o+'/spire/advanced_products/', o, spire=True)
                for i in range(0, len(spaxel)):
                    extract_noise(home+indir+'/'+o+'/spire/advanced_products/cube/', o, cube=spaxel[i], spire=True)
            else:
                if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/'+o+'_spire_corrected_noise.txt'):
                    os.remove(home+indir+'/'+o+'/spire/advanced_products/'+o+'_spire_corrected_noise.txt')
                spaxel = ['SLWA1','SLWA2','SLWA3','SLWB1','SLWB2','SLWB3','SLWB4','SLWC1','SLWC2','SLWC3','SLWC4','SLWC5','SLWD1','SLWD2','SLWD3','SLWD4','SLWE1','SLWE2','SLWE3',\
                          'SSWA1','SSWA2','SSWA3','SSWA4','SSWB1','SSWB2','SSWB3','SSWB4','SSWB5','SSWC1','SSWC2','SSWC3','SSWC4','SSWC5','SSWC6','SSWD1','SSWD2','SSWD3','SSWD4',\
                          'SSWD6','SSWD7','SSWE1','SSWE2','SSWE3','SSWE4','SSWE5','SSWE6','SSWF1','SSWF2','SSWF3','SSWF5','SSWG1','SSWG2','SSWG3','SSWG4']
                for i in range(0, len(spaxel)):
                    if os.path.exists(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_'+str(spaxel[i])+'_noise.txt'):
                        os.remove(home+indir+'/'+o+'/spire/advanced_products/cube/'+o+'_'+str(spaxel[i])+'_noise.txt')

    print min(ra_std), max(ra_std), np.mean(ra_std)
    print min(dec_std), max(dec_std), np.mean(dec_std)
    diff = np.array(diff)
    print min(diff), max(diff), np.mean(diff)
    if err == 0:
        print 'Passed the object test!'
        return ra_std, dec_std# True

def fits_com(indir):
    """
    Check the completeness of the FITS files that should be included in the release

    Usage:
        fits_com('/data/FWD_bettyjo/FWD_archive_slim/')
    """
    import os
    import glob
    from astropy.io import fits

    home = os.path.expanduser('~')

    # pre-define OBSID info

    objlist = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','Ced110-IRS4','DGTau','EC82','Elias29','FUOri','FUOri','GSS30-IRS1','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
               'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
               'HD97048','HD98922','HH46','HH100','HTLup','IRAM04191','IRAS03245','IRAS03245','IRAS03301','IRAS03301','IRAS12496','IRAS12496','IRAS15398','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3','L1455-IRS3',\
               'L1489','L1527','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RNO91','RULup','RYLup','SCra','SR21',\
               'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1057Cyg','V1331Cyg','V1331Cyg','V1515Cyg','V1515Cyg','V1735Cyg','V1735Cyg','VLA1623','VLA1623','WL12']

    jitter_exclude = ['IRAM04191','IRS46','L1014','L1455-IRS3','RCrA-IRS5A','RCrA-IRS7C','Serpens-SMM4','EC82','HD98922','HD245906','HD203024','HTLup','HD142666','HD35187']

    obsid = [['AB_Aur','1342217842','1342217843','0'],\
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
    err = 0
    for i in range(0, len(obsid)):
        for j in range(1,3):
            if obsid[i][j] != '0':
                # check FITS files for PACS 1d and cube
                # nojitter first
                # 1d - 9Spx
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_central9Spaxels_PointSourceCorrected_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_blue_central9Spaxels_PointSourceCorrected_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_central9Spaxels_PointSourceCorrected_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_red_central9Spaxels_PointSourceCorrected_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                # 1d - 3x3NO
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                # 1d - 3x3YES
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7_nojitter.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                    err += 1
                # cube - finalcube
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_finalcubes_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_blue_finalcubes_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_finalcubes_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_red_finalcubes_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                # cube - noda
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnoda_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_blue_rebinnedcubesnoda_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnoda_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_red_rebinnedcubesnoda_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                # cube - nodb
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnodb_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_blue_rebinnedcubesnodb_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnodb_*sf7_nojitter.fits')) == 0:
                    print '%s missing OBSID_%s_red_rebinnedcubesnodb_os8_sf7_nojitter.fits' % (objlist[i], obsid[i][j])
                    err += 1
                # print objlist[i]
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_central9Spaxels_PointSourceCorrected_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_central9Spaxels_PointSourceCorrected_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_finalcubes_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_finalcubes_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnoda_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnoda_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnodb_*sf7_nojitter.fits')[0])[0].header['DATE']
                # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnodb_*sf7_nojitter.fits')[0])[0].header['DATE']

                # check jitter FITS
                if objlist[i] in jitter_exclude == False:
                    # 1d - 9Spx
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_central9Spaxels_PointSourceCorrected_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_blue_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_central9Spaxels_PointSourceCorrected_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_red_central9Spaxels_PointSourceCorrected_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    # 1d - 3x3NO
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    # 1d - 3x3YES
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_%s_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_00_os8sf7.fits' % (objlist[i], obsid[i][j], obsid[i][0])
                        err += 1
                    # cube - finalcube
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_finalcubes_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_blue_finalcubes_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_finalcubes_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_red_finalcubes_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1
                    # cube - noda
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnoda_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_blue_rebinnedcubesnoda_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnoda_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_red_rebinnedcubesnoda_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1
                    # cube - nodb
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnodb_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_blue_rebinnedcubesnodb_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1
                    if len(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnodb_*sf7.fits')) == 0:
                        print '%s missing OBSID_%s_red_rebinnedcubesnodb_os8_sf7.fits' % (objlist[i], obsid[i][j])
                        err += 1

                    # print objlist[i]
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_central9Spaxels_PointSourceCorrected_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_central9Spaxels_PointSourceCorrected_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3NO_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_blue_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_'+obsid[i][0]+'_red_centralSpaxel_PointSourceCorrected_Corrected3x3YES_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_finalcubes_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_finalcubes_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnoda_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnoda_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_blue_rebinnedcubesnodb_*sf7.fits')[0])[0].header['DATE']
                    # print fits.open(glob.glob(home+indir+'/'+objlist[i]+'/pacs/data/fits/OBSID_'+obsid[i][j]+'_red_rebinnedcubesnodb_*sf7.fits')[0])[0].header['DATE']

        if obsid[i][3] != '0':
            # check FITS files for SPIRE 1d and cube
            if (os.path.exists(home+indir+'/'+objlist[i]+'/spire/data/fits/'+obsid[i][0]+'_spire_corrected.fits') or \
                os.path.exists(home+indir+'/'+objlist[i]+'/spire/data/fits/'+obsid[i][0].lower()+'_spire_corrected.fits')) == False:
                # spire-1d spectra of IRS46, HH100, and V1735Cyg are mis-pointed
                if (objlist[i] != 'IRS46') and (objlist[i] != 'HH100') and (objlist[i] != 'V1735Cyg'):
                    print '%s missing %s_spire_corrected.fits' % (objlist[i], obsid[i][0].lower())
                    err += 1
            if os.path.exists(home+indir+'/'+objlist[i]+'/spire/data/fits/'+obsid[i][3]+'_spectrum_extended_HR_aNB_15.fits') == False:
                print '%s missing %s_spectrum_extended_HR_aNB_15.fits' % (objlist[i], obsid[i][3])
                err += 1
    if err == 0:
        print 'passed the FITS test!'
        return True

def strong_line(indir):
    """
    Return the number of strong lines (SNR >= 10) for A given source. (Return the total number of strong lines found
    in 1d and all spaxels)
    Usage:
        strong_line('/FWD_archive/CDF_archive/BHR71/')
    """
    from astropy.io import ascii
    import numpy as np
    import os
    home = os.path.expanduser('~')

    # Child function to find the files that match the given pattern
    def find(pattern, path):
        import os, fnmatch
        result = []
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result

    # Search for the fitting results tables under the primary directory
    # only search in the the object directory
    path2fit = find('*_lines.txt', indir)
    # Header of the all cube fitting results
    # ====================================================================================================
    # Object,   Line,           LabWL(um),      ObsWL(um),      Sig_Cen(um),    Str(W/cm2),     Sig_str(W/cm2)
    # FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), Noise(W/cm2/um),SNR,            E_u(K),         A(s-1)        
    # g,        RA(deg),        Dec(deg),       Blend,          Validity
    # ====================================================================================================
    num_strong = 0
    # Read the data
    for path in path2fit:
        data = ascii.read(path)
        header = data.colnames
        data = data[(np.isnan(data['SNR'])!=True) & (data['Validity']==1)]  # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
        num_strong = num_strong + len(data[data['SNR'] >= 10.0])

    return num_strong

# print strong_line('/test/BHR71/pacs/advanced_products/')

def unc_test(filepath,plotdir,png=True):
    """
    Show the uncertainty relation of the measurement from baseline flucuation and the measurement from the fitting routine.
    Usage:
        unc_test('/FWD_archive/CDF_archive/CDF_archive_pacs_1d_lines.txt', '/analysis/')
    """
    from astropy.io import ascii
    import numpy as np
    # to avoid X server error
    import matplotlib as mpl
    mpl.use('Agg')
    #
    import matplotlib.pyplot as plt
    import os
    import astropy.constants as const
    from scipy.interpolate import interp1d
    home = os.path.expanduser('~')

    plotname = os.path.splitext(os.path.basename(filepath))[0]

    filepath = home + filepath
    data = ascii.read(filepath)
    if 'Str(W/cm2/as2)' in data.columns:
        unit = '/as2'
    else:
        unit = ''
    # Header of the all cube fitting results
    # ====================================================================================================
    # Object,   Line,           LabWL(um),      ObsWL(um),      Sig_Cen(um),    Str(W/cm2),     Sig_str(W/cm2)
    # FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), Noise(W/cm2/um),SNR,            E_u(K),         A(s-1)        
    # g,        RA(deg),        Dec(deg),       Blend,          Validity
    # ====================================================================================================
    header = data.colnames
    data = data[(np.isnan(data['SNR'])!=True) & (data['Validity']==1)]  # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
    snr = abs(data['SNR'][np.argsort(data['ObsWL(um)'])])

    snr_flux = (data['Str(W/cm2'+unit+')']/data['Sig_str(W/cm2'+unit+')'])[np.argsort(data['ObsWL(um)'])]
    wl = data['ObsWL(um)'][np.argsort(data['ObsWL(um)'])]

    sig_str = data['Sig_str(W/cm2'+unit+')'][np.argsort(data['ObsWL(um)'])]
    noise = data['Noise(W/cm2/um'+unit+')'][np.argsort(data['ObsWL(um)'])]
    fwhm = data['FWHM(um)'][np.argsort(data['ObsWL(um)'])]


    print 'Plotting scatter/histogram...'
    from matplotlib.ticker import NullFormatter

    nullfmt   = NullFormatter()         # no labels

    # a = sig_str[(snr >= 3.0) & (snr < 10.0)]
    # dum = data[np.argsort(data['ObsWL(um)'])]
    # print dum[(snr >= 3.0) & (snr < 10.0)][a <=0]

    # data
    x = np.log10( noise[(snr >= 3.0) & (snr < 10.0)] * fwhm[(snr >= 3.0) & (snr < 10.0)] )
    y = np.log10( sig_str[(snr >= 3.0) & (snr < 10.0)] )

    xx = np.log10( noise[(snr >= 10.0)] * fwhm[(snr >= 10.0)] )
    yy = np.log10( sig_str[(snr >= 10.0)] )

    line = np.hstack((np.log10(sig_str[snr >= 3.0]), np.log10(noise[snr >= 3.0]*fwhm[snr >= 3.0])))

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(1,figsize=(12,12))

    axScatter = fig.add_axes(rect_scatter)
    axHistx = fig.add_axes(rect_histx)
    axHisty = fig.add_axes(rect_histy)

    # set labels for scatter plot
    axScatter.set_xlabel(r'$\rm{log(Noise\times FWHM)\,(W\,cm^{-2})}$', fontsize=22)
    axScatter.set_ylabel(r'$\rm{log(\sigma_{fit})\,(W\,cm^{-2})}$', fontsize=22)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # no tick labels
    axHistx.yaxis.set_ticks([])
    axHisty.xaxis.set_ticks([])

    # the scatter plot:
    low = axScatter.scatter(x, y, c='g', s=27, lw=0.5)
    high = axScatter.scatter(xx, yy, c='r', s=27, lw=0.5)
    slope, = axScatter.plot(line, line, 'b-', linewidth=1.5)
    lg = plt.legend([low, high, slope],[r'$\rm{3<SNR<10}$',r'$\rm{10<SNR}$',r'$\rm{Equality}$'],\
                    loc='best', bbox_to_anchor=[0.35,1],bbox_transform=axScatter.transAxes, numpoints=1, scatterpoints=1, fontsize=18)
    
    # now determine nice limits by hand:
    binwidth = 0.05
    xymax = np.nanmax( [np.max(x), np.max(y)] )
    xymin = np.nanmin( [np.min(x), np.min(y)] )
    ulim = ( int(xymax/binwidth) + 1) * binwidth + 0.5
    llim = ( int(xymin/binwidth) + 1) * binwidth - 0.5

    axScatter.set_xlim( (llim, ulim) )
    axScatter.set_ylim( (llim, ulim) )

    bins = np.arange(llim, ulim + binwidth , binwidth)
    axHistx.hist(x, bins=bins, histtype='step', edgecolor='Green',hatch='////')
    axHisty.hist(y, bins=bins, histtype='step', edgecolor='Green',hatch='////', orientation='horizontal')

    axHistx.hist(xx, bins=bins, histtype='step', edgecolor='Red', hatch='\\\\\\\\',)
    axHisty.hist(yy, bins=bins, histtype='step', edgecolor='Red', hatch='\\\\\\\\', orientation='horizontal')

    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    # Tweak the axes thickness
    axes = [axHistx,axHisty]
    for ax in axes:
        ax.tick_params('both',labelsize=14,width=1.5,which='major')
        ax.tick_params('both',labelsize=14,width=1.5,which='minor')
        [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]
    ax = axScatter
    ax.tick_params('both',labelsize=14,width=1.5,which='major')
    ax.tick_params('both',labelsize=14,width=1.5,which='minor')
    [ax.spines[axis].set_linewidth(1.5) for axis in ['top','bottom','left','right']]

    if png == False:
        plt.savefig(home+plotdir+plotname+'_scatter_hist.pdf',format='pdf',dpi=300,bbox_inches='tight')
    else:
        plt.savefig(home+plotdir+plotname+'_scatter_hist.png',format='png',dpi=300,bbox_inches='tight')
    plt.cla()
    plt.clf()

    print 'Finished the uncertainty plots, check them!'

# unc_test('/test/fixedwidth/CDF_archive_pacs_1d_lines.txt', '/test/')

def fitting_check(indir,outdir):
    """
    Print the statistic of the fitting results separated by PACS and SPIRE.
    Usage:
        fitting_check('/FWD_archive/CDF_archive/', '/analysis/')
    """
    from astropy.io import ascii
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    home = os.path.expanduser('~')

    foo = open(home+outdir+'stat.txt','w')

    # Child function to find the files that match the given pattern
    def find(pattern, path):
        import os, fnmatch
        result = []
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result

    # Search for the fitting results tables under the primary directory
    # If the object name is specified, then only search in the the object directory
    # objdir = [x[0] for x in os.walk(home+indir)][1:]
    objdir = os.walk(home+indir).next()[1]
    pacspath = []
    spirepath = []
    num_pacs = 0.0
    num_spire = 0.0
    for o in objdir:
        # PACS search
        path_dum = find('*_lines.txt', home+indir+'/'+o+'/pacs/')
        if len(path_dum) > 0:
            pacspath.extend(path_dum)
            num_pacs  += 1
        # SPIRE search
        path_dum = find('*_lines.txt', home+indir+'/'+o+'/spire/')
        if len(path_dum) > 0:
            spirepath.extend(path_dum)
            num_spire += 1

    num_fit = 0.0
    num_line = 0.0
    num1 = 0.0
    num2 = 0.0
    num3 = 0.0
    num4 = 0.0
    num5 = 0.0
    num6 = 0.0
    num7 = 0.0
    num8 = 0.0
    num9 = 0.0
    num10 = 0.0
    num11 = 0.0
    num_line2 = 0.0
    num_dg = 0.0
    num_dg_line = 0.0

    num_test = 0.0

    # PACS statistic
    for path in pacspath:
        data = ascii.read(path)
        # Header of the 1-D fitting results
        # =========================================================================================
        # Line,     LabWL(um),      ObsWL(um),      Sig_Cen(um),    Str(W/cm2),     Sig_str(W/cm2)
        # FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,            E_u(K),         A(s-1)        
        # g,        RA(deg),        Dec(deg),       Blend,          Validity
        # =========================================================================================
        header = data.colnames
        # If there is a missing segment in the spectra, the SNR will return NaN.
        num1 = len(data[np.isnan(data['SNR'])==True]) + num1
        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
        data = data[np.isnan(data['SNR'])!=True]
        num_fit = len(data['SNR']) + num_fit
        num_line = len(data[data['SNR'] >= 3]) + num_line
        num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
        num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
        num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
        num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
        num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
        num7 = len(data[data['Sig_str(W/cm2)'] == 0.0]) + num7
        num8 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) + num8
        num9 = len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
        num10 = len(data[(data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & (data['Sig_str(W/cm2)'] != 0.0) & \
                         (data['Validity'] == 0)]) + num10
        num11 = len(data[(data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & (data['Sig_str(W/cm2)'] != 0.0) & \
                         (data['Validity'] == 0) & (data['SNR'] >= 3)]) + num11

        num_line2 = len(data[(data['SNR'] >=3) & (data['Validity'] ==1)]) + num_line2
        num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
        num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

        # Print out the detail information of the line that has zero in line strength uncertainty.
        if len(data[(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]) != 0:
            print path
            print data['Line'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['ObsWL(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['Sig_Cen(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['Str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['Sig_str(W/cm2)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]
            print data['Sig_FWHM(um)'][(data['Sig_str(W/cm2)'] == 0.0) & (data['SNR'] >= 3)]

    # Print out the statistic of the pacs fitting results
    foo.write('<PACS>\n')
    foo.write('\t Number of object: %d \n' % num_pacs)
    foo.write('\t %d lines fitted, %.2f lines fitted per object\n' % (num_fit,num_fit/num_pacs))
    foo.write('\t %d detections w/ anomalous, %.2f detections per object.\n' % (num_line,num_line/num_pacs))
    foo.write('\t %d lines fitted with blend Gaussian, %d lines detections among them.\n' % (num_dg,num_dg_line))
    foo.write('\t <<Anomaly>>\n')
    foo.write('\t \t SNR anomalies due to the missing spectra: %d \n' % num1)
    foo.write('\t \t Zeros in line centroid uncertainty: %d and %d with detections.\n' % (num2,num3))
    foo.write('\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.\n' % (num4,num5,num6))
    foo.write('\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.\n' % (num7,num8,num9))
    foo.write('\t \t Validity = 0 only due to blending with neighbors: %d, %d with detections.\n' % (num10, num11))
    foo.write('\t %d detections without anomalous, and %.2f lines per object.\n' % (num_line2,num_line2/num_pacs))

    foo.write('====================================================================================================================\n')
    num_fit = 0.0
    num_line = 0.0
    num1 = 0.0
    num2 = 0.0
    num3 = 0.0
    num4 = 0.0
    num5 = 0.0
    num6 = 0.0
    num7 = 0.0
    num8 = 0.0
    num9 = 0.0
    num10 = 0.0
    num11 = 0.0
    num_line2 = 0.0
    num_dg = 0.0
    num_dg_line = 0.0

    num_test = 0.0

    # SPIRE statistic
    for path in spirepath:
        data = ascii.read(path)
        if 'Str(W/cm2/as2)' in data.columns:
            unit = '/as2'
        else:
            unit = ''
        # Header of the 1-D fitting results
        # =========================================================================================
        # Line,     LabWL(um),      ObsWL(um),      Sig_Cen(um),    Str(W/cm2),     Sig_str(W/cm2)
        # FWHM(um), Sig_FWHM(um),   Base(W/cm2/um), SNR,            E_u(K),         A(s-1)        
        # g,        RA(deg),        Dec(deg),       Blend,          Validity
        # =========================================================================================
        header = data.colnames
        # If there is a missing segment in the spectra, the SNR will return NaN.
        num1 = len(data[np.isnan(data['SNR'])==True]) + num1
        # Temperory procedure to exclude the missing segment in the spectrum resulting in the NaN in SNR
        data = data[np.isnan(data['SNR'])!=True]
        num_fit = len(data['SNR']) + num_fit
        num_line = len(data[data['SNR'] >= 3]) + num_line
        num2 = len(data[data['Sig_Cen(um)'] == -999.0]) + num2
        num3 = len(data[(data['Sig_Cen(um)'] == -999.0) & (data['SNR'] >= 3)]) + num3
        num4 = len(data[data['Sig_FWHM(um)'] == -999.0]) + num4
        num5 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3)]) + num5
        num6 = len(data[(data['Sig_FWHM(um)'] == -999.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num6
        num7 = len(data[data['Sig_str(W/cm2'+unit+')'] == 0.0]) + num7
        num8 = len(data[(data['Sig_str(W/cm2'+unit+')'] == 0.0) & (data['SNR'] >= 3)]) + num8
        num9 = len(data[(data['Sig_str(W/cm2'+unit+')'] == 0.0) & (data['SNR'] >= 3) & (data['Blend'] == 'DoubleGaussian')]) + num9
        num10 = len(data[(data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & (data['Sig_str(W/cm2'+unit+')'] != 0.0) & \
                         (data['Validity'] == 0)]) + num10
        num11 = len(data[(data['Sig_Cen(um)'] != -999.0) & (data['Sig_FWHM(um)'] != -999.0) & (data['Sig_str(W/cm2'+unit+')'] != 0.0) & \
                         (data['Validity'] == 0) & (data['SNR'] >= 3)]) + num11

        num_line2 = len(data[(data['SNR'] >=3) & (data['Validity'] == 1)]) + num_line2
        num_dg = len(data[data['Blend'] == 'DoubleGaussian']) + num_dg
        num_dg_line = len(data[(data['Blend'] == 'DoubleGaussian') & (data['SNR'] >= 3)]) + num_dg_line

    # Print out the statistic of the pacs fitting results
    foo.write('<SPIRE>\n')
    foo.write('\t Number of object: %d \n' % num_spire)
    foo.write('\t %d lines fitted, %.2f lines fitted per object.\n' % (num_fit,num_fit/num_spire))
    foo.write('\t %d detections, %.2f detections per object.\n' % (num_line,num_line/num_spire))
    foo.write('\t %d lines fitted with blend Gaussian, %d lines detections among them.\n' % (num_dg,num_dg_line))
    foo.write('\t <<Anomaly>>\n')
    foo.write('\t \t SNR anomalies due to the missing spectra: %d \n' % num1)
    foo.write('\t \t Zeros in line centroid uncertainty: %d and %d with detections.\n' % (num2,num3))
    foo.write('\t \t Zeros in FWHM uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.\n' % (num4,num5,num6))
    foo.write('\t \t Zeros in line strength uncertainty: %d, %d with detections, and %d with detections and blend Gaussian.\n' % (num7,num8,num9))
    foo.write('\t \t Validity = 0 only due to blending with neighbors: %d, %d with detections.\n' % (num10, num11))
    foo.write('\t %d detections without anomalous, and %.2f lines per object.\n' % (num_line2,num_line2/num_spire))

    foo.close()

    print 'Finished the statistic of the fitting results, go check them at %s !' % home+outdir+'stat.txt'

def pointing_test(indir):
    import numpy as np
    from astropy.io import ascii

    pacsobj = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
               'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
               'HD97048','HD98922','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
               'L1489','L1527','L1551-IRS5','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RULup','RYLup','SCra','SR21',\
               'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']


    spireobj = ['B1-a','B1-c','B335','BHR71','Ced110-IRS4','FUOri','GSS30-IRS1','HH46','HH100','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','L1014',\
               'L1157','L1455-IRS3','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO91','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg',\
               'V1515Cyg','V1735Cyg','VLA1623','WL12']

    objlist = [pacsobj, spireobj]

    objlist = ['ABAur','AS205','B1-a','B1-c','B335','BHR71','Ced110-IRS4','DGTau','EC82','Elias29','FUOri','GSS30-IRS1','HD100453','HD100546','HD104237','HD135344B','HD139614',\
               'HD141569','HD142527','HD142666','HD144432','HD144668','HD150193','HD163296','HD169142','HD179218','HD203024','HD245906','HD35187','HD36112','HD38120','HD50138',\
               'HD97048','HD98922','HH46','HH100','HTLup','IRAM04191','IRAS03245','IRAS03301','IRAS12496','IRAS15398','IRS46','IRS48','IRS63','L1014','L1157','L1448-MM','L1455-IRS3',\
               'L1489','L1527','L1551-IRS5','L483','L723-MM','RCrA-IRS5A','RCrA-IRS7B','RCrA-IRS7C','RNO90','RNO91','RULup','RYLup','SCra','SR21',\
               'Serpens-SMM3','Serpens-SMM4','TMC1','TMC1A','TMR1','V1057Cyg','V1331Cyg','V1515Cyg','V1735Cyg','VLA1623','WL12']

    # read the PACS coordinates
    # only check the central spaxel
    for obj in pacsobj:
        coord_dum = np.genfromtxt(indir+'/'+obj+'/pacs/data/cube/'+obj+'_pacs_pixel13_os8_sf7_coord.txt')

def cdf_test(indir,outdir):
    """
    The main test script for checking the performance of the archive and fitting results.
    Usage:
        cdf_test('/FWD_archive/CDF_archive', '/FWD_archive/')
    """
    import numpy as np
    import os
    home = os.path.expanduser('~')

    # Check the existence of the outdir
    if os.path.exists(home+outdir) == False:
        os.makedirs(home+outdir)

    # Object completeness test
    obj_com(indir)

    # FITS files completeness test
    fits_com(indir)

    # Strong lines test
    objdir = np.array(os.walk(home+indir).next()[1])
    objdir = np.sort(objdir[objdir != 'contour'])
    for o in objdir:
        # test pacs fitting
        if os.path.exists(home+indir+'/'+o+'/pacs/advanced_products') == True:
            print o, 'PACS:  ', strong_line(home+indir+'/'+o+'/pacs/advanced_products/')
        # test spire fitting
        if os.path.exists(home+indir+'/'+o+'/spire/advanced_products') == True:
            print o, 'SPIRE: ', strong_line(home+indir+'/'+o+'/spire/advanced_products/')

    # Uncertainty relation plots
    unc_test(indir+'/CDF_archive_pacs_1d_lines.txt', outdir)
    unc_test(indir+'/CDF_archive_pacs_cube_lines.txt', outdir)
    unc_test(indir+'/CDF_archive_spire_1d_lines.txt', outdir)
    unc_test(indir+'/CDF_archive_spire_cube_lines.txt', outdir)

    # Stats of the fitting results
    fitting_check(indir, outdir)




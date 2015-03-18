def oh2o():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/oh2o_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/oh2o_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/oh2o_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()

def ph2o():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/ph2o_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/ph2o_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/ph2o_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()
def co():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/co_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/co_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/co_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()

def oh():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        hps = raw_input('Do you want to use hyperfine structure lines? (yes or no)')
        if hps == 'no':
            level = np.genfromtxt(home+'/data/oh_level.txt', dtype='str')
            ref = np.genfromtxt(home+'/data/oh_ref.txt','str')
            for i in range(0,len(level[0:])):
                    for j in range(0, len(ref[0:])):
                            if ref[j,1] == level[i,0]:
                                    ref[j,0] = level[i,2]
                                    ref[j,1] = level[i,3]
                            if ref[j,2] == level[i,0]:
                                    ref[j,2] = level[i,3]
            c = 2.998e8
            ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
            ref = ref[np.argsort(ref[:,4].astype(float))]
            ref_sort = ref
            dummy = np.copy(ref[:,0])
            ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
            ref_sort[:,5] = dummy
            ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
            slt_trans = ref_sort[ind,:]
            print slt_trans
            print len(slt_trans[0,:,0])
            foo = open(home+'/data/oh_ref_sort.txt','w')
            np.savetxt(foo,ref_sort, fmt='%s')
            foo.close()
        else:
            if hps != 'yes':
                print 'I cannot recognize your answer. I take that as a yes.'
            level = np.genfromtxt(home+'/data/oh_hps_level.txt', dtype='str')
            ref = np.genfromtxt(home+'/data/oh_hps_ref.txt', dtype='str')
            

def co13():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/co13_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/co13_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/co13_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()

def hco():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/hco_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/hco_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/hco_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()
def ci():
        import numpy as np
        import os
        home = os.path.expanduser('~')
        band = raw_input('Select the band:')
        if band == 'pacs':
                upper, lower = 200, 54
        if band == 'spire':
                upper, lower = 671, 200
        level = np.genfromtxt(home+'/data/ci_level.txt', dtype='str')
        ref = np.genfromtxt(home+'/data/ci_ref.txt','str')
        for i in range(0,len(level[0:])):
                for j in range(0, len(ref[0:])):
                        if ref[j,1] == level[i,0]:
                                ref[j,0] = level[i,2]
                                ref[j,1] = level[i,3]
                        if ref[j,2] == level[i,0]:
                                ref[j,2] = level[i,3]
        c = 2.998e8
        ref[:,4] = c/ref[:,4].astype(float)/1e9*1e6
        ref = ref[np.argsort(ref[:,4].astype(float))]
        ref_sort = ref
        dummy = np.copy(ref[:,0])
        ref_sort[:,0],ref_sort[:,1],ref_sort[:,2],ref_sort[:,3],ref_sort[:,4] = ref[:,1],ref[:,2],ref[:,4],ref[:,3],ref[:,5]
        ref_sort[:,5] = dummy
        ind = np.where((ref_sort[:,2].astype(float) >= lower) & (ref_sort[:,2].astype(float) <= upper))
        slt_trans = ref_sort[ind,:]
        print slt_trans
        print len(slt_trans[0,:,0])
        foo = open(home+'/data/ci_ref_sort.txt','w')
        np.savetxt(foo,ref_sort, fmt='%s')
        foo.close()


species = raw_input('What molecule:')
if species == 'oh2o':
	oh2o()
elif species == 'ph2o':
	ph2o()
elif species == 'co':
	co()
elif species == 'oh':
    oh()
elif species == '13co':
    co13()
elif species == 'hco':
    hco()
elif species == 'ci':
    ci()

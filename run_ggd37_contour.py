# # LH
# # stitch all fitting tables for each pixel into a single table
# max_pixel = 300
# foo = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/LH/ggd37_irs_LH_lines.txt', 'w')
# for i in range(1,max_pixel+1):
#     foo_dum = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/LH/ggd37_irs_pixel'+str(i)+'_lines.txt', 'r').readlines()
#     foo_dum[0] = foo_dum[0].replace('Pixel_No.', 'Pixel')
#     if i == 1:
#         foo.write('{:>12s}  '.format('Object')+foo_dum[0][:-1]+'  {:>12s}\n'.format('Pixel_No.'))
#     for line in foo_dum[1:]:
#         foo.write('{:>12s}  '.format('GGD37')+line[:-1]+'  {:>12d}\n'.format(i))
#
# foo.close()
#
# # SH
# # stitch all fitting tables for each pixel into a single table
# max_pixel = 1260
# foo = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SH/ggd37_irs_SH_lines.txt', 'w')
# for i in range(1,max_pixel+1):
#     foo_dum = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SH/ggd37_irs_pixel'+str(i)+'_lines.txt', 'r').readlines()
#     foo_dum[0] = foo_dum[0].replace('Pixel_No.', 'Pixel')
#     if i == 1:
#         foo.write('{:>12s}  '.format('Object')+foo_dum[0][:-1]+'  {:>12s}\n'.format('Pixel_No.'))
#     for line in foo_dum[1:]:
#         foo.write('{:>12s}  '.format('GGD37')+line[:-1]+'  {:>12d}\n'.format(i))
#
# foo.close()
#
# SL
# stitch all fitting tables for each pixel into a single table
max_pixel = 6448
foo1 = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SL/ggd37_irs_SL1_lines.txt', 'w')
foo2 = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SL/ggd37_irs_SL2_lines.txt', 'w')

foo_dum = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SL/ggd37_irs_pixel1_lines.txt', 'r').readlines()
foo_dum[0] = foo_dum[0].replace('Pixel_No.', 'Pixel')
foo1.write('{:>12s}  '.format('Object')+foo_dum[0][:-1]+'  {:>12s}\n'.format('Pixel_No.'))
foo2.write('{:>12s}  '.format('Object')+foo_dum[0][:-1]+'  {:>12s}\n'.format('Pixel_No.'))

for i in range(1,max_pixel+1):
    foo_dum = open('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/SL/ggd37_irs_pixel'+str(i)+'_lines.txt', 'r').readlines()
    foo_dum[0] = foo_dum[0].replace('Pixel_No.', 'Pixel')
    if 'FeII53' in foo_dum[1]:
        for line in foo_dum[1:]:
            foo1.write('{:>12s}  '.format('GGD37')+line[:-1]+'  {:>12d}\n'.format(i))
    else:
        for line in foo_dum[1:]:
            foo2.write('{:>12s}  '.format('GGD37')+line[:-1]+'  {:>12d}\n'.format(i))

foo1.close()
foo2.close()


# get line list from the fitting table
from astropy.io import ascii
from CDF_contour import CDF_contour
modules = ['SH','LH','SL1','SL2']
# modules = ['SL1', 'SL2']

for m in modules:
    linelist = list(set(ascii.read('/Volumes/SD-Mac/research/GGD37/ggd37_fitting/'+m[:2]+'/ggd37_irs_'+m+'_lines.txt')['Line']))
    print(linelist)

    for line in linelist:
        if m == 'SL2':
            pix_cen = 27
        else:
            pix_cen = 1
        print('processing ', line)
        coord_rebin, z_list = CDF_contour(line, 'GGD37',
                                          '/Volumes/SD-Mac/research/GGD37/ggd37_fitting/'+m[:2]+'/ggd37_irs_'+m+'_lines.txt',
                                          pix_cen=pix_cen, cont=False,
                                          output='/Volumes/SD-Mac/research/GGD37/ggd37_fitting/'+m[:2]+'/'+line+'.fits',
                                          ggd37=False)

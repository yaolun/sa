def phot_filter(fil_name, filter_dir):
    import numpy as np
    from astropy.io import ascii
    from pprint import pprint
    import os

    filename = filter_dir+'allfilters.dat'

    foo = open(filename, 'r')
    index = []
    filter_name = []
    for i, line in enumerate(foo):
        if line[0] == '#':
            index.append(i)
            filter_name.append(line.split('#')[1].lstrip().rstrip())
    index = np.array(index)
    filter_name = np.array(filter_name)

    ind, = np.where(filter_name == fil_name)

    while len(ind) == 0:
        if fil_name != 'ls':
            print('requested filter not found in database!')
            fil_name = input('Please enter the filter name (or ls for listing the filters in database): ')
        if fil_name == 'ls':
            pprint(filter_name)
            fil_name = input('Please enter the filter name (or ls for listing the filters in database): ')
        ind, =np.where(filter_name == fil_name)
    if fil_name != filter_name[-1]:
        phot_filter = ascii.read(filename, data_start=index[ind]-ind, data_end=index[ind+1]-ind-1, names=['wave','transmission'],header_start=None)
    else:
        phot_filter = ascii.read(filename, data_start=index[ind]-ind, names=['wave','transmission'],header_start=None)

    return phot_filter

def CDF_fitting_stitch(indir, outdir):
	pacs_cube = open(indir+'CDF_archive_pacs_cube_lines.txt','r').readlines()

	foo = open(outdir+'CDF_archive_all_lines.txt','w')
	foo.write('%s \n' % (pacs_cube[0].rstrip()) )

	for line in open(indir+'CDF_archive_pacs_1d_lines.txt','r').readlines()[1:]:
		foo.write('%s %21s %s \n' % ( line[0:373].rstrip(), 'c', line[373:].rstrip()) )
	for line in open(indir+'CDF_archive_pacs_cube_lines.txt','r').readlines()[1:]:
		foo.write('%s \n' % line.rstrip() )
	for line in open(indir+'CDF_archive_spire_1d_lines.txt','r').readlines()[1:]:
		foo.write('%s %21s %s \n' % ( line[0:373].rstrip(), 'c', line[373:].rstrip()) )
	for line in open(indir+'CDF_archive_spire_cube_lines.txt','r').readlines()[1:]:
		foo.write('%s \n' % line.rstrip() )

	foo.close()

indir = '/Users/yaolun/data/CDF_archive/'
outdir = indir
CDF_fitting_stitch(indir, outdir)
#!/usr/bin/python
# This is HIPE script

# some parameter setup
obsid1 = 1342212230
obsid2 = 1342212231

outdir = '/Users/yaolun/bhr71/data/HSA/'

obs1 = getObservation(obsid=obsid1,useHsa=True)
obs2 = getObservation(obsid=obsid2,useHsa=True)

# get calTree and print out the version of calTree in the data and on the disk
# disk
caltree=getCalTree(obs=obs1)

# BHR71 blue rebinned cube
bhr71_rebinnedcube_blue1 = obs1.refs["level2"].product.refs["HPS3DRB"].product.refs[0].product
bhr71_rebinnedcube_blue2 = obs2.refs["level2"].product.refs["HPS3DRB"].product.refs[0].product
# BHR71 red rebinned cube
bhr71_rebinnedcube_red1 = obs1.refs["level2"].product.refs["HPS3DRR"].product.refs[0].product
bhr71_rebinnedcube_red2 = obs2.refs["level2"].product.refs["HPS3DRR"].product.refs[0].product

simpleFitsWriter(product=bhr71_rebinnedcube_blue1, file=outdir+'bhr71_rebinnedcube_blue1.fits')
simpleFitsWriter(product=bhr71_rebinnedcube_red1, file=outdir+'bhr71_rebinnedcube_red1.fits')
simpleFitsWriter(product=bhr71_rebinnedcube_blue2, file=outdir+'bhr71_rebinnedcube_blue2.fits')
simpleFitsWriter(product=bhr71_rebinnedcube_red2, file=outdir+'bhr71_rebinnedcube_red2.fits')
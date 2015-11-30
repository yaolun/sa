#
# Advanced example script to show how to display images in PlotXY.
#
# This script shows how to produce a "publication-ready" figure of a Herschel map.
# Run in HIPE 8 or newer.
#
# Author: Pasquale Panuzzo, CEA Saclay Irfu/SAp
# pasquale.panuzzo@cea.fr
#
# Version: 9 June 2011
#
# Import some useful classes
from java.awt import Color
from java.lang import Math
from herschel.share.fltdyn.math import SexagesimalFormatter
from herschel.share.fltdyn.math.SexagesimalFormatter import Mode
# Get a public SPIRE observation
myobs=getObservation(1342224922, useHsa=True) # 1342224922
# Extract the PSW (70 µm band) map
map=myobs.level2_5.refs["HPPUNIMAPB"].product
# Extract the image data
image=map.image
# mask the NaN value
mask = IS_FINITE(image)
q = image.where(~mask)
image[q] = 1e-50
# image=imageLog10(image=map.image)

# Create the layer with the image
layIma=LayerImage(image)
# Create a PlotXY object
myPlot=PlotXY()
# Add the image layer to the plot
myPlot.addLayer(layIma)
#### Coordinates ####
# The image will be plotted in a reference system with origin in the lower-left
# corner of the image pixel [0,0] and with pixel size of 1x1.
#
# It is possible to change the position of the image respect to the PlotXY axes
# and change the pixel size, so that the PlotXY axes will show the
# Right Ascension and the Declination coordinates.
#
# The LayerImage provides methods to set the position of the image and the
# pixel size with a system similar to the FITS WCS.
#
# Please note:
# 1) PlotXY CANNOT rotate the image, so if you need to plot a map that
# is not aligned with the North on RA and Dec axis, you will need to rotate it
# before plotting it. In this script we assume that the map is aligned with
# the North.
#
# 2) The PlotXY axis system assumes that coordinates are a linear
# transformation of pixel coordinates. So coordinates on the plotted axes
# are not fully correct in some projections. They are correct only for small
# angles near the projection reference.
# Extract the WCS of the map and put some WCS info into variables
wcs=map.wcs
crpix1=wcs.crpix1
crpix2=wcs.crpix2
crval1=wcs.crval1
crval2=wcs.crval2
cdelt1=wcs.cdelt1
cdelt2=wcs.cdelt2
naxis1=wcs.naxis1
naxis2=wcs.naxis2
# cos(Dec)
cosd=Math.cos(Math.toRadians(crval2))
# Set the origin and the scale of the axes so that they coincide with the WCS.
myPlot[0].xcdelt=cdelt1/cosd # note the cos(Dec)!!!
myPlot[0].ycdelt=cdelt2
myPlot[0].xcrpix=crpix1
myPlot[0].ycrpix=crpix2
myPlot[0].xcrval=crval1
myPlot[0].ycrval=crval2
# Change the axis type so that we have ticks in degrees/hours, min, sec
# and the RA growing toward the left.
myPlot.xaxis.type=Axis.RIGHT_ASCENSION
myPlot.yaxis.type=Axis.DECLINATION
myPlot.xaxis.titleText="Right Ascension (J2000)"
myPlot.yaxis.titleText="Declination (J2000)"
# Set the axes ranges so that the image fills completely the plotting area
xrange=[crval1-(crpix1-0.5)*cdelt1/cosd*0.2,\
 crval1-(crpix1-naxis1-0.5)*cdelt1/cosd*0.2]
myPlot.xrange=xrange
yrange=[crval2-(crpix2-0.5)*cdelt2*0.2,\
 crval2-(crpix2-naxis2-0.5)*cdelt2*0.2]
myPlot.yrange=yrange
# Adjust ticks to be nicer
myPlot.xaxis.tick.autoAdjustNumber=0
myPlot.xaxis.tick.number=5
myPlot.xaxis.tick.minorNumber=3
myPlot.xaxis.getAuxAxis(0).tick.autoAdjustNumber=0
myPlot.xaxis.getAuxAxis(0).tick.number=5
myPlot.xaxis.getAuxAxis(0).tick.minorNumber=3

myPlot.yaxis.tick.autoAdjustNumber=0
myPlot.yaxis.tick.number=5
myPlot.yaxis.tick.minorNumber=3
myPlot.yaxis.getAuxAxis(0).tick.autoAdjustNumber=0
myPlot.yaxis.getAuxAxis(0).tick.number=5
myPlot.yaxis.getAuxAxis(0).tick.minorNumber=3

# Change the size of the plotting area so that proportions are as on the sky
myPlot.setPlotSize(4.0,4.0*(naxis2*-cdelt2)/(naxis1*cdelt1))
#### Colours and intensity manipulation ####
# Set the colour table to have a grey image and the intensity table to have
# sources in black and empty sky in white
myPlot[0].colorTable="Ramp"
myPlot[0].intensityTable="Negative"
# Set the intensity range
highCut=0.2
lowCut=0.0001
myPlot[0].highCut=highCut
myPlot[0].lowCut=lowCut
#### Coordinate grid ####
# We can draw a coordinate grid on the image. PlotXY doesn't provide a built-in
# way to generate a coordinate grid, so we need to compute in the script the
# positions of a number of meridians and parallels and draw them as LayerXY.
# Parallels every 30", meridians every 5s
deltaDec=30.0/3600.0
deltaRa=1.0*15.0/3600.0*5.0
# Compute the nearest parallel and meridian to the projection center
decCenter=(Math.round(crval2/deltaDec))*deltaDec
raCenter=(Math.round(crval1/deltaRa))*deltaRa
# Estimate how many parallels and meridians shall be drawn on each side
ndec=Integer(Math.round((yrange[1]-yrange[0])/deltaDec)).intValue()
ndec=1+ndec/2
nra=Integer(Math.round((xrange[0]-xrange[1])/deltaRa)).intValue()
nra=1+nra/2
# Draw parallels
dd=10
nn=(2*nra)*dd+1
for i in range(2*ndec+1):
	# Coordinates of parallels in the sky coordinates
	raPara=raCenter+(Float1d.range(nn)-nra*dd)*deltaRa/dd
	decPara=Float1d(nn)+decCenter+(i-ndec)*deltaDec
	# Coordinates of parallels in the pixels coordinates
	xpixPara=Float1d(nn)
	ypixPara=Float1d(nn)
	for j in range(nn):
		ypixPara[j],xpixPara[j]=wcs.getPixelCoordinates(raPara[j],decPara[j])
	pass
	# Coordinates of parallels in the plot axes coordinates
	xplotPara=(xpixPara-crpix1+1.0)*cdelt1/cosd+crval1
	yplotPara=(ypixPara-crpix2+1.0)*cdelt2+crval2
	layPara=LayerXY(xplotPara,yplotPara,color=Color.green,stroke=0.5)
	myPlot.addLayer(layPara)
pass
# Draw meridians
nn=(2*ndec)*dd+1
for i in range(2*nra+1):
	raMeri=Float1d(nn)+raCenter+(i-nra)*deltaRa
	decMeri=decCenter+(Float1d.range(nn)-ndec*dd)*deltaDec/dd
	xpixMeri=Float1d(nn)
	ypixMeri=Float1d(nn)
	for j in range(nn):
		ypixMeri[j],xpixMeri[j]=wcs.getPixelCoordinates(raMeri[j],decMeri[j])
	pass
	xplotMeri=(xpixMeri-crpix1+1.0)*cdelt1/cosd+crval1
	yplotMeri=(ypixMeri-crpix2+1.0)*cdelt2+crval2
	layMeri=LayerXY(xplotMeri,yplotMeri,color=Color.green,stroke=0.5)
	myPlot.addLayer(layMeri)
pass
# We can note now that meridians and parallels don't cross the axes exactly
# at the ticks positions. That's because axes are linear with respect to pixels,
# while sky coordinates are not (with the exception of some projections).
#
# We want now put ticks to coincide with meridians and parallels. To do this
# we need to compute the plot coordinates where meridians and parallels cross
# the plot axes; we will impose these positions as tick locations and we will
# set the correct labels.
# Setting up formatters for Right Ascension and Declination
format_ra = SexagesimalFormatter(Mode.RA_HMS_LOWER)
format_ra.decimals = 0
format_dec = SexagesimalFormatter(Mode.DEC_DMS_SYMBOL)
format_dec.decimals = 0
# Compute again the location of parallels and find where they cross the Y axes
nn=(2*nra)*dd+1
ycrossPara0=Double1d(2*ndec+1,Float.NaN)
ycrossPara1=Double1d(2*ndec+1,Float.NaN)
decCross=String1d(2*ndec+1)
for i in range(2*ndec+1):
	raPara=raCenter+(Float1d.range(nn)-nra*dd)*deltaRa/dd
	decPara=Float1d(nn)+decCenter+(i-ndec)*deltaDec
	decCross[i]=format_dec.formatDegrees(decCenter+(i-ndec)*deltaDec)
	xpixPara=Float1d(nn)
	ypixPara=Float1d(nn)
	for j in range(nn):
		ypixPara[j],xpixPara[j]=wcs.getPixelCoordinates(raPara[j],decPara[j])
	pass
	xplotPara=(xpixPara-crpix1+1.0)*cdelt1/cosd+crval1
	yplotPara=(ypixPara-crpix2+1.0)*cdelt2+crval2
	xplotPara0=xplotPara-xrange[0]
	xplotPara1=xplotPara-xrange[1]
	for j in range(nn-1):
		if xplotPara0[j]*xplotPara0[j+1] <= 0.0:
			ycrossPara0[i]=(yplotPara[j]*xplotPara0[j+1]-yplotPara[j
			+1]*xplotPara0[j])/ \
			(xplotPara0[j+1]-xplotPara0[j])
		if xplotPara1[j]*xplotPara1[j+1] <= 0.0:
			ycrossPara1[i]=(yplotPara[j]*xplotPara1[j+1]-yplotPara[j
			+1]*xplotPara1[j])/ \
			(xplotPara1[j+1]-xplotPara1[j])
pass
iii0=ycrossPara0.where(IS_FINITE(ycrossPara0))
iii1=ycrossPara0.where(IS_FINITE(ycrossPara1))
myPlot.yaxis.tick.setFixedValues(ycrossPara0[iii0])
myPlot.yaxis.tick.label.fixedStrings=decCross[iii0].toArray()
myPlot.yaxis.getAuxAxis(0).tick.setFixedValues(ycrossPara1[iii1])
# Remove minor ticks
myPlot.yaxis.tick.minorNumber=0
myPlot.yaxis.getAuxAxis(0).tick.minorNumber=0
# Compute the location of meridians and find where they cross the X axes
nn=(2*ndec)*dd+1
xcrossMeri0=Double1d(2*nra+1,Double.NaN)
xcrossMeri1=Double1d(2*nra+1,Double.NaN)
raCross=String1d(2*nra+1)
for i in range(2*nra+1):
	raMeri=Float1d(nn)+raCenter+(i-nra)*deltaRa
	decMeri=decCenter+(Float1d.range(nn)-ndec*dd)*deltaDec/dd
	raCross[i]=format_ra.formatDegrees(raCenter+(i-nra)*deltaRa)
	xpixMeri=Float1d(nn)
	ypixMeri=Float1d(nn)
	for j in range(nn):
		ypixMeri[j],xpixMeri[j]=wcs.getPixelCoordinates(raMeri[j],decMeri[j])
	pass
	xplotMeri=(xpixMeri-crpix1+1.0)*cdelt1/cosd+crval1
	yplotMeri=(ypixMeri-crpix2+1.0)*cdelt2+crval2
	yplotMeri0=(yplotMeri-yrange[0])
	yplotMeri1=(yplotMeri-yrange[1])
	for j in range(nn-1):
		if yplotMeri0[j]*yplotMeri0[j+1] <= 0.0:
			xcrossMeri0[i]=(xplotMeri[j]*yplotMeri0[j+1]-xplotMeri[j
			+1]*yplotMeri0[j])/ \
			(yplotMeri0[j+1]-yplotMeri0[j])
		if yplotMeri1[j]*yplotMeri1[j+1] <= 0.0:
			xcrossMeri1[i]=(xplotMeri[j]*yplotMeri1[j+1]-xplotMeri[j
			+1]*yplotMeri1[j])/ \
			(yplotMeri1[j+1]-yplotMeri1[j])
pass
iii0=xcrossMeri0.where(IS_FINITE(xcrossMeri0))
iii1=xcrossMeri0.where(IS_FINITE(xcrossMeri1))
myPlot.xaxis.tick.setFixedValues(xcrossMeri0[iii0])
myPlot.xaxis.tick.label.fixedStrings=raCross[iii0].toArray()
myPlot.xaxis.getAuxAxis(0).tick.setFixedValues(xcrossMeri1[iii1])
myPlot.xaxis.tick.minorNumber=0
myPlot.xaxis.getAuxAxis(0).tick.minorNumber=0
#### Contours #####
# We want to draw level contours. We don't have (yet) a specialized layer for
# contours, so we will need to plot each contour segment.
# Generate the contours
# contours = automaticContour(image=map,levels=4,min=0.05,max=0.2,distribution=0)
# Plot the contours as
# keys=contours.keySet()
# for key in keys:
# 	if key.startswith("Contour"):
# 	cont=contours[key]
# 	keysc=cont.keySet()
# 	for keyc in keysc:
# 	x=(cont[keyc].data[:,1]-crpix1+1.0)*cdelt1/cosd+crval1
# 	y=(cont[keyc].data[:,0]-crpix2+1.0)*cdelt2+crval2
# 	myPlot.addLayer(LayerXY(x,y,color=Color.green))
# pass
#### Annotations #####
# Let's plot a circle around the observed source. This will be done with a
# classical LayerXY
# Radius of the circle (60")
radius=60.0/3600.
phase=Float1d.range(101)*2.0*Math.PI/100.0
# Position of the source as entered in HSPOT
raNom=map.meta["raNominal"].value
decNom=map.meta["decNominal"].value
# Coordinates of the circle
xx=radius*COS(phase)/cosd+raNom
yy=radius*SIN(phase)+decNom
# Add the circle
# myPlot.addLayer(LayerXY(xx,yy,color=Color.blue,stroke=2))
# And let's put an annotation with the name of the observed source
# ann=Annotation(raNom-0.015,decNom+0.015,map.meta["object"].value.upper())
# ann.fontSize=12
# ann.color=Color.blue
# myPlot.addAnnotation(ann)
# We want also to put a line to show the angle scale
xx2=Float1d([0,-17./3600/cosd])+myPlot.xrange[1]-0.02
yy2=Float1d(2)+myPlot.yrange[0]+0.01
myPlot.addLayer(LayerXY(xx2,yy2,color=Color.blue,stroke=2))
ann=Annotation(myPlot.xrange[1]-0.021,myPlot.yrange[0]+0.012,"17"+u"\u02BA")
ann.fontSize=12
ann.color=Color.blue
myPlot.addAnnotation(ann)
# Finally add the wavelength
# annText="%3.0f"%map.meta["wavelength"].value
# annText=annText+'$\mu m $'
# ann=Annotation(myPlot.xrange[0]+0.25,myPlot.yrange[0]+0.032,annText)
# ann.fontSize=12
# ann.color=Color.blue
# myPlot.addAnnotation(ann)
#### Colour bar ####
# We now want to create a colour bar on the right side of the plot.
# The colour bar is just another image in a subplot.
# Create an overlay layout so that we can plot the colour bar
layout=PlotOverlayLayout(marginRight=1.0)
myPlot.setLayout(layout)
# The colour bar is just a subplot showing an image,
# create the SubPlot where we put the colour bar
# the numbers here define the position of the bar respect to the main plot
spBar = SubPlot(SubPlotBoundsConstraints(0.0, 1.02, 0.0, -0.07))
# Here we construct the colour bar image data
lenBar=2561
barIma=Float2d(lenBar,1)
barIma[:,0]=(Float1d.range(lenBar)*(highCut-lowCut)/lenBar)+lowCut
# Construct a LayerImage with it
layBar = LayerImage(barIma)
layBar.colorTable=myPlot[0].colorTable
layBar.intensityTable=myPlot[0].intensityTable
# layBar.highCut=Math.log10(highCut)
# layBar.lowCut=Math.log10(lowCut)
# layBar.ycdelt=(Math.log10(highCut)-Math.log10(lowCut))/lenBar
# layBar.ycrval=Math.log10(lowCut)
layBar.highCut=highCut
layBar.lowCut=lowCut
layBar.ycdelt=(highCut-lowCut)/lenBar
layBar.ycrval=lowCut
layBar.ycrpix=+0.5
# layBar.yrange=[Math.log10(lowCut),Math.log10(highCut)]
layBar.yrange=[lowCut,highCut]
layBar.xrange=[0,1]
# Add the layer to the subplot
spBar.addLayer(layBar)
# Put the colour bar on the plot
myPlot.addSubPlot(spBar)
# Adjust the axes characteristics
# The X axis
xaxis=spBar.baseLayer.xaxis
xaxis.titleText=""
xaxis.tick.label.visible=0
xaxis.tick.autoAdjustNumber=0
xaxis.tick.minorNumber=0
xaxis.tick.height=0.0
xaxis.getAuxAxis(0).tick.height=0.0
# The Y axis
yaxis=spBar.baseLayer.yaxis
yaxis.titleText=""
yaxis.tick.label.visible=0
yaxis.tick.autoAdjustNumber=0
yaxis.tick.number=5
yaxis.tick.minorNumber=0
yaxis.tick.height=0.04
yaxis.getAuxAxis(0).tick.label.visible=1
yaxis.getAuxAxis(0).titleText="Flux (Jy beam$^{-1}$)"
yaxis.getAuxAxis(0).title.visible=1
yaxis.getAuxAxis(0).tick.height=0.04
yaxis.getAuxAxis(0).tick.number=5
yaxis.getAuxAxis(0).tick.minorNumber=0

##### Save ####
# Finally, save the figure to PDF
myPlot.saveAsPDF(map.meta["object"].value.upper()+".pdf")
# End of the example.
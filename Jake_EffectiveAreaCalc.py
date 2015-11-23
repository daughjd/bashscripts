import numpy, pylab, dashi, tables
from matplotlib.font_manager import fontManager, FontProperties
font=FontProperties(size='x-small')

NFiles={}
NFiles['IC86 Level3']=99.
NFiles['IC79 Level3']=99.

f={}

f['IC86 Level3']=tables.openFile('/data/uwa/jfeintzeig/IC86L3/8778/Level3-NEW/Merged8778.hd5')
f['IC79 Level3']=tables.openFile('/data/uwa/jfeintzeig/IC86L3/IC79_6308/MergedIC79L3.hd5')
#f=tables.openFile('/data/uwa/jfeintzeig/StartingTracks/MonteCarlo/7314/CutVariables/MergedNugen.hd5')

c={}
c[0]='blue'
c[30]='red'
c[60]='black'
c[90]='green'
c[120]='magenta'
c[150]='cyan'

Zenith={}
Energy={}
OneWeight={}
Weight={}
content={}

for lab in ['IC86 Level3','IC79 Level3']:
	content[lab]={}
	NEvents=NFiles[lab]*f[lab].root.I3MCWeightDict.cols.NEvents[0]

	if 'IC79' in lab:
		Zenith[lab]=f[lab].root.Primary.cols.zenith[:]
	if 'IC86' in lab:
		Zenith[lab]=f[lab].root.MCPrimaryNeutrino.cols.zenith[:]
	Energy[lab]=f[lab].root.I3MCWeightDict.cols.PrimaryNeutrinoEnergy[:]
	OneWeight[lab]=f[lab].root.I3MCWeightDict.cols.OneWeight[:]

	MinE=2
	#MaxE=int(max(numpy.log10(Energy[lab])))+1
	MaxE=7
	NBins=(MaxE-MinE)*4

	#for ZenPair in [[0,30],[30,60],[60,90],[90,120],[120,150],[150,180]]:
	for ZenPair in [[0,90],[90,180]]:
		# zenith range in degrees
		ZenMin=ZenPair[0]
		ZenMax=ZenPair[1]
		print ZenMin, ZenMax

		SolidAngle=2*numpy.pi*(numpy.cos(numpy.radians(ZenMin))-numpy.cos(numpy.radians(ZenMax)))
		Weight[lab]=OneWeight[lab]*10**(-4)*1/NEvents*1/SolidAngle
		ind=(Zenith[lab]>numpy.radians(ZenMin))&(Zenith[lab]<numpy.radians(ZenMax))

		PlotEnergy=Energy[lab][ind]
		PlotWeight=Weight[lab][ind]
		content[lab][ZenMin],edges=numpy.histogram(numpy.log10(PlotEnergy),weights=PlotWeight,range=[MinE,MaxE],bins=NBins)

		lin_edges=10**edges
		widths=lin_edges[1:]-lin_edges[:-1]
		content[lab][ZenMin]/=widths
		centers=(edges[:-1]+edges[1:])/2

		if lab=='IC79 Level3':
			ls='dashed'
		else:
			ls='solid'
		if ZenMin==0:
			label='Downgoing'
		if ZenMin==90:
			label='Upgoing'
		pylab.semilogy(centers,content[lab][ZenMin],label=lab+' '+label,linestyle=ls,color=c[ZenMin])

pylab.legend(loc='lower right',prop=font)
pylab.xlabel('log(Neutrino Energy (GeV))')
pylab.ylabel('Effective Area (m^2)')
pylab.savefig('EffArea.png')
pylab.show()
pylab.clf()

for ZenMin in [0,90]:
	ratio=content['IC86 Level3'][ZenMin]/content['IC79 Level3'][ZenMin]
	if ZenMin==0:
		label='Downgoing'
	if ZenMin==90:
		label='Upgoing'
	pylab.plot(centers,ratio,label='%s, Avg=%5.3f' % (label,numpy.mean(ratio[numpy.isfinite(ratio)])),color=c[ZenMin])

pylab.axhline(1.0,linestyle='dotted',color='black')
pylab.ylabel('IC86/IC79 Effective Area Ratio')
pylab.xlabel('log(Neutrino Energy (GeV))')
pylab.ylim(0,1.5)
pylab.legend(loc='lower right')
pylab.savefig('EffAreaRatio.png')
pylab.show()
pylab.clf()

for file in f:
	f[file].close()

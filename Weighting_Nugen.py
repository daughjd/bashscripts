import tables,pylab,numpy,pickle
import matplotlib
#from scipy.optimize import curve_fit
matplotlib.use("Agg")

pylab.rcParams['font.size'] = 14.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']='large'
pylab.rcParams['xtick.labelsize']='large'
pylab.rcParams['lines.markeredgewidth']=1.0
pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Computer Modern Roman')


dat = tables.openFile('paper_tables/FinalSample_nugen_numu_IC86.2013.010090.Table.hdf5')

energy = dat.root.I3MCWeightDict.col('PrimaryNeutrinoEnergy')
oweight = dat.root.I3MCWeightDict.col('OneWeight')
zenith = dat.root.PrimaryNu.col('zenith')
atmo = dat.root.AtmoWeight.col('value')

dat.close()

nevents_generated = 1000. * 200000.
e2_flux_normalization = 10**-8

pylab.figure()
pylab.hist(numpy.log10(energy),bins=60,histtype='step',weights=e2_flux_normalization * (oweight / nevents_generated) * energy**-2, ec='k',label=r'E$^{-2}$')
#pylab.hist(numpy.log10(energy),bins=60,histtype='step',weights=e2_flux_normalization * (oweight / nevents_generated) * energy**-3, ec='b',label=r'E$^{-3}$',log=True)
#pylab.hist(numpy.log10(energy),bins=60,histtype='step',weights=e2_flux_normalization * (oweight / nevents_generated) * energy**-4, ec='g',label=r'E$^{-4}$',log=True)
pylab.xlabel(r"Log$_{10}$[Energy (GeV)]")
pylab.ylabel("Neutrino Events")
pylab.title("MC Nugen Muon Neutrino Expectation")
pylab.legend()
pylab.savefig("MCNugenMuonExpectation")


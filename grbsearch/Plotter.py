import grbllh, numpy, grblist
import matplotlib
matplotlib.use("Agg")

import pylab
import pickle


pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']=18.0
pylab.rcParams['xtick.labelsize']=18.0
pylab.rcParams['lines.markeredgewidth']=1.0

#pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Times')


grbs = pickle.load(open("data_objects/grbarray.pkl",'r'))
null_tsd = pickle.load(open('tsd.pickle','r'))


pylab.figure()
pylab.hist(numpy.log10(grbs.arrays.t_100),bins=numpy.linspace(-1.5,2.6,20),histtype="step",ec='k')
pylab.xlabel(r"GRB Duration (Log$_{10}$T$_{100}$)")
pylab.axis([-1.5,3.,0,22])
pylab.ylabel("GRBs")
pylab.title("NH GRB Duration Distriubtion for IC86-2")
pylab.savefig("plots/GRB_DurationDistribution")


pylab.figure()
pylab.hist(null_tsd.T_vals,bins=numpy.linspace(0,11,60),histtype="step",ec='k',log=True)
pylab.xlabel(r"Test Statistic $T$")
#pylab.axis([])
pylab.ylabel(r"Log$_{10}$(N Trials)")
pylab.text(8,4e3,"Trials: 1e5",fontsize=16)
pylab.text(8,1e3,"Plotted: 23820",fontsize=16)
pylab.title("Background Trial Test Statistic Distribution")
pylab.savefig("plots/NullTestStatisticDistribution")


#######################################################################
# Generate GRB Ensembles representative of actual 2012-2013 NH Sample #
#######################################################################



import grbllh, grblist
import datetime,astrodate
import arrays
import pickle
import tables
from vars_class import Vars
from arrays import Arrays
from icecube import hdfdataset
import fakeps,time,os
import numpy as np
import pylab
import subprocess
from shortgrb import CDFShortGRB

from scipy.interpolate import UnivariateSpline

spline_zs = np.linspace(0,8,10000)

shortspline = UnivariateSpline(CDFShortGRB(spline_zs),spline_zs,s=1e-10)

real_grbs = pickle.load(open("data_objects/grbarray.pkl",'r'))              ### GRB List between May 15th 2012 -> April 30th 2013
realz=real_grbs.arrays.z
realz = realz[(realz!=2.15)*(realz!=0.5)]    ### Measured Redshifts


area_spline_dec0 = pickle.load(open("data_objects/nugenarea_spline_dec0.pkl","r"))
area_spline_dec16 = pickle.load(open("data_objects/nugenarea_spline_dec16.pkl","r"))
area_spline_dec30 = pickle.load(open("data_objects/nugenarea_spline_dec30.pkl","r"))
area_spline_dec45 = pickle.load(open("data_objects/nugenarea_spline_dec45.pkl","r"))
area_spline_dec60 = pickle.load(open("data_objects/nugenarea_spline_dec60.pkl","r"))
area_spline_dec75 = pickle.load(open("data_objects/nugenarea_spline_dec75.pkl","r"))

subphoto_fluxspline = pickle.load(open('data_objects/subphotofluxspline.pkl','r'))

energy_bins = np.linspace(20,1000,10000.)    ### WARNING: Area splines go negative below 20 GeV
dE = np.diff(energy_bins)[0]

### Long-duration GRB Redshift Distribution ###
inverse_cdf_spline = pickle.load(open("data_objects/InverseRedShiftCDF.pkl",'r'))		### Spline generated from paper detailing Swift GRBS
gen_inverse_cdf_spline = pickle.load(open("data_objects/GeneratedInverseRedShiftCDF.pkl",'r'))  ### Spline generated from L_GRB list direct from Swift website


def swift_z_draw(distro=0):
	if distro==0:
		return inverse_cdf_spline(np.random.uniform(0.0012,1))   ### Cutoff at z=0.025
	else:
		return gen_inverse_cdf_spline(np.random.uniform(0.0072,1))    ### Cutoff at z=0.025
def short_grbz_zdraw():
	return shortspline(np.random.uniform(0.00008))

def nexpected(grb):
        d_dec = np.rad2deg(grb.dec)
        effarea = area_spline_dec16   ### Default Effective Area
        if d_dec > -5 and d_dec <= 10:
                effarea = area_spline_dec0
        elif d_dec > 25 and d_dec <= 40:
                effarea = area_spline_dec30
        elif d_dec > 40 and d_dec <= 55:
                effarea = area_spline_dec45
        elif d_dec > 55 and d_dec <= 70:
                effarea = area_spline_dec60
        elif d_dec > 75:
                effarea = area_spline_dec75
        DL_0 = 475.34   ### Luminosity Distance for Reference Flux
        DL = float(subprocess.check_output(["./cosmology.py",str(grb.z)]).split()[2])
        energy_bins = np.linspace(20,1000,10000.)
        energy_bins_adj = energy_bins * 1.1 / (1+grb.z)
        energy_bins_adj = energy_bins_adj[energy_bins_adj > 20.0]
	shortmodifier = 1.0
	if grb.short:
		shortmodifier = 0.1
        events = shortmodifier * subphoto_fluxspline(energy_bins[energy_bins_adj > 20.0])*(DL/DL_0)**-2 * effarea(energy_bins_adj) * dE
	return events.sum()


	
class simgrb(object):
	def __init__(self,ra,dec,z,nexp,short_bool=False):
		self.ra = ra
		self.dec = dec
		self.z = z
		self.nexp = nexp
		self.short = short_bool
class grblist(object):
	def __init__(self,grblist=[]):
		self.grblist=grblist
	def dec(self):
		declinations=[]
		for i in self.grblist:
			declinations.append(i.dec)
		return np.array(declinations)
	def ra(self):
                ras=[]
                for i in self.grblist:
                        ras.append(i.ra)
                return np.array(ras)
	def z(self):
                zs=[]
                for i in self.grblist:
                        zs.append(i.z)
                return np.array(zs)
	def nexp(self):
		nexps=[]
		for i in self.grblist:
			nexps.append(i.nexp)
		return np.array(nexps)

def random_grb_set(isotropic=True,ngrbs=len(real_grbs)):
	simlist=[]
	for event in range(ngrbs):
		shortflag=False
		if isotropic:
			decdraw = np.arccos(np.random.uniform(-1,0.087)) - np.deg2rad(90.)
		else:
			dex = grbs.arrays.dec
			decdraw = dex[np.random.random_integers(0,len(dex))]
		if len(simlist) < 15:
			zdraw = realz[len(simlist)]
		elif len(simlist) < 34:
			zdraw = short_grbz_zdraw()
			shortflag = True
		else:
			zdraw = swift_z_draw()
		simmy  = simgrb(np.random.uniform(0,2*np.pi),decdraw,zdraw,1.0,shortflag)
		simmy_wevents = simgrb(simmy.ra,simmy.dec,simmy.z,nexpected(simmy),shortflag)
		simlist.append(simmy_wevents)
	return simlist

simulated_grbs = grblist(random_grb_set())

event_summation = []

for instance in range(0,1):
	simulated_grbs = grblist(random_grb_set())
	event_summation.append(simulated_grbs.nexp().sum())
	print "Sim Number: ",instance+1

event_summation = np.array(event_summation)

pickle.dump(event_summation,open("grb_event_expectation.pkl",'w'))

surpassing_sens_pct = 1.0*(event_summation > 3.07).sum() / len(event_summation) * 100.
surpassing_disco_pct = 1.0*(event_summation > 6.8).sum() / len(event_summation) * 100.
'''
pylab.figure()
hist=pylab.hist(event_summation,bins=np.linspace(0,7.0,18),histtype='step',ec='k',lw=2,normed=True)
pylab.vlines(3.07,0,10.,linestyles='dashed',label="Sensitivity",colors='b')
pylab.text(0.9*3.07,0.55*hist[0].max(),str(surpassing_sens_pct)+"% of Trials Greater",verticalalignment='center',rotation='vertical',fontsize=14)
pylab.vlines(6.8,0,10.,linestyles='dashed',label="Discovery",colors='g')
pylab.text(0.95*6.8,0.55*hist[0].max(),str(surpassing_disco_pct)+"% of Trials Greater",verticalalignment='center',rotation='vertical',fontsize=14)
pylab.axis([0,7,0,1.1*hist[0].max()])
pylab.xlabel("Total Events")
pylab.ylabel("Trials")
pylab.legend(loc="upper right")
pylab.title("Sub-Photo GRB Event Expectation Trials")
pylab.savefig("SubPhotoEventExpectations")
'''

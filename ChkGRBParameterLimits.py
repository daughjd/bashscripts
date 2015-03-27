import ROOT
#from root_numpy import root2array, root2rec, tree2rec
from ROOT import TFile

import pylab
import pickle

import numpy as np
import glob,matplotlib
matplotlib.use("Agg")

pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']=18.0
pylab.rcParams['xtick.labelsize']=18.0
pylab.rcParams['lines.markeredgewidth']=1.0

pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Computer Modern Roman')


low_en_bins = pickle.load(open("./pickles/effarea_low_energy_bins.pkl",'r'))
high_en_bins = pickle.load(open("./pickles/effarea_high_energy_bins.pkl",'r'))


## GENIE (4-190 GeV) & Nugen (190-1000) in Units of m^2##
genie_avg_area = pickle.load(open("./pickles/numu_effarea_avg.pkl",'r'))
genie_dec0_area = pickle.load(open("./pickles/numu_effarea_dec0.pkl",'r'))
genie_dec16_area = pickle.load(open("./pickles/numu_effarea_dec16.pkl",'r'))
genie_dec30_area = pickle.load(open("./pickles/numu_effarea_dec30.pkl",'r'))
genie_dec45_area = pickle.load(open("./pickles/numu_effarea_dec45.pkl",'r'))
genie_dec60_area = pickle.load(open("./pickles/numu_effarea_dec60.pkl",'r'))
genie_dec75_area = pickle.load(open("./pickles/numu_effarea_dec75.pkl",'r'))

nugen_avg_area = pickle.load(open("./pickles/nugmu_effarea_avg.pkl",'r'))
nugen_dec0_area = pickle.load(open("./pickles/nugmu_effarea_dec0.pkl",'r'))
nugen_dec16_area = pickle.load(open("./pickles/nugmu_effarea_dec16.pkl",'r'))
nugen_dec30_area = pickle.load(open("./pickles/nugmu_effarea_dec30.pkl",'r'))
nugen_dec45_area = pickle.load(open("./pickles/nugmu_effarea_dec45.pkl",'r'))
nugen_dec60_area = pickle.load(open("./pickles/nugmu_effarea_dec60.pkl",'r'))
nugen_dec75_area = pickle.load(open("./pickles/nugmu_effarea_dec75.pkl",'r'))


def FluxByEnergy(energy_array,ejet,Gamma,which_string):
	emaxp = Gamma*1.4*10000/4
	hadcoolbreak_pionic = pow(Gamma,5)*(3.7037037037037037*pow(10.0,50))/ejet ### pionic break dependency -> E_j^-1 Gamma^7 Theta_j_^2 (~Gamma^-2) t_j (10-100s) t_v^2 (0.1s)
	radcoolbreak_pionic = 33.333333*Gamma  ### (e_e + e_B)^-1 * Gamma 

	emaxk = Gamma*1.4*10000/2
	hadcoolbreak_kaonic = pow(Gamma,5)*(2.4691358024691357*pow(10.0,51))/ejet
	radcoolbreak_kaonic = 6666.66666*Gamma
	total_flux_array = []
	pionic_flux_array = []
	kaonic_flux_array = []
	if hadcoolbreak_pionic <= radcoolbreak_pionic:
        	ebreak1_pionic = hadcoolbreak_pionic
	        ebreak2_pionic = radcoolbreak_pionic
	if hadcoolbreak_pionic > radcoolbreak_pionic:
	        ebreak1_pionic = radcoolbreak_pionic
	        ebreak2_pionic = hadcoolbreak_pionic

	if hadcoolbreak_kaonic <= radcoolbreak_kaonic:
	        ebreak1_kaonic = hadcoolbreak_kaonic
	        ebreak2_kaonic = radcoolbreak_kaonic
	if hadcoolbreak_kaonic > radcoolbreak_kaonic:
	        ebreak1_kaonic = radcoolbreak_kaonic
	        ebreak2_kaonic = hadcoolbreak_kaonic
 
	## For Gamma=3,Ejet=3*10**51 @ 30 GeV
	nom_fluence_pionic = 0.05 * (ejet/(3.0*10**51)) * (Gamma/3.0)**2
	nom_fluence_kaonic = (5.0*10**-5) * (ejet/(3.0*10**51)) * (Gamma/3.0)**2
	for energy in energy_array:
		pionic_contribution=0.0
		kaonic_contribution=0.0
		# pionic
		if energy > emaxp:
			pionic_contribution = nom_fluence_pionic * 900.0 * ebreak1_pionic * ebreak2_pionic * (energy)**-8
		elif energy < ebreak1_pionic:
			pionic_contribution = nom_fluence_pionic * 900.0 * (energy)**-2
		elif energy < ebreak2_pionic:
			pionic_contribution = nom_fluence_pionic * 900.0 * ebreak1_pionic * (energy)**-3
		else:
			pionic_contribution = nom_fluence_pionic * 900.0 * ebreak1_pionic * ebreak2_pionic * (energy)**-4
			
		#kaonic
		if energy > emaxk:
			kaonic_contribution = nom_fluence_kaonic * 40000.0 * ebreak1_kaonic * ebreak2_kaonic * (energy)**-8
		elif energy < ebreak1_kaonic:
                        kaonic_contribution = nom_fluence_kaonic * 40000.0 * (energy)**-2
                elif energy < ebreak2_kaonic:
                        kaonic_contribution = nom_fluence_kaonic * 40000.0 * ebreak1_kaonic * (energy)**-3
                else:
                        kaonic_contribution = nom_fluence_kaonic * 40000.0 * ebreak1_kaonic * ebreak2_kaonic * (energy)**-4
		
		total_flux_array.append(pionic_contribution + kaonic_contribution)
		pionic_flux_array.append(pionic_contribution)
		kaonic_flux_array.append(kaonic_contribution)
	total_flux_array = np.array(total_flux_array)
	pionic_flux_array = np.array(pionic_flux_array)
	kaonic_flux_array = np.array(kaonic_flux_array)
	if which_string == "pionic":
		return pionic_flux_array
	if which_string == "kaonic":
		return kaonic_flux_array
	if which_string == None or which_string == "total":
		return total_flux_array
	if which_string == "breaks":
		return [ebreak1_pionic,ebreak2_pionic,ebreak1_kaonic,ebreak2_kaonic]
class ChkGRB_ParameterLimits(object):
	def __init__(self,gamma_b,ejet,dlist_ul):
		self.gamma_boost = gamma_b
		self.ejet = ejet
		self.declist = [0,16,30,45,60,75]
		self.ul = {}
		for i,dec in enumerate(self.declist):
			self.ul[dec] = dlist_ul[i]
	def get_visible_dist(self):
		return 1

class wrapper(object):
	def __init__(self, gamma_b_list=None, ejet_list=None, dlist_ul_list=None):
		if not gamma_b_list == None:
			self.gamma_b_list = gamma_b_list.copy()
			self.ejet_list = ejet_list.copy()
			self.dlist_ul_list = dlist_ul_list.copy()
		else:
			self.gamma_b_list = np.array([])
			self.ejet_list = np.array([])
			self.dlist_ul_list = None

	def get_by_ejet(self,ejet_value):
		tmplist = []
		for i,ejet in self.ejet_list:
			if ejet == ejet_value:
				tmplist.append(i)
		return wrapper(self.gamma_b_list[tmplist], self.ejet_list[tmplist], self.dlist_ul_list[tmplist])
	def get_limit_by_ejet_gamma(self,ejet_value,gamma_value):
		tmp_ejet = []
		tmp_gamma = []
		for i,ejet in enumerate(self.ejet_list):
			if ejet == ejet_value:
				tmp_ejet.append(i)
		for i,gamma in enumerate(self.gamma_b_list):
			if gamma == gamma_value:
				tmp_gamma.append(i)
		n = np.intersect1d(tmp_ejet,tmp_gamma)
                if len(n) < 1:
			return 0
		else:
			n=n[0]
			return ChkGRB_ParameterLimits(self.gamma_b_list[n], self.ejet_list[n], self.dlist_ul_list[n]).ul[16]
		
	def get_entry(self,n):
		return ChkGRB_ParameterLimits(self.gamma_b_list[n], self.ejet_list[n], self.dlist_ul_list[n])

	def add_object(self,gamma_b,ejet,dlist_ul):
		self.gamma_b_list = np.append(self.gamma_b_list,gamma_b)
		self.ejet_list = np.append(self.ejet_list,ejet)
		if self.dlist_ul_list == None:
			self.dlist_ul_list = dlist_ul.copy()
		else:
			self.dlist_ul_list = np.concatenate((self.dlist_ul_list, dlist_ul),axis=0)

	def plot_values(self):
		pass

#### Event Upper Limits (10^-3 days [86.4 s] Width) [0,16,30,45,60,75] Declination####
sens_runs = glob.glob("./root_tables/*.root")

gamma_train = np.linspace(0.5,10,20)
#oldjetty = ["1.00e+50","2.00e+50","3.00e+50","4.00e+50","5.00e+50","6.00e+50","7.00e+50","8.00e+50","9.00e+50","1.00e+51","2.00e+51","3.00e+51","4.00e+51","5.00e+51","6.00e+51","7.00e+51","8.00e+51","9.00e+51","1.00e+52",]
jetty = ["3.16e+50","4.03e+50","5.13e+50","6.54e+50","8.34e+50","1.06e+51","1.35e+51","1.73e+51","2.20e+51","2.80e+51","3.57e+51","4.55e+51","5.80e+51","7.39e+51","9.41e+51","1.20e+52","1.53e+52","1.95e+52","2.48e+52","3.16e+52"]

x=ROOT.Double(0)
y=ROOT.Double(0)
upper_limits = wrapper()

for jet in jetty:
	for g_unit in gamma_train:
		globule = glob.glob("./root_tables/IC86II_LowEnNorthernSkyTimeDep_ChkGRBSens_Gamma_"+str(g_unit)+"_Ejet_"+jet+"_Dec16.root")
		if len(globule) > 2:
			print "We got multiples!"
		if len(globule) < 1:
			pass
		else:
			fd = TFile(globule[0])
			tc=fd.Get("canComp")
			tg=tc.GetPrimitive("Graph")
			tg.GetPoint(0,x,y)
			upper_limits.add_object(g_unit,float(jet),np.array([[0,y,0,0,0,0]]))
		

def VisibilityDistance(ejet,gamma,numu_frac=1/3.):
	event_ul = upper_limits.get_limit_by_ejet_gamma(ejet,gamma)
	expected_lowen_events = (FluxByEnergy(low_en_bins,ejet,gamma,"total")*genie_avg_area*10000.0*(np.diff(low_en_bins)[0])).sum()
	expected_highen_events = (FluxByEnergy(high_en_bins,ejet,gamma,"total")*nugen_avg_area*10000.0*(np.diff(high_en_bins)[0])).sum()
	tot_numu_expected_events = numu_frac*expected_lowen_events + expected_highen_events ## Events expected at 10 Mpc
	visible_dist = 10*(event_ul / tot_numu_expected_events)**-0.5
	return visible_dist

def VolumetricRateLimit(ejet,gamma,frak=1/3.):
	vis_distance = VisibilityDistance(ejet,gamma,frak)
	search_SA = 2*np.pi*(1-np.cos(np.deg2rad(95.0)))
	v90 = (1/3.0)*search_SA*vis_distance**3
	return 2.4 / (v90 * 28475585.0/31536000.0)


x=[]
y=[]
z=[]

range_rover = np.arange(0,len(upper_limits.ejet_list))

for indexical in range_rover:
	param_set = upper_limits.get_entry(indexical)
	x.append(param_set.ejet)
	y.append(param_set.gamma_boost)
	z.append(param_set.ul[16])

ux = np.unique(x)
uy = np.unique(y)

ulgrid = []
visgrid = []
rategrid = []
for y0 in uy:
	line = []
	visline = []
	rateline = []
	for x0 in ux:
		if upper_limits.get_limit_by_ejet_gamma(x0,y0) > 0:
			line.append(upper_limits.get_limit_by_ejet_gamma(x0,y0))
			visline.append(VisibilityDistance(x0,y0))
			rateline.append(VolumetricRateLimit(x0,y0))
		else:
			line.append(4.6)
			visline.append(0.1)
			rateline.append(40)
	ulgrid.append(line)
	visgrid.append(visline)
	rategrid.append(rateline)

visgrid = np.array(visgrid)
ulgrid = np.array(ulgrid)
ux = np.log10(ux)

xedges=np.append(min(ux)-np.diff(ux)[0]/2,ux[:-1] + np.diff(ux)/2.)
xedges=np.append(xedges,xedges[-1]+np.diff(ux)[-1])
yedges=np.append(min(uy)-np.diff(uy)[0]/2,uy[:-1] + np.diff(uy)/2.)
yedges=np.append(yedges,yedges[-1]+np.diff(uy)[-1])



X,Y = np.meshgrid(xedges,yedges)
pylab.figure()
pylab.pcolormesh(X,Y,ulgrid)
pylab.axis([xedges.min(),xedges.max(),yedges.min(),yedges.max()])
pylab.title("Choked GRB Event Upper Limit (90\% C.L.)")
pylab.xlabel(r"$\log(E_{jet})$ erg")
pylab.ylabel(r"$\Gamma_{b}$")
pylab.colorbar(label="Events Required")
pylab.savefig("UpperLimit_2DHisto")

pylab.figure()
pylab.pcolormesh(X,Y,np.log10(visgrid))
pylab.axis([xedges.min(),xedges.max(),yedges.min(),yedges.max()])
pylab.title("Choked GRB Visible Distance Limit (90\% C.L.)")
pylab.xlabel(r"$\log(E_{jet})$ erg")
pylab.ylabel(r"$\Gamma_{b}$")
pylab.colorbar(label="Distance (Log10[Mpc])")
pylab.savefig("LogDistanceLimit_2DHisto")

levels=[1.0,5.0,10.0,20.0,50.0,100.0]
strs=['1 Mpc','5 Mpc','10 Mpc','20 Mpc','50 Mpc','100 Mpc']
fmt={}

for l,s in zip(levels,strs):
	fmt[l] = s


pylab.rcParams['lines.markersize']=16.0
pylab.figure()
pylab.pcolormesh(X,Y,visgrid)
pylab.axis([xedges.min(),xedges.max(),yedges.min(),yedges.max()])
pylab.plot(51.5,3.0,'y*')
pylab.title("Choked GRB Visible Distance Limit (90\% C.L.)")
pylab.xlabel(r"$\log(E_{jet})$ erg")
pylab.ylabel(r"$\Gamma_{b}$")
pylab.colorbar(label="Distance (Mpc)")
tours = pylab.contour(visgrid,levels,linewidths=2.0*np.ones(6),colors=['w','w','w','w','k','k'],extent=[xedges.min(),xedges.max(),yedges.min(),yedges.max()])
clabel_locs = [(52,1.5), (52,2.6), (52,3.3), (52,5), (51.75,7),(52.1,8)]
pylab.clabel(tours,inline=1,fmt=fmt,fontsize=16,manual=clabel_locs)
pylab.savefig("DistanceLimit_2DHisto_WithContours")

pylab.figure()
pylab.pcolormesh(X,Y,np.log10(rategrid))
pylab.axis([xedges.min(),xedges.max(),yedges.min(),yedges.max()])
pylab.plot(51.5,3.0,'y*')
pylab.title("Choked GRB Volumetric Rate Limit (90\% C.L.)")
pylab.xlabel(r"$\log(E_{jet})$ erg")
pylab.ylabel(r"$\Gamma_{b}$")
pylab.colorbar(label=r"Choked GRB Rate (Mpc$^{-3} \cdot$ yr$^{-1}$)")
#pylab.show()
#tours = pylab.contour(visgrid,levels,linewidths=2.0*np.ones(6),colors=['w','w','w','w','k','k'],extent=[xedges.min(),xedges.max(),yedges.min(),yedges.max()])
#clabel_locs = [(52,1.5), (52,2.6), (52,3.3), (52,5), (51.75,7),(52.1,8)]
#pylab.clabel(tours,inline=1,fmt=fmt,fontsize=16,manual=clabel_locs)
pylab.savefig("RateLimit_2DHisto")

levels=[0,-1,-2,-3,-3.321,-4.1549,-5]
styles=['solid','solid','solid','solid','dashed','dashdot','solid']
pylab.figure()
pylab.pcolormesh(X[2:],Y[2:],np.log10(rategrid[2:]))
pylab.axis([xedges.min(),xedges.max(),yedges[2:].min(),yedges.max()])
pylab.plot(51.5,3.0,'y*')
pylab.title("Choked GRB Volumetric Rate Limit (90\% C.L.)")
pylab.xlabel(r"$\log(E_{jet})$ erg")
pylab.ylabel(r"$\Gamma_{b}$")
pylab.colorbar(label=r"Choked GRB Rate (Mpc$^{-3} \cdot$ yr$^{-1}$)")
tours = pylab.contour(np.log10(rategrid),levels,linewidths=2.0*np.ones(6),colors=['k','k','k','k','k','k','w'],extent=[xedges.min(),xedges.max(),yedges.min(),yedges.max()],linestyles=styles)
clabel_locs = [(52,1.5), (52,2.4), (52,2.9), (52,3.6),(51.3,5),(51.75,6),(52.1,8)]
pylab.clabel(tours,inline=1,fontsize=16,manual=clabel_locs,fmt='%1.1f')
pylab.savefig("RateLimit_2DHisto_wContours")






### Flux ###
energy_array = np.linspace(0,100000,10000)
fluxy_p = FluxByEnergy(energy_array,3.0*10**51,3.0,"pionic")
fluxy_k = FluxByEnergy(energy_array,3.0*10**51,3.0,"kaonic")
total_fluxy = FluxByEnergy(energy_array,3.0*10**51,3.0,"total")
e_breaks = FluxByEnergy(energy_array,3.0*10**51,3.0,"breaks")

pylab.figure()
pylab.loglog(energy_array,fluxy_p,'k--',lw=2,label="Pionic")
pylab.loglog(energy_array,fluxy_k,'k-.',lw=2,label="Kaonic")
pylab.loglog(energy_array,total_fluxy,'k-',lw=2,label="Total")
pylab.axis([0,10**4,1e-9,1e2])
pylab.legend(loc="upper right")
pylab.xlabel(r"$Log_{10}(E_{\nu})$ (GeV)")
pylab.ylabel(r"Flux Log$_{10}(\frac{dN}{dE})$ GeV$^{-1}$ cm$^{-2}$")
pylab.title("Choked GRB Neutrino Flux")
pylab.savefig("FluxPlot_Canonical_RespectiveFluenceNorm")

#pylab.show()

pylab.figure()
pylab.loglog(energy_array,energy_array**2*fluxy_p,'k--',lw=2,label="Pionic")
pylab.loglog(energy_array,energy_array**2*fluxy_k,'k-.',lw=2,label="Kaonic")
pylab.loglog(energy_array,energy_array**2*total_fluxy,'k-',lw=2,label="Total")
pylab.vlines(e_breaks[0],1e-3,1e2,colors='b',linestyle='dashed',lw=2,label=r"$E_{\pi}^{1}")
pylab.vlines(e_breaks[1],1e-3,1e2,colors='b',linestyle='dashed',lw=2,label=r"$E_{\pi}^{2}")
pylab.vlines(e_breaks[2],1e-3,1e2,colors='g',linestyle='dashed',lw=2,label=r"$E_{\/K}^{1}")
#pylab.vlines(e_breaks[3],1e-3,1e2,colors='g',linestyle='dashed',lw=2,label=r"$E_{\K}^{2}")
pylab.axis([0,10**4,1e-2,1e2])
pylab.title("Choked GRB Flux w/ Break Energies")
pylab.legend(loc="upper right")
pylab.xlabel("Log10(Energy) (GeV)")
pylab.ylabel("Flux Log10(E^2 dN/dE) GeV^-1 cm^-2")
pylab.savefig("FluxPlot_Canonical_E2_RespectiveFluenceNorm")

pylab.figure()
pylab.loglog(energy_array,energy_array**2*total_fluxy,'k-',lw=2,label=r"\Gamma = 3")
pylab.loglog(energy_array,energy_array**2*FluxByEnergy(energy_array,3.0*10**51,1.5,"total"),'g-',lw=2,label=r"\Gamma = 1.5")
pylab.loglog(energy_array,energy_array**2*FluxByEnergy(energy_array,3.0*10**51,2.0,"total"),'b-',lw=2,label=r"\Gamma = 2")
pylab.loglog(energy_array,energy_array**2*FluxByEnergy(energy_array,3.0*10**51,5.0,"total"),'r-',lw=2,label=r"\Gamma = 5")
pylab.loglog(energy_array,energy_array**2*FluxByEnergy(energy_array,3.0*10**51,7.0,"total"),'m-',lw=2,label=r"\Gamma = 7")
pylab.loglog(energy_array,energy_array**2*FluxByEnergy(energy_array,3.0*10**51,10.0,"total"),'c-',lw=2,label=r"\Gamma = 10")
pylab.axis([10,10**5,1e-3,1000.])
pylab.xlabel("Log10(Energy) (GeV)")
pylab.ylabel("Flux Log10(E^2 dN/dE) GeV^-1 cm^-2")
pylab.savefig("FluxPlot_GammaComparison_E2")


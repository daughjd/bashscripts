import ROOT
#from root_numpy import root2array, root2rec, tree2rec
import pylab,numpy,pickle
import matplotlib

pylab.rcParams['font.size'] = 14.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']='large'
pylab.rcParams['xtick.labelsize']='large'
pylab.rcParams['lines.markeredgewidth']=1.0
pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Computer Modern Roman')


log_sigma_days = numpy.array([-5,-4,-3,-2,-1,-0.52287874528033762,0,1])

### NEW GENIE 1460 Included ###
dec0_e3_foldedspectrum = (1072.916206382002,0)
dec16_e3_foldedspectrum = (1545.0315486757047,0)
dec30_e3_foldedspectrum = (1803.4879220886971,0)
dec45_e3_foldedspectrum = (1955.9670994116407,0)
dec60_e3_foldedspectrum = (2117.1599069802728,0)
dec75_e3_foldedspectrum = (2228.3197855702933,0)
sa_avg_foldedspectrum = (1654.0807981564465,0)

sys_adjustment = 0.89559693491089454

### Int(EffaE-3)  (JF,RH)###
#samp2_e3_foldedspectrum_sum = (1759.219287256351,0) ## 100 GeV flux equal to 1.0 GeV^-1 cm^-2 s^-1
#samp2_e35_foldedspectrum_sum = (2925.5560058208703,0) ## 
#samp2_e25_foldedspectrum_sum = (1320.5883336274608,0) ## 



sens_e3_dec0_meansrc_events = numpy.array([6.4656,6.70643,6.7344,7.38432,10.4106,13.2816,16.2928,28.1549])
sens_e3_dec16_meansrc_events = numpy.array([6.4384,6.62176,6.79315,7.4096,10.5558,13.0896,16.5709,30.3184])
sens_e3_dec30_meansrc_events = numpy.array([7.632,7.32,7.54048,8.00864,10.68,12.6272,16.0406,27.1056])
sens_e3_dec45_meansrc_events = numpy.array([6.86976,6.87104,7.09792,8.60768,11.3456,12.983,16.1408,27.0288])
sens_e3_dec60_meansrc_events = numpy.array([6.77216,6.54144,7.29088,8.584,11.0262,13.2019,15.5658,24.368])
sens_e3_dec75_meansrc_events = numpy.array([5.6608,5.64512,5.95296,7.37824,10.8947,12.7984,15.9766,28.8221])


ul_e3_dec0_meansrc_events = numpy.array([7.5456,8.09952,9.06432,11.376,17.5674,22.2304,29.9581,60.232])
ul_e3_dec16_meansrc_events = numpy.array([7.77754,8.51104,9.67872,11.8336,18.1984,23.208,30.528,64.568])
ul_e3_dec30_meansrc_events = numpy.array([8.95392,9.34349,10.2138,12.5501,18.1462,22.568,29.6342,59.744])
ul_e3_dec45_meansrc_events = numpy.array([8.45888,8.73325,9.74496,12.8112,19.0477,22.5107,29.5024,59.3357])
ul_e3_dec60_meansrc_events = numpy.array([8.17261,8.74912,10.1846,13.3968,19.3747,23.0784,30.0032,57.7504])
ul_e3_dec75_meansrc_events = numpy.array([7.30272,7.66144,8.52512,11.688,19.0272,24.0032,31.9216,64.608])

ilow_en_bins = pickle.load(open("./pickles/effarea_low_energy_bins.pkl",'r'))
high_en_bins = pickle.load(open("./pickles/effarea_high_energy_bins.pkl",'r'))


genie_avg_area = pickle.load(open("./pickles/g1460_numu_effarea_avg.pkl",'r'))
genie_dec0_area = pickle.load(open("./pickles/g1460_numu_effarea_dec0.pkl",'r'))
genie_dec16_area = pickle.load(open("./pickles/g1460_numu_effarea_dec16.pkl",'r'))
genie_dec30_area = pickle.load(open("./pickles/g1460_numu_effarea_dec30.pkl",'r'))
genie_dec45_area = pickle.load(open("./pickles/g1460_numu_effarea_dec45.pkl",'r'))
genie_dec60_area = pickle.load(open("./pickles/g1460_numu_effarea_dec60.pkl",'r'))
genie_dec75_area = pickle.load(open("./pickles/g1460_numu_effarea_dec75.pkl",'r'))

nugen_avg_area = pickle.load(open("./pickles/g1460_nugmu_effarea_avg.pkl",'r'))
nugen_dec0_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec0.pkl",'r'))
nugen_dec16_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec16.pkl",'r'))
nugen_dec30_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec30.pkl",'r'))
nugen_dec45_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec45.pkl",'r'))
nugen_dec60_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec60.pkl",'r'))
nugen_dec75_area = pickle.load(open("./pickles/g1460_nugmu_effarea_dec75.pkl",'r'))

sa0 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(95.))) - (1-numpy.cos(numpy.deg2rad(80.))))
sa16 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(80.))) - (1-numpy.cos(numpy.deg2rad(65.))))
sa30 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(65.))) - (1-numpy.cos(numpy.deg2rad(50.))))
sa45 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(50.))) - (1-numpy.cos(numpy.deg2rad(35.))))
sa60 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(35.))) - (1-numpy.cos(numpy.deg2rad(20.))))
sa75 = 2*numpy.pi*(1-numpy.cos(numpy.deg2rad(20.)))
saTotal = 2*numpy.pi*(1-numpy.cos(numpy.deg2rad(95.)))

sky_frac = [0.23989563791056959, 0.22901050354066707, 0.20251868181221927, 0.16222554659621455, 0.11087700847006936, 0.055472621670260208]



fluxnorm_dec16_e3 = ul_e3_dec16_meansrc_events/dec16_e3_foldedspectrum[0]
fluxnorm_dec0_e3 = ul_e3_dec0_meansrc_events/dec0_e3_foldedspectrum[0]
fluxnorm_dec30_e3 = ul_e3_dec30_meansrc_events/dec30_e3_foldedspectrum[0]
fluxnorm_dec45_e3 = ul_e3_dec45_meansrc_events/dec45_e3_foldedspectrum[0]
fluxnorm_dec60_e3 = ul_e3_dec60_meansrc_events/dec60_e3_foldedspectrum[0]
fluxnorm_dec75_e3 = ul_e3_dec75_meansrc_events/dec75_e3_foldedspectrum[0]

uls = [ul_e3_dec0_meansrc_events,ul_e3_dec16_meansrc_events,ul_e3_dec30_meansrc_events,ul_e3_dec45_meansrc_events,ul_e3_dec60_meansrc_events,ul_e3_dec75_meansrc_events]
event_ul_avg_list = [uls[i]*sky_frac[i] for i in range(len(sky_frac))]

event_ul_avg = numpy.array([0.,0.,0.,0.,0.,0.,0.,0.])
for listy in event_ul_avg_list:
        event_ul_avg+=listy


fluxnorm_sa_avg_e3 = event_ul_avg / sa_avg_foldedspectrum[0]

#fluxnorm_0 = sens_bdt0_e3_meansrc_events/samp2_e3_foldedspectrum_sum[0]
#fluxnorm_0_disco = disco_bdt0_e3_meansrc_events/samp2_e3_foldedspectrum_sum[0]

#fluxnorm_0_25 = sens_bdt0_e25_meansrc_events/samp2_e25_foldedspectrum_sum[0]
#fluxnorm_0_35 = sens_bdt0_e35_meansrc_events/samp2_e35_foldedspectrum_sum[0]



pylab.figure()
pylab.plot(log_sigma_days,event_ul_avg,'k-',lw=2,label="Averaged")
pylab.plot(log_sigma_days,ul_e3_dec0_meansrc_events,'k--',lw=2,label=r"$\delta=0^{\circ}$")
pylab.plot(log_sigma_days,ul_e3_dec30_meansrc_events,'k-.',lw=2,label=r"$\delta=30^{\circ}$")
pylab.plot(log_sigma_days,ul_e3_dec60_meansrc_events,'k:',lw=2,label=r"$\delta=60^{\circ}$")
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad (Days)$')
pylab.ylabel("NSrc Events")
pylab.axis([-5,1,3,60])
pylab.grid()
pylab.legend(loc="upper left")
matplotlib.pyplot.gcf().subplots_adjust(right=.85)
pylab.title(r"Upper Limit $E^{-3}$ 90% C.L.")
pylab.savefig("LowEnTransient_NEventUpperLimit_E3_G1460_MultiDec")

fig1=pylab.figure()
pylab.plot(log_sigma_days,event_ul_avg,'k-',lw=2)
#pylab.plot(0.77011529478710161,13.5279,"w*",ms=20.0,label="Most Significant Flare")
#pylab.plot(log_sigma_days,disco_bdt0_e3_meansrc_events,'k-',lw=2,label="Discovery Potential")
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad$ (Days)')
pylab.ylabel("NSrc Events")
pylab.axis([-5,1,0,62])
pylab.grid()
pylab.legend(loc="upper left")
matplotlib.pyplot.gcf().subplots_adjust(right=.85)
pylab.title(r"Upper Limit $E^{-3}$ 90$\%$ C.L.")
pylab.savefig("LowEnTransient_NEventUpperLimit_E3_G1460_Avg.pdf")

figgy=pylab.figure()
ax = figgy.add_subplot(111)
pylab.plot(log_sigma_days,fluxnorm_sa_avg_e3,'k-',lw=2,label=r"$E^{-3.0}$")
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad$ (Days)')
pylab.ylabel(r"$\frac{dN}{dE}$ @ 100 GeV ($10^{-2}$GeV$^{-1}$ cm$^{-2}$)")
pylab.axis([-5,1,0.00,0.037483054073961818])
pylab.yticks([0.0060456538828970677,0.012091307765794135,0.018136961648691202,0.024182615531588271,0.030228269414485337,0.036273923297382403],["0.6","1.21","1.81","2.42","3.02","3.63"])
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
matplotlib.pyplot.gcf().subplots_adjust(right=.85)
pylab.grid()
#pylab.legend(loc="upper left")
pylab.title(r"Time-Integrated Flux Upper Limit $E^{-3}$")
pylab.savefig("LowEnTransient_FluxUpperLimit_E3_G1460_Avg.pdf")

figgy=pylab.figure()
ax = figgy.add_subplot(111)
pylab.plot(log_sigma_days,event_ul_avg,'k-',lw=2,label=r"$E^{-3.0}$")
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad$ (Days)')
pylab.ylabel("NSrc Events")
pylab.axis([-5,1,0.00,62])
pylab.yticks([  0.,  10.,  20.,  30.,  40.,  50.,  60.])
pylab.grid()
ax2 = ax.twinx()
ax2.set_ylim(0,0.037483054073961818)
ax2.set_xlim(-5,1)
ax2.set_yticks([0.0060456538828970677,0.012091307765794135,0.018136961648691202,0.024182615531588271,0.030228269414485337,0.036273923297382403])
ax2.set_yticklabels(["0.6","1.21","1.81","2.42","3.02","3.63"])
ax2.set_ylabel(r"$\frac{dN}{dE}$ @ 100 GeV ($10^{-2}$GeV$^{-1}$ cm$^{-2}$)")
matplotlib.pyplot.gcf().subplots_adjust(right=.85)
#pylab.legend(loc="upper left")
pylab.title(r"Time-Integrated Flux Upper Limit $E^{-3}$")
pylab.savefig("LowEnTransient_FluxUpperLimit_E3_G1460_Avg_DoubleY.pdf")



figgy=pylab.figure()
ax = figgy.add_subplot(111)
pylab.plot(log_sigma_days,fluxnorm_dec0_e3,'k--',lw=2,label=r"$\delta = 0^{\circ}$")
pylab.plot(log_sigma_days,fluxnorm_dec16_e3,'k-',lw=2,label=r"$\delta = 16^{\circ}$")
pylab.plot(log_sigma_days,fluxnorm_dec30_e3,'k-.',lw=2,label=r"$\delta = 30^{\circ}$")
pylab.plot(log_sigma_days,fluxnorm_dec60_e3,'k:',lw=2,label=r"$\delta = 60^{\circ}$")
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad$ (Days)')
pylab.ylabel(r"$\frac{dN}{dE}$ @ 100 GeV ($10^{-2}$GeV$^{-1}$ cm$^{-2}$)")
pylab.axis([-5,1,0.00,0.058])
pylab.yticks([0.00 ,  0.00828571,  0.01657143,  0.02485714,  0.03314286, 0.04142857, 0.04971429, 0.058],["0.0","0.83","1.7","2.5","3.3","4.1","5.0","5.8"])
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
matplotlib.pyplot.gcf().subplots_adjust(right=.85)
pylab.grid()
pylab.legend(loc="upper left")
pylab.title(r"Time-Integrated Flux Upper Limit $E^{-3}$")
pylab.savefig("LowEnTransient_FluxUpperLimit_E3_G1460_MultiDec.pdf")


'''
pylab.figure(figsize=(10,8))
pylab.plot(log_sigma_days,fluxnorm_0,'b-',lw=2,label='Sensitivity (90% C.L.)')
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad (Days)$')
pylab.ylabel(r"$\frac{dN}{dE}$ [$GeV^{-1} cm^{-2} s^{-1}$] @ 100 GeV Pivot Energy")
#pylab.axis([-5,1,5e3,5e4])
pylab.yticks([0.001,0.005,0.01,0.015,0.02,0.025],["$1e-3$","$5.0e-3$","$1.0e-2$","1.5e-2","2.0e-2","2.5e-2"])
pylab.grid()
pylab.legend(loc="upper left")
pylab.title(r"Flux Sensitivity (MergedSim) $E^{-3}$")
pylab.savefig("LowEnTransient_FluenceSensitivity_E3_MergedSim_FinalCut")

pylab.figure(figsize=(10,8)) 
pylab.plot(log_sigma_days,fluxnorm_0,'b-',lw=2,label='Sensitivity (90% C.L.)')
pylab.plot(log_sigma_days,fluxnorm_0_disco,'k-',lw=2,label='Discovery Potential')
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad (Days)$') 
pylab.ylabel(r"$\frac{dN}{dE}$ [$GeV^{-1} cm^{-2} s^{-1}$] @ 100 GeV Pivot Energy") 
#pylab.axis([-5,1,5e3,5e4])
pylab.yticks([0.001,0.005,0.01,0.015,0.02,0.025],["$1e-3$","$5.0e-3$","$1.0e-2$","1.5e-2","2.0e-2","2.5e-2"])
pylab.grid() 
pylab.legend(loc="upper left")
pylab.title(r"Flux Sensitivity (MergedSim) $E^{-3}$") 
pylab.savefig("LowEnTransient_FluenceSensitivityAndDisco_E3_MergedSim_FinalCut") 


pylab.figure(figsize=(10,8)) 
pylab.plot(log_sigma_days,merged_samp1_e2_meansrc_events,'g-',lw=2,label='Sample 1') 
pylab.plot(log_sigma_days,merged_samp2_e2_meansrc_events,'b-',lw=2,label='Sample 2')  
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad (Days)$') 
pylab.ylabel("NSrc Events") 
pylab.axis([-6,1,3,15])
pylab.grid() 
pylab.title("Sensitivity (MergedSim)") 
pylab.legend(loc='upper left') 
pylab.savefig("LowEnTransient_DiscoPotential_E2_MergedSim_SampleComparison") 


pylab.figure(figsize=(10,8))
pylab.plot(log_sigma_days,nugen_samp1_e2_meansrc_events,'g--',lw=2,label='Sample 1 (Nugen)')
pylab.plot(log_sigma_days,nugen_samp2_e2_meansrc_events,'b--',lw=2,label='Sample 2 (Nugen)') 
pylab.plot(log_sigma_days,merged_samp1_e2_meansrc_events,'g-',lw=2,label='Sample 1 (MergedSim)') 
pylab.plot(log_sigma_days,merged_samp2_e2_meansrc_events,'b-',lw=2,label='Sample 2 (MergedSim)') 
pylab.xlabel(r'$Log_{10}(\sigma_{\omega})\quad (Days)$')
pylab.ylabel("NSrc Events")
pylab.axis([-6,1,3,15])
pylab.grid()
pylab.title("Sensitivity")
pylab.legend(loc='upper left')
pylab.savefig("LowEnTransient_DiscoPotential_E2_NugenANDMerged_SampleComparison")
'''



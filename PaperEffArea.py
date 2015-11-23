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



dat = tables.openFile('paper_tables/Merged_IC86.2012_genie_ic.1460.AllFiles.FINAL.hdf5')

energy = dat.root.I3MCWeightDict.col('PrimaryNeutrinoEnergy')
oweight = dat.root.I3MCWeightDict.col('OneWeight')
zenith = dat.root.PrimaryNu.col('zenith')
les = dat.root.L4Bool_LES.col('value')
hes = dat.root.L4Bool_HES.col('value')
bdt_les = dat.root.LES_L6_BDTScore.col('value')
bdt_hes = dat.root.HES_L6_BDTScore.col('value')
avg_distq = dat.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
spline_rlogl = dat.root.SplineMPEModFitParams.col('rlogl')
finite_x = dat.root.SplineMPEMod_Contained.col('x')
splinemod_zen = dat.root.SplineMPEMod.col('zenith')
dhd = dat.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
atmo = dat.root.AtmoWeight.col('value')
para1 = dat.root.SplineMPEModParaboloidFitParams.col('err1')
para2 = dat.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para = numpy.sqrt(para1**2 + para2**2) / numpy.sqrt(2)
les_nch_clean = dat.root.L4_Variables_LES.col('NCh_Clean')
les_hlc = dat.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean = dat.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc = dat.root.L4_Variables_HES.col('Nch_HLC')


dat.close()


dat = tables.openFile('third_time/Level5_nugen_numu_IC86.2013.010090.AllSubs.L6Out_Third.hdf5')

energy_nugmu = dat.root.I3MCWeightDict.col('PrimaryNeutrinoEnergy')
oweight_nugmu = dat.root.I3MCWeightDict.col('OneWeight')
zenith_nugmu = dat.root.PrimaryNu.col('zenith')
les_nugmu = dat.root.L4Bool_LES.col('value')
hes_nugmu = dat.root.L4Bool_HES.col('value')
bdt_les_nugmu = dat.root.LES_L6_BDTScore.col('value')
bdt_hes_nugmu = dat.root.HES_L6_BDTScore.col('value')
avg_distq_nugmu = dat.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
spline_rlogl_nugmu = dat.root.SplineMPEModFitParams.col('rlogl')
finite_x_nugmu = dat.root.SplineMPEMod_Contained.col('x')
splinemod_zen_nugmu = dat.root.SplineMPEMod.col('zenith')
dhd_nugmu = dat.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
atmo_nugmu = dat.root.AtmoWeight.col('value')
para1_nugmu = dat.root.SplineMPEModParaboloidFitParams.col('err1')
para2_nugmu = dat.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_nugmu = numpy.sqrt(para1_nugmu**2 + para2_nugmu**2) / numpy.sqrt(2)
les_nch_clean_nugmu = dat.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_nugmu = dat.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_nugmu = dat.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_nugmu = dat.root.L4_Variables_HES.col('Nch_HLC')


dat.close()


merged_nch = les_nch_clean + hes_nch_clean
merged_nch_nugmu = les_nch_clean_nugmu + hes_nch_clean_nugmu


samp1_pull_coeffs_numu = pickle.load(open('samp1_pull_coefficients_numu.pkl','r'))
samp1_pull_coeffs_nugmu = pickle.load(open('samp1_pull_coefficients_nugmu.pkl','r'))
samp2_pull_coeffs_numu = pickle.load(open('samp2_pull_coefficients_numu.pkl','r'))
samp2_pull_coeffs_nugmu = pickle.load(open('samp2_pull_coefficients_nugmu.pkl','r'))

def pullfit_numu(x,p):
        val = p[0]+p[1]*x+numpy.power(p[2]*x,2)+numpy.power(p[3]*x,3)+numpy.power(p[4]*x,4)
        return val

def ParaFixNuMu(sigma,nch,coeffs):
  corrected = sigma * pullfit_numu(nch,coeffs)
  return corrected

corrected_sig_para = ParaFixNuMu(sig_para,merged_nch,samp2_pull_coeffs_numu)
corrected_sig_para_nugmu = ParaFixNuMu(sig_para_nugmu,merged_nch_nugmu,samp2_pull_coeffs_numu)

#dat = tables.openFile('genie_ic.1304_1404.numu.FinalLevel_HES.hdf5')

#energy_hes = dat.root.I3MCWeightDict.col('PrimaryNeutrinoEnergy')
#oweight_hes = dat.root.I3MCWeightDict.col('OneWeight')
#zenith_hes = dat.root.PrimaryNu.col('zenith')
les_preselect = (les==1)*(numpy.isfinite(avg_distq))*(numpy.isfinite(finite_x))
hes_preselect = (hes==1)*(numpy.isfinite(avg_distq))*(numpy.isfinite(dhd))*(numpy.isfinite(spline_rlogl))

les_preselect_nugmu = (les_nugmu==1)*(numpy.isfinite(avg_distq_nugmu))*(numpy.isfinite(finite_x_nugmu))
hes_preselect_nugmu = (hes_nugmu==1)*(numpy.isfinite(avg_distq_nugmu))*(numpy.isfinite(dhd_nugmu))*(numpy.isfinite(spline_rlogl_nugmu))


samp1_final_level = ((les_preselect==1)*(spline_rlogl<7.5) + (hes_preselect==1)*(bdt_hes>-0.01))*(numpy.cos(splinemod_zen)<0.087)*((180/numpy.pi)*corrected_sig_para<45.0)

samp2_final_level = ((les_preselect==1)*(bdt_les>0.0) + (hes_preselect==1)*(bdt_hes>-0.01))*(numpy.cos(splinemod_zen)<0.087)*((180/numpy.pi)*corrected_sig_para<45.0)

samp2_final_level_nugmu = ((les_preselect_nugmu==1)*(bdt_les_nugmu>0.0) + (hes_preselect_nugmu==1)*(bdt_hes_nugmu>-0.01))*(numpy.cos(splinemod_zen_nugmu)<0.087)*((180/numpy.pi)*corrected_sig_para_nugmu<45.0)*(energy_nugmu>190.0)


zenith_samp1 = zenith[(samp1_final_level==1)]
zenith_samp2 = zenith[(samp2_final_level==1)]

energy_samp1 = energy[(samp1_final_level==1)]
energy_samp2 = energy[(samp2_final_level==1)]

oweight_samp1 = oweight[(samp1_final_level==1)]
oweight_samp2 = oweight[(samp2_final_level==1)]

energy_samp2_nugmu = energy_nugmu[(samp2_final_level_nugmu==1)]
oweight_samp2_nugmu = oweight_nugmu[(samp2_final_level_nugmu==1)]
zenith_samp2_nugmu = zenith[(samp2_final_level_nugmu==1)]

systematic_error_adj = (1-0.01*numpy.sqrt(10.0**2 + 3.0**2))

#dat.close()
fudgefactor = 1.57 ### Match NUGEN EffArea to GENIE

nevents = 300003*4000.
nevents_nugmu = 200000*1000. * fudgefactor

solidangle = 2*numpy.pi*(1-numpy.cos(95*(numpy.pi/180)))

jfenergy=numpy.array([1.7629,1.8888,2.0147,2.1405,2.2664,2.3922,2.5181,2.6440,2.7698,2.8957,3.0216,3.1474,3.2733,3.3991])
jfenergy=10**jfenergy
jf_upgoingarea=numpy.array([0.000017,0.000085,0.000429,0.001559,0.004536,0.010765,0.023325,0.047096,0.086437,0.160619,0.276211,0.468503,0.767271,1.293379]) # 90 -> 180 Zenith$
jf_downgoingarea=numpy.array([0.000000,0.000001,0.000012,0.000056,0.000214,0.000286,0.001049,0.002474,0.004479,0.009235,0.015264,0.028105,0.052265,0.092717])
ANTARES_energy = numpy.array([2.1,2.2,2.3,2.4,2.5,2.6])
ANTARES_ps_area=numpy.array([3.5e-5,5.5e-5,1e-4,1.9e-4,3.2e-4,5e-4])   ### -90 -> -35 Declination


upgoing_samp1 = zenith_samp1 >= 1.48
upgoing_samp2 = zenith_samp2 >= 1.48
upgoing_samp2_nugmu = zenith_samp2_nugmu >= 1.48


total_samp1 = zenith_samp1 >= 0.0
total_samp2 = zenith_samp2 >= 0.0


decfilter_0_numu = (numpy.rad2deg(zenith_samp2) > 85.)*(numpy.rad2deg(zenith_samp2) < 100.)
decfilter_16_numu = (numpy.rad2deg(zenith_samp2) > 100.)*(numpy.rad2deg(zenith_samp2) < 115.)
decfilter_30_numu = (numpy.rad2deg(zenith_samp2) > 115.)*(numpy.rad2deg(zenith_samp2) < 130.)
decfilter_45_numu = (numpy.rad2deg(zenith_samp2) > 130.)*(numpy.rad2deg(zenith_samp2) < 145.)
decfilter_60_numu = (numpy.rad2deg(zenith_samp2) > 145.)*(numpy.rad2deg(zenith_samp2) < 160.)
decfilter_75_numu = (numpy.rad2deg(zenith_samp2) > 160.)*(numpy.rad2deg(zenith_samp2) < 180.)

decfilter_0_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 85.)*(numpy.rad2deg(zenith_samp2_nugmu) < 100.)
decfilter_16_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 100.)*(numpy.rad2deg(zenith_samp2_nugmu) < 115.)
decfilter_30_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 115.)*(numpy.rad2deg(zenith_samp2_nugmu) < 130.)
decfilter_45_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 130.)*(numpy.rad2deg(zenith_samp2_nugmu) < 145.)
decfilter_60_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 145.)*(numpy.rad2deg(zenith_samp2_nugmu) < 160.)
decfilter_75_nugmu = (numpy.rad2deg(zenith_samp2_nugmu) > 160.)*(numpy.rad2deg(zenith_samp2_nugmu) < 180.)
### Dec Bands ###
sa0 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(95.))) - (1-numpy.cos(numpy.deg2rad(80.))))
sa16 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(80.))) - (1-numpy.cos(numpy.deg2rad(65.))))
sa30 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(65.))) - (1-numpy.cos(numpy.deg2rad(50.))))
sa45 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(50.))) - (1-numpy.cos(numpy.deg2rad(35.))))
sa60 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(35.))) - (1-numpy.cos(numpy.deg2rad(20.))))
sa75 = 2*numpy.pi*(1-numpy.cos(numpy.deg2rad(20.)))

binny = numpy.linspace(4.0,190,31)
binny_nugmu = numpy.linspace(190,1000,26)
dE=(binny[1]-binny[0])
dE_nugmu = (binny_nugmu[1]-binny_nugmu[0])

samp2_effarea_dec0 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_0_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec16 = pylab.hist(energy_samp2[decfilter_16_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_16_numu]/10000./sa16/nevents/dE,log=True,lw=2)
samp2_effarea_dec30 = pylab.hist(energy_samp2[decfilter_30_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_30_numu]/10000./sa30/nevents/dE,log=True,lw=2)
samp2_effarea_dec45 = pylab.hist(energy_samp2[decfilter_45_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_45_numu]/10000./sa45/nevents/dE,log=True,lw=2)
samp2_effarea_dec60 = pylab.hist(energy_samp2[decfilter_60_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_60_numu]/10000./sa60/nevents/dE,log=True,lw=2)
samp2_effarea_dec75 = pylab.hist(energy_samp2[decfilter_75_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_75_numu]/10000./sa75/nevents/dE,log=True,lw=2)

samp2_effarea_dec0_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_0_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_0_nugmu]/10000./sa0/nevents_nugmu/dE_nugmu,log=True,lw=2)
samp2_effarea_dec16_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_16_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_16_nugmu]/10000./sa16/nevents_nugmu/dE_nugmu,log=True,lw=2)
samp2_effarea_dec30_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_30_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_30_nugmu]/10000./sa30/nevents_nugmu/dE_nugmu,log=True,lw=2)
samp2_effarea_dec45_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_45_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_45_nugmu]/10000./sa45/nevents_nugmu/dE_nugmu,log=True,lw=2)
samp2_effarea_dec60_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_60_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_60_nugmu]/10000./sa60/nevents_nugmu/dE_nugmu,log=True,lw=2)
samp2_effarea_dec75_nugmu = pylab.hist(energy_samp2_nugmu[decfilter_75_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[decfilter_75_nugmu]/10000./sa75/nevents_nugmu/dE_nugmu,log=True,lw=2)


#samp1_effarea=pylab.hist(energy_samp1[upgoing_samp1],bins=binny,histtype='step',weights=oweight_samp1[upgoing_samp1]/10000./solidangle/nevents/dE,log=True,lw=2)
samp2_effarea=pylab.hist(energy_samp2[upgoing_samp2],bins=binny,histtype='step',weights=oweight_samp2[upgoing_samp2]/10000./solidangle/nevents/dE,log=True,lw=2)
samp2_effarea_nugmu=pylab.hist(energy_samp2_nugmu[upgoing_samp2_nugmu],bins=binny_nugmu,histtype='step',weights=oweight_samp2_nugmu[upgoing_samp2_nugmu]/10000./solidangle/nevents_nugmu/dE_nugmu,log=True,lw=2)

#totarea = leseffarea[0]+heseffarea[0]

plotbinny=binny[:-1]+dE/2
plotbinny_nugmu=binny_nugmu[:-1]+dE_nugmu/2

pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='g',label="IC86-2 Low-En Transient (GENIE)",ls='-')
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,190,10**-5,10**-1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area (GENIE)",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_DiffBin")

combined_plot_binny = numpy.hstack([plotbinny,plotbinny_nugmu])
combined_effa = numpy.hstack([samp2_effarea[0],samp2_effarea_nugmu[0]])

pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='g',label="IC86-2 Low-En Transient (GENIE)",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_nugmu[0],lw=2,c='k',label="IC86-2 Low-En Transient (Nugen)",ls='-')
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,400,10**-5,1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_Nugen")

pylab.figure()
pylab.semilogy(combined_plot_binny,systematic_error_adj*combined_effa,lw=2,c='k',label="Low-Energy Transient",ls='-')
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='k', ls='--',label="IceCube 86-String PS")
pylab.semilogy(10**ANTARES_energy,ANTARES_ps_area,lw=2,c='k',ls='-.',label=r"ANTARES $-90^{\circ} < \delta < -45^{\circ}$")
#pylab.plot(combined_plot_binny,areafit(combined_plot_binny,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]),'r--')
pylab.axis([4.0,400,10**-5,1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$\nu_{\mu}$ $A_{eff}$ (m$^2$)',fontsize=16)
#pylab.title("Muon Neutrino Effective Area")
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
#pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
pylab.grid()
pylab.savefig("PaperEffArea_ANTARES_StdPS.pdf")


reffluxvalue= 1.0 ### dN/dE (100 GeV) GeV^-1 cm^-2 s^-1

e3flux=plotbinny**-3
e3flux*=reffluxvalue/e3flux[15]
e25flux=plotbinny**-2.5
e25flux*=reffluxvalue/e25flux[15]
e35flux=plotbinny**-3.5
e35flux*=reffluxvalue/e35flux[15]


e3extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-3
e3extendedflux *= reffluxvalue/e3extendedflux[15] 
e25extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-2.5
e25extendedflux *= reffluxvalue/e25extendedflux[15]
e35extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-3.5
e35extendedflux *= reffluxvalue/e35extendedflux[15]



e3foldedflux=e3flux*samp2_effarea[0]*10000*dE
e35foldedflux=e35flux*samp2_effarea[0]*10000*dE
e25foldedflux=e25flux*samp2_effarea[0]*10000*dE

e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec16[0]*dE,samp2_effarea_dec16_nugmu[0]*dE_nugmu]))*10000
e35foldedfluence=e35extendedflux*(numpy.hstack([samp2_effarea[0]*dE,samp2_effarea_nugmu[0]*dE_nugmu]))*10000
e25foldedfluence=e25extendedflux*(numpy.hstack([samp2_effarea[0]*dE,samp2_effarea_nugmu[0]*dE_nugmu]))*10000



def areafit_nugmu(x,a,b,c,d,e):
        val = a+b*x+numpy.power(c*x,2)+numpy.power(d*x,3)+numpy.power(e*x,4)
        return val

#area_coeff_nugmu,covar_numu = curve_fit(areafit_nugmu, plotbinny_nugmu, 1.1819698556266118*samp2_effarea_nugmu[0], p0=(4,0.1,0.025,0.005,0.0001), sigma=None)

energizer=numpy.linspace(190,1000,10000)

pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='g',label="IC86-2 Low-En Transient (GENIE)",ls='-')
#pylab.semilogy(energizer,areafit_numu(energizer,area_coeff_numu[0],area_coeff_numu[1],area_coeff_numu[2],area_coeff_numu[3],area_coeff_numu[4]))
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,190,10**-5,10**-1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area (G1460)",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_DiffBin_WithFit")


pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='g',label="IC86-2 Low-En Transient (GENIE)",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_nugmu[0],lw=2,c='k',label="IC86-2 Low-En Transient (Nugen)",ls='-')
#pylab.semilogy(energizer,areafit_nugmu(energizer,area_coeff_nugmu[0],area_coeff_nugmu[1],area_coeff_nugmu[2],area_coeff_nugmu[3],area_coeff_nugmu[4]),'k--')
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,400,10**-5,1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_Nugen_WithFit")

pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='k',label="Low-En Transient Nominal",ls='--')
pylab.semilogy(plotbinny,systematic_error_adj*samp2_effarea[0],lw=2,c='k',label="Low-En Transient SysError Adjusted",ls='-')
#pylab.semilogy(energizer,areafit_numu(energizer,area_coeff_numu[0],area_coeff_numu[1],area_coeff_numu[2],area_coeff_numu[3],area_coeff_numu[4]))
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
#pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,190,10**-5,10**-2])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_WithSystematicAdjustedArea")

pylab.figure(figsize=(10,8))
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea_dec0[0],lw=2,c='g',label=r"$-5^{\circ} < \delta < 10^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec0_nugmu[0],lw=2,c='g',ls='--')
pylab.semilogy(plotbinny,samp2_effarea_dec16[0],lw=2,c='k',label=r"$10^{\circ} < \delta < 25^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec16_nugmu[0],lw=2,c='k',ls='--')
pylab.semilogy(plotbinny,samp2_effarea_dec30[0],lw=2,c='b',label=r"$25^{\circ} < \delta < 40^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec30_nugmu[0],lw=2,c='b',ls='--')
pylab.semilogy(plotbinny,samp2_effarea_dec45[0],lw=2,c='r',label=r"$40^{\circ} < \delta < 55^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec45_nugmu[0],lw=2,c='r',ls='--')
pylab.semilogy(plotbinny,samp2_effarea_dec60[0],lw=2,c='c',label=r"$55^{\circ} < \delta < 70^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec60_nugmu[0],lw=2,c='c',ls='--')
pylab.semilogy(plotbinny,samp2_effarea_dec75[0],lw=2,c='purple',label=r"$70^{\circ} < \delta < 90^{\circ}$",ls='-')
pylab.semilogy(plotbinny_nugmu,samp2_effarea_dec75_nugmu[0],lw=2,c='purple',ls='--')
pylab.axis([4.0,400,10**-5,1])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_NewGENIE_Nugen_DeclinationBands")

pickle.dump(samp2_effarea[0],open("g1460_numu_effarea_avg.pkl","w"))
pickle.dump(samp2_effarea_dec0[0],open("g1460_numu_effarea_dec0.pkl","w"))
pickle.dump(samp2_effarea_dec16[0],open("g1460_numu_effarea_dec16.pkl","w"))
pickle.dump(samp2_effarea_dec30[0],open("g1460_numu_effarea_dec30.pkl","w"))
pickle.dump(samp2_effarea_dec45[0],open("g1460_numu_effarea_dec45.pkl","w"))
pickle.dump(samp2_effarea_dec60[0],open("g1460_numu_effarea_dec60.pkl","w"))
pickle.dump(samp2_effarea_dec75[0],open("g1460_numu_effarea_dec75.pkl","w"))

pickle.dump(samp2_effarea_nugmu[0],open("g1460_nugmu_effarea_avg.pkl","w"))
pickle.dump(samp2_effarea_dec0_nugmu[0],open("g1460_nugmu_effarea_dec0.pkl","w"))
pickle.dump(samp2_effarea_dec16_nugmu[0],open("g1460_nugmu_effarea_dec16.pkl","w"))
pickle.dump(samp2_effarea_dec30_nugmu[0],open("g1460_nugmu_effarea_dec30.pkl","w"))
pickle.dump(samp2_effarea_dec45_nugmu[0],open("g1460_nugmu_effarea_dec45.pkl","w"))
pickle.dump(samp2_effarea_dec60_nugmu[0],open("g1460_nugmu_effarea_dec60.pkl","w"))
pickle.dump(samp2_effarea_dec75_nugmu[0],open("g1460_nugmu_effarea_dec75.pkl","w"))
systematic_error_adj
### Systematic Error Adjusted ###
pickle.dump(systematic_error_adj*samp2_effarea[0],open("g1460_numu_sysadj_effarea_avg.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec0[0],open("g1460_numu_sysadj_effarea_dec0.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec16[0],open("g1460_numu_sysadj_effarea_dec16.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec30[0],open("g1460_numu_sysadj_effarea_dec30.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec45[0],open("g1460_numu_sysadj_effarea_dec45.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec60[0],open("g1460_numu_sysadj_effarea_dec60.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec75[0],open("g1460_numu_sysadj_effarea_dec75.pkl","w"))

pickle.dump(systematic_error_adj*samp2_effarea_nugmu[0],open("g1460_nugmu_sysadj_effarea_avg.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec0_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec0.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec16_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec16.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec30_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec30.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec45_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec45.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec60_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec60.pkl","w"))
pickle.dump(systematic_error_adj*samp2_effarea_dec75_nugmu[0],open("g1460_nugmu_sysadj_effarea_dec75.pkl","w"))


dec0_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec0[0]*dE,samp2_effarea_dec0_nugmu[0]*dE_nugmu]))*10000
dec30_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec30[0]*dE,samp2_effarea_dec30_nugmu[0]*dE_nugmu]))*10000
dec45_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec45[0]*dE,samp2_effarea_dec45_nugmu[0]*dE_nugmu]))*10000
dec60_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec60[0]*dE,samp2_effarea_dec60_nugmu[0]*dE_nugmu]))*10000
dec75_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea_dec75[0]*dE,samp2_effarea_dec75_nugmu[0]*dE_nugmu]))*10000

sa_avg_e3foldedfluence=e3extendedflux*(numpy.hstack([samp2_effarea[0]*dE,samp2_effarea_nugmu[0]*dE_nugmu]))*10000

plotvalues_subphotospheric_fluence = numpy.array([4.1,3.5,3.09,2.9,2.7,2.65,2.6,2.55,2.46,2.41,2.39,2.37,2.35,2.34,2.325,2.315,2.315,2.315,2.315,2.315,2.315,2.315,2.315,2.32,2.33,2.331,2.34,2.345,2.35,2.37,2.38,2.4,2.47,2.55,2.65,2.75,2.85,2.95,3.03,3.15,3.25,3.35,3.45,3.55,3.64,3.75,3.85,3.95,4.05,4.12,4.14,4.18,4.22,4.26,4.3])

plotvalues_subphotospheric_fluence = -1*plotvalues_subphotospheric_fluence

subphoto_flux = 624.150934*10000.*(10**plotvalues_subphotospheric_fluence) / (combined_plot_binny**2)  ### Convert to GeV^-1 * m^-2!

pickle.dump(subphoto_flux,open('subphoto_flux.pkl',"w"))

def areafit(x,a,b,c,d,e):
        val = a+b*x+numpy.power(c*x,2)+numpy.power(d*x,3)+numpy.power(e*x,4)
        return val

coeffs = [-5.03210156e-05,1.68568118e-06,1.73773173e-04,-2.84286033e-04,-2.86752597e-04]






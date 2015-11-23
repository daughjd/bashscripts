import tables,pylab,numpy,pickle
from matplotlib import gridspec
#from scipy.optimize import curve_fit


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


dat = tables.openFile('third_time/genie_ic.1304_1404.combo_numu_numubar.L6Out_Third.hdf5')

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

samp2_final_level_nugmu = ((les_preselect_nugmu==1)*(bdt_les_nugmu>0.0) + (hes_preselect_nugmu==1)*(bdt_hes_nugmu>-0.01))*(numpy.cos(splinemod_zen_nugmu)<0.087)*((180/numpy.pi)*corrected_sig_para_nugmu<45.0)


zenith_samp1 = zenith[(samp1_final_level==1)]
zenith_samp2 = zenith[(samp2_final_level==1)]

energy_samp1 = energy[(samp1_final_level==1)]
energy_samp2 = energy[(samp2_final_level==1)]

oweight_samp1 = oweight[(samp1_final_level==1)]
oweight_samp2 = oweight[(samp2_final_level==1)]

energy_samp2_nugmu = energy_nugmu[(samp2_final_level_nugmu==1)]
oweight_samp2_nugmu = oweight_nugmu[(samp2_final_level_nugmu==1)]
zenith_samp2_nugmu = zenith_nugmu[(samp2_final_level_nugmu==1)]

systematic_error_adj = (1-0.01*numpy.sqrt(10.0**2 + 3.0**2))

#dat.close()
fudgefactor = 1.57 ### Match NUGEN EffArea to GENIE

nevents = 300003*4000.
nevents_nugmu = 300000*3999.

solidangle = 2*numpy.pi*(1-numpy.cos(95*(numpy.pi/180)))

jfenergy=numpy.array([1.7629,1.8888,2.0147,2.1405,2.2664,2.3922,2.5181,2.6440,2.7698,2.8957,3.0216,3.1474,3.2733,3.3991])
jfenergy=10**jfenergy
jf_upgoingarea=numpy.array([0.000017,0.000085,0.000429,0.001559,0.004536,0.010765,0.023325,0.047096,0.086437,0.160619,0.276211,0.468503,0.767271,1.293379])
jf_downgoingarea=numpy.array([0.000000,0.000001,0.000012,0.000056,0.000214,0.000286,0.001049,0.002474,0.004479,0.009235,0.015264,0.028105,0.052265,0.092717])

upgoing_samp1 = zenith_samp1 >= 1.48
upgoing_samp2 = zenith_samp2 >= 1.48
upgoing_samp2_nugmu = zenith_samp2_nugmu >= 1.48


total_samp1 = zenith_samp1 >= 0.0
total_samp2 = zenith_samp2 >= 0.0


binny = numpy.linspace(4.0,190,41)
binny_nugmu = numpy.linspace(190,400,21)
dE=(binny[1]-binny[0])
dE_nugmu = (binny_nugmu[1]-binny_nugmu[0])

#samp1_effarea=pylab.hist(energy_samp1[upgoing_samp1],bins=binny,histtype='step',weights=oweight_samp1[upgoing_samp1]/10000./solidangle/nevents/dE,log=True,lw=2)
samp2_effarea=pylab.hist(energy_samp2[upgoing_samp2],bins=binny,histtype='step',weights=oweight_samp2[upgoing_samp2]/10000./solidangle/nevents/dE,log=True,lw=2)
samp2_effarea_nugmu=pylab.hist(energy_samp2_nugmu[upgoing_samp2_nugmu],bins=binny,histtype='step',weights=oweight_samp2_nugmu[upgoing_samp2_nugmu]/10000./solidangle/nevents_nugmu/dE,log=True,lw=2)

#totarea = leseffarea[0]+heseffarea[0]
decfilter_0_numu = (numpy.rad2deg(zenith_samp2) > 85.)*(numpy.rad2deg(zenith_samp2) < 100.)
decfilter_16_numu = (numpy.rad2deg(zenith_samp2) > 100.)*(numpy.rad2deg(zenith_samp2) < 115.)
decfilter_30_numu = (numpy.rad2deg(zenith_samp2) > 115.)*(numpy.rad2deg(zenith_samp2) < 130.)
decfilter_45_numu = (numpy.rad2deg(zenith_samp2) > 130.)*(numpy.rad2deg(zenith_samp2) < 145.)
decfilter_60_numu = (numpy.rad2deg(zenith_samp2) > 145.)*(numpy.rad2deg(zenith_samp2) < 160.)
decfilter_75_numu = (numpy.rad2deg(zenith_samp2) > 160.)*(numpy.rad2deg(zenith_samp2) < 180.)

decfilter_0_nugmu = (numpy.rad2deg(zenith_samp2) > 85.)*(numpy.rad2deg(zenith_samp2) < 100.)
decfilter_16_nugmu = (numpy.rad2deg(zenith_samp2) > 100.)*(numpy.rad2deg(zenith_samp2) < 115.)
decfilter_30_nugmu = (numpy.rad2deg(zenith_samp2) > 115.)*(numpy.rad2deg(zenith_samp2) < 130.)
decfilter_45_nugmu = (numpy.rad2deg(zenith_samp2) > 130.)*(numpy.rad2deg(zenith_samp2) < 145.)
decfilter_60_nugmu = (numpy.rad2deg(zenith_samp2) > 145.)*(numpy.rad2deg(zenith_samp2) < 160.)
decfilter_75_nugmu = (numpy.rad2deg(zenith_samp2) > 160.)*(numpy.rad2deg(zenith_samp2) < 180.)
### Dec Bands ###
sa0 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(95.))) - (1-numpy.cos(numpy.deg2rad(80.))))
sa16 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(80.))) - (1-numpy.cos(numpy.deg2rad(65.))))
sa30 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(65.))) - (1-numpy.cos(numpy.deg2rad(50.))))
sa45 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(50.))) - (1-numpy.cos(numpy.deg2rad(35.))))
sa60 = 2*numpy.pi*((1-numpy.cos(numpy.deg2rad(35.))) - (1-numpy.cos(numpy.deg2rad(20.))))
sa75 = 2*numpy.pi*(1-numpy.cos(numpy.deg2rad(20.)))


samp2_effarea_dec0 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_0_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec16 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_16_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec30 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_30_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec45 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_45_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec60 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_60_numu]/10000./sa0/nevents/dE,log=True,lw=2)
samp2_effarea_dec75 = pylab.hist(energy_samp2[decfilter_0_numu],bins=binny,histtype='step',weights=oweight_samp2[decfilter_75_numu]/10000./sa0/nevents/dE,log=True,lw=2)


plotbinny=binny[:-1]+dE/2
plotbinny_nugmu=binny_nugmu[:-1]+dE_nugmu/2

reffluxvalue= 1.0 ### dN/dE (100 GeV) GeV^-1 cm^-2 s^-1

e3flux=plotbinny**-3
e3flux*=reffluxvalue/e3flux[31]
e25flux=plotbinny**-2.5
e25flux*=reffluxvalue/e25flux[31]
e35flux=plotbinny**-3.5
e35flux*=reffluxvalue/e35flux[31]


e3extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-3
e25extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-2.5
e35extendedflux = numpy.hstack([plotbinny,plotbinny_nugmu])**-3.5



e3foldedflux=e3flux*samp2_effarea[0]*10000*dE
e35foldedflux=e35flux*samp2_effarea[0]*10000*dE
e25foldedflux=e25flux*samp2_effarea[0]*10000*dE


#area_coeff_nugmu,covar_numu = curve_fit(areafit_nugmu, plotbinny_nugmu, 1.1819698556266118*samp2_effarea_nugmu[0], p0=(4,0.1,0.025,0.005,0.0001), sigma=None)

energizer=numpy.linspace(190,1000,10000)

pylab.figure(figsize=(10,12))
gs = gridspec.GridSpec(2, 1 ,height_ratios=[2, 1])
pylab.subplot(gs[0])
#pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 2 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='k',label="GENIE 1460",ls='-')
pylab.semilogy(plotbinny,samp2_effarea_nugmu[0],lw=2,c='k',label="GENIE 1304/1404",ls='--')
#pylab.semilogy(energizer,areafit_numu(energizer,area_coeff_numu[0],area_coeff_numu[1],area_coeff_numu[2],area_coeff_numu[3],area_coeff_numu[4]))
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
#pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,190,10**-5,10**-2])
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area",fontsize=20)
area_ratio = samp2_effarea[0]/samp2_effarea_nugmu[0]
pylab.subplot(gs[1])
pylab.plot(plotbinny,area_ratio,'ko')
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$A_{Old}/A_{New} (m^2)$',fontsize=16)
pylab.hlines(1.0,4,190,colors='k',linestyles='dashed')
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_GENIE_Comparison")

pickle.dump(samp2_effarea[0],(open("g1460_numu_effarea_avg.pkl","w"))



import tables,pylab,numpy,pickle

dat = tables.openFile('third_time/genie_ic.1304_1404.combo_numu_numubar.L6Out_Third.hdf5')

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

merged_nch = les_nch_clean + hes_nch_clean

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
#dat = tables.openFile('genie_ic.1304_1404.numu.FinalLevel_HES.hdf5')

#energy_hes = dat.root.I3MCWeightDict.col('PrimaryNeutrinoEnergy')
#oweight_hes = dat.root.I3MCWeightDict.col('OneWeight')
#zenith_hes = dat.root.PrimaryNu.col('zenith')
les_preselect = (les==1)*(numpy.isfinite(avg_distq))*(numpy.isfinite(finite_x))
hes_preselect = (hes==1)*(numpy.isfinite(avg_distq))*(numpy.isfinite(dhd))*(numpy.isfinite(spline_rlogl))

samp1_final_level = (les_preselect==1)*(spline_rlogl<7.5) + (hes_preselect==1)*(bdt_hes>-0.01)

samp2_final_level = (les_preselect==1)*(bdt_les>0.0) + (hes_preselect==1)*(bdt_hes>-0.01)


zenith_samp1 = zenith[(samp1_final_level==1)]
zenith_samp2 = zenith[(samp2_final_level==1)]

energy_samp1 = energy[(samp1_final_level==1)]
energy_samp2 = energy[(samp2_final_level==1)]

oweight_samp1 = oweight[(samp1_final_level==1)]
oweight_samp2 = oweight[(samp2_final_level==1)]



#dat.close()

nevents = 300000*2000.

#solidangle = 2*numpy.pi*(1-numpy.cos(85*(numpy.pi/180)))
solidangle = 2.17*numpy.pi

jfenergy=numpy.array([1.7629,1.8888,2.0147,2.1405,2.2664,2.3922,2.5181,2.6440,2.7698,2.8957,3.0216,3.1474,3.2733,3.3991])
jfenergy=10**jfenergy
jf_upgoingarea=numpy.array([0.000017,0.000085,0.000429,0.001559,0.004536,0.010765,0.023325,0.047096,0.086437,0.160619,0.276211,0.468503,0.767271,1.293379])
jf_downgoingarea=numpy.array([0.000000,0.000001,0.000012,0.000056,0.000214,0.000286,0.001049,0.002474,0.004479,0.009235,0.015264,0.028105,0.052265,0.092717])

upgoing_samp1 = zenith_samp1 >= 1.48
upgoing_samp2 = zenith_samp2 >= 1.48

total_samp1 = zenith_samp1 >= 0.0
total_samp2 = zenith_samp2 >= 0.0


binny = numpy.linspace(4.0,190,21)

samp1_effarea=pylab.hist(energy_samp1[upgoing_samp1],bins=binny,histtype='step',weights=oweight_samp1[upgoing_samp1]/10000./energy_samp1[upgoing_samp1]/solidangle/nevents,log=True,lw=2)
samp2_effarea=pylab.hist(energy_samp2[upgoing_samp2],bins=binny,histtype='step',weights=oweight_samp2[upgoing_samp2]/10000./energy_samp2[upgoing_samp2]/solidangle/nevents,log=True,lw=2)

#totarea = leseffarea[0]+heseffarea[0]

plotbinny=binny[:-1]+(binny[1]-binny[0])/2

pylab.figure(figsize=(10,8))
pylab.semilogy(plotbinny,samp1_effarea[0],lw=2,c='b',label="Sample 1 (GENIE)",ls='-')
pylab.semilogy(plotbinny,samp2_effarea[0],lw=2,c='g',label="Sample 2 (GENIE)",ls='-')
#pylab.semilogy(plotbinny,totarea,lw=2,c='k',label="Combined")
pylab.semilogy(jfenergy,jf_upgoingarea,lw=2, c='r', ls='--',label="IC86-1 PS (Nugen)")
pylab.axis([4.0,190,10**-7,10**-2])
pylab.xlabel(r'$E_{\nu}(GeV)$',fontsize=16)
pylab.ylabel(r'$Effective Area (m^2)$',fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.legend(loc='upper left')
pylab.title(r"Low-En Transient Analysis Effective Area (GENIE)",fontsize=20)
pylab.grid()
pylab.savefig("LowEnTransient_EffArea_GENIE")


e3foldedflux=(binny[1]-binny[0])*(plotbinny**-3)*samp2_effarea[0]*10000

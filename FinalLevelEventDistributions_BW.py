pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']=18.0
pylab.rcParams['xtick.labelsize']=18.0
pylab.rcParams['lines.markeredgewidth']=1.0

pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Computer Modern Roman')

corsika_adj = 0.75
combined_corsk = corsika_adj * combined_corsk

nue_adj = 0.7
atmo_nue = nue_adj * atmo_nue

### Errors ###
#terr=pylab.hist(merged_nch_nutau[nch_real_final_level_nutau==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_nutau[nch_real_final_level_nutau==1])**2)
eerr=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_nue[nch_real_final_level_nue==1])**2)
merr=pylab.hist(merged_nch_numu[nch_real_final_level_numu==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_numu[nch_real_final_level_numu==1])**2)
cerr=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*combined_corsk[nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(bdt_les_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*C9255EventWeight[nch_real_final_level_c9255==1])**2)
##############

pylab.figure()
pylab.hist(merged_nch_data[nch_real_final_level_data==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[nch_real_final_level_data==1])
#chist=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(merged_nch_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[nch_real_final_level_c9255==1])
muhist=pylab.hist(merged_nch_combined_numu[nch_real_final_level_combined_numu==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[nch_real_final_level_combined_numu==1])
ehist=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'$\nu_{e}$',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[nch_real_final_level_nue==1])
#pylab.hist(merged_nch_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]#+chist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(10,61,26)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([10,60,0.0,0.06])
pylab.legend()
pylab.grid()
#pylab.title("NCh Distribution")
pylab.xlabel('NCh')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_NchDistribution_Rates_NormedNugen_BW_NoCorsika_NuEAdj")

eerr=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,60,11),log=False,normed=False,weights=(1000*atmo_nue[nch_real_final_level_nue==1])**2)
merr=pylab.hist(merged_nch_numu[nch_real_final_level_numu==1],bins=numpy.linspace(10,60,11),log=False,normed=False,weights=(1000*atmo_numu[nch_real_final_level_numu==1])**2)
cerr=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,60,11),log=False,normed=False,weights=(1000*combined_corsk[nch_real_final_level_combined_corsk==1])**2)


pylab.figure()
pylab.hist(merged_nch_data[nch_real_final_level_data==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[nch_real_final_level_data==1])
chist=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(merged_nch_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[nch_real_final_level_c9255==1])
muhist=pylab.hist(merged_nch_combined_numu[nch_real_final_level_combined_numu==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[nch_real_final_level_combined_numu==1])
ehist=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'$\nu_{e}$',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[nch_real_final_level_nue==1])
#pylab.hist(merged_nch_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]#+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(10,60,11)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([10,60,0.0,0.12])
pylab.legend()
pylab.grid()
#pylab.title("NCh Distribution")
pylab.xlabel('NCh')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_NchDistribution_Rates_NormedNugen_BW_CoarseBinning")


eerr=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,60,11),log=False,normed=False,cumulative=True,weights=(1000*atmo_nue[nch_real_final_level_nue==1])**2)
merr=pylab.hist(merged_nch_numu[nch_real_final_level_numu==1],bins=numpy.linspace(10,60,11),log=False,normed=False,cumulative=True,weights=(1000*atmo_numu[nch_real_final_level_numu==1])**2)
cerr=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,60,11),log=False,normed=False,cumulative=True,weights=(1000*combined_corsk[nch_real_final_level_combined_corsk==1])**2)

pylab.figure()
pylab.hist(merged_nch_data[nch_real_final_level_data==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,cumulative=True,weights=1000*dataweight[nch_real_final_level_data==1])
chist=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,cumulative=True,weights=1000*combined_corsk[nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(merged_nch_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,cumulative=True,weights=1000*C9255EventWeight[nch_real_final_level_c9255==1])
muhist=pylab.hist(merged_nch_combined_numu[nch_real_final_level_combined_numu==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='k',normed=False,cumulative=True,ls='dashed',weights=1000*normed_combined_atmo_numu[nch_real_final_level_combined_numu==1])
ehist=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label=r'$\nu_{e}$',ec='k',normed=False,cumulative=True,ls='dashdot',weights=1000*atmo_nue[nch_real_final_level_nue==1])
#pylab.hist(merged_nch_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(10,60,11),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,cumulative=True,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]#+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(10,60,11)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.axis([10,60,0.0,0.12])
pylab.legend(loc='upper left')
pylab.grid()
#pylab.title("NCh Distribution")
pylab.xlabel('NCh')
pylab.ylabel('Cumulative')
pylab.savefig("FinalSample_NchDistribution_Rates_NormedNugen_BW_CoarseBinning_Cumulative")


pylab.figure())
muhist=pylab.hist(merged_nch_numu[nch_real_final_level_numu==1],bins=numpy.linspace(10,61,52),histtype='step',lw=2,log=False,label='(G)NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[nch_real_final_level_numu==1])
hemuhist=pylab.hist(merged_nch_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)],bins=numpy.linspace(10,61,52),histtype='step',lw=2,log=False,label=r'(N)NuMu',ec='g',normed=False,ls='solid',weights=1000*atmo_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)])
#pylab.hist(merged_nch_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(10,61,52),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
#errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
#pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([10,60,0.0,0.04])
pylab.legend()
pylab.grid()
pylab.title("NCh Distribution")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("FinalNuMuSample_NchDistribution_Rates")


pylab.figure(figsize=(10,8))
muhist=pylab.hist(primary_energy_numu[nch_real_final_level_numu==1],bins=numpy.linspace(5,400,20),histtype='step',lw=2,log=False,label='(G)NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[nch_real_final_level_numu==1])
hemuhist=pylab.hist(primary_energy_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)],bins=numpy.linspace(5,400,20),histtype='step',lw=2,log=False,label=r'(N)NuMu',ec='g',normed=False,ls='solid',weights=1000*atmo_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)])
adjhemuhist=pylab.hist(primary_energy_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)],bins=numpy.linspace(5,400,20),histtype='step',lw=2,log=False,label=r'(N)NuMu Normalized',ec='g',normed=False,ls='dashed',weights=0.625*1000*atmo_nugmu[(nch_real_final_level_nugmu==1)*(primary_energy_nugmu>5)])
#pylab.hist(primary_energy_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(5,400,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
#errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
#pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([5,400,0.0,0.1])
pylab.legend()
pylab.grid()
pylab.title("Simulation Energy Distribution")
pylab.xlabel('Energy (GeV)')
pylab.ylabel('mHz')
pylab.savefig("FinalNuMuSample_EnergyDistribution_Rates")



### Reco Zenith ###
eerr=pylab.hist(numpy.cos(splinemod_zen_nue[real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(numpy.cos(splinemod_zen_combined_numu[real_final_level_combined_numu==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])**2)
cerr=pylab.hist(numpy.cos(splinemod_zen_combined_corsk[real_final_level_combined_corsk==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*combined_corsk[real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(numpy.cos(splinemod_zen_c9255[real_final_level_c9255==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*C9255EventWeight[real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(numpy.cos(splinemod_zen_data[real_final_level_data==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
chist=pylab.hist(numpy.cos(splinemod_zen_combined_corsk[real_final_level_combined_corsk==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(numpy.cos(splinemod_zen_c9255[real_final_level_c9255==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1])
muhist=pylab.hist(numpy.cos(splinemod_zen_combined_numu[real_final_level_combined_numu==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])
ehist=pylab.hist(numpy.cos(splinemod_zen_nue[real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Nue',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(numpy.cos(splinemod_zen_nuge[real_final_level_nuge==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-1,0.087,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1,0.087,0.0,0.08])
pylab.legend(loc='upper left')
pylab.grid()
#pylab.title("Zenith Distribution")
pylab.xlabel(r'Cos($\Theta$)')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_SplineMPEZenith_Rates_NormedNugen_BW")


### CORSIKA MC INFO ###
pylab.figure()
pylab.hist(numpy.cos(zenith_c9622[real_final_level_c9622==1]),bins=numpy.linspace(0,1,10),histtype='step',lw=2,log=False,label=r'atm $\mu$ (9622)',ec='k',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1],ls="solid")
pylab.hist(numpy.cos(zenith_c9255[real_final_level_c9255==1]),bins=numpy.linspace(0,1,10),histtype='step',lw=2,log=False,label=r'atm $\mu$ (9255)',ec='k',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1],ls="dotted")
pylab.title("CORSIKA Primary Zenith")
pylab.xlabel(r'Cos($\Theta$)')
pylab.legend(loc="upper left")
pylab.ylabel('mHz')
pylab.grid()
pylab.savefig("CorsikaPrimaryZenith_FinalLevel")

pylab.figure()
pylab.hist(numpy.log10(energy_c9622[real_final_level_c9622==1]),bins=numpy.linspace(2,8,20),histtype='step',lw=2,log=False,label=r'atm $\mu$ (9622)',ec='k',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1],ls="solid")
pylab.hist(numpy.log10(energy_c9255[real_final_level_c9255==1]),bins=numpy.linspace(2,8,20),histtype='step',lw=2,log=False,label=r'atm $\mu$ (9255)',ec='k',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1],ls="dotted")
pylab.title("CORSIKA Primary Energy")
pylab.xlabel('Energy (GeV)')
pylab.axis([2,8,0,0.02])
pylab.legend(loc="upper right")
pylab.ylabel('mHz')
pylab.grid()
pylab.savefig("CorsikaPrimaryEnergy_FinalLevel")





### Reco Azimuth ###




eerr=pylab.hist(bdt_hes_nue[(hes_nue==1)*(real_final_level_nue==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nue[(hes_nue==1)*(real_final_level_nue==1)])**2)
merr=pylab.hist(bdt_hes_numu[(hes_numu==1)*(real_final_level_numu==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_numu[(hes_numu==1)*(real_final_level_numu==1)])**2)
#cerr=pylab.hist(bdt_hes_c7437[(hes_c7437==1)*(real_final_level_c7437==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*C7437EventWeight[(hes_c7437==1)*(real_final_level_c7437==1)])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(bdt_hes_data[(hes_data==1)*(real_final_level_data==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(hes_data==1)*(real_final_level_data==1)])
#thist=pylab.hist(bdt_hes_nutau[(hes_nutau)*(real_final_level_nutau==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[(hes_nutau)*(real_final_level_nutau==1])
#pylab.hist(bdt_hes_nugmu[(hes_)*(real_final_level_nugmu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[(hes_)*(real_final_level_nugmu==1])
#chist=pylab.hist(bdt_hes_c7437[(hes_c7437==1)*(real_final_level_c7437==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[(hes_c7437==1)*(real_final_level_c7437==1)])
muhist=pylab.hist(bdt_hes_nugmu[(hes_nugmu==1)*(real_final_level_nugmu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='b',normed=False,ls='solid',weights=1000*atmo_nugmu[(hes_nugmu==1)*(real_final_level_nugmu==1)])
ehist=pylab.hist(bdt_hes_nue[(hes_nue==1)*(real_final_level_nue==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_e$',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[(hes_nue==1)*(real_final_level_nue==1)])
#pylab.hist(bdt_hes_nuge[(hes_)*(real_final_level_nuge==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[(hes_)*(real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(-1,1,40)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-0.2,1.0,0.0,0.2])
pylab.legend()
pylab.grid()
pylab.title("Final Level LES BDT Score Distribution")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("HES_PostCutBDTScoreDist_NormalizedRates_G1460")


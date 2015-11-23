pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=18.0
pylab.rcParams['axes.titlesize']=20.0
pylab.rcParams['ytick.labelsize']=18.0
pylab.rcParams['xtick.labelsize']=18.0
pylab.rcParams['lines.markeredgewidth']=1.0

pylab.rc ('text', usetex=True)
pylab.rc ('font', family='serif')
pylab.rc ('font', serif='Computer Modern Roman')

### Errors ###
#terr=pylab.hist(merged_nch_nutau[nch_real_final_level_nutau==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_nutau[nch_real_final_level_nutau==1])**2)
eerr=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_nue[nch_real_final_level_nue==1])**2)
merr=pylab.hist(merged_nch_numu[nch_real_final_level_numu==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*atmo_numu[nch_real_final_level_numu==1])**2)
cerr=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*combined_corsk[nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(bdt_les_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,61,26),log=False,normed=False,weights=(1000*C9255EventWeight[nch_real_final_level_c9255==1])**2)
##############

pylab.figure()
pylab.hist(merged_nch_data[nch_real_final_level_data==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[nch_real_final_level_data==1])
chist=pylab.hist(merged_nch_combined_corsk[nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='r',normed=False,weights=1000*combined_corsk[nch_real_final_level_combined_corsk==1])
#hechist=pylab.hist(merged_nch_c9255[nch_real_final_level_c9255==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[nch_real_final_level_c9255==1])
muhist=pylab.hist(merged_nch_combined_numu[nch_real_final_level_combined_numu==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='b',normed=False,ls='solid',weights=1000*normed_combined_atmo_numu[nch_real_final_level_combined_numu==1])
ehist=pylab.hist(merged_nch_nue[nch_real_final_level_nue==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label=r'$\nu_{e}$',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[nch_real_final_level_nue==1])
#pylab.hist(merged_nch_nuge[nch_real_final_level_nuge==1],bins=numpy.linspace(10,61,26),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(10,61,26)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([10,60,0.0,0.06])
pylab.legend()
pylab.grid()
pylab.title("NCh Distribution")
pylab.xlabel('NCh')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_NchDistribution_Rates_NormedNugen")

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
chist=pylab.hist(numpy.cos(splinemod_zen_combined_corsk[real_final_level_combined_corsk==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='r',normed=False,weights=1000*combined_corsk[real_final_level_combined_corsk==1])
#hechist=pylab.hist(numpy.cos(splinemod_zen_c9255[real_final_level_c9255==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1])
muhist=pylab.hist(numpy.cos(splinemod_zen_combined_numu[real_final_level_combined_numu==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])
ehist=pylab.hist(numpy.cos(splinemod_zen_nue[real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(numpy.cos(splinemod_zen_nuge[real_final_level_nuge==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-1,0.087,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1,0.087,0.0,0.08])
pylab.legend(loc='upper left')
pylab.grid()
pylab.title("Zenith Distribution")
pylab.xlabel(r'Cos($\Theta$)')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_SplineMPEZenith_Rates_NormedNugen")

### Reco Azimuth ###
eerr=pylab.hist(57.2957795*splinemod_azi_nue[real_final_level_nue==1],bins=numpy.linspace(0,360,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(57.2957795*splinemod_azi_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(0,360,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])**2)
cerr=pylab.hist(57.2957795*splinemod_azi_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(0,360,20),log=False,normed=False,weights=(1000*combined_corsk[real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(57.2957795*splinemod_azi_c9255[real_final_level_c9255==1],bins=numpy.linspace(0,360,20),log=False,normed=False,weights=(1000*C9255EventWeight[real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(57.2957795*splinemod_azi_data[real_final_level_data==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
chist=pylab.hist(57.2957795*splinemod_azi_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='r',normed=False,weights=1000*combined_corsk[real_final_level_combined_corsk==1])
#hechist=pylab.hist(57.2957795*splinemod_azi_c9255[real_final_level_c9255==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1])
muhist=pylab.hist(57.2957795*splinemod_azi_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])
ehist=pylab.hist(57.2957795*splinemod_azi_nue[real_final_level_nue==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(57.2957795*splinemod_azi_nuge[real_final_level_nuge==1],bins=numpy.linspace(0,360,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,360,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,360,0.0,0.08])
pylab.legend(loc='upper left')
pylab.grid()
pylab.title("Azimuth Distribution")
pylab.xlabel(r'$\Phi$)')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_SplineMPEAzimuth_Rates_NormedNugen")


### Reco RLogL ###
eerr=pylab.hist(spline_rlogl_nue[real_final_level_nue==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(spline_rlogl_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])**2)
cerr=pylab.hist(spline_rlogl_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*combined_corsk[real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(spline_rlogl_c9255[real_final_level_c9255==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*C9255EventWeight[real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(spline_rlogl_data[real_final_level_data==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
chist=pylab.hist(spline_rlogl_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='r',normed=False,weights=1000*combined_corsk[real_final_level_combined_corsk==1])
#hechist=pylab.hist(spline_rlogl_c9255[real_final_level_c9255==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1])
muhist=pylab.hist(spline_rlogl_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])
ehist=pylab.hist(spline_rlogl_nue[real_final_level_nue==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(spline_rlogl_nuge[real_final_level_nuge==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(5,12,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([5,12,0.0,0.15])
pylab.legend(loc='upper right')
pylab.grid()
pylab.title("Spline MPE RLogL")
pylab.xlabel('RLogL')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_SplineMPERLogL_Rates_NormedNugen")

### Finite Z ###
eerr=pylab.hist(finite_z_nue[real_final_level_nue==1],bins=numpy.linspace(-700,-150,30),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(finite_z_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(-700,-150,30),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])**2)
cerr=pylab.hist(finite_z_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(-700,-150,30),log=False,normed=False,weights=(1000*combined_corsk[real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(finite_z_c9255[real_final_level_c9255==1],bins=numpy.linspace(-700,-150,30),log=False,normed=False,weights=(1000*C9255EventWeight[real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(finite_z_data[real_final_level_data==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
chist=pylab.hist(finite_z_combined_corsk[real_final_level_combined_corsk==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='r',normed=False,weights=1000*combined_corsk[real_final_level_combined_corsk==1])
#hechist=pylab.hist(finite_z_c9255[real_final_level_c9255==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[real_final_level_c9255==1])
muhist=pylab.hist(finite_z_combined_numu[real_final_level_combined_numu==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*normed_combined_atmo_numu[real_final_level_combined_numu==1])
ehist=pylab.hist(finite_z_nue[real_final_level_nue==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(finite_z_nuge[real_final_level_nuge==1],bins=numpy.linspace(-700,-150,30),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-700,-150,30)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-700,-150,0.0,0.08])
pylab.legend(loc='upper right')
pylab.grid()
pylab.title("Finite Reco Z")
pylab.xlabel('Z(m)')
pylab.ylabel('mHz')
pylab.savefig("FinalSample_FiniteZ_Rates_NormedNugen")

### Final Level Resolution ###

energy_bins = numpy.linspace(30,300,20)
finalsample_median_resolution,plotbins_energy=EnergyBinMedianRes(combined_primary_energy_numu[real_final_level_combined_numu],combined_spline_primary_diff_numu[real_final_level_combined_numu],energy_bins)
pylab.figure()
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_median_resolution)),'k-',lw=2)
pylab.axis([30,300,0.0,20.0])
pylab.xlabel("Neutrino Energy (GeV)")
pylab.ylabel("Angular Error ($^{\circ}$)")
pylab.title("Median Muon Neutrino Resolution")
pylab.grid()
pylab.savefig("MedianResolution_CombinedSimulationFinalSample")

energy_bins = numpy.linspace(30,300,20)
finalsample_contained_resolution,plotbins_energy=EnergyBinOneSigmaContainedRes(combined_primary_energy_numu[real_final_level_combined_numu],combined_spline_primary_diff_numu[real_final_level_combined_numu],energy_bins)
pylab.figure()
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_contained_resolution)),'k-',lw=2)
pylab.axis([30,300,0.0,30.0])
pylab.xlabel("Neutrino Energy (GeV)")
pylab.ylabel("Angular Error ($^{\circ}$)")
pylab.title("Muon Neutrino Resolution (68$\%$ Containment)")
pylab.grid()
pylab.savefig("ContainedResolution_CombinedSimulationFinalSample")

pylab.figure()
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_contained_resolution)),'k--',lw=2,label='68$\%$')
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_median_resolution)),'k-',lw=2,label="Median")
pylab.axis([30,300,0.0,30.0])
pylab.xlabel("Neutrino Energy (GeV)")
pylab.ylabel("Angular Error ($^{\circ}$)")
#pylab.title("Muon Neutrino Resolution")
pylab.legend(loc='upper right')
pylab.grid()
pylab.savefig("CombinedResolution_CombinedSimulationFinalSample")

finalsample_mediankinematic_error,plotbins_energy = EnergyBinMedianRes(combined_primary_energy_numu[real_final_level_combined_numu_cc],combined_primary_daughter_diff_numu[real_final_level_combined_numu_cc],energy_bins)

pylab.figure()
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_contained_resolution)),'k--',lw=2,label='68$\%$')
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_median_resolution)),'k-',lw=2,label="Median")
pylab.plot(plotbins_energy,numpy.rad2deg(numpy.array(finalsample_mediankinematic_error)),'k-.',lw=2,label="Kinematic Angle (Median)")
pylab.axis([30,300,0.0,30.0])
pylab.xlabel("Neutrino Energy (GeV)")
pylab.ylabel("Angular Error ($^{\circ}$)")
pylab.title("Muon Neutrino Resolution")
pylab.legend(loc='upper right')
pylab.grid()
pylab.savefig("CombinedResolution_CombinedSimulationFinalSample_WKinematic_CCOnly")



### Include CC and NC EVENTS?!? ###


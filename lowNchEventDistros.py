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
#terr=pylab.hist(merged_nch_nutau[low_nch_real_final_level_nutau==1],bins=numpy.linspace(10,20,11),log=False,normed=False,weights=(1000*atmo_nutau[low_nch_real_final_level_nutau==1])**2)
eerr=pylab.hist(merged_nch_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(10,20,11),log=False,normed=False,weights=(1000*atmo_nue[low_nch_real_final_level_nue==1])**2)
merr=pylab.hist(merged_nch_numu[low_nch_real_final_level_numu==1],bins=numpy.linspace(10,20,11),log=False,normed=False,weights=(1000*atmo_numu[low_nch_real_final_level_numu==1])**2)
cerr=pylab.hist(merged_nch_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,20,11),log=False,normed=False,weights=(1000*combined_corsk[low_nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(bdt_les_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(10,20,11),log=False,normed=False,weights=(1000*C9255EventWeight[low_nch_real_final_level_c9255==1])**2)
##############

pylab.figure()
pylab.hist(merged_nch_data[low_nch_real_final_level_data==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[low_nch_real_final_level_data==1])
chist=pylab.hist(merged_nch_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[low_nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(merged_nch_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[low_nch_real_final_level_c9255==1])
muhist=pylab.hist(merged_nch_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])
ehist=pylab.hist(merged_nch_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label=r'$\nu_{e}$',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[low_nch_real_final_level_nue==1])
#pylab.hist(merged_nch_nuge[low_nch_real_final_level_nuge==1],bins=numpy.linspace(10,20,11),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[low_nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(10,20,11)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([10,20,0.0,0.15])
pylab.legend(loc='upper left')
pylab.grid()
pylab.title(r"NCh Distribution (Nch$<=20$)")
pylab.xlabel('NCh')
pylab.ylabel('mHz')
pylab.savefig("mc_investigations/FinalSample_NchDistribution_Rates_NormedNugen_LowNChOnly")


### Reco Zenith ###
eerr=pylab.hist(numpy.cos(splinemod_zen_nue[low_nch_real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nue[low_nch_real_final_level_nue==1])**2)
merr=pylab.hist(numpy.cos(splinemod_zen_combined_numu[low_nch_real_final_level_combined_numu==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])**2)
cerr=pylab.hist(numpy.cos(splinemod_zen_combined_corsk[low_nch_real_final_level_combined_corsk==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*combined_corsk[low_nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(numpy.cos(splinemod_zen_c9255[low_nch_real_final_level_c9255==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*C9255EventWeight[low_nch_real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(numpy.cos(splinemod_zen_data[low_nch_real_final_level_data==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[low_nch_real_final_level_data==1])
chist=pylab.hist(numpy.cos(splinemod_zen_combined_corsk[low_nch_real_final_level_combined_corsk==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[low_nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(numpy.cos(splinemod_zen_c9255[low_nch_real_final_level_c9255==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[low_nch_real_final_level_c9255==1])
muhist=pylab.hist(numpy.cos(splinemod_zen_combined_numu[low_nch_real_final_level_combined_numu==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])
ehist=pylab.hist(numpy.cos(splinemod_zen_nue[low_nch_real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Nue',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[low_nch_real_final_level_nue==1])
#pylab.hist(numpy.cos(splinemod_zen_nuge[low_nch_real_final_level_nuge==1]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[low_nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-1,0.087,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1,0.087,0.0,0.03])
pylab.legend(loc='upper left')
pylab.grid()
pylab.title(r"Zenith Distribution (Nch$<=20$)")
pylab.xlabel(r'Cos($\Theta$)')
pylab.ylabel('mHz')
pylab.savefig("mc_investigations/FinalSample_SplineMPEZenith_Rates_NormedNugen_LowNch")


### Reco Azimuth ###
eerr=pylab.hist(splinemod_azi_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(0,6.2831821490221742,20),log=False,normed=False,weights=(1000*atmo_nue[low_nch_real_final_level_nue==1])**2)
merr=pylab.hist(splinemod_azi_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(0,6.2831821490221742,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])**2)
cerr=pylab.hist(splinemod_azi_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(0,6.2831821490221742,20),log=False,normed=False,weights=(1000*combined_corsk[low_nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(splinemod_azi_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(0,6.2831821490221742,20),log=False,normed=False,weights=(1000*C9255EventWeight[low_nch_real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(splinemod_azi_data[low_nch_real_final_level_data==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[low_nch_real_final_level_data==1])
chist=pylab.hist(splinemod_azi_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[low_nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(splinemod_azi_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[low_nch_real_final_level_c9255==1])
muhist=pylab.hist(splinemod_azi_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label='NuMu',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])
ehist=pylab.hist(splinemod_azi_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label='Nue',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[low_nch_real_final_level_nue==1])
#pylab.hist(splinemod_azi_nuge[low_nch_real_final_level_nuge==1],bins=numpy.linspace(0,6.2831821490221742,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[low_nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,6.2831821490221742,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,6.2831821490221742,0.0,0.03])
pylab.legend(loc='upper left')
pylab.grid()
pylab.title(r"Azimuth Distribution (Nch$<=20$)")
pylab.xlabel(r'$\phi$')
pylab.ylabel('mHz')
pylab.savefig("mc_investigations/FinalSample_SplineMPEAzimuth_Rates_NormedNugen_LowNch")

### Error Estimation ###
eerr=pylab.hist(corrected_samp2_sig_para_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(0,0.78539816339744828,20),log=False,normed=False,weights=(1000*atmo_nue[low_nch_real_final_level_nue==1])**2)
merr=pylab.hist(corrected_samp2_sig_para_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(0,0.78539816339744828,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])**2)
cerr=pylab.hist(corrected_samp2_sig_para_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(0,0.78539816339744828,20),log=False,normed=False,weights=(1000*combined_corsk[low_nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(corrected_samp2_sig_para_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(0,0.78539816339744828,20),log=False,normed=False,weights=(1000*C9255EventWeight[low_nch_real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(corrected_samp2_sig_para_data[low_nch_real_final_level_data==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[low_nch_real_final_level_data==1])
chist=pylab.hist(corrected_samp2_sig_para_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[low_nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(corrected_samp2_sig_para_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[low_nch_real_final_level_c9255==1])
muhist=pylab.hist(corrected_samp2_sig_para_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label='NuMu',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])
ehist=pylab.hist(corrected_samp2_sig_para_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label='Nue',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[low_nch_real_final_level_nue==1])
#pylab.hist(corrected_samp2_sig_para_nuge[low_nch_real_final_level_nuge==1],bins=numpy.linspace(0,0.78539816339744828,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[low_nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,0.78539816339744828,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,0.78539816339744828,0.0,0.03])
pylab.legend(loc='upper right')
pylab.grid()
pylab.title(r"Paraboloid Sigma (Nch$<=20$)")
pylab.xlabel(r'$\sigma_{est}$')
pylab.ylabel('mHz')
pylab.savefig("mc_investigations/FinalSample_ParaboloidErrorEst_Rates_NormedNugen_LowNch")

### Spline RLogL ###
eerr=pylab.hist(spline_rlogl_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*atmo_nue[low_nch_real_final_level_nue==1])**2)
merr=pylab.hist(spline_rlogl_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])**2)
cerr=pylab.hist(spline_rlogl_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*combined_corsk[low_nch_real_final_level_combined_corsk==1])**2)
#hecerr=pylab.hist(spline_rlogl_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(5,12,20),log=False,normed=False,weights=(1000*C9255EventWeight[low_nch_real_final_level_c9255==1])**2)

pylab.figure()
pylab.hist(spline_rlogl_data[low_nch_real_final_level_data==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[low_nch_real_final_level_data==1])
chist=pylab.hist(spline_rlogl_combined_corsk[low_nch_real_final_level_combined_corsk==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label=r'atm $\mu$',ec='k',normed=False,weights=1000*combined_corsk[low_nch_real_final_level_combined_corsk==1],ls="dotted")
#hechist=pylab.hist(spline_rlogl_c9255[low_nch_real_final_level_c9255==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[low_nch_real_final_level_c9255==1])
muhist=pylab.hist(spline_rlogl_combined_numu[low_nch_real_final_level_combined_numu==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='NuMu',ec='k',normed=False,ls='dashed',weights=1000*normed_combined_atmo_numu[low_nch_real_final_level_combined_numu==1])
ehist=pylab.hist(spline_rlogl_nue[low_nch_real_final_level_nue==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='Nue',ec='k',normed=False,ls='dashdot',weights=1000*atmo_nue[low_nch_real_final_level_nue==1])
#pylab.hist(spline_rlogl_nuge[low_nch_real_final_level_nuge==1],bins=numpy.linspace(5,12,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[low_nch_real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(5,12,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="k",ms=8,label="MC Total",ls='none')
#pylab.plot(binzo,mctot,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([5,12,0.0,0.2])
pylab.legend(loc='upper right')
pylab.grid()
pylab.title(r"SplineMPE RLogL (Nch$<=20$)")
pylab.xlabel(r'rLogL')
pylab.ylabel('mHz')
pylab.savefig("mc_investigations/FinalSample_SplineRLogL_Rates_NormedNugen_LowNch")


### Pulse Series Disperion ###
pylab.figure()
pylab.hist(low_samplepulsearrays_data,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='k',label="Data",log=True,normed=True)
pylab.hist(low_samplepulsearrays_nugmu,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='b',label="Nugen",log=True,normed=True)
pylab.hist(low_samplepulsearrays_numu,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='g',label="GENIE",log=True,normed=True)
pylab.xlabel("Pulse Time (W.R.T. Median Pulse)")
pylab.ylabel("Log(Counts)")
pylab.title("Low NCh Pulse Dispersion")
pylab.grid()
pylab.legend()
pylab.savefig("PulseDispersion_LowNch")

pylab.figure()
pylab.hist(high_samplepulsearrays_data,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='k',label="Data",log=True,normed=True)
pylab.hist(high_samplepulsearrays_nugmu,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='b',label="Nugen",log=True,normed=True)
pylab.hist(high_samplepulsearrays_numu,bins=numpy.linspace(-10000,10000,100),histtype='step',ec='g',label="GENIE",log=True,normed=True)
pylab.xlabel("Pulse Time (W.R.T. Median Pulse)")
pylab.ylabel("Log(Counts)")
pylab.title("High NCh (>20) Pulse Dispersion")
pylab.grid()
pylab.legend()
pylab.savefig("PulseDispersion_HighNch")





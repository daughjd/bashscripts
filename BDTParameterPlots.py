pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=16.0
pylab.rcParams['axes.titlesize']=18.0
pylab.rcParams['ytick.labelsize']='large'
pylab.rcParams['xtick.labelsize']='large'
pylab.rcParams['lines.markeredgewidth']=1.0

### Adding Spline-ParaSpline Difference ###

### Errors ###
#terr=pylab.hist(57.3*spline_paraspline_diff_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,25,30),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(57.3*spline_paraspline_diff_nue[real_final_level_nue==1],bins=numpy.linspace(0,25,30),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(57.3*spline_paraspline_diff_numu[real_final_level_numu==1],bins=numpy.linspace(0,25,30),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(57.3*spline_paraspline_diff_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,25,30),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(57.3*spline_paraspline_diff_data[real_final_level_data==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#thist=pylab.hist(57.3*spline_paraspline_diff_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(57.3*spline_paraspline_diff_nugmu[real_final_level_nugmu==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[real_final_level_nugmu==1])
chist=pylab.hist(57.3*spline_paraspline_diff_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
muhist=pylab.hist(57.3*spline_paraspline_diff_numu[real_final_level_numu==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[real_final_level_numu==1])
ehist=pylab.hist(57.3*spline_paraspline_diff_nue[real_final_level_nue==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#sig1muhist=pylab.hist(57.3*spline_paraspline_diff_numu[(real_final_level_numu==1)*(badly_recoed_numu==0)],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3}$',ec='b',normed=False,ls='dashed',weights=0.00006*(primary_energy_numu[(real_final_level_numu==1)*(badly_recoed_numu==0)])**-0.5)
#pylab.hist(57.3*spline_paraspline_diff_nuge[real_final_level_nuge==1],bins=numpy.linspace(0,25,30),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,25,30)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
#pylab.axis([0,120,0.0,0.4])
pylab.legend()
pylab.grid()
pylab.title("Space Angle (SplineMod-SplinePara)")
pylab.xlabel('Space Angle (Deg)')
pylab.ylabel('mHz')
pylab.savefig("LES_Spline_ParaSpline_Diff_FinalLevelDist_G1460")


### Errors ###
#terr=pylab.hist(vetotrackcharge_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,6,36),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(vetotrackcharge_nue[real_final_level_nue==1],bins=numpy.linspace(0,6,36),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(vetotrackcharge_numu[real_final_level_numu==1],bins=numpy.linspace(0,6,36),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(vetotrackcharge_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,6,36),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(vetotrackcharge_data[real_final_level_data==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#thist=pylab.hist(vetotrackcharge_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(vetotrackcharge_nugmu[real_final_level_nugmu==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[real_final_level_nugmu==1])
chist=pylab.hist(vetotrackcharge_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
muhist=pylab.hist(vetotrackcharge_numu[real_final_level_numu==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[real_final_level_numu==1])
ehist=pylab.hist(vetotrackcharge_nue[real_final_level_nue==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(vetotrackcharge_nuge[real_final_level_nuge==1],bins=numpy.linspace(0,6,36),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,6,36)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,6,0.0,0.4])
pylab.legend()
pylab.grid()
pylab.title("L5VetoTrack Charge")
pylab.xlabel('pe')
pylab.ylabel('mHz')
pylab.savefig("LES_VetoTrackCharge_FinalLevelDist")

### Errors ###
#terr=pylab.hist(avg_distq_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(avg_distq_nue[real_final_level_nue==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(avg_distq_numu[real_final_level_numu==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(avg_distq_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(avg_distq_data[numpy.isfinite(avg_distq_data)*(real_final_level_data==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#pylab.hist(avg_distq_nugmu[numpy.isfinite(avg_distq_nugmu)*(real_final_level_nugmu==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
muhist=pylab.hist(avg_distq_numu[numpy.isfinite(avg_distq_numu)*(real_final_level_numu==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#thist=pylab.hist(avg_distq_nutau[numpy.isfinite(avg_distq_nutau)*(real_final_level_nutau==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(avg_distq_c7437[numpy.isfinite(avg_distq_c7437)*(real_final_level_c7437==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[real_final_level_c7437==1])
chist=pylab.hist(avg_distq_c9622[numpy.isfinite(avg_distq_c9622)*(real_final_level_c9622==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(avg_distq_nue[numpy.isfinite(avg_distq_nue)*(real_final_level_nue==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(avg_distq_nuge[numpy.isfinite(avg_distq_nuge)*(real_final_level_nuge==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,200,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,200,0.0,.7])
pylab.legend()
pylab.grid()
pylab.title('Avg DomDist Q')
pylab.xlabel('Distance (m)')
pylab.ylabel('Normed Counts')
pylab.savefig('LES_AvgDomDistQ_PreBDTCutL5_Rates_FinalLevelDist.png')

### Errors ###
#terr=pylab.hist(dhd_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(dhd_nue[real_final_level_nue==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(dhd_numu[real_final_level_numu==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(dhd_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############


pylab.figure(figsize=(10,8))
pylab.hist(dhd_data[real_final_level_data==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
muhist=pylab.hist(dhd_numu[real_final_level_numu==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#pylab.hist(dhd_c7437[real_final_level_c7437==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
chist=pylab.hist(dhd_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(dhd_nue[real_final_level_nue==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#thist=pylab.hist(dhd_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(dhd_nuge[real_final_level_nuge==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
#pylab.hist(dhd_nugmu[real_final_level_nugmu==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,20,21)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend(loc='upper right')
#pylab.vlines(x=3,ymin=10**-4,ymax=1,linestyle='dashed',lw=2)
#pylab.fill_between([0,3],10**-4,1,alpha=0.2)
#pylab.axis([0,10,10**-2,1])
pylab.grid()
pylab.title('DirectHitsD (SPE6)')
pylab.xlabel('DirectHits')
pylab.ylabel('Normed Counts')
pylab.savefig("LES_DirectHitsD_PreBDTCutL5_Rates_FinalLevelDist")

### Errors ###
#terr=pylab.hist(finite_z_nutau[real_final_level_nutau==1],bins=numpy.linspace(-600,200,41),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(finite_z_nue[real_final_level_nue==1],bins=numpy.linspace(-600,200,41),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(finite_z_numu[real_final_level_numu==1],bins=numpy.linspace(-600,200,41),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(finite_z_c9622[real_final_level_c9622==1],bins=numpy.linspace(-600,200,41),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(finite_z_data[real_final_level_data==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
muhist=pylab.hist(finite_z_numu[real_final_level_numu==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#pylab.hist(finite_z_c7437[real_final_level_c7437==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
chist=pylab.hist(finite_z_c9622[real_final_level_c9622==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(finite_z_nue[real_final_level_nue==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#thist=pylab.hist(finite_z_nutau[real_final_level_nutau==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(finite_z_nuge[real_final_level_nuge==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
#pylab.hist(finite_z_nugmu[real_final_level_nugmu==1],bins=numpy.linspace(-600,200,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-600,200,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend(loc='upper right')
pylab.axis([-600,200,0,0.2])
pylab.grid()
pylab.title('FiniteReco Z')
pylab.xlabel('Z (m)')
pylab.ylabel('Normed Counts')
pylab.savefig('LES_BDTParam_FiniteReco_Z_FinalLevelDist')

### Errors ###
#terr=pylab.hist(fr_R_nutau[real_final_level_nutau==1],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(fr_R_nue[real_final_level_nue==1],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(fr_R_numu[real_final_level_numu==1],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(fr_R_c9622[real_final_level_c9622==1],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(fr_R_data[numpy.isfinite(fr_R_data)*(real_final_level_data==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#pylab.hist(fr_R_nugmu[numpy.isfinite(fr_R_nugmu)*(real_final_level_nugmu==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
muhist=pylab.hist(fr_R_numu[numpy.isfinite(fr_R_numu)*(real_final_level_numu==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#thist=pylab.hist(fr_R_nutau[numpy.isfinite(fr_R_nutau)*(real_final_level_nutau==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(fr_R_c7437[numpy.isfinite(fr_R_c7437)*(real_final_level_c7437==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[real_final_level_c7437==1])
chist=pylab.hist(fr_R_c9622[numpy.isfinite(fr_R_c9622)*(real_final_level_c9622==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(fr_R_nue[numpy.isfinite(fr_R_nue)*(real_final_level_nue==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(fr_R_nuge[numpy.isfinite(fr_R_nuge)*(real_final_level_nuge==1)],bins=numpy.linspace(0,400,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,400,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.title('FiniteReco R (S36)')
pylab.xlabel('R (m)')
pylab.ylabel('mHz')
pylab.savefig('LES_BDTParam_FiniteReco_R_FinalLevelDist')



### Errors ###
#terr=pylab.hist(spline_rlogl_nutau[les_preselect_nutau==1],bins=numpy.linspace(0.0,25.0,41),log=False,normed=False,weights=(1000*atmo_nutau[les_preselect_nutau==1])**2)
eerr=pylab.hist(spline_rlogl_nue[les_preselect_nue==1],bins=numpy.linspace(0.0,25.0,41),log=False,normed=False,weights=(1000*atmo_nue[les_preselect_nue==1])**2)
merr=pylab.hist(spline_rlogl_numu[les_preselect_numu==1],bins=numpy.linspace(0.0,25.0,41),log=False,normed=False,weights=(1000*atmo_numu[les_preselect_numu==1])**2)
cerr=pylab.hist(spline_rlogl_c9622[les_preselect_c9622==1],bins=numpy.linspace(0.0,25.0,41),log=False,normed=False,weights=(1000*H3A_C9622EventWeight[les_preselect_c9622==1])**2)
hecerr=pylab.hist(spline_rlogl_c9255[les_preselect_c9255==1],bins=numpy.linspace(0.0,25.0,41),log=False,normed=False,weights=(1000*C9255EventWeight[les_preselect_c9255==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(spline_rlogl_data[numpy.isfinite(spline_rlogl_data)*(les_preselect_data==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[les_preselect_data==1])
#pylab.hist(spline_rlogl_nugmu[numpy.isfinite(spline_rlogl_nugmu)*(les_preselect_nugmu==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[les_preselect_nugmu==1],ls='dashed')
muhist=pylab.hist(spline_rlogl_numu[numpy.isfinite(spline_rlogl_numu)*(les_preselect_numu==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[les_preselect_numu==1])
#thist=pylab.hist(spline_rlogl_nutau[numpy.isfinite(spline_rlogl_nutau)*(les_preselect_nutau==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[les_preselect_nutau==1])
#pylab.hist(spline_rlogl_c7437[numpy.isfinite(spline_rlogl_c7437)*(les_preselect_c7437==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[les_preselect_c7437==1])
chist=pylab.hist(spline_rlogl_c9622[numpy.isfinite(spline_rlogl_c9622)*(les_preselect_c9622==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*H3A_C9622EventWeight[les_preselect_c9622==1])
hechist=pylab.hist(spline_rlogl_c9255[numpy.isfinite(spline_rlogl_c9255)*(les_preselect_c9255==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[les_preselect_c9255==1])
ehist=pylab.hist(spline_rlogl_nue[numpy.isfinite(spline_rlogl_nue)*(les_preselect_nue==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[les_preselect_nue==1])
#pylab.hist(spline_rlogl_nuge[numpy.isfinite(spline_rlogl_nuge)*(les_preselect_nuge==1)],bins=numpy.linspace(0.0,25.0,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[les_preselect_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]+hechist[0]
errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
binzo = numpy.linspace(0.0,25.0,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0.0,25.0,0.0,0.5])
pylab.legend()
pylab.grid()
pylab.title('SplineMPEMod RLogL')
pylab.xlabel('rllh')
pylab.ylabel('mHz')
pylab.savefig('LES_BDTParam_SplineMPEMod_RLogL_PreCutDist_WC9255')

### Errors ###
#terr=pylab.hist(dhd_nutau[les_preselect_nutau==1],bins=numpy.linspace(0,20.0,21),log=False,normed=False,weights=(1000*atmo_nutau[les_preselect_nutau==1])**2)
eerr=pylab.hist(dhd_nue[les_preselect_nue==1],bins=numpy.linspace(0,20.0,21),log=False,normed=False,weights=(1000*atmo_nue[les_preselect_nue==1])**2)
merr=pylab.hist(dhd_numu[les_preselect_numu==1],bins=numpy.linspace(0,20.0,21),log=False,normed=False,weights=(1000*atmo_numu[les_preselect_numu==1])**2)
cerr=pylab.hist(dhd_c9622[les_preselect_c9622==1],bins=numpy.linspace(0,20.0,21),log=False,normed=False,weights=(1000*H3A_C9622EventWeight[les_preselect_c9622==1])**2)
hecerr=pylab.hist(dhd_c9255[les_preselect_c9255==1],bins=numpy.linspace(0,20.0,21),log=False,normed=False,weights=(1000*C9255EventWeight[les_preselect_c9255==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(dhd_data[numpy.isfinite(dhd_data)*(les_preselect_data==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[les_preselect_data==1])
#pylab.hist(dhd_nugmu[numpy.isfinite(dhd_nugmu)*(les_preselect_nugmu==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[les_preselect_nugmu==1],ls='dashed')
muhist=pylab.hist(dhd_numu[numpy.isfinite(dhd_numu)*(les_preselect_numu==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[les_preselect_numu==1])
#thist=pylab.hist(dhd_nutau[numpy.isfinite(dhd_nutau)*(les_preselect_nutau==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[les_preselect_nutau==1])
#pylab.hist(dhd_c7437[numpy.isfinite(dhd_c7437)*(les_preselect_c7437==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[les_preselect_c7437==1])
chist=pylab.hist(dhd_c9622[numpy.isfinite(dhd_c9622)*(les_preselect_c9622==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*H3A_C9622EventWeight[les_preselect_c9622==1])
hechist=pylab.hist(dhd_c9255[numpy.isfinite(dhd_c9255)*(les_preselect_c9255==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[les_preselect_c9255==1])
ehist=pylab.hist(dhd_nue[numpy.isfinite(dhd_nue)*(les_preselect_nue==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[les_preselect_nue==1])
#pylab.hist(dhd_nuge[numpy.isfinite(dhd_nuge)*(les_preselect_nuge==1)],bins=numpy.linspace(0,20.0,21),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[les_preselect_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,20.0,21)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0,20.0,0.0,0.4])
pylab.legend()
pylab.grid()
pylab.title('SplineMPE DirectHits D')
pylab.xlabel('DHD')
pylab.ylabel('mHz')
pylab.savefig('LES_BDTParam_SplineMPE_DHD_FinalLevelDist_WC9255')


### Errors ###
#terr=pylab.hist(spline_plogl_nutau[les_preselect_nutau==1],bins=numpy.linspace(2.0,15,41),log=False,normed=False,weights=(1000*atmo_nutau[les_preselect_nutau==1])**2)
eerr=pylab.hist(spline_plogl_nue[les_preselect_nue==1],bins=numpy.linspace(2.0,15,41),log=False,normed=False,weights=(1000*atmo_nue[les_preselect_nue==1])**2)
merr=pylab.hist(spline_plogl_numu[les_preselect_numu==1],bins=numpy.linspace(2.0,15,41),log=False,normed=False,weights=(1000*atmo_numu[les_preselect_numu==1])**2)
cerr=pylab.hist(spline_plogl_c9622[les_preselect_c9622==1],bins=numpy.linspace(2.0,15,41),log=False,normed=False,weights=(1000*C9622EventWeight[les_preselect_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(spline_plogl_data[numpy.isfinite(spline_plogl_data)*(les_preselect_data==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[les_preselect_data==1])
#pylab.hist(spline_plogl_nugmu[numpy.isfinite(spline_plogl_nugmu)*(les_preselect_nugmu==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[les_preselect_nugmu==1],ls='dashed')
muhist=pylab.hist(spline_plogl_numu[numpy.isfinite(spline_plogl_numu)*(les_preselect_numu==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[les_preselect_numu==1])
#thist=pylab.hist(spline_plogl_nutau[numpy.isfinite(spline_plogl_nutau)*(les_preselect_nutau==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[les_preselect_nutau==1])
#pylab.hist(spline_plogl_c7437[numpy.isfinite(spline_plogl_c7437)*(les_preselect_c7437==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[les_preselect_c7437==1])
chist=pylab.hist(spline_plogl_c9622[numpy.isfinite(spline_plogl_c9622)*(les_preselect_c9622==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[les_preselect_c9622==1])
ehist=pylab.hist(spline_plogl_nue[numpy.isfinite(spline_plogl_nue)*(les_preselect_nue==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[les_preselect_nue==1])
#pylab.hist(spline_plogl_nuge[numpy.isfinite(spline_plogl_nuge)*(les_preselect_nuge==1)],bins=numpy.linspace(2.0,15,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[les_preselect_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(2.0,15,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([2.0,15,0.0,0.5])
pylab.legend()
pylab.grid()
pylab.title('SplineMPEMod PLogL')
pylab.xlabel('rllh')
pylab.ylabel('mHz')
pylab.savefig('LES_BDTParam_SplineMPEMod_PLogL_PreBDTCutL5_Rates')




### Errors ###
#terr=pylab.hist(57.3*splinemod_azi_nutau[real_final_level_nutau==1],bins=numpy.linspace(0.0,360.0,20),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(57.3*splinemod_azi_nue[real_final_level_nue==1],bins=numpy.linspace(0.0,360.0,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(57.3*splinemod_azi_numu[real_final_level_numu==1],bins=numpy.linspace(0.0,360.0,20),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(57.3*splinemod_azi_c9622[real_final_level_c9622==1],bins=numpy.linspace(0.0,360.0,20),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(57.3*splinemod_azi_data[(real_final_level_data==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#pylab.hist(57.3*splinemod_azi_nugmu[numpy.isfinite(57.3*splinemod_azi_nugmu)*(real_final_level_nugmu==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
muhist=pylab.hist(57.3*splinemod_azi_numu[(real_final_level_numu==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#thist=pylab.hist(57.3*splinemod_azi_nutau[(real_final_level_nutau==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(57.3*splinemod_azi_c7437[numpy.isfinite(57.3*splinemod_azi_c7437)*(real_final_level_c7437==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[real_final_level_c7437==1])
chist=pylab.hist(57.3*splinemod_azi_c9622[(real_final_level_c9622==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(57.3*splinemod_azi_nue[(real_final_level_nue==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(57.3*splinemod_azi_nuge[(real_final_level_nuge==1)],bins=numpy.linspace(0.0,360.0,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0.0,360.0,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0.0,360,0.0,0.07])
pylab.legend()
pylab.grid()
pylab.title('FinalLevel SplineMPE Azimuth Distribution')
pylab.xlabel(r'$\phi$')
pylab.ylabel('mHz')
pylab.savefig('LES_SplineMPE_AzimuthDist_FinalLevelDist')

### Errors ###
#terr=pylab.hist(numpy.cos(splinemod_zen_nutau[real_final_level_nutau==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nutau[real_final_level_nutau==1])**2)
eerr=pylab.hist(numpy.cos(splinemod_zen_nue[real_final_level_nue==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nue[real_final_level_nue==1])**2)
merr=pylab.hist(numpy.cos(splinemod_zen_numu[real_final_level_numu==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_numu[real_final_level_numu==1])**2)
cerr=pylab.hist(numpy.cos(splinemod_zen_c9622[real_final_level_c9622==1]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*C9622EventWeight[real_final_level_c9622==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(numpy.cos(splinemod_zen_data[(real_final_level_data==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[real_final_level_data==1])
#pylab.hist(numpy.cos(splinemod_zen_nugmu[numpy.isfinite(numpy.cos(splinemod_zen_nugmu)*(real_final_level_nugmu==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[real_final_level_nugmu==1],ls='dashed')
muhist=pylab.hist(numpy.cos(splinemod_zen_numu[(real_final_level_numu==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[real_final_level_numu==1])
#thist=pylab.hist(numpy.cos(splinemod_zen_nutau[(real_final_level_nutau==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[real_final_level_nutau==1])
#pylab.hist(numpy.cos(splinemod_zen_c7437[numpy.isfinite(numpy.cos(splinemod_zen_c7437)*(real_final_level_c7437==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[real_final_level_c7437==1])
chist=pylab.hist(numpy.cos(splinemod_zen_c9622[(real_final_level_c9622==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[real_final_level_c9622==1])
ehist=pylab.hist(numpy.cos(splinemod_zen_nue[(real_final_level_nue==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nue[real_final_level_nue==1])
#pylab.hist(numpy.cos(splinemod_zen_nuge[(real_final_level_nuge==1)]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[real_final_level_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-1,0.087,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1,0.087,0.0,0.07])
pylab.legend()
pylab.grid()
pylab.title('FinalLevel SplineMPE Zenith Distribution')
pylab.xlabel(r'$Cos(\Theta_{Zenith})$')
pylab.ylabel('mHz')
pylab.savefig('LES_SplineMPE_ZenithDist_FinalLevelDist')


                                                      

pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=16.0
pylab.rcParams['axes.titlesize']=18.0
pylab.rcParams['ytick.labelsize']='large'
pylab.rcParams['xtick.labelsize']='large'
pylab.rcParams['lines.markeredgewidth']=1.0

### LES ###

### Errors ###
#terr=pylab.hist(bdt_les_nutau[les_preselect_nutau==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nutau[les_preselect_nutau==1])**2)
eerr=pylab.hist(bdt_les_nue[les_preselect_nue==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nue[les_preselect_nue==1])**2)
merr=pylab.hist(bdt_les_numu[les_preselect_numu==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_numu[les_preselect_numu==1])**2)
cerr=pylab.hist(bdt_les_c9622[les_preselect_c9622==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*H3A_C9622EventWeight[les_preselect_c9622==1])**2)
hecerr=pylab.hist(bdt_les_c9255[les_preselect_c9255==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*C9255EventWeight[les_preselect_c9255==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(bdt_les_data[les_preselect_data==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[les_preselect_data==1])
#thist=pylab.hist(bdt_les_nutau[les_preselect_nutau==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[les_preselect_nutau==1])
#pylab.hist(bdt_les_nugmu[les_preselect_nugmu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[les_preselect_nugmu==1])
chist=pylab.hist(bdt_les_c9622[les_preselect_c9622==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*H3A_C9622EventWeight[les_preselect_c9622==1])
hechist=pylab.hist(bdt_les_c9255[les_preselect_c9255==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[les_preselect_c9255==1])
muhist=pylab.hist(bdt_les_numu[les_preselect_numu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[les_preselect_numu==1])
#sig1muhist=pylab.hist(bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3}$',ec='b',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-0.5)
#sig2muhist=pylab.hist(bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.5}$',ec='orange',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-1.0)
#sig3muhist=pylab.hist(bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.75}$',ec='green',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-1.25)
ehist=pylab.hist(bdt_les_nue[les_preselect_nue==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[les_preselect_nue==1])
#pylab.hist(bdt_les_nuge[les_preselect_nuge==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[les_preselect_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]+hechist[0]
errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
binzo = numpy.linspace(-1,1,40)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1.0,1.0,0.0,0.3])
pylab.legend()
pylab.grid()
pylab.title("LES BDT Score Distribution (Nch>15)")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("LES_BDTScoreDist_NormalizedRates_WC9255_And_H3AC9622_L6_NchCut_G1460")

### Errors ###
#terr=pylab.hist(alt_bdt_les_nutau[les_preselect_nutau==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nutau[les_preselect_nutau==1])**2)
eerr=pylab.hist(alt_bdt_les_nue[les_preselect_nue==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nue[les_preselect_nue==1])**2)
merr=pylab.hist(alt_bdt_les_numu[les_preselect_numu==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_numu[les_preselect_numu==1])**2)
cerr=pylab.hist(alt_bdt_les_c9622[les_preselect_c9622==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*H3A_C9622EventWeight[les_preselect_c9622==1])**2)
hecerr=pylab.hist(bdt_les_c9255[les_preselect_c9255==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*C9255EventWeight[les_preselect_c9255==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(alt_bdt_les_data[les_preselect_data==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[les_preselect_data==1])
#thist=pylab.hist(alt_bdt_les_nutau[les_preselect_nutau==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[les_preselect_nutau==1])
#pylab.hist(alt_bdt_les_nugmu[les_preselect_nugmu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[les_preselect_nugmu==1])
chist=pylab.hist(alt_bdt_les_c9622[les_preselect_c9622==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*H3A_C9622EventWeight[les_preselect_c9622==1])
hechist=pylab.hist(bdt_les_c9255[les_preselect_c9255==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[les_preselect_c9255==1])
muhist=pylab.hist(alt_bdt_les_numu[les_preselect_numu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[les_preselect_numu==1])
#sig1muhist=pylab.hist(alt_bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3}$',ec='b',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-0.5)
#sig2muhist=pylab.hist(alt_bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.5}$',ec='orange',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-1.0)
#sig3muhist=pylab.hist(alt_bdt_les_numu[(les_preselect_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.75}$',ec='green',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[(les_preselect_numu==1)*(good_angular_numu==1)])**-1.25)
ehist=pylab.hist(alt_bdt_les_nue[les_preselect_nue==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[les_preselect_nue==1])
#pylab.hist(alt_bdt_les_nuge[les_preselect_nuge==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[les_preselect_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]+hechist[0]
errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
binzo = numpy.linspace(-1,1,40)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1.0,1.0,0.0,0.3])
pylab.legend()
pylab.grid()
pylab.title("LES BDT Score Distribution")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("LES_BDT_NoRecoParam_ScoreDist_WC9255_And_H3AC9622_Rates_L6")


### HES ###

eerr=pylab.hist(bdt_hes_nuge[hes_preselect_nuge==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nuge[hes_preselect_nuge==1])**2)
merr=pylab.hist(bdt_hes_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nugmu[hes_preselect_nugmu==1])**2)
cerr=pylab.hist(bdt_hes_c7437[hes_preselect_c7437==1],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*C7437EventWeight[hes_preselect_c7437==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(bdt_hes_data[hes_preselect_data==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[hes_preselect_data==1])
#thist=pylab.hist(bdt_hes_nutau[hes_preselect_nutau==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[hes_preselect_nutau==1])
muhist=pylab.hist(bdt_hes_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[hes_preselect_nugmu==1])
chist=pylab.hist(bdt_hes_c7437[hes_preselect_c7437==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[hes_preselect_c7437==1])
#muhist=pylab.hist(bdt_hes_numu[hes_preselect_numu==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[hes_preselect_numu==1])
#ehist=pylab.hist(bdt_hes_nue[hes_preselect_nue==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Nue',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[hes_preselect_nue==1])
ehist=pylab.hist(bdt_hes_nuge[hes_preselect_nuge==1],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[hes_preselect_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(-1,1,40)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-1.0,1.0,0.0,.4])
pylab.legend()
pylab.grid()
pylab.title("HES BDT Score Distribution")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("HES_BDTScoreDist_Rates_L6")


### Errors ###
#terr=pylab.hist(bdt_les_nutau[(les_nutau==1)*(real_final_level_nutau==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nutau[(les_)*(real_final_level_nutau==1])**2)
eerr=pylab.hist(bdt_les_nue[(les_nue==1)*(real_final_level_nue==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_nue[(les_nue==1)*(real_final_level_nue==1)])**2)
merr=pylab.hist(bdt_les_numu[(les_numu==1)*(real_final_level_numu==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*atmo_numu[(les_numu==1)*(real_final_level_numu==1)])**2)
cerr=pylab.hist(bdt_les_c9622[(les_c9622==1)*(real_final_level_c9622==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*H3A_C9622EventWeight[(les_c9622==1)*(real_final_level_c9622==1)])**2)
hecerr=pylab.hist(bdt_les_c9255[(les_c9255==1)*(real_final_level_c9255==1)],bins=numpy.linspace(-1,1,40),log=False,normed=False,weights=(1000*C9255EventWeight[(les_c9255==1)*(real_final_level_c9255==1)])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(bdt_les_data[(les_data==1)*(real_final_level_data==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(les_data==1)*(real_final_level_data==1)])
#thist=pylab.hist(bdt_les_nutau[(les_nutau)*(real_final_level_nutau==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,ls='solid',weights=1000*atmo_nutau[(les_nutau)*(real_final_level_nutau==1])
#pylab.hist(bdt_les_nugmu[(les_)*(real_final_level_nugmu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[(les_)*(real_final_level_nugmu==1])
chist=pylab.hist(bdt_les_c9622[(les_c9622==1)*(real_final_level_c9622==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*H3A_C9622EventWeight[(les_c9622==1)*(real_final_level_c9622==1)])
hechist=pylab.hist(bdt_les_c9255[(les_c9255==1)*(real_final_level_c9255==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='C9255',ec='magenta',normed=False,weights=1000*C9255EventWeight[(les_c9255==1)*(real_final_level_c9255==1)])
muhist=pylab.hist(bdt_les_numu[(les_numu==1)*(real_final_level_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu}$',ec='b',normed=False,ls='solid',weights=1000*atmo_numu[(les_numu==1)*(real_final_level_numu==1)])
#sig1muhist=pylab.hist(bdt_les_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3}$',ec='b',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)])**-0.5)
#sig2muhist=pylab.hist(bdt_les_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.5}$',ec='orange',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)])**-1.0)
#sig3muhist=pylab.hist(bdt_les_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_{\mu} E^{-3.75}$',ec='green',normed=False,ls='dashed',weights=0.0002*(primary_energy_numu[((les_)*(real_final_level_numu==1)*(good_angular_numu==1)])**-1.25)
ehist=pylab.hist(bdt_les_nue[(les_nue==1)*(real_final_level_nue==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label=r'$\nu_e$',ec='g',normed=False,ls='solid',weights=1000*atmo_nue[(les_nue==1)*(real_final_level_nue==1)])
#pylab.hist(bdt_les_nuge[(les_)*(real_final_level_nuge==1)],bins=numpy.linspace(-1,1,40),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[(les_)*(real_final_level_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]+hechist[0]
errors = (eerr[0]+merr[0]+cerr[0]+hecerr[0])**0.5
binzo = numpy.linspace(-1,1,40)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([-0.2,1.0,0.0,0.2])
pylab.legend()
pylab.grid()
pylab.title("Final Level LES BDT Score Distribution")
pylab.xlabel('Score')
pylab.ylabel('mHz')
pylab.savefig("LES_PostCutBDTScoreDist_NormalizedRates_WC9255_And_H3AC9622_L6_NchCut_G1460")



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


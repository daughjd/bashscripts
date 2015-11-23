pylab.rcParams['font.size'] = 12.0
pylab.rcParams['axes.labelsize']=16.0
pylab.rcParams['axes.titlesize']=18.0
pylab.rcParams['ytick.labelsize']='large'
pylab.rcParams['xtick.labelsize']='large'
pylab.rcParams['lines.markeredgewidth']=1.0

### Errors ###
eerr=pylab.hist(ldird_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0,650,31),log=False,normed=False,weights=(1000*atmo_nuge[hes_preselect_nuge==1])**2)
merr=pylab.hist(ldird_numu[hes_preselect_numu==1],bins=numpy.linspace(0,650,31),log=False,normed=False,weights=(1000*atmo_numu[hes_preselect_numu==1])**2)
cerr=pylab.hist(ldird_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0,650,31),log=False,normed=False,weights=(1000*C7437EventWeight[hes_preselect_c7437==1])**2)
##############

pylab.figure(figsize=(10,8))
pylab.hist(ldird_data[hes_preselect_data==1],bins=numpy.linspace(0,650,31),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[hes_preselect_data==1])
muhist=pylab.hist(ldird_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0,650,31),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,ls='dashed',weights=1000*atmo_nugmu[hes_preselect_nugmu==1])
chist=pylab.hist(ldird_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0,650,31),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[hes_preselect_c7437==1])
ehist=pylab.hist(ldird_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0,650,31),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,ls='dashed',weights=1000*atmo_nuge[hes_preselect_nuge==1])
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,650,31)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.title("LDirD SplineMPE")
pylab.xlabel('Length (m)')
pylab.ylabel('mHz')
pylab.savefig("HES_LDirD_L5PreSelected_Rates")

### Errors ###
eerr=pylab.hist(avg_distq_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*atmo_nuge[hes_preselect_nuge==1])**2)
merr=pylab.hist(avg_distq_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*atmo_nugmu[hes_preselect_nugmu==1])**2)
cerr=pylab.hist(avg_distq_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0,200,20),log=False,normed=False,weights=(1000*C7437EventWeight[hes_preselect_c7437==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(avg_distq_data[numpy.isfinite(avg_distq_data)*(hes_preselect_data==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[hes_preselect_data==1])
muhist=pylab.hist(avg_distq_nugmu[numpy.isfinite(avg_distq_nugmu)*(hes_preselect_nugmu==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[hes_preselect_nugmu==1],ls='dashed')
chist=pylab.hist(avg_distq_c7437[numpy.isfinite(avg_distq_c7437)*(hes_preselect_c7437==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[hes_preselect_c7437==1])
ehist=pylab.hist(avg_distq_nuge[numpy.isfinite(avg_distq_nuge)*(hes_preselect_nuge==1)],bins=numpy.linspace(0,200,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[hes_preselect_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,200,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.title('Avg DomDist Q')
pylab.xlabel('Distance (m)')
pylab.ylabel('Normed Counts')
pylab.savefig('HES_AvgDomDistQ_L5PreSelected_Rates.png')

### Errors ###
eerr=pylab.hist(dhd_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*atmo_nuge[hes_preselect_nuge==1])**2)
merr=pylab.hist(dhd_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*atmo_nugmu[hes_preselect_nugmu==1])**2)
cerr=pylab.hist(dhd_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0,20,21),log=True,normed=False,weights=(1000*C7437EventWeight[hes_preselect_c7437==1])**2)
#############


pylab.figure(figsize=(10,8))
pylab.hist(dhd_data[hes_preselect_data==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[hes_preselect_data==1])
#muhist=pylab.hist(dhd_numu[hes_preselect_numu==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_numu[hes_preselect_numu==1])
chist=pylab.hist(dhd_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[hes_preselect_c7437==1])
ehist=pylab.hist(dhd_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[hes_preselect_nuge==1],ls='dashed')
muhist=pylab.hist(dhd_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0,20,21),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[hes_preselect_nugmu==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0,20,21)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend(loc='upper right')
#pylab.vlines(x=3,ymin=10**-4,ymax=1,linestyle='dashed',lw=2)
#pylab.fill_between([0,3],10**-4,1,alpha=0.2)
pylab.grid()
pylab.title('DirectHitsD (SplineMPEMod)')
pylab.xlabel('DirectHits')
pylab.ylabel('Normed Counts')
pylab.savefig("HES_DirectHitsD_L5PreSelected_Rates")

### Errors ###
eerr=pylab.hist(57.3*lf_spline_diff_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0.0,30,20),log=False,normed=False,weights=(1000*atmo_nuge[hes_preselect_nuge==1])**2)
merr=pylab.hist(57.3*lf_spline_diff_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0.0,30,20),log=False,normed=False,weights=(1000*atmo_nugmu[hes_preselect_nugmu==1])**2)
cerr=pylab.hist(57.3*lf_spline_diff_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0.0,30,20),log=False,normed=False,weights=(1000*C7437EventWeight[hes_preselect_c7437==1])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(57.3*lf_spline_diff_data[hes_preselect_data==1],bins=numpy.linspace(0.0,30,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[hes_preselect_data==1])
muhist=pylab.hist(57.3*lf_spline_diff_nugmu[hes_preselect_nugmu==1],bins=numpy.linspace(0.0,30,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[hes_preselect_nugmu==1],ls='dashed')
chist=pylab.hist(57.3*lf_spline_diff_c7437[hes_preselect_c7437==1],bins=numpy.linspace(0.0,30,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[hes_preselect_c7437==1])
ehist=pylab.hist(57.3*lf_spline_diff_nuge[hes_preselect_nuge==1],bins=numpy.linspace(0.0,30,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[hes_preselect_nuge==1],ls='dashed')
mctot=ehist[0]+muhist[0]+chist[0]
errors = (eerr[0]+merr[0]+cerr[0])**0.5
binzo = numpy.linspace(0.0,30,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.axis([0.0,30,0.00,0.24])
pylab.legend(loc='upper right')
pylab.grid()
pylab.title('Space Angle LF-SplineMPE')
pylab.xlabel('Space Angle (Degrees)')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_LF_SplineMPE_Diff_L5PreSelected_Rates')

### Errors ###
#terr=pylab.hist(fr_R_nutau[(real_final_level_nutau)*(hes_==1)],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*atmo_nutau[(real_final_level_nutau)*(hes_==1)])**2)
eerr=pylab.hist(fr_R_nuge[(real_final_level_nuge)*(hes_nuge==1)],bins=numpy.linspace(0,600,41),log=False,normed=False,weights=(1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])**2)
merr=pylab.hist(fr_R_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)],bins=numpy.linspace(0,600,41),log=False,normed=False,weights=(1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])**2)
#cerr=pylab.hist(fr_R_c9622[(real_final_level_c9622)*(hes_==1)],bins=numpy.linspace(0,600,41),log=False,normed=False,weights=(1000*C9622EventWeight[(real_final_level_c9622)*(hes_==1)])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(fr_R_data[numpy.isfinite(fr_R_data)*((real_final_level_data)*(hes_data==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(real_final_level_data)*(hes_data==1)])
#pylab.hist(fr_R_nugmu[numpy.isfinite(fr_R_nugmu)*((real_final_level_nugmu)*(hes_==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_==1)],ls='dashed')
muhist=pylab.hist(fr_R_nugmu[numpy.isfinite(fr_R_nugmu)*((real_final_level_nugmu)*(hes_nugmu==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])
#thist=pylab.hist(fr_R_nutau[numpy.isfinite(fr_R_nutau)*((real_final_level_nutau)*(hes_==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[(real_final_level_nutau)*(hes_==1)])
#pylab.hist(fr_R_c7437[numpy.isfinite(fr_R_c7437)*((real_final_level_c7437)*(hes_==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[(real_final_level_c7437)*(hes_==1)])
#chist=pylab.hist(fr_R_c9622[numpy.isfinite(fr_R_c9622)*((real_final_level_c9622)*(hes_==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[(real_final_level_c9622)*(hes_==1)])
ehist=pylab.hist(fr_R_nuge[numpy.isfinite(fr_R_nuge)*((real_final_level_nuge)*(hes_nuge==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])
#pylab.hist(fr_R_nuge[numpy.isfinite(fr_R_nuge)*((real_final_level_nuge)*(hes_==1))],bins=numpy.linspace(0,600,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_==1)],ls='dashed')
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(0,600,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend(loc="upper left")
pylab.grid()
pylab.axis([0,600,0.0,0.02])
pylab.title('FiniteReco R (S36)')
pylab.xlabel('R (m)')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_FiniteReco_R_L5PreSelected_Rates')

### Errors ###
#terr=pylab.hist(finite_z_nutau[(real_final_level_nutau)*(hes_==1)],bins=numpy.linspace(0,400,41),log=False,normed=False,weights=(1000*atmo_nutau[(real_final_level_nutau)*(hes_==1)])**2)
eerr=pylab.hist(finite_z_nuge[(real_final_level_nuge)*(hes_nuge==1)],bins=numpy.linspace(-600,0,41),log=False,normed=False,weights=(1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])**2)
merr=pylab.hist(finite_z_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)],bins=numpy.linspace(-600,0,41),log=False,normed=False,weights=(1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])**2)
#cerr=pylab.hist(finite_z_c9622[(real_final_level_c9622)*(hes_==1)],bins=numpy.linspace(-600,0,41),log=False,normed=False,weights=(1000*C9622EventWeight[(real_final_level_c9622)*(hes_==1)])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(finite_z_data[numpy.isfinite(finite_z_data)*((real_final_level_data)*(hes_data==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(real_final_level_data)*(hes_data==1)])
#pylab.hist(finite_z_nugmu[numpy.isfinite(finite_z_nugmu)*((real_final_level_nugmu)*(hes_==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_==1)],ls='dashed')
muhist=pylab.hist(finite_z_nugmu[numpy.isfinite(finite_z_nugmu)*((real_final_level_nugmu)*(hes_nugmu==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='NuMu',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])
#thist=pylab.hist(finite_z_nutau[numpy.isfinite(finite_z_nutau)*((real_final_level_nutau)*(hes_==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='NuTau',ec='orange',normed=False,weights=1000*atmo_nutau[(real_final_level_nutau)*(hes_==1)])
#pylab.hist(finite_z_c7437[numpy.isfinite(finite_z_c7437)*((real_final_level_c7437)*(hes_==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*atmo_c7437[(real_final_level_c7437)*(hes_==1)])
#chist=pylab.hist(finite_z_c9622[numpy.isfinite(finite_z_c9622)*((real_final_level_c9622)*(hes_==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='C9622',ec='r',normed=False,weights=1000*C9622EventWeight[(real_final_level_c9622)*(hes_==1)])
ehist=pylab.hist(finite_z_nuge[numpy.isfinite(finite_z_nuge)*((real_final_level_nuge)*(hes_nuge==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='NuE',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])
#pylab.hist(finite_z_nuge[numpy.isfinite(finite_z_nuge)*((real_final_level_nuge)*(hes_==1))],bins=numpy.linspace(-600,0,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_==1)],ls='dashed')
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(-600,0,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend(loc="upper right")
pylab.grid()
pylab.axis([-600,0,0.0,0.02])
pylab.title('FiniteReco Z')
pylab.xlabel('Z (m)')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_FiniteReco_Z_L5PreSelected_Rates')



### Errors ###
eerr=pylab.hist(spline_rlogl_nuge[(real_final_level_nuge)*(hes_nuge==1)],bins=numpy.linspace(0.0,24,41),log=False,normed=False,weights=(1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])**2)
merr=pylab.hist(spline_rlogl_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)],bins=numpy.linspace(0.0,24,41),log=False,normed=False,weights=(1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])**2)
#cerr=pylab.hist(spline_rlogl_c7437[(real_final_level_c7437)*(hes_c7437==1)],bins=numpy.linspace(0.0,24,41),log=False,normed=False,weights=(1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(spline_rlogl_data[numpy.isfinite(spline_rlogl_data)*((real_final_level_data)*(hes_data==1))],bins=numpy.linspace(0.0,24,41),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(real_final_level_data)*(hes_data==1)])
muhist=pylab.hist(spline_rlogl_nugmu[numpy.isfinite(spline_rlogl_nugmu)*((real_final_level_nugmu)*(hes_nugmu==1))],bins=numpy.linspace(0.0,24,41),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])
#chist=pylab.hist(spline_rlogl_c7437[numpy.isfinite(spline_rlogl_c7437)*((real_final_level_c7437)*(hes_c7437==1))],bins=numpy.linspace(0.0,24,41),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])
ehist=pylab.hist(spline_rlogl_nuge[numpy.isfinite(spline_rlogl_nuge)*((real_final_level_nuge)*(hes_nuge==1))],bins=numpy.linspace(0.0,24,41),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(0.0,24,41)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.title('SplineMPEMod RLogL')
pylab.xlabel('rllh')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_SplineMPEMod_RLogL_FinalLevel_Rates')

### Errors ###
eerr=pylab.hist(57.3*splinemod_azi_nuge[(real_final_level_nuge)*(hes_nuge==1)],bins=numpy.linspace(0.,360.,20),log=False,normed=False,weights=(1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])**2)
merr=pylab.hist(57.3*splinemod_azi_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)],bins=numpy.linspace(0.,360.,20),log=False,normed=False,weights=(1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])**2)
#cerr=pylab.hist(57.3*splinemod_azi_c7437[(real_final_level_c7437)*(hes_c7437==1)],bins=numpy.linspace(0.,360.,20),log=False,normed=False,weights=(1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(57.3*splinemod_azi_data[numpy.isfinite(57.3*splinemod_azi_data)*((real_final_level_data)*(hes_data==1))],bins=numpy.linspace(0.,360.,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(real_final_level_data)*(hes_data==1)])
muhist=pylab.hist(57.3*splinemod_azi_nugmu[numpy.isfinite(57.3*splinemod_azi_nugmu)*((real_final_level_nugmu)*(hes_nugmu==1))],bins=numpy.linspace(0.,360.,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])
#chist=pylab.hist(57.3*splinemod_azi_c7437[numpy.isfinite(57.3*splinemod_azi_c7437)*((real_final_level_c7437)*(hes_c7437==1))],bins=numpy.linspace(0.,360.,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])
ehist=pylab.hist(57.3*splinemod_azi_nuge[numpy.isfinite(57.3*splinemod_azi_nuge)*((real_final_level_nuge)*(hes_nuge==1))],bins=numpy.linspace(0.,360.,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(0.,360.,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.axis([0.0,360,0.0,0.025])
pylab.title('FinalLevel SplineMPE Azimuth Distribution')
pylab.xlabel(r'$\phi$')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_SplineMPEMod_AzimuthDist_FinalLevel_Rates')

### Errors ###
eerr=pylab.hist(numpy.cos(splinemod_zen_nuge[(real_final_level_nuge)*(hes_nuge==1)]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])**2)
merr=pylab.hist(numpy.cos(splinemod_zen_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])**2)
#cerr=pylab.hist(numpy.cos(splinemod_zen_c7437[(real_final_level_c7437)*(hes_c7437==1)]),bins=numpy.linspace(-1,0.087,20),log=False,normed=False,weights=(1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])**2)
#############

pylab.figure(figsize=(10,8))
pylab.hist(numpy.cos(splinemod_zen_data[numpy.isfinite(splinemod_zen_data)*((real_final_level_data)*(hes_data==1))]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='Data',ec='k',normed=False,weights=1000*dataweight[(real_final_level_data)*(hes_data==1)])
muhist=pylab.hist(numpy.cos(splinemod_zen_nugmu[numpy.isfinite(splinemod_zen_nugmu)*((real_final_level_nugmu)*(hes_nugmu==1))]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',normed=False,weights=1000*atmo_nugmu[(real_final_level_nugmu)*(hes_nugmu==1)])
#chist=pylab.hist(numpy.cos(splinemod_zen_c7437[numpy.isfinite(splinemod_zen_c7437)*((real_final_level_c7437)*(hes_c7437==1))]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='C7437',ec='r',normed=False,weights=1000*C7437EventWeight[(real_final_level_c7437)*(hes_c7437==1)])
ehist=pylab.hist(numpy.cos(splinemod_zen_nuge[numpy.isfinite(splinemod_zen_nuge)*((real_final_level_nuge)*(hes_nuge==1))]),bins=numpy.linspace(-1,0.087,20),histtype='step',lw=2,log=False,label='NuE(N)',ec='g',normed=False,weights=1000*atmo_nuge[(real_final_level_nuge)*(hes_nuge==1)])
mctot=ehist[0]+muhist[0]
errors = (eerr[0]+merr[0])**0.5
binzo = numpy.linspace(-1,0.087,20)[:-1]
binzo = binzo+((binzo[1]-binzo[0])/2.)
pylab.errorbar(binzo,mctot,yerr=errors,marker='o',color="purple",ms=8,label="MC Total",ls='none')
pylab.legend()
pylab.grid()
pylab.axis([-1.0,0.087,0.0,0.025])
pylab.title('FinalLevel SplineMPE Zenith Distribution')
pylab.xlabel(r'$Cos(\Theta_{Zenith})$')
pylab.ylabel('mHz')
pylab.savefig('HES_BDTParam_SplineMPEMod_ZenithDist_FinalLevel_Rates')


### Final Sample 1 (Spline MPE RLogL Straight Cut) ###
pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(sig_para_numu<(50/57.3))]
y1=57.3*sig_para_numu[(sig_para_numu<(50/57.3))*(merged_nch_numu<78)*(samp1_final_level_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(80,68),weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(sig_para_numu<(50/57.3))])
extent = [10,78,0,50]
#pylab.plot(xval,yval,'w--',lw=2.0)
pylab.xlabel('Nch')
pylab.ylabel('Paraboloid Sigma (Deg)')
pylab.title('Nch vs. Paraboloid (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_Paraboloid_NuMu')


### Long Range ###
pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp1_final_level_numu==1)*(merged_nch_numu<142)*(para_pull_numu<13)]
y1=para_pull_numu[(merged_nch_numu<142)*(samp1_final_level_numu==1)*(para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,132))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<142)*(para_pull_numu<12)])
extent = [10,142,0,13]
pylab.plot(samp1_pull_nch_numu,samp1_median_pull_numu,'wo',ms=8,label='Median Pull')
pylab.plot(samp1_pull_nch_numu,pullfit_numu(samp1_pull_nch_numu,samp1_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_ParaboloidPull_NuMu_UnWeighted_MuonError_ExtendedNch')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(para_pull_numu<13)]
y1=para_pull_numu[(merged_nch_numu<78)*(samp1_final_level_numu==1)*(para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp1_pull_nch_numu,samp1_median_pull_numu,'wo',ms=8,label='Median Pull')
pylab.plot(samp1_pull_nch_numu,pullfit_numu(samp1_pull_nch_numu,samp1_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_ParaboloidPull_NuMu_UnWeighted_MuonError')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(corrected_samp1_para_pull_numu<13)]
y1=corrected_samp1_para_pull_numu[(merged_nch_numu<78)*(samp1_final_level_numu==1)*(corrected_samp1_para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp1_pull_nch_numu,corrected_samp1_median_pull_numu,'wo',ms=8,label='Median Pull')
#pylab.plot(samp1_pull_nch_numu,pullfit_numu(samp1_pull_nch_numu,samp1_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Corrected Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_Corrected_ParaboloidPull_NuMu_UnWeighted_MuonError')


pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(real_final_level_numu==1)*(merged_nch_numu<78)]
y1=primary_energy_numu[(merged_nch_numu<78)*(real_final_level_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(60,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,4,190]
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'Primary Neutrino Energy (GeV)')
pylab.title('Nch vs. Primary Energy (Genie NuMu)')
pylab.imshow(numpy.log10(H),extent=extent,interpolation='nearest',aspect='auto',origin='lower')
cbar=pylab.colorbar()
#cbar.set_label(r'$Log_{10}(Counts)$', rotation=270)
pylab.savefig('FinalLevel_Nch_Vs_Energy_NuMu')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(real_final_level_numu==1)*(merged_nch_numu<78)*(cc_numu==1)]
y1=primary_energy_numu[(merged_nch_numu<78)*(real_final_level_numu==1)*(cc_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(60,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,4,190]
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'Primary Neutrino Energy (GeV)')
pylab.title('Nch vs. Primary Energy (Genie NuMu)')
pylab.imshow(numpy.log10(H),extent=extent,interpolation='nearest',aspect='auto',origin='lower')
cbar=pylab.colorbar()
#cbar.set_label(r'$Log_{10}(Counts)$', rotation=270)
pylab.savefig('FinalLevel_Nch_Vs_Energy_NuMu_CCOnly')


pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(real_final_level_numu==1)*(merged_nch_numu<78)*(les_numu==1)]
y1=primary_energy_numu[(merged_nch_numu<78)*(real_final_level_numu==1)*(les_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(60,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,4,190]
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'Primary Neutrino Energy (GeV)')
pylab.title('Nch vs. Primary Energy LES Only (Genie NuMu)')
pylab.imshow(numpy.log10(H),extent=extent,interpolation='nearest',aspect='auto',origin='lower')
cbar=pylab.colorbar()
#cbar.set_label(r'$Log_{10}(Counts)$', rotation=270)
pylab.savefig('FinalLevel_Nch_Vs_Energy_NuMu_LESOnly')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(real_final_level_numu==1)*(merged_nch_numu<78)*(les_numu==0)]
y1=primary_energy_numu[(merged_nch_numu<78)*(real_final_level_numu==1)*(les_numu==0)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(60,68))#,weights=atmo_numu[(samp1_final_level_numu==1)*(merged_nch_numu<78)*(samp1_para_pull_numu<12)])
extent = [10,78,4,190]
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'Primary Neutrino Energy (GeV)')
pylab.title('Nch vs. Primary Energy HES Only (Genie NuMu)')
pylab.imshow(numpy.log10(H),extent=extent,interpolation='nearest',aspect='auto',origin='lower')
cbar=pylab.colorbar()
#cbar.set_label(r'$Log_{10}(Counts)$', rotation=270)
pylab.savefig('FinalLevel_Nch_Vs_Energy_NuMu_HESOnly')


### Final Sample 2 (BDT Trained on Well-Reconstructed) ###
pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(sig_para_numu<(50/57.3))]
y1=57.3*sig_para_numu[(sig_para_numu<(50/57.3))*(merged_nch_numu<78)*(samp2_final_level_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(80,68),weights=atmo_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(sig_para_numu<(50/57.3))])
extent = [10,78,0,50]
#pylab.plot(xval,yval,'w--',lw=2.0)
pylab.xlabel('Nch')
pylab.ylabel('Paraboloid Sigma (Deg)')
pylab.title('Nch vs. Paraboloid (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_Paraboloid_NuMu')


### Long Range ###
pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp2_final_level_numu==1)*(merged_nch_numu<142)*(para_pull_numu<13)]
y1=para_pull_numu[(merged_nch_numu<142)*(samp2_final_level_numu==1)*(para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,132))#,weights=atmo_numu[(samp2_final_level_numu==1)*(merged_nch_numu<142)*(para_pull_numu<12)])
extent = [10,142,0,13]
pylab.plot(samp2_pull_nch_numu,samp2_median_pull_numu,'wo',ms=8,label='Median Pull')
pylab.plot(samp2_pull_nch_numu,pullfit_numu(samp2_pull_nch_numu,samp2_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_NuMu_UnWeighted_MuonError_ExtendedNch')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(para_pull_numu<13)]
y1=para_pull_numu[(merged_nch_numu<78)*(samp2_final_level_numu==1)*(para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu,samp2_median_pull_numu,'wo',ms=8,label='Median Pull')
pylab.plot(samp2_pull_nch_numu,pullfit_numu(samp2_pull_nch_numu,samp2_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_NuMu_UnWeighted_MuonError')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(corrected_samp2_para_pull_numu<13)]
y1=corrected_samp2_para_pull_numu[(merged_nch_numu<78)*(samp2_final_level_numu==1)*(corrected_samp2_para_pull_numu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(samp2_final_level_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu,corrected_samp2_median_pull_numu,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_numu,pullfit_numu(samp2_pull_nch_numu,samp2_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_Corrected_ParaboloidPull_NuMu_UnWeighted_MuonError')


################ NUMU NUGEN ######################

pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(samp1_final_level_nugmu==1)*(merged_nch_nugmu<78)*(sig_para_nugmu<(50/57.3))]
y1=57.3*sig_para_nugmu[(sig_para_nugmu<(50/57.3))*(merged_nch_nugmu<78)*(samp1_final_level_nugmu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(80,68),weights=atmo_nugmu[(samp1_final_level_nugmu==1)*(merged_nch_nugmu<78)*(sig_para_nugmu<(50/57.3))])
extent = [10,78,0,50]
#pylab.plot(xval,yval,'w--',lw=2.0)
pylab.xlabel('Nch')
pylab.ylabel('Paraboloid Sigma (Deg)')
pylab.title('Nch vs. Paraboloid (Genie NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_Paraboloid_NuMuNugen')


### Long Range ###
pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(samp1_final_level_nugmu==1)*(merged_nch_nugmu<142)*(para_pull_nugmu<13)]
y1=para_pull_nugmu[(merged_nch_nugmu<142)*(samp1_final_level_nugmu==1)*(para_pull_nugmu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,132))#,weights=atmo_nugmu[(samp1_final_level_nugmu==1)*(merged_nch_nugmu<142)*(para_pull_nugmu<12)])
extent = [10,142,0,13]
pylab.plot(samp1_pull_nch_nugmu,samp1_median_pull_nugmu,'wo',ms=8,label='Median Pull')
pylab.plot(samp1_pull_nch_nugmu,pullfit_nugmu(samp1_pull_nch_nugmu,samp1_pull_coeffs_nugmu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_ParaboloidPull_NuMuNugen_UnWeighted_MuonError_ExtendedNch')

pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(para_pull_nugmu<13)]
y1=para_pull_nugmu[(merged_nch_nugmu<78)*(good_angular_nugmu==1)*(para_pull_nugmu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp1_para_pull_nugmu<12)])
extent = [10,78,0,13]
pylab.plot(samp1_pull_nch_nugmu,samp1_median_pull_nugmu,'wo',ms=8,label='Median Pull')
pylab.plot(samp1_pull_nch_nugmu,pullfit_nugmu(samp1_pull_nch_nugmu,samp1_pull_coeffs_nugmu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_ParaboloidPull_NuMuNugen_UnWeighted_MuonError')


pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(corrected_samp1_para_pull_nugmu<13)]
y1=corrected_samp1_para_pull_nugmu[(merged_nch_nugmu<78)*(good_angular_nugmu==1)*(corrected_samp1_para_pull_nugmu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp1_para_pull_nugmu<12)])
extent = [10,78,0,13]
pylab.plot(samp1_pull_nch_nugmu,corrected_samp1_median_pull_nugmu,'wo',ms=8,label='Median Pull')
#pylab.plot(samp1_pull_nch_nugmu,pullfit_nugmu(samp1_pull_nch_nugmu,samp1_pull_coeffs_nugmu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Corrected Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample1_Nch_Vs_Corrected_ParaboloidPull_NuMuNugen_UnWeighted_MuonError')

### Final Sample 2 (BDT Trained on Well-Reconstructed) ###
pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(samp2_final_level_nugmu==1)*(merged_nch_nugmu<78)*(sig_para_nugmu<(50/57.3))]
y1=57.3*sig_para_nugmu[(sig_para_nugmu<(50/57.3))*(merged_nch_nugmu<78)*(samp2_final_level_nugmu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(80,68),weights=atmo_nugmu[(samp2_final_level_nugmu==1)*(merged_nch_nugmu<78)*(sig_para_nugmu<(50/57.3))])
extent = [10,78,0,50]
#pylab.plot(xval,yval,'w--',lw=2.0)
pylab.xlabel('Nch')
pylab.ylabel('Paraboloid Sigma (Deg)')
pylab.title('Nch vs. Paraboloid (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_Paraboloid_NuMuNugen')


### Long Range!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 ###
pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(samp2_final_level_nugmu==1)*(merged_nch_nugmu<142)*(para_pull_nugmu<13)]
y1=para_pull_nugmu[(merged_nch_nugmu<142)*(samp2_final_level_nugmu==1)*(para_pull_nugmu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,132))#,weights=atmo_nugmu[(samp2_final_level_nugmu==1)*(merged_nch_nugmu<142)*(para_pull_nugmu<12)])
extent = [10,142,0,13]
pylab.plot(samp2_pull_nch_nugmu,samp2_median_pull_nugmu,'wo',ms=8,label='Median Pull')
pylab.plot(samp2_pull_nch_nugmu,pullfit_nugmu(samp2_pull_nch_nugmu,samp2_pull_coeffs_nugmu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_NuMuNugen_UnWeighted_MuonError_ExtendedNch')

pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(merged_nch_nugmu<78)*(para_pull_nugmu<13)*(hes_nugmu==1)]
y1=para_pull_nugmu[(merged_nch_nugmu<78)*(para_pull_nugmu<13)*(hes_nugmu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp2_para_pull_nugmu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_nugmu,tailcut_samp2_median_pull_nugmu_hes,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_nugmu,samp2_peak_pull_nugmu,'bo',ms=8,label='Peak Pull')
pylab.plot(samp2_pull_nch_nugmu,pullfit_nugmu(samp2_pull_nch_nugmu,samp2_pull_coeffs_nugmu_hes),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('HES Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_HESNuMuNugen_UnWeighted_NeutrinoError_Tailcut')

pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(merged_nch_nugmu<142)*(para_pull_nugmu<13)*(hes_nugmu==1)]
y1=para_pull_nugmu[(merged_nch_nugmu<142)*(para_pull_nugmu<13)*(hes_nugmu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,132))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp2_para_pull_nugmu<12)])
extent = [10,142,0,13]
pylab.plot(samp2_pull_nch_nugmu,tailcut_samp2_median_pull_nugmu_hes,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_nugmu,samp2_peak_pull_nugmu,'bo',ms=8,label='Peak Pull')
pylab.plot(samp2_pull_nch_nugmu,pullfit_nugmu(samp2_pull_nch_nugmu,samp2_pull_coeffs_nugmu_hes),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('HES Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_HESNuMuNugen_UnWeighted_MuonError_LONGRANGE_TailCutMedian')


pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(merged_nch_nugmu<78)*(para_pull_nugmu<13)*(hes_nugmu==1)]
y1=corrected_samp2_para_pull_nugmu[(merged_nch_nugmu<78)*(para_pull_nugmu<13)*(hes_nugmu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp2_para_pull_nugmu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_nugmu,corrected_samp2_OneSig2D_pull_nugmu_hes,'wo',ms=8,label=r'2D Gaussian 1$\sigma$ Pull')
#pylab.plot(samp2_pull_nch_nugmu,samp2_peak_pull_nugmu,'bo',ms=8,label='Peak Pull')
#pylab.plot(samp2_pull_nch_nugmu,pullfit_nugmu(samp2_pull_nch_nugmu,samp2_pull_coeffs_nugmu_hes),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('HES Nch vs. Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_HESNuMuNugen_UnWeighted_NeutrinoError_CORRECTED')


pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(hes_numu==1)]
y1=para_pull_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(hes_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(good_angular_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu_hes,samp2_median_pull_numu_hes,'wo',ms=8,label='Median Pull')
pylab.plot(samp2_pull_nch_numu_hes,samp2_peak_pull_numu_hes,'bo',ms=8,label='Peak Pull')
#pylab.plot(samp2_pull_nch_numu,pullfit_numu(samp2_pull_nch_numu,samp2_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('HES Nch vs. Paraboloid Pull (GENIE NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_HESNuMuGENIE_UnWeighted_MuonError')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(les_numu==1)]
y1=para_pull_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(les_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(good_angular_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu_les,tailcut_samp2_median_pull_numu_les,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_numu_les,samp2_peak_pull_numu_les,'bo',ms=8,label='Peak Pull')
pylab.plot(samp2_pull_nch_numu_les,pullfit_numu(samp2_pull_nch_numu_les,samp2_pull_coeffs_numu_les),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('LES Nch vs. Paraboloid Pull (GENIE NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_LESNuMuGENIE_UnWeighted_NeutrinoError_Tailcut')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(les_numu==1)]
y1=corrected_samp2_para_pull_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(les_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(good_angular_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu_les,corrected_samp2_tailcut_pull_numu_les,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_numu_les,samp2_peak_pull_numu_les,'bo',ms=8,label='Peak Pull')
#pylab.plot(samp2_pull_nch_numu_les,pullfit_numu(samp2_pull_nch_numu_les,samp2_pull_coeffs_numu_les),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('LES Nch vs. Paraboloid Pull (GENIE NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_LESNuMuGENIE_UnWeighted_NeutrinoError_CORRECTED_Tailcut')


pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(merged_nch_numu<142)*(para_pull_numu<13)*(les_numu==1)]
y1=para_pull_numu[(merged_nch_numu<142)*(para_pull_numu<13)*(les_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(good_angular_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,142,0,13]
pylab.plot(samp2_pull_nch_numu_les,onesig2d_samp2_median_pull_numu_les,'wo',ms=8,label=r'2D Gaussian 1$\sigma$ Pull')
#pylab.plot(samp2_pull_nch_numu_les,samp2_peak_pull_numu_les,'bo',ms=8,label='Peak Pull')
pylab.plot(samp2_pull_nch_numu_les,pullfit_numu(samp2_pull_nch_numu_les,samp2_pull_coeffs_numu_les),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('LES Nch vs. Paraboloid Pull (GENIE NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_LESNuMuGENIE_UnWeighted_NeutrinoError_LongRange')

pylab.figure(figsize=(10,8))
x1=merged_nch_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(hes_numu==1)]
y1=para_pull_numu[(merged_nch_numu<78)*(para_pull_numu<13)*(hes_numu==1)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_numu[(good_angular_numu==1)*(merged_nch_numu<78)*(samp2_para_pull_numu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_numu,samp2_median_pull_numu,'wo',ms=8,label='Median Pull')
pylab.plot(samp2_pull_nch_numu,samp2_peak_pull_numu,'bo',ms=8,label='Peak Pull')
#pylab.plot(samp2_pull_nch_numu,pullfit_numu(samp2_pull_nch_numu,samp2_pull_coeffs_numu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{\nu} / \sigma$)')
pylab.title('HES Nch vs. Paraboloid Pull (GENIE NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_ParaboloidPull_HESNuMuGENIE_UnWeighted_MuonError')

pylab.figure(figsize=(10,8))
x1=merged_nch_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(corrected_samp2_para_pull_nugmu<13)]
y1=corrected_samp2_para_pull_nugmu[(merged_nch_nugmu<78)*(good_angular_nugmu==1)*(corrected_samp2_para_pull_nugmu<13)]
H,xedges,yedges=numpy.histogram2d(y1,x1,bins=(40,68))#,weights=atmo_nugmu[(good_angular_nugmu==1)*(merged_nch_nugmu<78)*(samp2_para_pull_nugmu<12)])
extent = [10,78,0,13]
pylab.plot(samp2_pull_nch_nugmu,corrected_samp2_median_pull_nugmu,'wo',ms=8,label='Median Pull')
#pylab.plot(samp2_pull_nch_nugmu,pullfit_nugmu(samp2_pull_nch_nugmu,samp2_pull_coeffs_nugmu),'w--',lw=2)
pylab.axhline(1,color='k',ls='dashed',linewidth=2.0)
pylab.xlabel('Nch')
pylab.legend()
pylab.ylabel(r'SplineMPE Paraboloid Pull ($\Delta \Psi_{muon} / \sigma$)')
pylab.title('Nch vs. Corrected Paraboloid Pull (Nugen NuMu)')
pylab.imshow(H,extent=extent,interpolation='nearest',aspect='auto',origin='lower')
pylab.savefig('Sample2_Nch_Vs_Corrected_ParaboloidPull_NuMuNugen_UnWeighted_MuonError')


### Energy Efficiencies ###
enhist_numu = pylab.hist(primary_energy_numu[(les_numu==1)*(cc_numu==1)],bins=numpy.linspace(4,190,40),histtype='step',lw=2,log=False,label='NuMu(N)',ec='b',weights=atmo_numu[(les_numu==1)*(cc_numu==1)])
samp1hist_numu = pylab.hist(primary_energy_numu[(cc_numu==1)*(samp1_final_level_numu==1)*(les_numu==1)],bins=numpy.linspace(4,190,40),histtype='step',lw=2,log=False,ec='b',ls='dashed',weights=atmo_numu[(samp1_final_level_numu==1)*(les_numu==1)*(cc_numu==1)])
samp2hist_numu = pylab.hist(primary_energy_numu[(cc_numu==1)*(samp2_final_level_numu==1)*(les_numu==1)],bins=numpy.linspace(4,190,40),histtype='step',lw=2,log=False,ec='b',ls='dashed',weights=atmo_numu[(samp2_final_level_numu==1)*(les_numu==1)*(cc_numu==1)])


cutrat1_numu=(1.0*samp1hist_numu[0]/enhist_numu[0])
cutrat1_numu[numpy.isnan(cutrat1_numu)==1]=0
cutrat2_numu=(1.0*samp2hist_numu[0]/enhist_numu[0])
cutrat2_numu[numpy.isnan(cutrat1_numu)==1]=0


genie_effa = numpy.array([5.73726197e-07,   3.24120773e-06,   8.54690148e-06,
         1.69890763e-05,   2.90093427e-05,   4.38059590e-05,
         5.91137176e-05,   8.00555758e-05,   9.82415859e-05,
         1.27731718e-04,   1.49674915e-04,   1.89958823e-04,
         2.26629158e-04,   2.70022005e-04,   3.03548251e-04,
         3.57449301e-04,   4.06208182e-04,   4.62600126e-04,
         5.47016694e-04,   6.05635029e-04,   6.70480602e-04,
         7.27587623e-04,   8.05674253e-04,   8.90001835e-04,
         9.34552254e-04,   1.06447711e-03,   1.16515091e-03,
         1.18622044e-03,   1.24991524e-03,   1.35586984e-03,
         1.55677488e-03,   1.67240598e-03,   1.70294141e-03,
         1.82606214e-03,   1.87869561e-03,   2.10777639e-03,
         2.22308725e-03,   2.16879972e-03,   2.45456387e-03,
         2.32325028e-03])
nugen_effa = numpy.array([0.0020747 ,  0.00215379,  0.0024384 ,  0.00261351,  0.00287766,
        0.00289464,  0.00306835,  0.00344779,  0.00341507,  0.00397283,
        0.00480283,  0.00513204,  0.0043073 ,  0.00510067,  0.00526595,
        0.0051937 ,  0.00525962,  0.00611304,  0.00665335,  0.00757188])

plotbinny_nugmu = numpy.array([195.25,  205.75,  216.25,  226.75,  237.25,  247.75,  258.25,
        268.75,  279.25,  289.75,  300.25,  310.75,  321.25,  331.75,
        342.25,  352.75,  363.25,  373.75,  384.25,  394.75])
plotbinny = numpy.array([6.325,   10.975,   15.625,   20.275,   24.925,   29.575,
         34.225,   38.875,   43.525,   48.175,   52.825,   57.475,
         62.125,   66.775,   71.425,   76.075,   80.725,   85.375,
         90.025,   94.675,   99.325,  103.975,  108.625,  113.275,
        117.925,  122.575,  127.225,  131.875,  136.525,  141.175,
        145.825,  150.475,  155.125,  159.775,  164.425,  169.075,
        173.725,  178.375,  183.025,  187.675])

area_coeff_nugmu = numpy.array([ -1.67663369e-04,   1.08457883e-05,   5.63771277e-08,  4.07278287e-04,  -2.34010444e-07])

def areafit_nugmu(x,a,b,c,d,e):
        val = a+b*x+numpy.power(c*x,2)+numpy.power(d*x,3)+numpy.power(e*x,4)
        return val


pylab.figure(figsize=(10,8))
binny = numpy.linspace(4,190,41)
pylab.plot(binny[:-1],cutrat1_numu,'b',lw=2,ls='steps',label='NuMu (CC) Sample 1')
pylab.plot(binny[:-1],cutrat2_numu,'g',lw=2,ls='steps',label='NuMu (CC) Sample 2')
pylab.grid()
pylab.xlabel('Energy (GeV)')
pylab.ylabel("Ratio")
pylab.title("LES Energy Efficiency on NuMu Tracks Sample 1 (GENIE)")
pylab.legend(loc='lower right')
pylab.axis([4,190,0.0,1.1])
pylab.savefig("FinalLevel_LESEnergyEfficiency_BothSamples_GenieNuMu_CCOnly")


his1=pylab.hist(primary_energy_numu[real_final_level_numu==1],bins=numpy.linspace(4,190,41),histtype='step',lw=2,label='GENIE',weights=primary_energy_numu[real_final_level_numu==1]**0,normed=True,log=True,ec='b')
his2=pylab.hist(primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)],bins=numpy.linspace(190,1000,31),weights=primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]**-0.5,histtype="step",lw=2,ec='b',ls='dashed',label="Nugen",normed=True)

plotbinny_nugmu = (his2[1])[:-1] + numpy.diff(his2[1])[0]/2
nugen_effa = areafit_nugmu(plotbinny_nugmu,area_coeff_nugmu[0],area_coeff_nugmu[1],area_coeff_nugmu[2],area_coeff_nugmu[3],area_coeff_nugmu[4])

histot = (genie_effa*his1[0]).sum() + (0.12699587233481074*nugen_effa*his2[0]).sum()

pylab.figure(figsize=(10,8))
pylab.semilogy(plotbinny_nugmu,0.12699587233481074*his2[0]*nugen_effa,'b--',label="Nugen",lw=2)
pylab.semilogy(plotbinny,his1[0]*genie_effa,'b-',label="GENIE",lw=2)
#pylab.vlines(ninetycontainment,1e-8,1e-5,colors='k',linestyles='dashed',label=r"$90\%$",lw=2)
#pylab.plot(b[:-1],a,'b_')
#pylab.hist(x=b[:-1]+numpy.diff(b)/2,bins=b,weights=a,histtype='step',label="Nugen",ls='dashed',ec="b",log=True,lw=2)
pylab.xlabel(r"Energy (GeV)")
pylab.ylabel(r"Flux (Arbitrary Units)")
pylab.legend(loc='upper right')
pylab.axis([4,1000,1e-7,5e-6])
pylab.grid()
pylab.title(r"Final Level Sample $E^{-2.5}$")
pylab.savefig("FinalLevel_NeutrinoEnergy_2_5")

his1=pylab.hist(primary_energy_numu[real_final_level_numu==1],bins=numpy.linspace(4,190,41),histtype='step',lw=2,label='GENIE',weights=primary_energy_numu[real_final_level_numu==1]**-0.5,normed=True,log=True,ec='b')
his2=pylab.hist(primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)*(primary_energy_nugmu<401)],bins=numpy.linspace(190,400,21),weights=primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)*(primary_energy_nugmu<401)]**-1.0,histtype="step",lw=2,ec='b',ls='dashed',label="Nugen",normed=True)

plotbinny_nugmu = (his2[1])[:-1] + numpy.diff(his2[1])[0]/2
nugen_effa = areafit_nugmu(plotbinny_nugmu,area_coeff_nugmu[0],area_coeff_nugmu[1],area_coeff_nugmu[2],area_coeff_nugmu[3],area_coeff_nugmu[4])

histot = (genie_effa*his1[0]).sum() + (0.02652560104008789*nugen_effa*his2[0]).sum()

contL = 29.575
contR = 331.75

pylab.figure(figsize=(10,8))
pylab.semilogy(plotbinny_nugmu,0.02652560104008789*his2[0]*nugen_effa,'g--',label="Nugen",lw=2)
pylab.semilogy(plotbinny,his1[0]*genie_effa,'g-',label="GENIE",lw=2)
pylab.vlines(contL,1e-8,5e-6,colors='k',linestyles='dashed',lw=2)
pylab.vlines(contR,1e-8,5e-6,colors='k',linestyles='dashed',lw=2)
pylab.fill_between([contL,contR],10**-8,5e-6,alpha=0.2)
#pylab.plot(b[:-1],a,'b_')
#pylab.hist(x=b[:-1]+numpy.diff(b)/2,bins=b,weights=a,histtype='step',label="Nugen",ls='dashed',ec="b",log=True,lw=2)
pylab.xlabel(r"Energy (GeV)")
pylab.ylabel(r"Flux (Arbitrary Units)")
pylab.legend(loc='upper right')
pylab.axis([4,400,1e-7,5e-6])
pylab.grid()
pylab.title(r"Final Level Sample $E^{-3}$")
pylab.savefig("FinalLevel_NeutrinoEnergy_3")

his1=pylab.hist(primary_energy_numu[real_final_level_numu==1],bins=numpy.linspace(4,190,41),histtype='step',lw=2,label='GENIE',weights=primary_energy_numu[real_final_level_numu==1]**-1.0,normed=True,log=True,ec='b')
his2=pylab.hist(primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)],bins=numpy.linspace(190,400,21),weights=primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]**-1.5,histtype="step",lw=2,ec='b',ls='dashed',label="Nugen",normed=True)

plotbinny_nugmu = (his2[1])[:-1] + numpy.diff(his2[1])[0]/2
nugen_effa = areafit_nugmu(plotbinny_nugmu,area_coeff_nugmu[0],area_coeff_nugmu[1],area_coeff_nugmu[2],area_coeff_nugmu[3],area_coeff_nugmu[4])

histot = (genie_effa*his1[0]).sum() + (0.006629808157329814*nugen_effa*his2[0]).sum()

contL = 10.975
contR = 258.25

pylab.figure(figsize=(10,8))
pylab.semilogy(plotbinny_nugmu,0.006629808157329814*his2[0]*nugen_effa,'r--',label="Nugen",lw=2)
pylab.semilogy(plotbinny,his1[0]*genie_effa,'r-',label="GENIE",lw=2)
pylab.vlines(contL,1e-8,5e-6,colors='k',linestyles='dashed',lw=2)
pylab.vlines(contR,1e-8,5e-6,colors='k',linestyles='dashed',lw=2)
pylab.fill_between([contL,contR],10**-8,5e-6,alpha=0.2)
#pylab.plot(b[:-1],a,'b_')
#pylab.hist(x=b[:-1]+numpy.diff(b)/2,bins=b,weights=a,histtype='step',label="Nugen",ls='dashed',ec="b",log=True,lw=2)
pylab.xlabel(r"Energy (GeV)")
pylab.ylabel(r"Flux (Arbitrary Units)")
pylab.legend(loc='upper right')
pylab.axis([4,400,2e-8,1e-6])
pylab.grid()
pylab.title(r"Final Level Sample $E^{-3.5}$")
pylab.savefig("FinalLevel_NeutrinoEnergy_3_5")


pylab.figure(figsize=(10,8))
pylab.hist(primary_energy_numu[real_final_level_numu==1],bins=numpy.linspace(4,190,40),histtype='step',lw=2,label='GENIE',weights=oweight_numu[real_final_level_numu==1]*primary_energy_numu[real_final_level_numu==1]**-0.5,normed=True,log=True,ec='g')
pylab.hist(primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)],bins=numpy.linspace(190,500,30),weights=(7.36*10**-12)*oweight_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]*primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]**-1,ec='g',histtype='step',lw=2,ls='dashed',label="Nugen")
pylab.xlabel(r"Energy (GeV)")
pylab.ylabel(r"Normed Counts")
pylab.legend(loc='upper right')
pylab.axis([4,500,1e-5,1])
pylab.grid()
pylab.title(r"Final Level Sample $E^{-3}$")
pylab.savefig("FinalLevel_NeutrinoEnergy_3")


pylab.figure(figsize=(10,8))
pylab.hist(primary_energy_numu[real_final_level_numu==1],bins=numpy.linspace(4,190,40),histtype='step',lw=2,label='GENIE',weights=oweight_numu[real_final_level_numu==1]*primary_energy_numu[real_final_level_numu==1]**-1,normed=True,log=True,ec='r')
pylab.hist(primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)],bins=numpy.linspace(190,500,30),histtype='step',lw=2,ls='dashed',label='Nugen',weights=(7.943658413674662e-11)*oweight_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]*primary_energy_nugmu[(real_final_level_nugmu==1)*(primary_energy_nugmu>189.99999)]**-1.5,normed=False,ec='r')
pylab.xlabel(r"Energy (GeV)")
pylab.ylabel(r"Normed Counts")
pylab.legend(loc='upper right')
pylab.axis([4,500,1e-5,1])
pylab.grid()
pylab.title(r"Final Level Sample $E^{-3.5}$")
pylab.savefig("FinalLevel_NeutrinoEnergy_3_5")



### Resolution ###
pylab.figure(figsize=(10,8))
pylab.plot(samp1_pull_nch_nugmu,57.3*samp1_median_res_nugmu,'g-',lw=2,label='Sample 1 (Straight Cut)')
pylab.plot(samp2_pull_nch_nugmu,57.3*samp2_median_res_nugmu,'b-',lw=2,label='Sample 2 (BDT)')
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,30])
pylab.grid()
pylab.title("Median Neutrino Resolution by Nch")
pylab.savefig("Nugen_OLDSELECTION_FinalLevel_NeutrinoResolutionByNch")

pylab.figure(figsize=(10,8))
pylab.plot(numpy.linspace(10,merged_nch_numu[(real_final_level_numu==1)].max(),merged_nch_numu[(real_final_level_numu==1)].max()-9),57.3*real_final_level_median_res_numu,'b-',lw=2,label="GENIE")
pylab.plot(numpy.linspace(10,merged_nch_nugmu[(real_final_level_nugmu==1)].max(),merged_nch_nugmu[(real_final_level_nugmu==1)].max()-9),57.3*real_final_level_median_res_nugmu,'g-',lw=2,label="Nugen")
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,45])
pylab.grid()
pylab.title("Median Neutrino Resolution by Nch")
pylab.savefig("FinalLevel_NeutrinoResolutionByNch")

pylab.figure(figsize=(10,8))
pylab.plot(numpy.linspace(10,merged_nch_nugmu[(real_final_level_nugmu==1)].max(),merged_nch_nugmu[(real_final_level_nugmu==1)].max()-9),57.3*real_final_level_median_res_nugmu,'b-',lw=2)
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,45])
pylab.grid()
pylab.title("Median Neutrino Resolution by Nch")
pylab.savefig("Nugen_FinalLevel_NeutrinoResolutionByNch")

pylab.figure(figsize=(10,8))
pylab.plot(numpy.linspace(10,merged_nch_numu[(les_numu==1)*(real_final_level_numu==1)].max(),merged_nch_numu[(les_numu==1)*(real_final_level_numu==1)].max()-9),57.3*les_real_final_level_median_res_numu,'b-',lw=2,label="GENIE")
pylab.plot(numpy.linspace(10,merged_nch_nugmu[(les_nugmu==1)*(real_final_level_nugmu==1)].max(),merged_nch_nugmu[(real_final_level_nugmu==1)*(les_nugmu==1)].max()-9),57.3*les_real_final_level_median_res_nugmu,'g-',lw=2,label="Nugen")
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,50])
pylab.grid()
pylab.title("LES Final Level Median Neutrino Resolution by Nch")
pylab.savefig("FinalLevel_LESNeutrinoResolutionByNch")

pylab.figure(figsize=(10,8))
pylab.plot(numpy.linspace(10,merged_nch_numu[(hes_numu==1)*(real_final_level_numu==1)].max(),merged_nch_numu[(hes_numu==1)*(real_final_level_numu==1)].max()-9),57.3*hes_real_final_level_median_res_numu,'b-',lw=2,label="GENIE")
pylab.plot(numpy.linspace(10,merged_nch_nugmu[(hes_nugmu==1)*(real_final_level_nugmu==1)].max(),merged_nch_nugmu[(real_final_level_nugmu==1)*(hes_nugmu==1)].max()-9),57.3*hes_real_final_level_median_res_nugmu,'g-',lw=2,label="Nugen")
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,15])
pylab.grid()
pylab.title("HES Final Level Median Neutrino Resolution by Nch")
pylab.savefig("FinalLevel_HESNeutrinoResolutionByNch")


pylab.figure(figsize=(10,8))
pylab.plot(energy_res_bins,57.3*real_final_level_median_res_energy_numu,'b-',lw=2,label="GENIE")
pylab.plot(energy_res_bins_nugmu,57.3*real_final_level_median_res_energy_nugmu,'g-',lw=2,label="Nugen")
#pylab.plot(energy_res_bins,57.3*real_final_level_median_paraboloid_numu,'b--',lw=2,label="GENIE (Paraboloid)")
#pylab.plot(energy_res_bins_nugmu,57.3*real_final_level_median_paraboloid_nugmu,'g--',lw=2,label="Nugen (Paraboloid)")
pylab.xlabel("Energy (GeV)")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([4,190,0,50])
pylab.grid()
pylab.title("Final Level Median Neutrino Resolution by Energy")
pylab.savefig("FinalLevel_NeutrinoResolutionByEnergy")

pylab.figure(figsize=(10,8))
pylab.plot(energy_res_bins,57.3*real_final_level_median_paraboloid_numu,'b-',lw=2,label="GENIE")
pylab.plot(energy_res_bins_nugmu,57.3*real_final_level_median_paraboloid_nugmu,'g-',lw=2,label="Nugen")
pylab.xlabel("Energy (GeV)")
pylab.ylabel(r"Paraboloid $\sigma$ Error Estimation $\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([4,190,0,50])
pylab.grid()
pylab.title("Final Level Median Neutrino Resolution by Energy")
pylab.savefig("FinalLevel_ParaboloidResolutionEstimationByEnergy")



pylab.figure(figsize=(10,8))
#pylab.plot(energy_res_bins,57.3*real_final_level_median_res_energy_numu,'b-',lw=2,label="GENIE")
pylab.plot(energy_res_bins_nugmu,57.3*real_final_level_median_res_energy_nugmu,'g-',lw=2,label="Nugen")
pylab.xlabel("Energy (GeV)")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([200,800,0,5.0])
pylab.grid()
pylab.title("Final Level Median Neutrino Resolution by Energy")
pylab.savefig("FinalLevel_NeutrinoResolutionByEnergy_NugenOnly")


pylab.figure(figsize=(10,8))
pylab.plot(numpy.linspace(10,merged_nch_numu[(real_final_level_numu==1)].max(),merged_nch_numu[(real_final_level_numu==1)].max()-9),57.3*real_final_level_median_res_numu,'b-',lw=2)
pylab.xlabel("Nch")
pylab.ylabel(r"$\Delta \Psi_{\nu}$")
pylab.legend(loc='upper right')
pylab.axis([10,100,0,45])
pylab.grid()
pylab.title("Median Neutrino Resolution by Nch")
pylab.savefig("GENIE_FinalLevel_NeutrinoResolutionByNch")


#### Sample 1 Plots ###
binzo = numpy.linspace(-1,0.087,20)
pylab.figure(figsize=(10,8))
pylab.hist(numpy.cos(splinemod_zen_data[samp1_final_level_data==1]),bins=binzo,histtype='step',label="Data")
pylab.grid()
pylab.xlabel(r"$\cos{\Theta_{SplineMPE}}$")
pylab.title('Final Level Zenith Distribution')
pylab.ylabel('Counts')
pylab.savefig('FinalLevel_ZenithDist_Sample1')

### Checking Kent Sigma Distribution ###
fudgey = 0.0003046174197867085
from numpy import pi,e,sqrt,power,cos,sin,sinh,exp

def Kent(r,sig,kludge=True,altk=False):
        if sig < 3:
                return Gaussian(r,sig)
        sig= sig*(pi/180)
        r = r*(pi/180)
        k = sig**-2
        if altk:
                k = numpy.log(1-0.39)/(cos(sig)-1)
        val = (k/(4*pi*sinh(k)))*exp(k*cos(r))
        if kludge:
                return val*fudgey
        return val

def Gaussian(r,sig):
        #sig = sig*(pi/180)
        #r = r*(pi/180)
        val = (1/(sig*sig*2*pi))*power(e,-(r**2)/(2*sig**2))
        return val

r_sep = numpy.linspace(0,50,6000)
sigma=20.5

pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_primary_diff_numu[(corrected_samp2_sig_para_numu*57.3>20)*(corrected_samp2_sig_para_numu*57.3<21)*(real_final_level_numu==1)],bins=numpy.linspace(0,50,25),normed=True,histtype="step")
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,sigma,altk=False),'-',c='purple',lw=2,label='Kent Dist',alpha=1.0)
#pylab.plot(r_sep,88.0*Gaussian(r_sep,sigma),'-',c='green',lw=2,label='Gaussian Dist',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\nu}$')
pylab.legend()
pylab.grid()
pylab.savefig("CorrectSigmaVsKentDistribution_20Deg_Nu")

r_sep = numpy.linspace(0,70,6000)
sigma=30.5

pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_daughter_diff_numu[(corrected_samp2_sig_para_numu*57.3>30)*(corrected_samp2_sig_para_numu*57.3<31)*(real_final_level_numu==1)],bins=numpy.linspace(0,70,30),normed=True,histtype="step")
pylab.plot(r_sep,1.66*88.0*Kent(r_sep,sigma,altk=False),'-',c='purple',lw=2,label='Kent Dist',alpha=1.0)
pylab.plot(r_sep,1.66*88.0*Gaussian(r_sep,sigma),'-',c='green',lw=2,label='Gaussian Dist',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\mu}$')
pylab.legend()
pylab.grid()
pylab.savefig("CorrectSigmaVsKentDistribution_30Deg_mu")

r_sep = numpy.linspace(0,80,6000)
sigma=40.5

pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_daughter_diff_numu[(cc_numu==1)*(corrected_samp2_sig_para_numu*57.3>40)*(corrected_samp2_sig_para_numu*57.3<41)*(real_final_level_numu==1)],bins=numpy.linspace(0,80,30),normed=True,histtype="step")
pylab.plot(r_sep,2.2*88.0*Kent(r_sep,sigma,altk=False),'-',c='purple',lw=2,label='Kent Dist',alpha=1.0)
pylab.plot(r_sep,2.2*88.0*Gaussian(r_sep,sigma),'-',c='green',lw=2,label='Gaussian Dist',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\mu}$')
pylab.legend()
pylab.grid()
pylab.savefig("CorrectSigmaVsKentDistribution_40Deg_mu")

r_sep = numpy.linspace(0,40,6000)
sigma=10.5

pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_primary_diff_numu[(corrected_samp2_sig_para_numu*57.3>10)*(corrected_samp2_sig_para_numu*57.3<11)*(real_final_level_numu==1)],bins=numpy.linspace(0,40,30),normed=True,histtype="step")
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,sigma,altk=False),'-',c='purple',lw=2,label='Kent Dist',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\nu}$')
pylab.grid()
pylab.savefig("CorrectSigmaVsKentDistribution_10Deg_Nu")


r_sep = numpy.linspace(0,60,6000)
sigma=3.


pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_primary_diff_numu[(corrected_samp2_sig_para_numu*57.3>2)*(corrected_samp2_sig_para_numu*57.3<4)*(real_final_level_numu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='purple',lw=2,weights=atmo_numu[(corrected_samp2_sig_para_numu*57.3>2)*(corrected_samp2_sig_para_numu*57.3<4)*(real_final_level_numu==1)])
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,sigma,altk=False),'-',c='purple',lw=2,label=r'Kent Dist $\sigma=3^{\circ}$',alpha=1.0)
his2=pylab.hist(57.3*spline_primary_diff_numu[(corrected_samp2_sig_para_numu*57.3>7)*(corrected_samp2_sig_para_numu*57.3<9)*(real_final_level_numu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='green',lw=2,weights=atmo_numu[(corrected_samp2_sig_para_numu*57.3>7)*(corrected_samp2_sig_para_numu*57.3<9)*(real_final_level_numu==1)])
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,8.0,altk=False),'-',c='green',lw=2,label=r'Kent Dist $\sigma=8^{\circ}$',alpha=1.0)
his3=pylab.hist(57.3*spline_primary_diff_numu[(corrected_samp2_sig_para_numu*57.3>14)*(corrected_samp2_sig_para_numu*57.3<16)*(real_final_level_numu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='k',lw=2,weights=atmo_numu[(corrected_samp2_sig_para_numu*57.3>14)*(corrected_samp2_sig_para_numu*57.3<16)*(real_final_level_numu==1)])
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,15.0,altk=False),'-',c='k',lw=2,label=r'Kent Dist $\sigma=15^{\circ}$',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\nu}$')
pylab.title("PSF Model vs. Actual Event Error NormalCorrection")
pylab.grid()
pylab.legend()
pylab.savefig("CorrectSigmaVsKentDistribution_Multi_Nu_SeparateStreamsOnly_Weighted")


r_sep = numpy.linspace(0,20,1000)
pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>.5)*(corrected_samp2_sig_para_nugmu*57.3<1.5)*(real_final_level_nugmu==1)*(hes_nugmu==1)],bins=numpy.linspace(0,20,40),normed=True,histtype="step",ec='purple',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,1.0,altk=False),'-',c='purple',lw=2,label=r'Kent Dist $\sigma=1^{\circ}$',alpha=1.0)
his2=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>2.5)*(corrected_samp2_sig_para_nugmu*57.3<3.5)*(real_final_level_nugmu==1)*(hes_nugmu==1)],bins=numpy.linspace(0,20,40),normed=True,histtype="step",ec='green',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,3.0,altk=False),'-',c='green',lw=2,label=r'Kent Dist $\sigma=3^{\circ}$',alpha=1.0)
his3=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>4.5)*(corrected_samp2_sig_para_nugmu*57.3<5.5)*(real_final_level_nugmu==1)*(hes_nugmu==1)],bins=numpy.linspace(0,20,40),normed=True,histtype="step",ec='k',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Kent(r_sep,5.0,altk=False),'-',c='k',lw=2,label=r'Kent Dist $\sigma=5^{\circ}$',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\nu}$')
pylab.title("PSF Model vs. Actual Event Error (HES Nugen) NormalCorrection")
pylab.grid()
pylab.legend()
pylab.savefig("CorrectSigmaVsKentDistribution_Multi_Nu_Nugen_HES_HighRes_SeparateStreamsOnly")


pylab.figure(figsize=(10,8))
his=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>2)*(corrected_samp2_sig_para_nugmu*57.3<4)*(real_final_level_nugmu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='purple',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Gaussian(r_sep,sigma),'-',c='purple',lw=2,label=r'Gaussian Dist $\sigma=3^{\circ}$',alpha=1.0)
his2=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>7)*(corrected_samp2_sig_para_nugmu*57.3<9)*(real_final_level_nugmu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='green',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Gaussian(r_sep,8.0),'-',c='green',lw=2,label=r'Gaussian Dist $\sigma=8^{\circ}$',alpha=1.0)
his3=pylab.hist(57.3*spline_primary_diff_nugmu[(corrected_samp2_sig_para_nugmu*57.3>14)*(corrected_samp2_sig_para_nugmu*57.3<16)*(real_final_level_nugmu==1)],bins=numpy.linspace(0,60,45),normed=True,histtype="step",ec='k',lw=2)
pylab.plot(r_sep,2*pi*(r_sep)*Gaussian(r_sep,15.0),'-',c='k',lw=2,label=r'Gaussian Dist $\sigma=15^{\circ}$',alpha=1.0)
pylab.xlabel(r'Angular Separation $\Delta\Psi_{\nu}$')
pylab.title("PSF Model vs. Actual Event Error (Nugen)")
pylab.grid()
pylab.legend()
pylab.savefig("CorrectSigmaVsKentDistribution_Multi_Nu_Nugen_Gaussian")


pylab.figure()
his=pylab.hist(numpy.log10(para_pull_numu[(samp2_final_level_numu==1)*(merged_nch_numu==15)]),bins=20,histtype='step')
pylab.vlines(numpy.median(para_pull_numu[(samp2_final_level_numu==1)*(merged_nch_numu==15)]),0,his[0].max(),colors='b',linestyles='dashed',lw=2.0,label="Median")
pylab.vlines(numpy.median(numpy.log10(para_pull_numu[(samp2_final_level_numu==1)*(merged_nch_numu==15)])),0,his[0].max(),colors='k',linestyles='dashed',lw=2.0,label="Median LOG")
pylab.legend()
pylab.savefig("CheckingPullDist")


pylab.figure()
pylab.hist(para_pull_numu[(les_numu==1)*(merged_nch_numu==15)*(para_pull_numu<50)],bins=20,histtype='step',ec='k',lw=2)
pylab.vlines(numpy.median(para_pull_numu[(les_numu==1)*(samp2_final_level_numu==1)*(merged_nch_numu==15)]),0,8000,colors='b',linestyles='dashed',lw=2.0,label="Median")
pylab.vlines(numpy.mean(para_pull_numu[(les_numu==1)*(samp2_final_level_numu==1)*(merged_nch_numu==15)]),0,8000,colors='g',linestyles='dashed',lw=2.0,label="Mean")
pylab.vlines(10**(scipy.stats.norm.fit(numpy.log10(para_pull_numu[(samp2_final_level_numu==1)*(merged_nch_numu==15)]))[0]),0,8000,colors='orange',linestyles='dashed',lw=2.0,label="Gaussian LOG")
pylab.ylabel("Counts")
pylab.xlabel("Paraboloid Pull")
pylab.legend()
pylab.title("Paraboloid Pull @ Nch = 15")
pylab.savefig("ExaminePull")



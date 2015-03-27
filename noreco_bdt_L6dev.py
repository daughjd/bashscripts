import sys,tables, numpy,pickle,pylab

from subprocess import call


### Acquire CORSIKA table info ###
#datums = tables.openFile('bdt_scored/Corsika_7437_0000-3999_Level3.L5Out.hdf5')
datums = tables.openFile('third_time/Corsika_7437_0000-3999_Level3.L6Out_Third_WNoRecoScore.hdf5')

dc12_c7437 = datums.root.FilterMask.col('DeepCoreFilter_12')
tldc12_c7437 = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_c7437 = datums.root.L4Bool_LES.col('value')
hes_c7437 = datums.root.L4Bool_HES.col('value')
bdt_les_c7437 = datums.root.LES_L6_BDTScore.col('value')
alt_bdt_les_c7437 = datums.root.NoReco_LES_L6_BDTScore.col('value')
bdt_hes_c7437 = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_c7437 = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_c7437 = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_c7437 = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_c7437 = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_c7437 = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_c7437 = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_c7437 = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_c7437 = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_c7437 = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_c7437 = datums.root.LineFit_DC.col('zenith')
linefit_azi_c7437 = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_c7437 = datums.root.FirstHLC.col('x')
#firsthlc_y_c7437 = datums.root.FirstHLC.col('y')
#firsthlc_z_c7437 = datums.root.FirstHLC.col('z')
dhd_c7437 = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_c7437 = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_c7437 = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_c7437 = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
#bdtscore_c7437 = datums.root.LES_L5_BDTScore.col('value')
finite_x_c7437 = datums.root.SplineMPEMod_Contained.col('x')
finite_y_c7437 = datums.root.SplineMPEMod_Contained.col('y')
finite_z_c7437 = datums.root.SplineMPEMod_Contained.col('z')
finite_length_c7437 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_c7437 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_c7437 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_c7437 = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_c7437 = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_c7437 = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_c7437 = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_c7437 = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_c7437 = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_c7437 = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_c7437 = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_c7437 = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_c7437 = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_c7437 = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_c7437 = numpy.sqrt(para1_c7437**2 + para2_c7437**2) / numpy.sqrt(2)

datums.close()

### 9622 ###
#datums = tables.openFile('third_time/bdt_scored/')
datums = tables.openFile('third_time/Level4_IC86.2011_corsika.009622.100K.L6Out_Third_WNoRecoScore.hdf5')

dc12_c9622 = datums.root.FilterMask.col('DeepCoreFilter_11')
#tldc12_c9622 = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_c9622 = datums.root.L4Bool_LES.col('value')
alt_bdt_les_c9622 = datums.root.NoReco_LES_L6_BDTScore.col('value')
#hes_c9622 = datums.root.L4Bool_HES.col('value')
bdt_les_c9622 = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_c9622 = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_c9622 = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_c9622 = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_c9622 = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_c9622 = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_c9622 = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_c9622 = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_c9622 = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_c9622 = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_c9622 = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_c9622 = datums.root.LineFit_DC.col('zenith')
linefit_azi_c9622 = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_c9622 = datums.root.FirstHLC.col('x')
#firsthlc_y_c9622 = datums.root.FirstHLC.col('y')
#firsthlc_z_c9622 = datums.root.FirstHLC.col('z')
dhd_c9622 = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_c9622 = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_c9622 = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_c9622 = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
#bdtscore_c9622 = datums.root.LES_L5_BDTScore.col('value')
finite_x_c9622 = datums.root.SplineMPEMod_Contained.col('x')
finite_y_c9622 = datums.root.SplineMPEMod_Contained.col('y')
finite_z_c9622 = datums.root.SplineMPEMod_Contained.col('z')
finite_length_c9622 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_c9622 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_c9622 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_c9622 = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_c9622 = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_c9622 = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_c9622 = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_c9622 = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_c9622 = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_c9622 = datums.root.L4_Variables_LES.col('Nch_HLC')
#hes_nch_clean_c9622 = datums.root.L4_Variables_HES.col('NCh_Clean')
#hes_hlc_c9622 = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_c9622 = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_c9622 = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_c9622 = numpy.sqrt(para1_c9622**2 + para2_c9622**2) / numpy.sqrt(2)

cweight_9622 = datums.root.CorsikaWeightMap.col('Weight')
polyweight_9622 = datums.root.CorsikaWeightMap.col('Polygonato')
h3aweight_9622 = datums.root.GaisserH3aWeight.col('value')
ts_9622 = datums.root.CorsikaWeightMap.col("TimeScale")

datums.close()

### 9255 ###
#datums = tables.openFile('third_time/bdt_scored/')
datums = tables.openFile('third_time/Level3_IC86.2011_corsika.009255.49805Files.L6Out_Third_WNoRecoScore.hdf5')

dc12_c9255 = datums.root.FilterMask.col('DeepCoreFilter_11')
#tldc12_c9255 = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_c9255 = datums.root.L4Bool_LES.col('value')
alt_bdt_les_c9255 = datums.root.NoReco_LES_L6_BDTScore.col('value')
#hes_c9255 = datums.root.L4Bool_HES.col('value')
bdt_les_c9255 = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_c9255 = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_c9255 = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_c9255 = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_c9255 = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_c9255 = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_c9255 = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_c9255 = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_c9255 = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_c9255 = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_c9255 = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_c9255 = datums.root.LineFit_DC.col('zenith')
linefit_azi_c9255 = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_c9255 = datums.root.FirstHLC.col('x')
#firsthlc_y_c9255 = datums.root.FirstHLC.col('y')
#firsthlc_z_c9255 = datums.root.FirstHLC.col('z')
dhd_c9255 = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_c9255 = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_c9255 = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_c9255 = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
#bdtscore_c9255 = datums.root.LES_L5_BDTScore.col('value')
finite_x_c9255 = datums.root.SplineMPEMod_Contained.col('x')
finite_y_c9255 = datums.root.SplineMPEMod_Contained.col('y')
finite_z_c9255 = datums.root.SplineMPEMod_Contained.col('z')
finite_length_c9255 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_c9255 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_c9255 = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_c9255 = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_c9255 = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_c9255 = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_c9255 = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_c9255 = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_c9255 = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_c9255 = datums.root.L4_Variables_LES.col('Nch_HLC')
#hes_nch_clean_c9255 = datums.root.L4_Variables_HES.col('NCh_Clean')
#hes_hlc_c9255 = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_c9255 = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_c9255 = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_c9255 = numpy.sqrt(para1_c9255**2 + para2_c9255**2) / numpy.sqrt(2)

cweight_9255 = datums.root.CorsikaWeightMap.col('Weight')
polyweight_9255 = datums.root.CorsikaWeightMap.col('Polygonato')
ts_9255 = datums.root.CorsikaWeightMap.col("TimeScale")

datums.close()

### Acquire Data Table Info ###

#datums = tables.openFile('third_time/bdt_scored/DCStream_IC86.2012_data_Run00120810_11_12_AllSubs.L5Out.hdf5')
datums  = tables.openFile('third_time/Level5_IC86.2012.CompleteDataSet.L6Out.330Days_Third_WNoRecoScore.hdf5')

dc12_data = datums.root.FilterMask.col('DeepCoreFilter_12')
tldc12_data = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_data = datums.root.L4Bool_LES.col('value')
alt_bdt_les_data = datums.root.NoReco_LES_L6_BDTScore.col('value')
hes_data = datums.root.L4Bool_HES.col('value')
bdt_les_data = datums.root.LES_L6_BDTScore.col('value')
event_data = datums.root.I3EventHeader.col('event')
bdt_hes_data = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_data = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_data = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_data = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_data = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_data = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_data = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_data = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_data = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_data = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_data = datums.root.LineFit_DC.col('zenith')
linefit_azi_data = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_data = datums.root.FirstHLC.col('x')
#firsthlc_y_data = datums.root.FirstHLC.col('y')
#firsthlc_z_data = datums.root.FirstHLC.col('z')
dhd_data = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_data = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_data = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_data = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
#bdtscore_data = datums.root.LES_L5_BDTScore.col('value')
finite_x_data = datums.root.SplineMPEMod_Contained.col('x')
finite_y_data = datums.root.SplineMPEMod_Contained.col('y')
finite_z_data = datums.root.SplineMPEMod_Contained.col('z')
finite_length_data = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_data = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_data = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_data = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_data = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_data = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_data = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_data = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_data = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_data = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_data = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_data = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_data = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_data = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_data = numpy.sqrt(para1_data**2 + para2_data**2) / numpy.sqrt(2)


datums.close()

'''
### Acquire NuTau Info ###
#datums = tables.openFile('third_time/bdt_scored/genie_ic.1504_1604.nutau_allruns.AtmoWeight.L5Out.hdf5')
datums = tables.openFile('third_time/genie_ic.1504_1604.combo_nutau_nutaubar.L6Out_Third.hdf5')

les_nutau = datums.root.L4Bool_LES.col('value')
hes_nutau = datums.root.L4Bool_HES.col('value')
les_l5bool_nutau = datums.root.L5Bool_LES.col('value')
hes_l5bool_nutau = datums.root.L5Bool_HES.col('value')
bdt_les_nutau = datums.root.LES_L6_BDTScore.col('value')
bdt_hes_nutau = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_nutau = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nutau = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_nutau = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_nutau = datums.root.SplineMPEModParaboloid.col('azimuth')
tight_spe6_zen_nutau = datums.root.SPEFit6_DC_Tight.col('zenith')
tight_spe6_azi_nutau = datums.root.SPEFit6_DC_Tight.col('azimuth')
noearly_spe2_zen_nutau = datums.root.SPEFit2_EarlyRemoved.col('zenith')
noearly_spe2_azi_nutau = datums.root.SPEFit2_EarlyRemoved.col('azimuth')
avg_distq_nutau = datums.root.SPEFit6_DCCharacteristics.col('avg_dom_dist_q_tot_dom')
lempty_nutau = datums.root.SPEFit6_DCCharacteristics.col('empty_hits_track_length')
separation_nutau = datums.root.SPEFit6_DCCharacteristics.col('track_hits_separation_length')
spe6_rlogl_nutau = datums.root.SPEFit6_DCFitParams.col('rlogl')
spe6_logl_nutau = datums.root.SPEFit6_DCFitParams.col('logl')
linefit_zen_nutau = datums.root.LineFit_DC.col('zenith')
linefit_azi_nutau = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_nutau = datums.root.FirstHLC.col('x')
#firsthlc_y_nutau = datums.root.FirstHLC.col('y')
#firsthlc_z_nutau = datums.root.FirstHLC.col('z')
dhd_nutau = datums.root.SPEFit6_DC_DirectHitsD.col('n_dir_pulses')
nlated_nutau = datums.root.SPEFit6_DC_DirectHitsD.col('n_late_pulses')
nearlyd_nutau = datums.root.SPEFit6_DC_DirectHitsD.col('n_early_pulses')
ldird_nutau = datums.root.SPEFit6_DC_DirectHitsD.col('dir_track_length')
#bdtscore_nutau = datums.root.LES_L5_BDTScore.col('value')
finite_x_nutau = datums.root.SPEFit6_DC_Finite.col('x')
finite_y_nutau = datums.root.SPEFit6_DC_Finite.col('y')
finite_z_nutau = datums.root.SPEFit6_DC_Finite.col('z')
finite_length_nutau = datums.root.SPEFit6_DC_FiniteCuts.col('length')
finite_startfract_nutau = datums.root.SPEFit6_DC_FiniteCuts.col('start_fraction')
finite_endfract_nutau = datums.root.SPEFit6_DC_FiniteCuts.col('end_fraction')
#fr_llh_inf_nutau = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_nutau = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_nutau = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_nutau = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_nutau = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_nutau = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_nutau = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_nutau = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_nutau = datums.root.L4_Variables_HES.col('Nch_HLC')
split_hes_nutau = datums.root.L4_Variables_HES.col('SplitEvent')

#para1_nutau = datums.root.SPEFit6Paraboloid_DCFitParams.col('err1')
#para2_nutau = datums.root.SPEFit6Paraboloid_DCFitParams.col('err2')
#sig_para_nutau = numpy.sqrt(para1_nutau**2 + para2_nutau**2) / numpy.sqrt(2)


intpartx_nutau = datums.root.InteractionParticle.col('x')
intparty_nutau = datums.root.InteractionParticle.col('y')
intpartz_nutau = datums.root.InteractionParticle.col('z')
int_zen_nutau = datums.root.InteractionParticle.col('zenith')
int_azi_nutau = datums.root.InteractionParticle.col('azimuth')
chkgrbeventw_nutau = datums.root.ChkGRBEventWeight.col('value')
primary_energy_nutau = datums.root.PrimaryNu.col('energy')
primary_zenith_nutau = datums.root.PrimaryNu.col('zenith')
primary_azi_nutau = datums.root.PrimaryNu.col('azimuth')
indc_nutau = datums.root.InDC.col('value')
intldc_nutau = datums.root.InExpDC.col('value')
atmo_nutau = datums.root.AtmoWeight.col('value')
oweight_nutau = datums.root.I3MCWeightDict.col('OneWeight')

datums.close()

'''
### Acquire NuMu Table Info ###

#datums = tables.openFile('third_time/bdt_scored/genie_ic.1304_1404_Combo_NuMu.L5.hdf5')
datums = tables.openFile('third_time/Genie_NuMu_Combo_AllSubs.L6Out_Third_WNoRecoScore.hdf5')

dc12_numu = datums.root.FilterMask.col('DeepCoreFilter_12')
tldc12_numu = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_numu = datums.root.L4Bool_LES.col('value')
hes_numu = datums.root.L4Bool_HES.col('value')
alt_bdt_les_numu = datums.root.NoReco_LES_L6_BDTScore.col('value')
bdt_les_numu = datums.root.LES_L6_BDTScore.col('value')
bdt_hes_numu = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_numu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_numu = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_numu = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_numu = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_numu = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_numu = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_numu = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_numu = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_numu = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_numu = datums.root.LineFit_DC.col('zenith')
linefit_azi_numu = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_numu = datums.root.FirstHLC.col('x')
#firsthlc_y_numu = datums.root.FirstHLC.col('y')
#firsthlc_z_numu = datums.root.FirstHLC.col('z')
dhd_numu = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_numu = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_numu = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_numu = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
finite_x_numu = datums.root.SplineMPEMod_Contained.col('x')
finite_y_numu = datums.root.SplineMPEMod_Contained.col('y')
finite_z_numu = datums.root.SplineMPEMod_Contained.col('z')
finite_length_numu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_numu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_numu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_numu = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_numu = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_numu = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_numu = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_numu = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_numu = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_numu = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_numu = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_numu = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_numu = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_numu = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_numu = numpy.sqrt(para1_numu**2 + para2_numu**2) / numpy.sqrt(2)

intpartx_numu = datums.root.InteractionParticle.col('x')
intparty_numu = datums.root.InteractionParticle.col('y')
intpartz_numu = datums.root.InteractionParticle.col('z')
int_zen_numu = datums.root.InteractionParticle.col('zenith')
int_azi_numu = datums.root.InteractionParticle.col('azimuth')
#chkgrbeventw_numu = datums.root.ChkGRBEventWeight.col('value')
primary_energy_numu = datums.root.PrimaryNu.col('energy')
primary_zenith_numu = datums.root.PrimaryNu.col('zenith')
primary_azi_numu = datums.root.PrimaryNu.col('azimuth')
indc_numu = datums.root.InDC.col('value')
intldc_numu = datums.root.InExpDC.col('value')
atmo_numu = datums.root.AtmoWeight.col('value')
oweight_numu = datums.root.I3MCWeightDict.col('OneWeight')
cc_numu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()

les_nch_clean_numu[les_nch_clean_numu==-1] = 0
hes_nch_clean_numu[hes_nch_clean_numu==-1] = 0
merged_nch_numu = les_nch_clean_numu + hes_nch_clean_numu


### Acquire NuMu Nugen ###

#datums = tables.openFile('third_time/bdt_scored/nugen_numu_IC86.2012_Filtered.008644.AtmoWeight.L4Out.L5Out.hdf5')
datums = tables.openFile('third_time/Level5_nugen_numu_IC86.2013.010090.AllSubs.L6Out_Third_WNoRecoScore.hdf5')

dc12_nugmu = datums.root.FilterMask.col('DeepCoreFilter_13')
tldc12_nugmu = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_13')
alt_bdt_les_nugmu = datums.root.NoReco_LES_L6_BDTScore.col('value')
les_nugmu = datums.root.L4Bool_LES.col('value')
hes_nugmu = datums.root.L4Bool_HES.col('value')
bdt_les_nugmu = datums.root.LES_L6_BDTScore.col('value')
bdt_hes_nugmu = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_nugmu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nugmu = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_nugmu = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_nugmu = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_nugmu = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_nugmu = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_nugmu = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_nugmu = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_nugmu = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_nugmu = datums.root.LineFit_DC.col('zenith')
linefit_azi_nugmu = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_nugmu = datums.root.FirstHLC.col('x')
#firsthlc_y_nugmu = datums.root.FirstHLC.col('y')
#firsthlc_z_nugmu = datums.root.FirstHLC.col('z')
dhd_nugmu = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_nugmu = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_nugmu = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_nugmu = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
finite_x_nugmu = datums.root.SplineMPEMod_Contained.col('x')
finite_y_nugmu = datums.root.SplineMPEMod_Contained.col('y')
finite_z_nugmu = datums.root.SplineMPEMod_Contained.col('z')
finite_length_nugmu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_nugmu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_nugmu = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_nugmu = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_nugmu = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_nugmu = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_nugmu = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_nugmu = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_nugmu = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_nugmu = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_nugmu = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_nugmu = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_nugmu = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_nugmu = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_nugmu = numpy.sqrt(para1_nugmu**2 + para2_nugmu**2) / numpy.sqrt(2)

intpartx_nugmu = datums.root.InteractionParticle.col('x')
intparty_nugmu = datums.root.InteractionParticle.col('y')
intpartz_nugmu = datums.root.InteractionParticle.col('z')
int_zen_nugmu = datums.root.InteractionParticle.col('zenith')
int_azi_nugmu = datums.root.InteractionParticle.col('azimuth')
#chkgrbeventw_nugmu = datums.root.ChkGRBEventWeight.col('value')
primary_energy_nugmu = datums.root.PrimaryNu.col('energy')
primary_zenith_nugmu = datums.root.PrimaryNu.col('zenith')
primary_azi_nugmu = datums.root.PrimaryNu.col('azimuth')
#indc_nugmu = datums.root.InDC.col('value')
#intldc_nugmu = datums.root.InExpDC.col('value')
atmo_nugmu = datums.root.AtmoWeight.col('value')
oweight_nugmu = datums.root.I3MCWeightDict.col('OneWeight')
cc_nugmu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()

les_nch_clean_nugmu[les_nch_clean_nugmu==-1] = 0
hes_nch_clean_nugmu[hes_nch_clean_nugmu==-1] = 0
merged_nch_nugmu = les_nch_clean_nugmu + hes_nch_clean_nugmu

### Acquire NuE Table Info ###    

#datums = tables.openFile('third_time/bdt_scored/genie_ic.1104_1204_Combo_Nue.L5.hdf5')
datums = tables.openFile('third_time/genie_ic.1104_1204.combo_nue_nuebar.L6Out_Third_WNoRecoScore.hdf5')

dc12_nue = datums.root.FilterMask.col('DeepCoreFilter_12')
tldc12_nue = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12')
les_nue = datums.root.L4Bool_LES.col('value')
hes_nue = datums.root.L4Bool_HES.col('value')
alt_bdt_les_nue = datums.root.NoReco_LES_L6_BDTScore.col('value')
bdt_les_nue = datums.root.LES_L6_BDTScore.col('value')
bdt_hes_nue = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_nue = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nue = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_nue = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_nue = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_nue = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_nue = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_nue = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_nue = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_nue = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_nue = datums.root.LineFit_DC.col('zenith')
linefit_azi_nue = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_nue = datums.root.FirstHLC.col('x')
#firsthlc_y_nue = datums.root.FirstHLC.col('y')
#firsthlc_z_nue = datums.root.FirstHLC.col('z')
dhd_nue = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_nue = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_nue = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_nue = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
finite_x_nue = datums.root.SplineMPEMod_Contained.col('x')
finite_y_nue = datums.root.SplineMPEMod_Contained.col('y')
finite_z_nue = datums.root.SplineMPEMod_Contained.col('z')
finite_length_nue = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_nue = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_nue = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_nue = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_nue = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_nue = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_nue = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_nue = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_nue = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_nue = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_nue = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_nue = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_nue = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_nue = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_nue = numpy.sqrt(para1_nue**2 + para2_nue**2) / numpy.sqrt(2)

intpartx_nue = datums.root.InteractionParticle.col('x')
intparty_nue = datums.root.InteractionParticle.col('y')
intpartz_nue = datums.root.InteractionParticle.col('z')
int_zen_nue = datums.root.InteractionParticle.col('zenith')
int_azi_nue = datums.root.InteractionParticle.col('azimuth')
#chkgrbeventw_nue = datums.root.ChkGRBEventWeight.col('value')
primary_energy_nue = datums.root.PrimaryNu.col('energy')
primary_zenith_nue = datums.root.PrimaryNu.col('zenith')
primary_azi_nue = datums.root.PrimaryNu.col('azimuth')
indc_nue = datums.root.InDC.col('value')
intldc_nue = datums.root.InExpDC.col('value')
atmo_nue = datums.root.AtmoWeight.col('value')
oweight_nue = datums.root.I3MCWeightDict.col('OneWeight')


datums.close()

### Grab Nugen Nue ###

#datums = tables.openFile('third_time/bdt_scored/nugen_nue_IC86.2012_Filtered.008654.AtmoWeight.L4Out.L5Out.hdf5')
datums = tables.openFile('third_time/Level5_nugen_nue_IC86.2013.010193.AllSubs.L6Out_Third_WNoRecoScore.hdf5')

dc12_nuge = datums.root.FilterMask.col('DeepCoreFilter_13')
tldc12_nuge = datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_13')
les_nuge = datums.root.L4Bool_LES.col('value')
hes_nuge = datums.root.L4Bool_HES.col('value')
alt_bdt_les_nuge = datums.root.NoReco_LES_L6_BDTScore.col('value')
bdt_les_nuge = datums.root.LES_L6_BDTScore.col('value')
bdt_hes_nuge = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_nuge = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nuge = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_nuge = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_nuge = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_nuge = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_nuge = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_nuge = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_nuge = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_nuge = datums.root.SplineMPEModFitParams.col('logl')
linefit_zen_nuge = datums.root.LineFit_DC.col('zenith')
linefit_azi_nuge = datums.root.LineFit_DC.col('azimuth')
#firsthlc_x_nuge = datums.root.FirstHLC.col('x')
#firsthlc_y_nuge = datums.root.FirstHLC.col('y')
#firsthlc_z_nuge = datums.root.FirstHLC.col('z')
dhd_nuge = datums.root.SplineMPEMod_DirectHitsD.col('n_dir_pulses')
nlated_nuge = datums.root.SplineMPEMod_DirectHitsD.col('n_late_pulses')
nearlyd_nuge = datums.root.SplineMPEMod_DirectHitsD.col('n_early_pulses')
ldird_nuge = datums.root.SplineMPEMod_DirectHitsD.col('dir_track_length')
#bdtscore_nuge = datums.root.LES_L5_BDTScore.col('value')
finite_x_nuge = datums.root.SplineMPEMod_Contained.col('x')
finite_y_nuge = datums.root.SplineMPEMod_Contained.col('y')
finite_z_nuge = datums.root.SplineMPEMod_Contained.col('z')
finite_length_nuge = datums.root.SplineMPEMod_Contained_FiniteCuts.col('length')
finite_startfract_nuge = datums.root.SplineMPEMod_Contained_FiniteCuts.col('start_fraction')
finite_endfract_nuge = datums.root.SplineMPEMod_Contained_FiniteCuts.col('end_fraction')
#fr_llh_inf_nuge = datums.root.StartStopProb_StartStopParams.col('llh_inf_track')
#fr_llh_start_nuge = datums.root.StartStopProb_StartStopParams.col('llh_starting_track')
#fr_llh_stop_nuge = datums.root.StartStopProb_StartStopParams.col('llh_stopping_track')
vetotrackcharge_nuge = datums.root.L5VetoTrackVetoCharge.col('value')
vetotrackhits_nuge = datums.root.L5VetoTrackVetoChannels.col('value')
les_nch_clean_nuge = datums.root.L4_Variables_LES.col('NCh_Clean')
les_hlc_nuge = datums.root.L4_Variables_LES.col('Nch_HLC')
hes_nch_clean_nuge = datums.root.L4_Variables_HES.col('NCh_Clean')
hes_hlc_nuge = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_nuge = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_nuge = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_nuge = numpy.sqrt(para1_nuge**2 + para2_nuge**2) / numpy.sqrt(2)

intpartx_nuge = datums.root.InteractionParticle.col('x')
intparty_nuge = datums.root.InteractionParticle.col('y')
intpartz_nuge = datums.root.InteractionParticle.col('z')
int_zen_nuge = datums.root.InteractionParticle.col('zenith')
int_azi_nuge = datums.root.InteractionParticle.col('azimuth')
#chkgrbeventw_nuge = datums.root.ChkGRBEventWeight.col('value')
oweight_nuge = datums.root.I3MCWeightDict.col('OneWeight')
primary_energy_nuge = datums.root.PrimaryNu.col('energy')
primary_zenith_nuge = datums.root.PrimaryNu.col('zenith')
primary_azi_nuge = datums.root.PrimaryNu.col('azimuth')
#indc_nuge = datums.root.InDC.col('value')
#intldc_nuge = datums.root.InExpDC.col('value')
atmo_nuge = datums.root.AtmoWeight.col('value')

datums.close()

##############################################

dataweight=1.0/(28475585.0)
dataweight=dataweight*numpy.ones(len(les_data))

C7437EventWeight=(10.03598*3991)**-1
C7437EventWeight=C7437EventWeight*numpy.ones(len(les_c7437))
C9622EventWeight=cweight_9622*polyweight_9622/ts_9622
C9622EventWeight=C9622EventWeight/100000.

H3A_C9622EventWeight=h3aweight_9622/100000.

C9255EventWeight=cweight_9255*polyweight_9255/ts_9255
C9255EventWeight=C9255EventWeight/49805.

atmo_numu = atmo_numu*10**9/4000.
atmo_nue = atmo_nue*10**9/2000.
atmo_nugmu = atmo_nugmu*10**9/1000.
atmo_nuge = atmo_nuge*10**9/1000.

les_nch_clean_nugmu[les_nch_clean_nugmu==-1] = 0
hes_nch_clean_nugmu[hes_nch_clean_nugmu==-1] = 0
merged_nch_nugmu = les_nch_clean_nugmu + hes_nch_clean_nugmu

les_nch_clean_numu[les_nch_clean_numu==-1] = 0
hes_nch_clean_numu[hes_nch_clean_numu==-1] = 0
merged_nch_numu = les_nch_clean_numu + hes_nch_clean_numu

les_nch_clean_c7437[les_nch_clean_c7437==-1] = 0
hes_nch_clean_c7437[hes_nch_clean_c7437==-1] = 0
merged_nch_c7437 = les_nch_clean_c7437 + hes_nch_clean_c7437

les_nch_clean_data[les_nch_clean_data==-1] = 0
hes_nch_clean_data[hes_nch_clean_data==-1] = 0
merged_nch_data = les_nch_clean_data + hes_nch_clean_data

les_nch_clean_nue[les_nch_clean_nue==-1] = 0
hes_nch_clean_nue[hes_nch_clean_nue==-1] = 0
merged_nch_nue = les_nch_clean_nue + hes_nch_clean_nue

les_nch_clean_nuge[les_nch_clean_nuge==-1] = 0
hes_nch_clean_nuge[hes_nch_clean_nuge==-1] = 0
merged_nch_nuge = les_nch_clean_nuge + hes_nch_clean_nuge



fr_R_c7437 = numpy.sqrt((finite_x_c7437-46)**2 + (finite_y_c7437+35)**2)
fr_R_nuge = numpy.sqrt((finite_x_nuge-46)**2 + (finite_y_nuge+35)**2)
fr_R_nue = numpy.sqrt((finite_x_nue-46)**2 + (finite_y_nue+35)**2)
fr_R_numu = numpy.sqrt((finite_x_numu-46)**2 + (finite_y_numu+35)**2)
fr_R_nugmu = numpy.sqrt((finite_x_nugmu-46)**2 + (finite_y_nugmu+35)**2)
fr_R_data = numpy.sqrt((finite_x_data-46)**2 + (finite_y_data+35)**2)
fr_R_c9622 = numpy.sqrt((finite_x_c9622-46)**2 + (finite_y_c9622+35)**2)
fr_R_c9255 = numpy.sqrt((finite_x_c9255-46)**2 + (finite_y_c9255+35)**2)

#firsthlc_R_c7437 = numpy.sqrt((firsthlc_x_c7437-46)**2 + (firsthlc_y_c7437+35)**2)
#firsthlc_R_nuge = numpy.sqrt((firsthlc_x_nuge-46)**2 + (firsthlc_y_nuge+35)**2)
#firsthlc_R_nue = numpy.sqrt((firsthlc_x_nue-46)**2 + (firsthlc_y_nue+35)**2)
#firsthlc_R_numu = numpy.sqrt((firsthlc_x_numu-46)**2 + (firsthlc_y_numu+35)**2)
#firsthlc_R_nugmu = numpy.sqrt((firsthlc_x_nugmu-46)**2 + (firsthlc_y_nugmu+35)**2)
#firsthlc_R_nutau = numpy.sqrt((firsthlc_x_nutau-46)**2 + (firsthlc_y_nutau+35)**2)
#firsthlc_R_d10 = numpy.sqrt((firsthlc_x_d10-46)**2 + (firsthlc_y_d10+35)**2)

spline_plogl_c7437 = spline_logl_c7437/(merged_nch_c7437 - 2.5)
spline_plogl_numu = spline_logl_numu/(merged_nch_numu - 2.5)
spline_plogl_nugmu = spline_logl_nugmu/(merged_nch_nugmu - 2.5)
spline_plogl_nue = spline_logl_nue/(merged_nch_nue - 2.5)
spline_plogl_nuge = spline_logl_nuge/(merged_nch_nuge - 2.5)
spline_plogl_data = spline_logl_data/(merged_nch_data - 2.5)
spline_plogl_c9622 = spline_logl_c9622/(les_nch_clean_c9622 - 2.5)
spline_plogl_c9255 = spline_logl_c9255/(les_nch_clean_c9255 - 2.5)

def AngleBetweenAngles(theta1,phi1,theta2,phi2):
   x1 = numpy.sin(theta1)*numpy.cos(phi1)
   y1 = numpy.sin(theta1)*numpy.sin(phi1)
   z1 = numpy.cos(theta1)
   x2 = numpy.sin(theta2)*numpy.cos(phi2)
   y2 = numpy.sin(theta2)*numpy.sin(phi2)
   z2 = numpy.cos(theta2)
   return numpy.arccos(x1*x2+y1*y2+z1*z2)

spline_paraspline_diff_numu = AngleBetweenAngles(spline_para_zen_numu,spline_para_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
spline_paraspline_diff_nugmu = AngleBetweenAngles(spline_para_zen_nugmu,spline_para_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)
spline_paraspline_diff_nue = AngleBetweenAngles(spline_para_zen_nue,spline_para_azi_nue,splinemod_zen_nue,splinemod_azi_nue)
spline_paraspline_diff_nuge = AngleBetweenAngles(spline_para_zen_nuge,spline_para_azi_nuge,splinemod_zen_nuge,splinemod_azi_nuge)
spline_paraspline_diff_data = AngleBetweenAngles(spline_para_zen_data,spline_para_azi_data,splinemod_zen_data,splinemod_azi_data)
spline_paraspline_diff_c7437 = AngleBetweenAngles(spline_para_zen_c7437,spline_para_azi_c7437,splinemod_zen_c7437,splinemod_azi_c7437)
spline_paraspline_diff_c9622 = AngleBetweenAngles(spline_para_zen_c9622,spline_para_azi_c9622,splinemod_zen_c9622,splinemod_azi_c9622)
spline_paraspline_diff_c9255 = AngleBetweenAngles(spline_para_zen_c9255,spline_para_azi_c9255,splinemod_zen_c9255,splinemod_azi_c9255)


lf_spline_diff_data = AngleBetweenAngles(linefit_zen_data,linefit_azi_data,splinemod_zen_data,splinemod_azi_data)
lf_spline_diff_c7437 = AngleBetweenAngles(linefit_zen_c7437,linefit_azi_c7437,splinemod_zen_c7437,splinemod_azi_c7437)
lf_spline_diff_c9622 = AngleBetweenAngles(linefit_zen_c9622,linefit_azi_c9622,splinemod_zen_c9622,splinemod_azi_c9622)
lf_spline_diff_c9255 = AngleBetweenAngles(linefit_zen_c9255,linefit_azi_c9255,splinemod_zen_c9255,splinemod_azi_c9255)
lf_spline_diff_numu = AngleBetweenAngles(linefit_zen_numu,linefit_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
lf_spline_diff_nugmu = AngleBetweenAngles(linefit_zen_nugmu,linefit_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)
lf_spline_diff_nuge = AngleBetweenAngles(linefit_zen_nuge,linefit_azi_nuge,splinemod_zen_nuge,splinemod_azi_nuge)
lf_spline_diff_nue = AngleBetweenAngles(linefit_zen_nue,linefit_azi_nue,splinemod_zen_nue,splinemod_azi_nue)

spline_daughter_diff_numu = AngleBetweenAngles(int_zen_numu,int_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
spline_daughter_diff_nugmu = AngleBetweenAngles(int_zen_nugmu,int_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)
spline_daughter_diff_nuge = AngleBetweenAngles(int_zen_nuge,int_azi_nuge,splinemod_zen_nuge,splinemod_azi_nuge)
spline_daughter_diff_nue = AngleBetweenAngles(int_zen_nue,int_azi_nue,splinemod_zen_nue,splinemod_azi_nue)

spline_primary_diff_numu = AngleBetweenAngles(primary_zenith_numu,primary_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
spline_primary_diff_nugmu = AngleBetweenAngles(primary_zenith_nugmu,primary_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)
spline_primary_diff_nuge = AngleBetweenAngles(primary_zenith_nuge,primary_azi_nuge,splinemod_zen_nuge,splinemod_azi_nuge)
spline_primary_diff_nue = AngleBetweenAngles(primary_zenith_nue,primary_azi_nue,splinemod_zen_nue,splinemod_azi_nue)

para_pull_numu = spline_primary_diff_numu/sig_para_numu
para_pull_nugmu = spline_primary_diff_nugmu/sig_para_nugmu

#noearly_spe_primary_diff_numu = AngleBetweenAngles(primary_zenith_numu,primary_azi_numu,noearly_spe2_zen_numu,noearly_spe2_azi_numu)
#noearly_spe_primary_diff_nugmu = AngleBetweenAngles(primary_zenith_nugmu,primary_azi_nugmu,noearly_spe2_zen_nugmu,noearly_spe2_azi_nugmu)
#noearly_spe_primary_diff_nuge = AngleBetweenAngles(primary_zenith_nuge,primary_azi_nuge,noearly_spe2_zen_nuge,noearly_spe2_azi_nuge)
#noearly_spe_primary_diff_nue = AngleBetweenAngles(primary_zenith_nue,primary_azi_nue,noearly_spe2_zen_nue,noearly_spe2_azi_nue)

#tight_spe_primary_diff_numu = AngleBetweenAngles(primary_zenith_numu,primary_azi_numu,tight_spe6_zen_numu,tight_spe6_azi_numu)
#tight_spe_primary_diff_nugmu = AngleBetweenAngles(primary_zenith_nugmu,primary_azi_nugmu,tight_spe6_zen_nugmu,tight_spe6_azi_nugmu)
#tight_spe_primary_diff_nuge = AngleBetweenAngles(primary_zenith_nuge,primary_azi_nuge,tight_spe6_zen_nuge,tight_spe6_azi_nuge)
#tight_spe_primary_diff_nue = AngleBetweenAngles(primary_zenith_nue,primary_azi_nue,tight_spe6_zen_nue,tight_spe6_azi_nue)

fr_vertex_error_nugmu = numpy.sqrt((finite_x_nugmu-intpartx_nugmu)**2 + (finite_y_nugmu-intparty_nugmu)**2 + (finite_z_nugmu-intpartz_nugmu)**2)
fr_vertex_error_nuge = numpy.sqrt((finite_x_nuge-intpartx_nuge)**2 + (finite_y_nuge-intparty_nuge)**2 + (finite_z_nuge-intpartz_nuge)**2)
fr_vertex_error_numu = numpy.sqrt((finite_x_numu-intpartx_numu)**2 + (finite_y_numu-intparty_numu)**2 + (finite_z_numu-intpartz_numu)**2)
fr_vertex_error_nue = numpy.sqrt((finite_x_nue-intpartx_nue)**2 + (finite_y_nue-intparty_nue)**2 + (finite_z_nue-intpartz_nue)**2)

les_preselect_numu = (les_numu==1)*(numpy.isfinite(avg_distq_numu))*(numpy.isfinite(finite_x_numu))#*(merged_nch_numu>15)
les_preselect_nugmu = (les_nugmu==1)*(numpy.isfinite(avg_distq_nugmu))*(numpy.isfinite(finite_x_nugmu))#*(merged_nch_nugmu>15)
les_preselect_nuge = (les_nuge==1)*(numpy.isfinite(avg_distq_nuge))*(numpy.isfinite(finite_x_nuge))#*(merged_nch_nuge>15)
les_preselect_nue = (les_nue==1)*(numpy.isfinite(avg_distq_nue))*(numpy.isfinite(finite_x_nue))#*(merged_nch_nue>15)
les_preselect_data = (les_data==1)*(numpy.isfinite(avg_distq_data))*(numpy.isfinite(finite_x_data))#*(merged_nch_data>15)
les_preselect_c7437 = (les_c7437==1)*(numpy.isfinite(avg_distq_c7437))*(numpy.isfinite(finite_x_c7437))#*(merged_nch_c7437>15)
les_preselect_c9622 = (les_c9622==1)*(numpy.isfinite(avg_distq_c9622))*(numpy.isfinite(finite_x_c9622))#*(les_nch_clean_c9622>15)
les_preselect_c9255 = (les_c9255==1)*(numpy.isfinite(avg_distq_c9255))*(numpy.isfinite(finite_x_c9255))#*(les_nch_clean_c9255>15)


hes_preselect_numu = (hes_numu==1)*(numpy.isfinite(avg_distq_numu))*(numpy.isfinite(dhd_numu))*(numpy.isfinite(spline_rlogl_numu))
hes_preselect_nugmu = (hes_nugmu==1)*(numpy.isfinite(avg_distq_nugmu))*(numpy.isfinite(dhd_nugmu))*(numpy.isfinite(spline_rlogl_nugmu))
hes_preselect_nuge = (hes_nuge==1)*(numpy.isfinite(avg_distq_nuge))*(numpy.isfinite(dhd_nuge))*(numpy.isfinite(spline_rlogl_nuge))
hes_preselect_nue = (hes_nue==1)*(numpy.isfinite(avg_distq_nue))*(numpy.isfinite(dhd_nue))*(numpy.isfinite(spline_rlogl_nue))
hes_preselect_data = (hes_data==1)*(numpy.isfinite(avg_distq_data))*(numpy.isfinite(dhd_data))*(numpy.isfinite(spline_rlogl_data))
hes_preselect_c7437 = (hes_c7437==1)*(numpy.isfinite(avg_distq_c7437))*(numpy.isfinite(dhd_c7437))*(numpy.isfinite(spline_rlogl_c7437))

good_angular_numu = (les_preselect_numu==1)*(spline_rlogl_numu<7.5)*(spline_paraspline_diff_numu<(25/57.3))
good_angular_nugmu = (les_preselect_nugmu==1)*(spline_rlogl_nugmu<7.5)*(spline_paraspline_diff_nugmu<(25/57.3))
good_angular_nue = (les_preselect_nue==1)*(spline_rlogl_nue<7.5)*(spline_paraspline_diff_nue<(25/57.3))
good_angular_nuge = (les_preselect_nuge==1)*(spline_rlogl_nuge<7.5)*(spline_paraspline_diff_nuge<(25/57.3))
good_angular_c7437 = (les_preselect_c7437==1)*(spline_rlogl_c7437<7.5)*(spline_paraspline_diff_c7437<(25/57.3))
good_angular_c9622 = (les_preselect_c9622==1)*(spline_rlogl_c9622<7.5)*(spline_paraspline_diff_c9622<(25/57.3))
good_angular_c9255 = (les_preselect_c9255==1)*(spline_rlogl_c9255<7.5)*(spline_paraspline_diff_c9255<(25/57.3))
good_angular_data = (les_preselect_data==1)*(spline_rlogl_data<7.5)*(spline_paraspline_diff_data<(25/57.3))


badly_recoed_numu = spline_daughter_diff_numu > 5/57.3
badly_recoed_nugmu = spline_daughter_diff_nugmu > 5/57.3

les_trained_numu = (les_preselect_numu==1)*(badly_recoed_numu==0)

#les_l6bool_c7437 = (les_l5bool_c7437==1)*(bdt_les_c7437 > 0.0)
#les_l6bool_c9622 = (les_l5bool_c9622==1)*(bdt_les_c9622 > 0.0)
#les_l6bool_data = (les_l5bool_data==1)*(bdt_les_data > 0.0)
#les_l6bool_nuge = (les_l5bool_nuge==1)*(bdt_les_nuge > 0.0)
#les_l6bool_nue = (les_l5bool_nue==1)*(bdt_les_nue > 0.0)
#les_l6bool_numu = (les_l5bool_numu==1)*(bdt_les_numu > 0.0)


#hes_l6bool_data = (hes_l5bool_data==1)*(bdt_hes_data > -0.025)
#hes_l6bool_nugmu = (hes_l5bool_nugmu==1)*(bdt_hes_nugmu > -0.025)
#hes_l6bool_nuge = (hes_l5bool_nuge==1)*(bdt_hes_nuge > -0.025)
#hes_l6bool_nue = (hes_l5bool_nue==1)*(bdt_hes_nue > -0.025)
#hes_l6bool_numu = (hes_l5bool_numu==1)*(bdt_hes_numu > -0.025)

samp1_final_level_numu = ((les_preselect_numu==1)*(spline_rlogl_numu<7.5) + (hes_preselect_numu==1)*(bdt_hes_numu>-0.01))*(numpy.cos(splinemod_zen_numu)<0.087)
samp1_final_level_nugmu = ((les_preselect_nugmu==1)*(spline_rlogl_nugmu<7.5) + (hes_preselect_nugmu==1)*(bdt_hes_nugmu>-0.01))*(numpy.cos(splinemod_zen_nugmu)<0.087)
samp1_final_level_data = ((les_preselect_data==1)*(spline_rlogl_data<7.5) + (hes_preselect_data==1)*(bdt_hes_data>-0.01))*(numpy.cos(splinemod_zen_data)<0.087)
samp1_final_level_nue = ((les_preselect_nue==1)*(spline_rlogl_nue<7.5) + (hes_preselect_nue==1)*(bdt_hes_nue>-0.01))*(numpy.cos(splinemod_zen_nue)<0.087)
samp1_final_level_nuge = ((les_preselect_nuge==1)*(spline_rlogl_nuge<7.5) + (hes_preselect_nuge==1)*(bdt_hes_nuge>-0.01))*(numpy.cos(splinemod_zen_nuge)<0.087)
samp1_final_level_c7437 = ((les_preselect_c7437==1)*(spline_rlogl_c7437<7.5) + (hes_preselect_c7437==1)*(bdt_hes_c7437>-0.01))*(numpy.cos(splinemod_zen_c7437)<0.087)
samp1_final_level_c9622 = (les_preselect_c9622==1)*(spline_rlogl_c9622<7.5)*(numpy.cos(splinemod_zen_c9622)<0.087)
samp1_final_level_c9255 = (les_preselect_c9255==1)*(spline_rlogl_c9255<7.5)*(numpy.cos(splinemod_zen_c9255)<0.087)


samp2_final_level_numu = (((les_preselect_numu==1)*(bdt_les_numu>0.0) + (hes_preselect_numu==1)*(bdt_hes_numu>-0.01))*(numpy.cos(splinemod_zen_numu)<0.087))*numpy.isfinite(fr_R_numu)
samp2_final_level_nugmu = (((les_preselect_nugmu==1)*(bdt_les_nugmu>0.0) + (hes_preselect_nugmu==1)*(bdt_hes_nugmu>-0.01))*(numpy.cos(splinemod_zen_nugmu)<0.087))*numpy.isfinite(fr_R_nugmu)
samp2_final_level_data = (((les_preselect_data==1)*(bdt_les_data>0.0) + (hes_preselect_data==1)*(bdt_hes_data>-0.01))*(numpy.cos(splinemod_zen_data)<0.087))*numpy.isfinite(fr_R_data)
samp2_final_level_nue = (((les_preselect_nue==1)*(bdt_les_nue>0.0) + (hes_preselect_nue==1)*(bdt_hes_nue>-0.01))*(numpy.cos(splinemod_zen_nue)<0.087))*numpy.isfinite(fr_R_nue)
samp2_final_level_nuge = (((les_preselect_nuge==1)*(bdt_les_nuge>0.0) + (hes_preselect_nuge==1)*(bdt_hes_nuge>-0.01))*(numpy.cos(splinemod_zen_nuge)<0.087))*numpy.isfinite(fr_R_nuge)
samp2_final_level_c7437 = (((les_preselect_c7437==1)*(bdt_les_c7437>0.0) + (hes_preselect_c7437==1)*(bdt_hes_c7437>-0.01))*(numpy.cos(splinemod_zen_c7437)<0.087))*numpy.isfinite(fr_R_c7437)
samp2_final_level_c9622 = ((les_preselect_c9622==1)*(bdt_les_c9622>0.0)*(numpy.cos(splinemod_zen_c9622)<0.087))*numpy.isfinite(fr_R_c9622)
samp2_final_level_c9255 = ((les_preselect_c9255==1)*(bdt_les_c9255>0.0)*(numpy.cos(splinemod_zen_c9255)<0.087))*numpy.isfinite(fr_R_c9255)


alt_samp2_final_level_numu = (((les_preselect_numu==1)*(alt_bdt_les_numu>0.0) + (hes_preselect_numu==1)*(bdt_hes_numu>-0.01))*(numpy.cos(splinemod_zen_numu)<0.087))*numpy.isfinite(fr_R_numu)
alt_samp2_final_level_nugmu = ((les_preselect_nugmu==1)*(alt_bdt_les_nugmu>0.0) + (hes_preselect_nugmu==1)*(bdt_hes_nugmu>-0.01))*(numpy.cos(splinemod_zen_nugmu)<0.087)
alt_samp2_final_level_data = ((les_preselect_data==1)*(alt_bdt_les_data>0.0) + (hes_preselect_data==1)*(bdt_hes_data>-0.01))*(numpy.cos(splinemod_zen_data)<0.087)
alt_samp2_final_level_nue = ((les_preselect_nue==1)*(alt_bdt_les_nue>0.0) + (hes_preselect_nue==1)*(bdt_hes_nue>-0.01))*(numpy.cos(splinemod_zen_nue)<0.087)
alt_samp2_final_level_nuge = ((les_preselect_nuge==1)*(alt_bdt_les_nuge>0.0) + (hes_preselect_nuge==1)*(bdt_hes_nuge>-0.01))*(numpy.cos(splinemod_zen_nuge)<0.087)
alt_samp2_final_level_c7437 = ((les_preselect_c7437==1)*(alt_bdt_les_c7437>0.0) + (hes_preselect_c7437==1)*(bdt_hes_c7437>-0.01))*(numpy.cos(splinemod_zen_c7437)<0.087)
alt_samp2_final_level_c9622 = (les_preselect_c9622==1)*(alt_bdt_les_c9622>0.0)*(numpy.cos(splinemod_zen_c9622)<0.087)
alt_samp2_final_level_c9255 = (les_preselect_c9255==1)*(alt_bdt_les_c9255>0.0)*(numpy.cos(splinemod_zen_c9255)<0.087)

import scipy.misc 
import scipy
scipy.factorial = scipy.misc.factorial
import scipy.signal
import scipy.stats

samp1_pull_coeffs_numu = pickle.load(open('samp1_pull_coefficients_numu.pkl','r'))
samp1_pull_coeffs_nugmu = pickle.load(open('samp1_pull_coefficients_nugmu.pkl','r'))
samp2_pull_coeffs_numu = pickle.load(open('samp2_pull_coefficients_numu.pkl','r'))
samp2_pull_coeffs_nugmu = pickle.load(open('samp2_pull_coefficients_nugmu.pkl','r'))

#samp2_pull_coeffs_numu_les = pickle.load(open('new_samp2_pull_coefficients_numu_les.pkl','r'))
#samp2_pull_coeffs_nugmu_hes = pickle.load(open('new_samp2_pull_coefficients_nugmu_hes.pkl','r'))

#samp2_pull_coeffs_numu_les = pickle.load(open('OneSig2dGaus_samp2_pull_coefficients_numu_les.pkl','r'))
#samp2_pull_coeffs_nugmu_hes = pickle.load(open('OneSig2dGaus_samp2_pull_coefficients_nugmu_hes.pkl','r'))

samp2_pull_coeffs_numu_les = pickle.load(open('Tailcut_samp2_pull_coefficients_numu_les.pkl','r'))
samp2_pull_coeffs_nugmu_hes = pickle.load(open('Tailcut_samp2_pull_coefficients_nugmu_hes.pkl','r'))



def pullfit_numu(x,p):
        val = p[0]+p[1]*x+numpy.power(p[2]*x,2)+numpy.power(p[3]*x,3)+numpy.power(p[4]*x,4)
        return val
def pullfit_nugmu(x,p):
        val = p[0]+p[1]*x+numpy.power(p[2]*x,2)+numpy.power(p[3]*x,3)+numpy.power(p[4]*x,4)
        return val


def NchBinMedianPull(nch,pull):
  channybins=numpy.linspace(10,nch.max(),nch.max()-9)
  median_pull=[]
  for k in channybins:
        median_pull.append(numpy.median(pull[nch==k]))
  return channybins,median_pull

def RameezPullMethod(nch,pull):
  channybins=numpy.linspace(10,nch.max(),nch.max()-9)
  peak_pull=[]
  for k in channybins:
	peak_pull.append(10**scipy.stats.norm.fit(numpy.log10(pull[nch==k]))[0])
  return channybins,peak_pull

def OneSigma2dGaussianPull(nch,pull):
  channybins=numpy.linspace(10,nch.max(),nch.max()-9)
  sigma_pull=[]
  for k in channybins:
	sigmavalue_2d=0.0
	dpull=0.05
	if len(pull[nch==k]) < 1:
		sigmavalue_2d=numpy.nan
	else:
		while 1.0*len(pull[(nch==k)*(pull<sigmavalue_2d)]) / len(pull[(nch==k)]) < 0.39347:
			sigmavalue_2d+=dpull
	sigma_pull.append(sigmavalue_2d)
  return channybins,sigma_pull

def NchBinMedianRes(nch,res):
  channybins=numpy.linspace(10,nch.max(),nch.max()-9)
  median_res=[]
  for k in channybins:
        median_res.append(numpy.median(res[nch==k]))
  return median_res

def EnergyBinMedianRes(energy,res,binzo):
  channybins=binzo
  de=numpy.diff(channybins)[0]
  median_res=[]
  for k in channybins:
        median_res.append(numpy.median(res[(energy>=k)*(energy<k+de)]))
  return median_res[:-1],channybins[:-1]+de/2.0

samp2_pull_nch_numu,samp2_median_pull_numu = NchBinMedianPull(merged_nch_numu[(samp2_final_level_numu==1)],para_pull_numu[samp2_final_level_numu==1])
samp2_pull_nch_numu,samp2_peak_pull_numu = RameezPullMethod(merged_nch_numu[(samp2_final_level_numu==1)],para_pull_numu[samp2_final_level_numu==1])

samp2_pull_nch_nugmu,samp2_median_pull_nugmu = NchBinMedianPull(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)],para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)])
samp2_pull_nch_nugmu,samp2_peak_pull_nugmu = RameezPullMethod(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)],para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)])

samp2_pull_nch_numu_hes,samp2_median_pull_numu_hes = NchBinMedianPull(merged_nch_numu[(samp2_final_level_numu==1)*(hes_numu==1)],para_pull_numu[(samp2_final_level_numu==1)*(hes_numu==1)])
samp2_pull_nch_numu_hes,samp2_peak_pull_numu_hes = RameezPullMethod(merged_nch_numu[(samp2_final_level_numu==1)*(hes_numu==1)],para_pull_numu[(samp2_final_level_numu==1)*(hes_numu==1)])

samp2_pull_nch_numu_les,samp2_median_pull_numu_les = NchBinMedianPull(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)],para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)])
samp2_pull_nch_numu_les,samp2_peak_pull_numu_les = RameezPullMethod(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)],para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)])

tailcut_samp2_pull_nch_numu_les,tailcut_samp2_median_pull_numu_les = NchBinMedianPull(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)*(para_pull_numu<25)],para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)*(para_pull_numu<25)])
tailcut_samp2_pull_nch_nugmu_hes,tailcut_samp2_median_pull_nugmu_hes = NchBinMedianPull(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)*(para_pull_nugmu<25)],para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)*(para_pull_nugmu<25)])



onesig2d_samp2_pull_nch_numu_les,onesig2d_samp2_median_pull_numu_les = OneSigma2dGaussianPull(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)],para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)])

onesig2d_samp2_pull_nch_nugmu_hes,onesig2d_samp2_median_pull_nugmu_hes = OneSigma2dGaussianPull(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)],para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)])



samp1_median_res_numu = numpy.array(NchBinMedianRes(merged_nch_numu[(samp1_final_level_numu==1)],spline_primary_diff_numu[(samp1_final_level_numu==1)]))
samp1_median_res_nugmu = numpy.array(NchBinMedianRes(merged_nch_nugmu[(samp1_final_level_nugmu==1)],spline_primary_diff_nugmu[(samp1_final_level_nugmu==1)]))

samp2_median_res_numu = numpy.array(NchBinMedianRes(merged_nch_numu[(samp2_final_level_numu==1)],spline_primary_diff_numu[(samp2_final_level_numu==1)]))
samp2_median_res_nugmu = numpy.array(NchBinMedianRes(merged_nch_nugmu[(samp2_final_level_nugmu==1)],spline_primary_diff_nugmu[(samp2_final_level_nugmu==1)]))



def ParaFixNuMu(sigma,nch,lescoeffs,hescoeffs,lesbool,hesbool):
  corrected = sigma * lesbool * (nch<140) * pullfit_numu(nch,lescoeffs) + sigma * hesbool * (nch<140) * pullfit_numu(nch,hescoeffs) + sigma * (nch>=140) * 2.0
  return corrected
def ParaFixNuGMu(sigma,nch,lescoeffs,hescoeffs):
  corrected = sigma * lesbool * (nch<140) * pullfit_numu(nch,lescoeffs) + sigma * hesbool * (nch<140) * pullfit_numu(nch,hescoeffs) + sigma * (nch>=140) * 2.0
  return corrected

def ParaFixC9622(sigma,nch,lescoeffs):
  corrected = sigma * (nch<140) * pullfit_numu(nch,lescoeffs) +  sigma * (nch>140) * 2.0
  return corrected


corrected_samp2_sig_para_nugmu = ParaFixNuMu(sig_para_nugmu,merged_nch_nugmu,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_nugmu,hes_nugmu)
corrected_samp2_sig_para_numu = ParaFixNuMu(sig_para_numu,merged_nch_numu,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_numu,hes_numu)

corrected_samp2_para_pull_numu = spline_primary_diff_numu/corrected_samp2_sig_para_numu
corrected_samp2_para_pull_nugmu = spline_primary_diff_nugmu/corrected_samp2_sig_para_nugmu

corrected_samp2_sig_para_data = ParaFixNuMu(sig_para_data,merged_nch_data,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_data,hes_data)
corrected_samp2_sig_para_c9622 = ParaFixC9622(sig_para_c9622,les_nch_clean_c9622,samp2_pull_coeffs_numu_les)
corrected_samp2_sig_para_c9255 = ParaFixC9622(sig_para_c9255,les_nch_clean_c9255,samp2_pull_coeffs_numu_les)
corrected_samp2_sig_para_c7437 = ParaFixNuMu(sig_para_c7437,merged_nch_c7437,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_c7437,hes_c7437)

corrected_samp2_sig_para_nue = ParaFixNuMu(sig_para_nue,merged_nch_nue,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_nue,hes_nue)
corrected_samp2_sig_para_nuge = ParaFixNuMu(sig_para_nuge,merged_nch_nuge,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_nuge,hes_nuge)

samp2_pull_nch_numu_les, corrected_samp2_OneSig2D_pull_numu_les = OneSigma2dGaussianPull(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)],corrected_samp2_para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)])
samp2_pull_nch_nugmu_hes, corrected_samp2_OneSig2D_pull_nugmu_hes = OneSigma2dGaussianPull(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)],corrected_samp2_para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)])

samp2_pull_nch_numu_les, corrected_samp2_tailcut_pull_numu_les = NchBinMedianPull(merged_nch_numu[(samp2_final_level_numu==1)*(les_numu==1)],corrected_samp2_para_pull_numu[(samp2_final_level_numu==1)*(les_numu==1)])
samp2_pull_nch_nugmu_hes, corrected_samp2_tailcut_pull_nugmu_hes = NchBinMedianPull(merged_nch_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)],corrected_samp2_para_pull_nugmu[(samp2_final_level_nugmu==1)*(hes_nugmu==1)])


'''
import pickle
pickle.dump(onesig2d_samp2_median_pull_numu_les,open("samp2_one2dsig_pull_numu_les.pkl",'w'))
pickle.dump(onesig2d_samp2_median_pull_nugmu_hes,open("samp2_one2dsig_pull_nugmu_hes.pkl",'w'))

pickle.dump(tailcut_samp2_median_pull_numu_les,open("samp2_tailcut_pull_numu_les.pkl",'w'))
pickle.dump(tailcut_samp2_median_pull_nugmu_hes,open("samp2_tailcut_pull_nugmu_hes.pkl",'w'))


'''

real_final_level_data = samp2_final_level_data*((180/numpy.pi)*corrected_samp2_sig_para_data < 45.0)*numpy.isfinite(fr_R_data)
real_final_level_c9622 = samp2_final_level_c9622*((180/numpy.pi)*corrected_samp2_sig_para_c9622 < 45.0)*numpy.isfinite(fr_R_c9622)
real_final_level_c9255 = samp2_final_level_c9255*((180/numpy.pi)*corrected_samp2_sig_para_c9255 < 45.0)*numpy.isfinite(fr_R_c9255)
real_final_level_numu = samp2_final_level_numu*((180/numpy.pi)*corrected_samp2_sig_para_numu < 45.0)*numpy.isfinite(fr_R_numu)
real_final_level_nugmu = samp2_final_level_nugmu*((180/numpy.pi)*corrected_samp2_sig_para_nugmu < 45.0)*numpy.isfinite(fr_R_nugmu)
real_final_level_nuge = samp2_final_level_nuge*((180/numpy.pi)*corrected_samp2_sig_para_nuge < 45.0)*numpy.isfinite(fr_R_nuge)
real_final_level_nue = samp2_final_level_nue*((180/numpy.pi)*corrected_samp2_sig_para_nue < 45.0)*numpy.isfinite(fr_R_nue)
real_final_level_c7437 = samp2_final_level_c7437*((180/numpy.pi)*corrected_samp2_sig_para_c7437 < 45.0)*numpy.isfinite(fr_R_c7437)

real_final_level_median_res_numu = numpy.array(NchBinMedianRes(merged_nch_numu[(real_final_level_numu==1)],spline_primary_diff_numu[(real_final_level_numu==1)]))
real_final_level_median_res_nugmu = numpy.array(NchBinMedianRes(merged_nch_nugmu[(real_final_level_nugmu==1)],spline_primary_diff_nugmu[(real_final_level_nugmu==1)]))

real_final_level_median_paraboloid_numu,energy_res_bins = numpy.array(EnergyBinMedianRes(primary_energy_numu[(real_final_level_numu==1)],corrected_samp2_sig_para_numu[(real_final_level_numu==1)],numpy.linspace(4,190,20)))
real_final_level_median_paraboloid_nugmu,energy_res_bins_nugmu = numpy.array(EnergyBinMedianRes(primary_energy_nugmu[(real_final_level_nugmu==1)],corrected_samp2_sig_para_nugmu[(real_final_level_nugmu==1)],numpy.linspace(10,1000,40)))

real_final_level_median_res_energy_numu,energy_res_bins = numpy.array(EnergyBinMedianRes(primary_energy_numu[(real_final_level_numu==1)],spline_primary_diff_numu[(real_final_level_numu==1)],numpy.linspace(4,190,20)))
real_final_level_median_res_energy_nugmu,energy_res_bins_nugmu = numpy.array(EnergyBinMedianRes(primary_energy_nugmu[(real_final_level_nugmu==1)],spline_primary_diff_nugmu[(real_final_level_nugmu==1)],numpy.linspace(10,1000,40)))


les_real_final_level_median_res_numu = numpy.array(NchBinMedianRes(merged_nch_numu[(les_numu==1)*(real_final_level_numu==1)],spline_primary_diff_numu[(les_numu==1)*(real_final_level_numu==1)]))
les_real_final_level_median_res_nugmu = numpy.array(NchBinMedianRes(merged_nch_nugmu[(les_nugmu==1)*(real_final_level_nugmu==1)],spline_primary_diff_nugmu[(les_nugmu==1)*(real_final_level_nugmu==1)]))
hes_real_final_level_median_res_numu = numpy.array(NchBinMedianRes(merged_nch_numu[(hes_numu==1)*(real_final_level_numu==1)],spline_primary_diff_numu[(hes_numu==1)*(real_final_level_numu==1)]))
hes_real_final_level_median_res_nugmu = numpy.array(NchBinMedianRes(merged_nch_nugmu[(hes_nugmu==1)*(real_final_level_nugmu==1)],spline_primary_diff_nugmu[(hes_nugmu==1)*(real_final_level_nugmu==1)]))



#from scipy.optimize import curvefit

#def pullfit_numu(x,a,b,c,d,e):
#        val = a+b*x+numpy.power(c*x,2)+numpy.power(d*x,3)+numpy.power(e*x,4)
#        return val
#def pullfit_nugmu(x,a,b,c,d,e):
#        val = a+b*x+numpy.power(c*x,2)+numpy.power(d*x,3)+numpy.power(e*x,4)
#        return val

#coeff_numu,covar_numu = curve_fit(pullfit_numu, nchbin_numu[nchbin_numu>15], med_pull_numu[nchbin_numu>15], p0=(4,0.1,0.025,0.005,0.0001), sigma=None)
#coeff_nugmu,covar_nugmu = curve_fit(pullfit_nugmu, nchbin_nugmu[nchbin_nugmu>15], med_pull_nugmu[nchbin_nugmu>15], p0=(4,0.1,0.025,0.005,0.0001), sigma=None)


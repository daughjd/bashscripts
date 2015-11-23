from __future__ import print_function
import sys,tables, numpy, pylab
from pybdt.ml import DataSet
from pybdt.util import mkdir, save



### Acquire CORSIKA table info ###
datums = tables.openFile('Corsika_7437_0000-3999_Level3.L6Out.hdf5')

les_c7437 = datums.root.L4Bool_LES.col('value')
hes_c7437 = datums.root.L4Bool_HES.col('value')
#bdt_les_c7437 = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_c7437 = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_c7437 = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_c7437 = datums.root.SplineMPEMod.col('azimuth')
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
sig_para_c7437 = numpy.sqrt(para1_c7437**2 + para2_c7437**2)

datums.close()

### Acquire Data Table Info ###

datums = tables.openFile('Level5_IC86.2012.CompleteFullData.L6Out.330Days.hdf5')

les_data = datums.root.L4Bool_LES.col('value')
hes_data = datums.root.L4Bool_HES.col('value')
#bdt_les_data = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_data = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_data = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_data = datums.root.SplineMPEMod.col('azimuth')
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
sig_para_data = numpy.sqrt(para1_data**2 + para2_data**2)


datums.close()


### Acquire NuMu Table Info ###
datums = tables.openFile('genie_ic.1304_1404.combo_numu_numubar.L5Only.L6Out.hdf5')

les_numu = datums.root.L4Bool_LES.col('value')
hes_numu = datums.root.L4Bool_HES.col('value')
#bdt_les_numu = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_numu = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_numu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_numu = datums.root.SplineMPEMod.col('azimuth')
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
#bdtscore_numu = datums.root.LES_L5_BDTScore.col('value')
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
#hes_nch_clean_numu = datums.root.L4_Variables_HES.col('NCh_Clean')
#hes_hlc_numu = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_numu = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_numu = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_numu = numpy.sqrt(para1_numu**2 + para2_numu**2)

intpartx_numu = datums.root.InteractionParticle.col('x')
intparty_numu = datums.root.InteractionParticle.col('y')
intpartz_numu = datums.root.InteractionParticle.col('z')
chkgrbeventw_numu = datums.root.ChkGRBEventWeight.col('value')
primary_energy_numu = datums.root.PrimaryNu.col('energy')
primary_zenith_numu = datums.root.PrimaryNu.col('zenith')
atmo_numu = datums.root.AtmoWeight.col('value')
cc_numu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()


### Acquire NuMu Nugen ###

datums = tables.openFile('Level5_nugen_numu_IC86.2013.AllSubs.010090.L6Out.hdf5')

les_nugmu = datums.root.L4Bool_LES.col('value')
hes_nugmu = datums.root.L4Bool_HES.col('value')
#bdt_les_nugmu = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_nugmu = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_nugmu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nugmu = datums.root.SplineMPEMod.col('azimuth')
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
#bdtscore_nugmu = datums.root.LES_L5_BDTScore.col('value')
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
#hes_nch_clean_nugmu = datums.root.L4_Variables_HES.col('NCh_Clean')
#hes_hlc_nugmu = datums.root.L4_Variables_HES.col('Nch_HLC')

para1_nugmu = datums.root.SplineMPEModParaboloidFitParams.col('err1')
para2_nugmu = datums.root.SplineMPEModParaboloidFitParams.col('err2')
sig_para_nugmu = numpy.sqrt(para1_nugmu**2 + para2_nugmu**2)

intpartx_nugmu = datums.root.InteractionParticle.col('x')
intparty_nugmu = datums.root.InteractionParticle.col('y')
intpartz_nugmu = datums.root.InteractionParticle.col('z')
primary_energy_nugmu = datums.root.PrimaryNu.col('energy')
primary_zenith_nugmu = datums.root.PrimaryNu.col('zenith')
#indc_nugmu = datums.root.InDC.col('value')
#intldc_nugmu = datums.root.InExpDC.col('value')
atmo_nugmu = datums.root.AtmoWeight.col('value')
cc_nugmu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()

##############################################

dataweight=1.0/(28475585.0)
dataweight=dataweight*numpy.ones(len(les_data))

C7437EventWeight=(10.03598*3991)**-1
C7437EventWeight=C7437EventWeight*numpy.ones(len(les_c7437))

def AngleBetweenAngles(theta1,phi1,theta2,phi2):
   x1 = numpy.sin(theta1)*numpy.cos(phi1)
   y1 = numpy.sin(theta1)*numpy.sin(phi1)
   z1 = numpy.cos(theta1)
   x2 = numpy.sin(theta2)*numpy.cos(phi2)
   y2 = numpy.sin(theta2)*numpy.sin(phi2)
   z2 = numpy.cos(theta2)
   return numpy.arccos(x1*x2+y1*y2+z1*z2)

lf_spline_diff_data = AngleBetweenAngles(linefit_zen_data,linefit_azi_data,splinemod_zen_data,splinemod_azi_data)
lf_spline_diff_c7437 = AngleBetweenAngles(linefit_zen_c7437,linefit_azi_c7437,splinemod_zen_c7437,splinemod_azi_c7437)
lf_spline_diff_numu = AngleBetweenAngles(linefit_zen_numu,linefit_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
lf_spline_diff_nugmu = AngleBetweenAngles(linefit_zen_nugmu,linefit_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)

### BDT TIME ###

hes_preselect_numu = (hes_numu==1)*(numpy.isfinite(avg_distq_numu))*(numpy.isfinite(dhd_numu))*(numpy.isfinite(spline_rlogl_numu))
hes_preselect_nugmu = (hes_nugmu==1)*(numpy.isfinite(avg_distq_nugmu))*(numpy.isfinite(dhd_nugmu))*(numpy.isfinite(spline_rlogl_nugmu))
hes_preselect_data = (hes_data==1)*(numpy.isfinite(avg_distq_data))*(numpy.isfinite(dhd_data))*(numpy.isfinite(spline_rlogl_data))
hes_preselect_c7437 = (hes_c7437==1)*(numpy.isfinite(avg_distq_c7437))*(numpy.isfinite(dhd_c7437))*(numpy.isfinite(spline_rlogl_c7437))

quality_nugmu = cc_nugmu
quality_numu = cc_numu
hes_preselect_numu = hes_preselect_numu*quality_numu
hes_preselect_nugmu = hes_preselect_nugmu*quality_nugmu


print ("Data rate after preselection: %f" % dataweight[hes_preselect_data].sum())

def make_data():
  a = lf_spline_diff_data[hes_preselect_data]
  b = ldird_data[hes_preselect_data]
  c = dhd_data[hes_preselect_data]
  d = avg_distq_data[hes_preselect_data]
  e = spline_rlogl_data[hes_preselect_data]
  weight = dataweight[hes_preselect_data]
  return dict (a=a, b=b, c=c, d=d, e=e, weight=weight)

def make_bgsim():
  a = lf_spline_diff_c7437[hes_preselect_c7437]
  b = ldird_c7437[hes_preselect_c7437]
  c = dhd_c7437[hes_preselect_c7437]
  d = avg_distq_c7437[hes_preselect_c7437]
  e = spline_rlogl_c7437[hes_preselect_c7437]
  weight = C7437EventWeight[hes_preselect_c7437]
  return dict (a=a, b=b, c=c, d=d, e=e, weight=weight)

def make_sig_sim():
  a = lf_spline_diff_nugmu[hes_preselect_nugmu]
  b = ldird_nugmu[hes_preselect_nugmu]
  c = dhd_nugmu[hes_preselect_nugmu]
  d = avg_distq_nugmu[hes_preselect_nugmu]
  e = spline_rlogl_nugmu[hes_preselect_nugmu]
  weight = primary_energy_nugmu[hes_preselect_nugmu]
  return dict (a=a, b=b, c=c, d=d, e=e, weight=weight)

### Make Training Samples ###
train_sig_sim = make_sig_sim ()
test_sig_sim = make_sig_sim ()
bg_sim = make_bgsim ()

### Datamerize ###
train_data = make_data ()
test_data = make_data ()

### Make them Sets ###
train_sig = DataSet (train_sig_sim)
train_data = DataSet (train_data)
test_sig = DataSet (test_sig_sim)
test_data = DataSet (test_data)

bg = DataSet (bg_sim)

# create a data directory, if one does not already exist
mkdir ('hesdata')

print ('Saving data in hesdata/ ...')

# save these files
save (train_sig, 'hesdata/train_sig.ds')
save (train_data, 'hesdata/train_data.ds')
save (test_sig, 'hesdata/test_sig.ds')
save (test_data, 'hesdata/test_data.ds')
save (bg, 'hesdata/bg.ds')


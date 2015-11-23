from __future__ import print_function
import sys,tables, numpy
from pybdt.ml import DataSet
from pybdt.util import mkdir, save



### Acquire CORSIKA table info ###
datums = tables.openFile('Level4_IC86.2011_corsika.009622.100K.L6Out.hdf5')

les_c9622 = datums.root.L4Bool_LES.col('value')
#hes_c9622 = datums.root.L4Bool_HES.col('value')
#bdt_les_c9622 = datums.root.LES_L6_BDTScore.col('value')
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
#linefit_zen_c9622 = datums.root.LineFit_DC.col('zenith')
#linefit_azi_c9622 = datums.root.LineFit_DC.col('azimuth')
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
sig_para_c9622 = numpy.sqrt(para1_c9622**2 + para2_c9622**2)

cweight_9622 = datums.root.CorsikaWeightMap.col('Weight')
polyweight_9622 = datums.root.CorsikaWeightMap.col('Polygonato')
ts_9622 = datums.root.CorsikaWeightMap.col("TimeScale")


datums.close()

### Acquire Data Table Info ###

datums = tables.openFile('Level5_IC86.2012.CompleteFullData.L6Out.330Days.hdf5')

les_data = datums.root.L4Bool_LES.col('value')
hes_data = datums.root.L4Bool_HES.col('value')
#bdt_les_data = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_data = datums.root.HES_L6_BDTScore.col('value')
splinemod_zen_data = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_data = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_data = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_data = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_data = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_data = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_data = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_data = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_data = datums.root.SplineMPEModFitParams.col('logl')
#linefit_zen_data = datums.root.LineFit_DC.col('zenith')
#linefit_azi_data = datums.root.LineFit_DC.col('azimuth')
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
spline_para_zen_numu = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_numu = datums.root.SplineMPEModParaboloid.col('azimuth')
avg_distq_numu = datums.root.SplineMPEMod_Characteristics.col('avg_dom_dist_q_tot_dom')
lempty_numu = datums.root.SplineMPEMod_Characteristics.col('empty_hits_track_length')
separation_numu = datums.root.SplineMPEMod_Characteristics.col('track_hits_separation_length')
spline_rlogl_numu = datums.root.SplineMPEModFitParams.col('rlogl')
spline_logl_numu = datums.root.SplineMPEModFitParams.col('logl')
#linefit_zen_numu = datums.root.LineFit_DC.col('zenith')
#linefit_azi_numu = datums.root.LineFit_DC.col('azimuth')
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
primary_azimuth_numu = datums.root.PrimaryNu.col('azimuth')
int_zen_numu = datums.root.InteractionParticle.col('zenith')
int_azi_numu = datums.root.InteractionParticle.col('azimuth')
atmo_numu = datums.root.AtmoWeight.col('value')
cc_numu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()


### Acquire NuMu Nugen ###

datums = tables.openFile('Level5_nugen_numu_IC86.2013.010090.AllSubs.L6Out.hdf5')

les_nugmu = datums.root.L4Bool_LES.col('value')
hes_nugmu = datums.root.L4Bool_HES.col('value')
#bdt_les_nugmu = datums.root.LES_L6_BDTScore.col('value')
#bdt_hes_nugmu = datums.root.HES_L6_BDTScore.col('value')
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
primary_azimuth_nugmu = datums.root.PrimaryNu.col('azimuth')
int_zen_nugmu = datums.root.InteractionParticle.col('zenith')
int_azi_nugmu = datums.root.InteractionParticle.col('azimuth')
#indc_nugmu = datums.root.InDC.col('value')
#intldc_nugmu = datums.root.InExpDC.col('value')
atmo_nugmu = datums.root.AtmoWeight.col('value')
cc_nugmu = numpy.abs(datums.root.InteractionParticle.col('pdg_encoding')) ==13

datums.close()

##############################################

dataweight=1.0/(28475585.0)
dataweight=dataweight*numpy.ones(len(les_data))

C9622EventWeight=cweight_9622*polyweight_9622/ts_9622
C9622EventWeight=C9622EventWeight/100000.


fr_R_c9622 = numpy.sqrt((finite_x_c9622-46)**2 + (finite_y_c9622+35)**2)
fr_R_numu = numpy.sqrt((finite_x_numu-46)**2 + (finite_y_numu+35)**2)
fr_R_nugmu = numpy.sqrt((finite_x_nugmu-46)**2 + (finite_y_nugmu+35)**2)
fr_R_data = numpy.sqrt((finite_x_data-46)**2 + (finite_y_data+35)**2)

def para_nch_2dcut(para,nch):
   bool_array=[]
   indexical=0
   for i in nch:
     boo = False
     if i <= 30:
       if para[indexical] <= -0.682*i+35.45:
          boo = True
     elif i <= 52:
       if para[indexical] <= -.454545*i+28.63:
          boo = True
     else:
       if para[indexical] <= 5:
          boo = True
     bool_array.append(boo)
     indexical+=1
   bool_array = numpy.array(bool_array)
   return bool_array

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
spline_paraspline_diff_data = AngleBetweenAngles(spline_para_zen_data,spline_para_azi_data,splinemod_zen_data,splinemod_azi_data)
spline_paraspline_diff_c9622 = AngleBetweenAngles(spline_para_zen_c9622,spline_para_azi_c9622,splinemod_zen_c9622,splinemod_azi_c9622)

spline_daughter_diff_numu = AngleBetweenAngles(int_zen_numu,int_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
spline_daughter_diff_nugmu = AngleBetweenAngles(int_zen_nugmu,int_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)

badly_recoed_numu = spline_daughter_diff_numu > 5/57.3
badly_recoed_nugmu = spline_daughter_diff_nugmu > 5/57.3

#para_nch_2dbool_numu = para_nch_2dcut(57.3*sig_para_numu,nch_clean_numu)
#para_nch_2dbool_nugmu = para_nch_2dcut(57.3*sig_para_nugmu,nch_clean_nugmu)
#para_nch_2dbool_c7437 = para_nch_2dcut(57.3*sig_para_c7437,nch_clean_c7437)
#para_nch_2dbool_data = para_nch_2dcut(57.3*sig_para_data,nch_clean_data)

### BDT TIME ###

les_preselect_data = numpy.isfinite(finite_z_data)*(les_data==1)*numpy.isfinite(avg_distq_data)
les_preselect_c9622 = numpy.isfinite(finite_z_c9622)*(les_c9622==1)*numpy.isfinite(avg_distq_c9622)
les_preselect_numu = numpy.isfinite(finite_z_numu)*(les_numu==1)*numpy.isfinite(avg_distq_numu)
les_preselect_nugmu = numpy.isfinite(finite_z_nugmu)*(les_nugmu==1)*numpy.isfinite(avg_distq_nugmu)

#good_angular_numu = (les_numu==1)*(spline_rlogl_numu<7.5)*(spline_paraspline_diff_numu<(25/57.3))
#good_angular_nugmu = (les_nugmu==1)*(spline_rlogl_nugmu<7.5)*(spline_paraspline_diff_nugmu<(25/57.3))
#good_angular_c9622 = (les_c9622==1)*(spline_rlogl_c9622<7.5)*(spline_paraspline_diff_c9622<(25/57.3))
#good_angular_data = (les_data==1)*(spline_rlogl_data<7.5)*(spline_paraspline_diff_data<(25/57.3))

good_angular_numu = (les_numu==1)*(badly_recoed_numu==0)
good_angular_nugmu = (les_nugmu==1)*(badly_recoed_nugmu==0)

quality_nugmu = cc_nugmu
quality_numu = cc_numu
les_preselect_numu = les_preselect_numu*quality_numu*good_angular_numu
les_preselect_nugmu = les_preselect_nugmu*quality_nugmu*good_angular_nugmu


print ("Data rate after preselection: %f" % dataweight[les_preselect_data].sum())

def make_data():
  a = fr_R_data[les_preselect_data]
  b = finite_z_data[les_preselect_data]
  c = dhd_data[les_preselect_data]
  d = avg_distq_data[les_preselect_data] 
  e = spline_rlogl_data[les_preselect_data]
  f = vetotrackcharge_data[les_preselect_data]
  g = spline_paraspline_diff_data[les_preselect_data] 
  weight = dataweight[les_preselect_data]
  return dict (a=a, b=b, c=c, d=d, e=e, f=f, g=g,weight=weight)

def make_bgsim():
  a = fr_R_c9622[les_preselect_c9622]
  b = finite_z_c9622[les_preselect_c9622]
  c = dhd_c9622[les_preselect_c9622]
  d = avg_distq_c9622[les_preselect_c9622]
  e = spline_rlogl_c9622[les_preselect_c9622]
  f = vetotrackcharge_c9622[les_preselect_c9622]
  g = spline_paraspline_diff_c9622[les_preselect_c9622]
  weight = C9622EventWeight[les_preselect_c9622]
  return dict (a=a, b=b, c=c, d=d, e=e, f=f, g=g,weight=weight)

def make_sig_sim():
  a = fr_R_numu[les_preselect_numu]
  b = finite_z_numu[les_preselect_numu]
  c = dhd_numu[les_preselect_numu]
  d = avg_distq_numu[les_preselect_numu]
  e = spline_rlogl_numu[les_preselect_numu]
  f = vetotrackcharge_numu[les_preselect_numu]
  g = spline_paraspline_diff_numu[les_preselect_numu]
  weight = primary_energy_numu[les_preselect_numu]**0
  return dict (a=a, b=b, c=c, d=d, e=e, f=f, g=g,weight=weight)

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
mkdir ('data')

print ('Saving data in data/ ...')

# save these files
save (train_sig, 'data/train_sig_thirds.ds')
save (train_data, 'data/train_data_thirds.ds')
save (test_sig, 'data/test_sig_thirds.ds')
save (test_data, 'data/test_data_thirds.ds')
save (bg, 'data/bg_thirds.ds')






import sys,tables, numpy,pickle,pylab

from subprocess import call


### Acquire Data Table Info ###

#datums = tables.openFile('third_time/bdt_scored/DCStream_IC86.2012_data_Run00120810_11_12_AllSubs.L5Out.hdf5')
datums  = tables.openFile('third_time/Level5_IC86.2012.CompleteDataSet.L6Out.330Days_Third.hdf5')

dc12_data = []
for i in datums.root.FilterMask.col('DeepCoreFilter_12'):
  dc12_data.append(i[0])
dc12_data = numpy.array(dc12_data)

tldc12_data = []
for i in datums.root.FilterMask.col('DeepCoreFilter_TwoLayerExp_12'):
  tldc12_data.append(i[0])
tldc12_data = numpy.array(tldc12_data)

muon12_data = []
for i in datums.root.FilterMask.col('MuonFilter_12'):
  muon12_data.append(i[0])
muon12_data = numpy.array(muon12_data)

moon12_data = []
for i in datums.root.FilterMask.col('MoonFilter_12'):
  moon12_data.append(i[0])
moon12_data = numpy.array(moon12_data)

fss12_data = []
for i in datums.root.FilterMask.col('FSSFilter_12'):
  fss12_data.append(i[0])
fss12_data = numpy.array(fss12_data)

gcf12_data = []
for i in datums.root.FilterMask.col('GCFilter_12'):
  gcf12_data.append(i[0])
gcf12_data = numpy.array(gcf12_data)

lowup12_data = []
for i in datums.root.FilterMask.col('LowUp_12'):
  lowup12_data.append(i[0])
lowup12_data = numpy.array(lowup12_data)

les_data = datums.root.L4Bool_LES.col('value')
hes_data = datums.root.L4Bool_HES.col('value')
bdt_les_data = datums.root.LES_L6_BDTScore.col('value')
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

dataweight=1.0/(28475585.0)
dataweight=dataweight*numpy.ones(len(les_data))

les_nch_clean_data[les_nch_clean_data==-1] = 0
hes_nch_clean_data[hes_nch_clean_data==-1] = 0
merged_nch_data = les_nch_clean_data + hes_nch_clean_data


fr_R_data = numpy.sqrt((finite_x_data-46)**2 + (finite_y_data+35)**2)
spline_plogl_data = spline_logl_data/(merged_nch_data - 2.5)

def AngleBetweenAngles(theta1,phi1,theta2,phi2):
   x1 = numpy.sin(theta1)*numpy.cos(phi1)
   y1 = numpy.sin(theta1)*numpy.sin(phi1)
   z1 = numpy.cos(theta1)
   x2 = numpy.sin(theta2)*numpy.cos(phi2)
   y2 = numpy.sin(theta2)*numpy.sin(phi2)
   z2 = numpy.cos(theta2)
   return numpy.arccos(x1*x2+y1*y2+z1*z2)

spline_paraspline_diff_data = AngleBetweenAngles(spline_para_zen_data,spline_para_azi_data,splinemod_zen_data,splinemod_azi_data)

lf_spline_diff_data = AngleBetweenAngles(linefit_zen_data,linefit_azi_data,splinemod_zen_data,splinemod_azi_data)

les_preselect_data = (les_data==1)*(numpy.isfinite(avg_distq_data))*(numpy.isfinite(finite_x_data))

hes_preselect_data = (hes_data==1)*(numpy.isfinite(avg_distq_data))*(numpy.isfinite(dhd_data))*(numpy.isfinite(spline_rlogl_data))

good_angular_data = (les_preselect_data==1)*(spline_rlogl_data<7.5)*(spline_paraspline_diff_data<(25/57.3))

samp1_final_level_data = ((les_preselect_data==1)*(spline_rlogl_data<7.5) + (hes_preselect_data==1)*(bdt_hes_data>-0.01))*(numpy.cos(splinemod_zen_data)<0.087)

samp2_final_level_data = ((les_preselect_data==1)*(bdt_les_data>0.0) + (hes_preselect_data==1)*(bdt_hes_data>-0.01))*(numpy.cos(splinemod_zen_data)<0.087)

import scipy.misc 
import scipy
scipy.factorial = scipy.misc.factorial
import scipy.signal
import scipy.stats

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



def ParaFixNuMu(sigma,nch,lescoeffs,hescoeffs,lesbool,hesbool):
  corrected = sigma * lesbool * (nch<140) * pullfit_numu(nch,lescoeffs) + sigma * hesbool * (nch<140) * pullfit_numu(nch,hescoeffs) + sigma * (nch>=140) * 2.0
  return corrected
def ParaFixNuGMu(sigma,nch,lescoeffs,hescoeffs):
  corrected = sigma * lesbool * (nch<140) * pullfit_numu(nch,lescoeffs) + sigma * hesbool * (nch<140) * pullfit_numu(nch,hescoeffs) + sigma * (nch>=140) * 2.0
  return corrected

def ParaFixC9622(sigma,nch,lescoeffs):
  corrected = sigma * (nch<140) * pullfit_numu(nch,lescoeffs) +  sigma * (nch>140) * 2.0
  return corrected


corrected_samp2_sig_para_data = ParaFixNuMu(sig_para_data,merged_nch_data,samp2_pull_coeffs_numu_les,samp2_pull_coeffs_nugmu_hes,les_data,hes_data)


real_final_level_data = samp2_final_level_data*((180/numpy.pi)*corrected_samp2_sig_para_data < 45.0)*numpy.isfinite(fr_R_data)


final_level_tlexclusive = (dc12_data==0)*(real_final_level_data==1)*(tldc12_data==1)

any_other_filter = (muon12_data==1)+(moon12_data==1)+(lowup12_data==1)


#!/usr/bin/env python

from I3Tray import *

from os.path import expandvars
from icecube import dataclasses, dataio, icetray, gulliver, lilliput, linefit,tableio,cramer_rao
from optparse import OptionParser
from icecube.common_variables import direct_hits
from icecube import MuonVariables
from icecube import rootwriter
from icecube import hdfwriter
from icecube import finiteReco
from icecube import photonics_service
from icecube import paraboloid
from icecube import millipede
from icecube.common_variables import direct_hits
from icecube.common_variables import track_characteristics
import pickle

import os
import sys

import numpy


import L4Segments

parser = OptionParser()

parser.add_option("-i", "--input", action="store",
        type="string", default="", dest="infile",
        help="Input i3 file(s)  (use comma separated list for multiple files)")

#parser.add_option("-g", "--gcd", action="store",
#        type="string", default="", dest="gcdfile",
#        help="GCD file for input i3 file")

#parser.add_option("-o", "--output", action="store",
#        type="string", default="automate", dest="outfile",
#        help="Output i3 file")

parser.add_option("--sample", action="store",
	type="int", default=1, dest="sample", help="Which final event sample?")

parser.add_option("-n", "--num", action="store",
        type="int", default=-1, dest="num",
        help="Number of frames to process"
	)
parser.add_option("--henugen",action="store_true",
		default=False,dest="he_nugen", help="Only Nugen > 190 GeV for combination with GENIE")

#parser.add_option("-t","--tableoutput", action="store_true",
#        default=False, dest="hdfout", help="Generate hdf5 file?")

(options,args) = parser.parse_args()

### Read in files ; generate output file names#

#infiles = [options.infile]
infiles = ['/nv/hp11/jdaughhetee3/data2/processed/official_L6/wombo_combo/Level5_nugen_numu_IC86.2013.010090.AllSubs.L6Out_HEOnly.i3','/nv/hp11/jdaughhetee3/data2/processed/official_L6/wombo_combo/Genie_NuMu_Combo_AllSubs.L6Out.i3']
infile = infiles[0]

tray = I3Tray()

def LES(frame):
        if "LES" in frame:
                return frame["LES"].value
        return False
def HES(frame):
        if "HES" in frame:
                return frame["HES"].value
        return False


tray.AddModule("I3Reader","i3reader")(
    ("FilenameList", infiles),
    ("SkipKeys", ['TWOfflinePulsesExpDCFidCleanedKeys','RecoChain_DipoleFit_DC_rusage','OfflinePulsesSLC','CascadeFilter_13','CascadeFilter_12']),
    )

from icecube.pybdtmodule import PyBDTModule
from pybdt.util import load
lesbdt_locale='/nv/hp11/jdaughhetee3/data2/software/bdt_files/les_final_level/third_les_l6_spline.bdt'
hesbdt_locale='/nv/hp11/jdaughhetee3/data2/software/bdt_files/hes_final_level/hes_l6_spline.bdt'

def AngleBetweenAngles(theta1,phi1,theta2,phi2):
   x1 = numpy.sin(theta1)*numpy.cos(phi1)
   y1 = numpy.sin(theta1)*numpy.sin(phi1)
   z1 = numpy.cos(theta1)
   x2 = numpy.sin(theta2)*numpy.cos(phi2)
   y2 = numpy.sin(theta2)*numpy.sin(phi2)
   z2 = numpy.cos(theta2)
   return numpy.arccos(x1*x2+y1*y2+z1*z2)

def lesvarsfunc (frame):
  if not frame['L4Bool_LES']:
        a = 0
        b = 0
        c = 0
        d = 0
        e = 0
        f = 0
        g = 0
        out = dict (a=a, b=b, c=c, d=d, e=e, f=f, g=g)
        return out
  a = numpy.sqrt((frame['SplineMPEMod_Contained'].pos.x-46)**2 + (frame['SplineMPEMod_Contained'].pos.y+35)**2)
  b = frame['SplineMPEMod_Contained'].pos.z
  c = frame['SplineMPEMod_DirectHitsD'].n_dir_pulses
  d = frame['SplineMPEMod_Characteristics'].avg_dom_dist_q_tot_dom
  e = frame['SplineMPEModFitParams'].rlogl
  f = frame['L5VetoTrackVetoCharge'].value
  g = AngleBetweenAngles(frame['SplineMPEModParaboloid'].dir.zenith,frame['SplineMPEModParaboloid'].dir.azimuth,frame['SplineMPEMod'].dir.zenith,frame['SplineMPEMod'].dir.azimuth)
  out = dict (a=a, b=b, c=c, d=d, e=e, f=f, g=g)
  return out

def hesvarsfunc (frame):
  if not frame['L4Bool_HES']:
        a = 0
        b = 0
        c = 0
        d = 0
        e = 0
        out = dict (a=a, b=b, c=c, d=d, e=e)
        return out
  a = AngleBetweenAngles(frame['LineFit_DC'].dir.zenith,frame['LineFit_DC'].dir.azimuth,frame['SplineMPEMod'].dir.zenith,frame['SplineMPEMod'].dir.azimuth)
  b = frame['SplineMPEMod_DirectHitsD'].dir_track_length
  c = frame['SplineMPEMod_DirectHitsD'].n_dir_pulses
  d = frame['SplineMPEMod_Characteristics'].avg_dom_dist_q_tot_dom
  e = frame['SplineMPEModFitParams'].rlogl
  out = dict (a=a, b=b, c=c, d=d, e=e)
  return out


tray.AddModule (PyBDTModule, 'lesl6bdt',
      BDTFilename=lesbdt_locale,
      varsfunc=lesvarsfunc,
      OutputName='LES_L6_BDTScore')
if not options.infile.__contains__('009622'):
  tray.AddModule (PyBDTModule, 'hesl6bdt',
      BDTFilename=hesbdt_locale,
      varsfunc=hesvarsfunc,
      OutputName='HES_L6_BDTScore')

### Filter Out ###

def Samp1Filter(frame):
	hes_bool=False
	les_bool=False
	final_bool=False
	if HES(frame):
		if frame['HES_L6_BDTScore'].value > -0.01:
			hes_bool = True
	if LES(frame):
		if (numpy.isfinite(frame['SplineMPEMod_Characteristics'].avg_dom_dist_q_tot_dom) and
		    numpy.isfinite(frame['SplineMPEMod_Contained'].pos.z) and
		    frame['SplineMPEModFitParams'].rlogl < 7.5):
			les_bool = True
	if hes_bool or les_bool:
		final_bool = True
        if numpy.cos(frame['SplineMPEMod'].dir.zenith) >= 0.087:
                final_bool = False
	frame['Final_Level_Bool'] = icetray.I3Bool(final_bool)
	return

def Samp2Filter(frame):
        hes_bool=False
        les_bool=False
        final_bool=False
        if HES(frame):
                if frame['HES_L6_BDTScore'].value > -0.01:
                        hes_bool = True
        if LES(frame):
                if (numpy.isfinite(frame['SplineMPEMod_Characteristics'].avg_dom_dist_q_tot_dom) and
                    numpy.isfinite(frame['SplineMPEMod_Contained'].pos.z) and
                    frame['LES_L6_BDTScore'].value > 0.0):
                        les_bool = True
        if hes_bool or les_bool:
                final_bool = True
	if numpy.cos(frame['SplineMPEMod'].dir.zenith) >= 0.087:
		final_bool = False
        frame['Final_Level_Bool'] = icetray.I3Bool(final_bool)
        return

def FinalFilter(frame):
	decision = False
	if 'Final_Level_Bool' in frame:
		return frame['Final_Level_Bool'].value
	return decision

if options.sample == 1:
	tray.AddModule(Samp1Filter,'sampy1')
	tray.AddModule(FinalFilter,'TheEnd')
if options.sample == 2:
	tray.AddModule(Samp2Filter,'sampy2')
	tray.AddModule(FinalFilter,'TheEnd2')

def Above190_Only(frame):
	boolio = False
	if 'PrimaryNu' in frame:
		if frame['PrimaryNu'].energy >189.999:
			boolio = True
	if not 'PrimaryNu' in frame:
		print "Uh oh..."
	return

if options.he_nugen:
	tray.AddModule(Above190_Only,'Energetics')


def UnifiedNch(frame):
	frame['FinalLevelNch'] = icetray.I3Int(len(frame['LET_RecoPulses'].apply(frame)))
	return
tray.AddModule(UnifiedNch,'unity')

### Add a Resolution Frame Object (Dummy Value) ###
samp1_coeffs = pickle.load(open('/nv/hp11/jdaughhetee3/data2/software/splines/samp1_pull_coefficients_numu.pkl','r'))
samp2_coeffs = pickle.load(open('/nv/hp11/jdaughhetee3/data2/software/splines/samp2_pull_coefficients_numu.pkl','r'))

def pullfit_numu(x,p):
        val = p[0]+p[1]*x+numpy.power(p[2]*x,2)+numpy.power(p[3]*x,3)+numpy.power(p[4]*x,4)
        return val

def RecoSigma(frame,sample):
	reco_zen = frame['SplineMPEMod'].dir.zenith
	chan = len(frame['LET_RecoPulses'].apply(frame))
	if sample == 1:
		ceffs=samp1_coeffs
	else:
		ceffs=samp2_coeffs
	para1 = frame['SplineMPEModParaboloidFitParams'].pbfErr1
	para2 = frame['SplineMPEModParaboloidFitParams'].pbfErr2
	para_sig = numpy.sqrt(para1**2 + para2**2) / numpy.sqrt(2)
	fixed_sig = para_sig * pullfit_numu(chan,ceffs)
	frame["ReScaled_Paraboloid_Sigma_SplineMPEMod"] = dataclasses.I3Double(fixed_sig)
	frame["OneWeight"] = dataclasses.I3Double(frame["I3MCWeightDict"]["OneWeight"])
	return True

tray.AddModule(RecoSigma,'sig',sample=options.sample)


KeysForTable=['PrimaryNu',
		      'FinalLevelNch',
                      'LET_RecoPulses',
                      'I3EventHeader',
                      'ChkGRBEventWeight',
                      'ChkGRBSpectrumWeight',
                      'SplineMPEMod',
                      'SplineMPEMod_DC_CramerRao',
		      'Final_Level_Bool',
                      'SplineMPEModParaboloid',
                      'SplineMPEMod_Contained',
                      'NchHLC',
                      'ReScaled_Paraboloid_Sigma_SplineMPEMod',
                      'AtmoWeight',
		      'OneWeight']

#taboutfile = options.infile.replace('.i3','.Sample'+str(options.sample)+'.root')
taboutfile = "/nv/hp11/jdaughhetee3/data2/processed/official_L6/wombo_combo/ZombieHybridNuMuSimulation_Sample2.root"
#if options.he_nugen:
#	taboutfile = options.infile.replace('.i3','.Sample'+str(options.sample)+'_HighEnergyOnly.root')
#taboutfile = infile.replace('.i3','.Sample'+str(options.sample)+'.root')

if taboutfile == options.infile:
	print "Table same name as input file!!!"
	exit()

tabler = rootwriter.I3ROOTTableService(taboutfile)
#tabler = hdfwriter.I3HDFTableService(taboutfile,1)

tray.AddModule(tableio.I3TableWriter,'writer1',tableservice = tabler, SubEventStreams = ['InIceSplit',],keys=KeysForTable)

tray.AddModule("TrashCan","Uncanny")

if options.num >0:
  tray.Execute(options.num)
else:
  tray.Execute()

tray.Finish()



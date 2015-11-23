#!/usr/bin/env python

#########################################################
# This script takes 2012 L2 files as input and will apply
# the L3 cuts with the option to filter.
#########################################################

from optparse import OptionParser
from os.path import expandvars
import time
import pickle

from I3Tray import *
from icecube import icetray, dataclasses, dataio,neutrinoflux
from icecube import tableio,hdfwriter,DeepCore_Filter,DomTools
from icecube.DeepCore_Filter import DOMS


from icecube.level3_filter_lowen import retro_level3_segments
import L4Segments


icetray.set_log_level(icetray.I3LogLevel.LOG_WARN)

#----------------------------------------------------------------------
# Allows the parser to take lists as inputs. Cleans out any white
# spaces that may be included in the names of the Official or
# Personal scripts to be run. Also removes the .py which is required
# to load a python module.

def list_callback(option, opt, value, parser):
    noWhiteSpace = value.replace(' ', '')
    noPyAppend   = noWhiteSpace.replace('.py', '')
    cleanList    = noPyAppend.split(',')
    setattr(parser.values, option.dest, cleanList)

def string_callback(option, opt, value, parser):
    lower_case   = value.lower()
    no_hyphens   = lower_case.replace('-', '')
    no_underscore= no_hyphens.replace('_', '')
    setattr(parser.values, option.dest, no_underscore)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# The inputfile can be added singly using the -i option
# or many input files can be processed by adding them
# as arguments.

usage = "%prog [options] <inputfiles>"

parser = OptionParser(usage=usage)

# i/o options
parser.add_option(
                  "-i", "--inputfile",
                  type      = "string",
                  action    = "store",
                  default   = '',
                  metavar   = "<input file>",
                  help      = "Name of the input file",
                  )
parser.add_option(
                  "-g", "--gcdfile",
                  type      = "string",
                  action    = "store",
                  default   = "None",
                  metavar   = "<geo file>",
                  help      = "Name of GCD file",
                  )
parser.add_option(
                  "-n", "--numevents",
                  type      = "int",
                  action    = "store",
                  default   = -1,
                  help      = "Number of events to process (default: all)",
                  )
parser.add_option(
                  "-o", "--outputfile",
                  type      = "string",
                  action    = "store",
                  default   = "TEST",
                  metavar   = "<output file(s) name>",
                  help      = "Name of the output file(s), i.e. .root and .i3.gz names",
                  )
parser.add_option(
                  "-b", "--branches",
                  type      = "int",
                  action    = "store",
                  default   = 1,
                  help      = "Process which filter branches (1=Both,2=StdDCFilter,3=ExpFidFilter)",
                  )
parser.add_option("-t","--tableoutput", action="store_true",
        default=False, dest="hdfout", help="Generate hdf5 file?")

parser.add_option("--filter", action="store_true",
	default=False, dest='filter', help="If true, only events passing specified branch will be written.")

(options, args) = parser.parse_args()

tray = I3Tray() 

icetray.load('libneutrinoflux')
icetray.load('libatmo-weights')
icetray.load('libDomTools')

### For BS processing ###
inputfile = options.inputfile
outputfile = '/data/user/daughjd/data/systematic_sets/10437/'+inputfile[76:-6]+'L4Out.i3'
#outputfile = options.outputfile

gcd = options.gcdfile

if inputfile == outputfile:
	print "Input same as Output Name!!! Exiting..."
	exit()

#gcdlist  = pickle.load(open('/nv/hp11/jdaughhetee3/runners/burn_sample_process/gcdlist.pkl','r'))
'''
if gcd == 'None':
	count=0
	for geo in gcdlist:
		if inputfile[:-21] == geo[:-14]:
			print "FOUND IT!"
			gcd = gcdlist.pop(count)
			print gcd
		count+=1	
	if gcd == 'None':
		print "Geometry File Not Found: Exiting!"
		print inputfile
		exit()
'''
files = []
files.append(gcd)
files.append(inputfile)

print files
tray.AddModule('I3Reader', 'reader',
               FileNameList = files)

#tray.AddModule("Rename",'saymyname',keys=['OfflinePulses','SplitInIcePulses'])

def AddSplitInIce(frame):
	frame["SplitInIcePulses"] = frame['MaskedOfflinePulses']
	return True

tray.AddModule(AddSplitInIce,'TransitionToTheFuture')

#### Plug other filter events into TwoLayer Veto Stream ####

def ExpFidFilterEventInput(frame):
    if frame.Has("FilterMask"):
	if (frame["FilterMask"].get("DeepCoreFilter_11").condition_passed or frame["FilterMask"].get("MuonFilter_11").condition_passed or frame["FilterMask"].get("CascadeFilter_11").condition_passed or frame["FilterMask"].get("MoonFilter_11").condition_passed or frame["FilterMask"].get("SDST_LowUp_11").condition_passed):
		return True
	return False

dlist = DOMS.DOMS("IC86EDC")
tldlist = DOMS.DOMS("IC86TwoLayVeto")

### Add DC Pulses if lacking ###

def NeedDCSRT(frame):
	if not "SRTTWOfflinePulsesDC" in frame:
		if  ExpFidFilterEventInput(frame):
			frame["SRTTWOfflinePulsesDC"] = frame["SRTOfflinePulses"]
	return True

tray.AddModule(NeedDCSRT,'ugh')

tray.AddModule("I3LCPulseCleaning","I3LCCleaning_Pulses_TLDCV",
                   If = lambda f: ExpFidFilterEventInput(f),
                   Input = "SplitInIcePulses",
                   OutputHLC = "InIcePulses_HLC",
                   OutputSLC = "InIcePulses_SLC"
                   )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>",'FiducialHLCRequirement',
                   If = lambda f: ExpFidFilterEventInput(f),
                   OmittedKeys=tldlist.DeepCoreFiducialDOMs,
                   SelectInverse = True,
                   InputResponse = "InIcePulses_HLC",
                   OutputResponse = 'HLC_DCFID',
                   OutputOMSelection = 'HLC_DCFIDSelection',
                   )

def FiducialTrig(frame):
	if "HLC_DCFID" in frame:
		if len(frame["HLC_DCFID"]) < 4:
			frame["DeepCoreFilter_TwoLayerExp_11"] = icetray.I3Bool(False)
	return True

tray.AddModule(FiducialTrig,'fiddy')

tray.AddModule("I3OMSelection<I3RecoPulseSeries>",'selectTwoLayICVetoDOMS',
                   OmittedKeys = tldlist.DeepCoreVetoDOMs,
                   OutputOMSelection = 'GoodSelection_TwoLayICVeto',
                   InputResponse = 'SRTTWOfflinePulsesDC',
                   OutputResponse = 'TwoLaySRTPulseICVeto',
                   SelectInverse = True,
                   If = lambda f: ExpFidFilterEventInput(f)
                   )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>",'selectTwoLayDCFidDOMs',
                   OmittedKeys= tldlist.DeepCoreFiducialDOMs,
                   SelectInverse = True,
                   InputResponse = 'SRTTWOfflinePulsesDC',
                   OutputResponse = 'TwoLaySRTPulseDCFid',
                   OutputOMSelection = 'GoodSelection_TwoLayDCFid',
                   If = lambda f: ExpFidFilterEventInput(f)
                   )

tray.AddModule("I3DeepCoreVeto<I3RecoPulse>",'twolaydeepcore_filter',
                   ChargeWeightCoG = False,
                   DecisionName = 'DeepCoreFilter_TwoLayerExp_11',
                   FirstHitOnly = True,
                   InputFiducialHitSeries = 'TwoLaySRTPulseDCFid',
                   InputVetoHitSeries = 'TwoLaySRTPulseICVeto',
                   If = lambda f: ExpFidFilterEventInput(f) and not f.Has("DeepCoreFilter_TwoLayerExp_11")
                   )

# ##############################
# 2012 DeepCore L3 tray segment
# ##############################

if (options.branches == 1) or (options.branches == 2):
  tray.AddSegment( retro_level3_segments.RetroDeepCoreCuts, "L3DeepCoreCuts",splituncleaned='SplitInIcePulses',year='11')

#########################################
# 2012 Expanded Fiducial DCFilter segment
#########################################

if (options.branches == 1) or (options.branches == 3):
  tray.AddSegment( retro_level3_segments.RetroExpFidDeepCoreCuts,'L3ExpDeepCoreCuts_2012',splituncleaned='SplitInIcePulses',year='11')


### Define Filter Module ###

def Filterizer(frame,branch=1):
	passbool=False
	stdbranch=False
	expbranch=False
	if frame.Has('IC2012_LE_L3'):
		stdbranch=frame['IC2012_LE_L3'].value
		#stdbranch=True
	if frame.Has('IC2012_ExpLE_L3'):
		expbranch=frame['IC2012_ExpLE_L3'].value
		#expbranch = True
	if branch == 1:
		passbool = expbranch or stdbranch
	if branch == 2:
		passbool = stdbranch
	if branch == 3:
		passbool = expbranch
	return passbool

if options.filter:
  tray.AddModule(Filterizer,'FilterToL3',branch = options.branches)

### Now Add L4 Processing! ###
def LES(frame):
        if "LES" in frame:
                return frame["LES"].value
        return False
def HES(frame):
        if "HES" in frame:
                return frame["HES"].value
        return False

tray.AddSegment(L4Segments.PreL4Processing,'PreL4')

### Add LES L4 Segment ###

tray.AddSegment(L4Segments.LESProcessing,'4Lorn')

### Add HES L4 Segment ###

tray.AddSegment(L4Segments.HESProcessing,'4Loko')

### Cleanup ###

tray.AddModule("Delete",'detritus',keys=['GoodSelection_EDCFid_NVetoSRT','GoodSelection_ICVeto_NVetoSRT','SRTTWOfflinePulsesExpDCFidCleanedKeys',
                                         'SRT300_InIcePulsesICVetoCleanedKeys','TwoLaySRTPulseDCFidCleanedKeys','TwoLaySRTPulseICVetoCleanedKeys',
					 'HLC_DCFIDCleanedKeys','HLC_DCFIDSelection','GoodSelection_TwoLayICVeto','GoodSelection_TwoLayDCFid',
					 'GoodSelection_TwoLayDCFid','TwoLaySRTPulseDCFid','TwoLaySRTPulseICVeto'])


# #########################################
# List of objects to write to table
# #########################################

FrameObjectsToKeep=['PrimaryNu',
                      'InteractionParticle',
                      'VolumeTriggerFlag',
                      'AtmoWeight',
                      'AtmoWeightNoOsc',
                      'CorsikaWeightMap',
                      'InIceSMTTriggered',
                      'DeepCoreSMTTriggered',
                      'InDC',
                      'ChkGRBEventWeight',
                      'ChkGRBSpectrumWeight',
                      'NchHLC',
                      'InExpDC',
                      'LineFit_DC',
                      'SPEFit6_DC',
                      'SPEFit6_DCFitParams',
                      'HES',
                      'LES',
                      'IC2012_LE_L3',
                      'IC2012_LE_L3_Vars',
                      'IC2012_ExpLE_L3',
                      'IC2012_ExpLE_L3_Vars',
                      "L4Bool_LES",
                      "L4Bool_HES",
                      'L4_Variables_LES',
                      'L4_Variables_HES',
                      'SPEFit6CramerRao_DCParams',
                      'ToI_DC',
                      'TwoLaySRTPulseDCFidCleanedKeys',
                      'TwoLaySRTPulseICVetoCleanedKeys',
                      'HLC_DCFIDCleanedKeys',
                      'HLC_DCFIDSelection',
                      'GoodSelection_TwoLayICVeto',
                      'GoodSelection_TwoLayDCFid',
                      'NString']

if options.hdfout:
	taboutfile = outputfile.replace('.i3','.hdf5')
	tabler = hdfwriter.I3HDFTableService(taboutfile,1)
	tray.AddModule(tableio.I3TableWriter,'writer1',tableservice = tabler, SubEventStreams = ['in_ice',],keys=FrameObjectsToKeep)


def L4Sieve(frame,stream='Both'):
        lesbool = False
        hesbool = False
        if "L4Bool_LES" in frame:
                lesbool = frame["L4Bool_LES"].value
        if "L4Bool_HES" in frame:
                hesbool = frame["L4Bool_HES"].value
        bothbool = lesbool or hesbool
        if stream == 'Both':
                return bothbool
        if stream == 'LES':
                return lesbool
        if stream == "HES":
                return hesbool

if options.filter:
  tray.AddModule(L4Sieve,'L4Filtering_LES',stream="Both")


### Add Atmospheric Weighting ###
def ParticleFetcher(frame):
        if "I3MCTree" in frame:
                plist=frame["I3MCTree"].in_ice
                primary = plist[0]
                interaction_particle = plist[1]
                frame["InteractionParticle"] = interaction_particle
                frame["PrimaryNu"] = primary
        return True

tray.AddModule(ParticleFetcher,'fetching')

load("libneutrinoflux")

tray.AddService("NeutrinoFluxFactory", "numu-flux-service-honda",
                ConventionalModelName = "honda2006_numu",
                InstallServiceAs = "honda2006_numu",
                )
tray.AddService("NeutrinoFluxFactory", "nue-flux-service-honda",
                ConventionalModelName = "honda2006_nue",
                InstallServiceAs = "honda2006_nue",
                PromptModelName = "naumov_rqpm_nue",
                )

load("libatmo-weights")

tray.AddModule("AtmoWeights", "hondaWeights",
                NuEFluxService = 'honda2006_nue',
                NuMuFluxService = 'honda2006_numu',
                DumpToFrame=False,
                If = lambda frame: frame["I3MCTree"][0].energy>10
                )

tray.AddService("NeutrinoFluxFactory", "numu-flux-service-bartol",
                InstallServiceAs="bartol_numu",
                )
tray.AddService("NeutrinoFluxFactory", "nue-flux-service-bartol",
                ConventionalModelName = "bartol_nue",
                InstallServiceAs = "bartol_nue",
                PromptModelName = "naumov_rqpm_nue",
                )
tray.AddModule("AtmoWeights", "bartolWeights",
                OutputName = "AtmoWeights_bartol",
                If = lambda frame: frame["I3MCTree"][0].energy>10,
                DumpToFrame = False
        )

##################################

tray.AddModule( 'I3Writer', 'EventWriter',
                Filename          = outputfile,
                Streams           = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ],
                DropOrphanStreams = [icetray.I3Frame.DAQ]
                )

tray.AddModule( 'TrashCan' , 'Done' )

if options.numevents < 0:
    tray.Execute()
else:
    tray.Execute(options.numevents)
tray.Finish ()


del tray

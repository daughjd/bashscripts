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
from icecube import icetray, dataclasses, dataio
from icecube import tableio,hdfwriter

from icecube.DeepCore_L3_2012 import level3_segments
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

### For BS processing ###
inputfile = options.inputfile
outputfile = '/net/user/daughjd/data/CorskL4/'+inputfile[71:-6]+'L4Out.i3'

gcd = options.gcdfile

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

tray.AddModule("Delete","Deli",keys=['TWOfflinePulsesDC'])
tray.AddModule("Rename",'saymyname',keys=['OfflinePulses','SplitInIcePulses'])
# ##############################
# 2012 DeepCore L3 tray segment
# ##############################

if (options.branches == 1) or (options.branches == 2):
  tray.AddSegment( level3_segments.DeepCoreCuts, "L3DeepCoreCuts",splituncleaned='SplitInIcePulses',year='11')

#########################################
# 2012 Expanded Fiducial DCFilter segment
#########################################

if (options.branches == 1) or (options.branches == 3):
  tray.AddSegment( level3_segments.ExpFidDeepCoreCuts,'L3ExpDeepCoreCuts_2012',splituncleaned='SplitInIcePulses',year='11')


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
                                         'SRT300_InIcePulsesICVetoCleanedKeys'])


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

#tray.AddModule(L4Sieve,'L4Filtering_LES',stream="Both")


#tray.AddModule( 'I3Writer', 'EventWriter',
#                Filename          = outputfile,
#                Streams           = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ],
#                DropOrphanStreams = [icetray.I3Frame.DAQ]
#                )

tray.AddModule( 'TrashCan' , 'Done' )

if options.numevents < 0:
    tray.Execute()
else:
    tray.Execute(options.numevents)
tray.Finish ()


del tray

#!/usr/bin/env python

#########################################################
# This script takes 2012 L2 files as input and will apply
# the L3 cuts with the option to filter.
#########################################################

from optparse import OptionParser
from os.path import expandvars
import time

from I3Tray import *
from icecube import icetray, dataclasses, dataio
from icecube import tensor_of_inertia, mcsummary, improvedLinefit,SeededRTCleaning
from icecube import tableio,hdfwriter

from icecube.DeepCore_L3_2012 import level3_segments
from icecube.std_processing import level2_globals,level2_HitCleaning,level2_DeepCoreReco

icetray.load("DomTools",False)
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
                  action    = "callback",
                  callback  = list_callback,
                  default   = [],
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

parser.add_option("--l2", action="store_true", default=False, dest="leveltwo", help="Perform DC L2 processing?")

parser.add_option("--filter", action="store_true",
	default=False, dest='filter', help="If true, only events passing specified branch will be written.")

(options, args) = parser.parse_args()

tray = I3Tray() 

files = []
files.append(options.gcdfile)
files.extend(options.inputfile)

print files

tray.AddModule('I3Reader', 'reader',
               FileNameList = files)


def PassedDCFilter(frame):
        if frame.Has("FilterMask"):
            if frame["FilterMask"].get("DeepCoreFilter_12").condition_passed:
                return True
            # end if()
        # end if()
        return False
    # end PassedDCFilter()

############ DO L2 ##############
if options.leveltwo:
  tray.AddModule("Rename","NamingConventionCorrector",keys=['SplitInIcePulses',level2_globals.masked_offline_pulses,'SplitUncleanedInIcePulses',level2_globals.masked_offline_pulses])

  tray.AddSegment(level2_HitCleaning.DeepCoreHitCleaning,"L2DCHitCleaning")
  
  tray.AddSegment(level2_DeepCoreReco.OfflineDeepCoreReco, 'RecoTime',
                  If = lambda f: PassedDCFilter(f),
                  Pulses = level2_globals.srt_tw_offline_pulses_dc,
                  suffix = level2_globals.deepcore_name
                 )

# ##############################
# 2012 DeepCore L3 tray segment
# ##############################

if (options.branches == 1) or (options.branches == 2):
  tray.AddSegment( level3_segments.DeepCoreCuts, "L3DeepCoreCuts",splituncleaned=level2_globals.masked_offline_pulses)

#########################################
# 2012 Expanded Fiducial DCFilter segment
#########################################

if (options.branches == 1) or (options.branches == 3):
  tray.AddSegment( level3_segments.ExpFidDeepCoreCuts,'L3ExpDeepCoreCuts_2012',splituncleaned=level2_globals.masked_offline_pulses)

# #########################################
# List of objects to write to table
# #########################################

FrameObjectsToKeep = [ "CascadeLast_DC",
                       "CorsikaWeightMap",
                       "DecayParticle",
		       "iLinefit_LE_L3",                      
                       "I3EventHeader",
                       "I3MCTree",
                       "AtmoWeight",
                       "AtmoWeightNoOsc",
                       "AtmoWeight_bartol",
                       "AtmoWeight_bartolNoOsc",
                       'VertexGuessInsideExpDC',
                       "I3MCWeightDict",
                       "I3TriggerHierarchy",
		       "IC2012_LE_L3_Vars",
		       "IC2012_LE_L3",
		       "IC2012_ExpLE_L3_Vars",
		       "IC2012_ExpLE_L3",
                       "IniceNu",
                       "InteractionVertex",
                       "InteractionParticle",
                       "FilterMask",
                       "LineFit_DC",
                       "NoiseEngine_bool",
                       "PrimaryNu",
                       "ToI_DC",
                       ]

### Define Filter Module ###

def Filterizer(frame,branch=1):
	passbool=False
	stdbranch=False
	expbranch=False
        if PassedDCFilter(frame):
                stdbranch=True
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

tray.AddModule( 'I3Writer', 'EventWriter',
                Filename          = options.outputfile,
                Streams           = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ],
                DropOrphanStreams = [icetray.I3Frame.DAQ]
                )



if options.hdfout:
	taboutfile = options.outputfile.replace('.i3','.hdf5')
	tabler = hdfwriter.I3HDFTableService(taboutfile,1)
	tray.AddModule(tableio.I3TableWriter,'writer1',tableservice = tabler, SubEventStreams = ['InIceSplit',],keys=FrameObjectsToKeep)

tray.AddModule( 'TrashCan' , 'Done' )

if options.numevents < 0:
    tray.Execute()
else:
    tray.Execute(options.numevents)
tray.Finish ()


del tray

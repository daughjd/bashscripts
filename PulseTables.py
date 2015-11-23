#!/usr/bin/env python

from I3Tray import *
from icecube import dataclasses, dataio, icetray, tableio
import os
import sys
from optparse import OptionParser

from icecube import rootwriter
from icecube import hdfwriter


parser = OptionParser()

parser.add_option("-o", "--output", action="store",
                        type="string", default="automate", dest="outfile",
                                help="Output i3 file")

parser.add_option("--filter", action="store_true",
                        default=False, dest="filter", help="Filter events that fail L4 cuts?")

parser.add_option("-n", "--num", action="store",
                        type="int", default=-1, dest="num",
                                help="Number of frames to process"
                                        )

parser.add_option("-t","--tableoutput", action="store_true",
                        default=False, dest="hdfout", help="Generate hdf5 file?")

(options,args) = parser.parse_args()


runfiles = sys.argv[1:-2]
outfile = options.outfile

print runfiles



tray = I3Tray()

tray.AddModule("I3Reader","i3reader")(
    ("FilenameList", runfiles),
    ("SkipKeys",[])
    )

## Add N Channel Cut to obtain events in region of disagreement ##
def NChannelCut(frame,nchanmax):
        boolio = False
        if 'LET_RecoPulses' in frame:
                if len(frame['LET_RecoPulses'].apply(frame)) <= nchanmax:
                        boolio = True
        if not 'LET_RecoPulses' in frame:
                print "Frame Object 'LET_RecoPulses' missing..."
        return boolio

tray.AddModule(NChannelCut,'ResolutionSlice',nchanmax=20)


### Still some non-unique GENIE events ###
import pickle
unique_genies = list(pickle.load(open("UniqueIDGENIEEvents.pkl",'r')))

def FishingExpedition(frame):
	keepbool = False
	event_id = frame["I3EventHeader"].event_id
	if event_id in unique_genies:
		keepbool = True
		unique_genies.remove(event_id)
	return keepbool
		
if runfiles[0].__contains__("genie"):
  tray.AddModule(FishingExpedition,'Catch')

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

taboutfile = outfile

#tabler = rootwriter.I3ROOTTableService(taboutfile)
tabler = hdfwriter.I3HDFTableService(taboutfile,1)


tray.AddModule(tableio.I3TableWriter,'writer1',tableservice = tabler, SubEventStreams = ['InIceSplit',],keys=KeysForTable)

tray.AddModule("TrashCan", "the can");


tray.Execute()
tray.Finish()
	

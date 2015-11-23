#!/usr/bin/env python

from I3Tray import *

from os.path import expandvars
from icecube import dataclasses, dataio, icetray
import os
import sys

from optparse import OptionParser

load("libicepick")
load("libicetray")
load("libdataclasses")
load("libdataio")

usage = "%prog [options] <inputfiles>"
parser = OptionParser(usage=usage)

parser.add_option(
                  "-i", "--inputfile",
                  type      = "string",
                  metavar   = "<input file>",
                  help      = "Name of the input file",
                  )
parser.add_option(
                  "-o", "--outputfile",
                  type      = "string",
                  metavar   = "<input file>",
                  help      = "Name of the input file",
                  )

(options, args) = parser.parse_args()

runfiles = options.inputfile

writefile = options.outputfile

print writefile

def TransientAnalysisEventGrab(frame):
        preciousfood = False
        if frame.Has('FilterMask'):
                if not frame['FilterMask'].get('DeepCoreFilter_12').condition_passed and frame['FilterMask'].get('DeepCoreFilter_TwoLayerExp_12').condition_passed:
                        preciousfood = True
        frame.Put('DCFilterStreamsEvent',icetray.I3Bool(preciousfood))
	return preciousfood

tray = I3Tray()

tray.AddModule("I3Reader","i3reader")(
    ("FilenameList", [runfiles]),
    ("SkipKeys",[])
    )

tray.AddModule(TransientAnalysisEventGrab,'grabule',Streams=[icetray.I3Frame.Physics])

tray.AddModule( "I3Writer", "EventWriter", filename=writefile,
                Streams=[icetray.I3Frame.Physics,icetray.I3Frame.DAQ],
                DropOrphanStreams=[icetray.I3Frame.DAQ]
                )


tray.AddModule("TrashCan", "the can");


tray.Execute(100)
tray.Finish()
	

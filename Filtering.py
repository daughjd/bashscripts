#!/usr/bin/env python

from I3Tray import *
from icecube import dataclasses, dataio, icetray
import os
import sys
from optparse import OptionParser

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

tray.AddModule("I3Writer", "i3writer",
                  CompressionLevel = 0,
                  DropOrphanStreams = [icetray.I3Frame.DAQ],
                  Streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                  Filename = outfile,
                  )




tray.AddModule("TrashCan", "the can");


tray.Execute()
tray.Finish()
	

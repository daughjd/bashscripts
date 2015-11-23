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

def CorrectedParaboloidCut(frame):
        boolio = False
        if 'ReScaled_Paraboloid_Sigma_SplineMPEMod' in frame:
                if frame['ReScaled_Paraboloid_Sigma_SplineMPEMod'].value < 0.6 and frame['ReScaled_Paraboloid_Sigma_SplineMPEMod'].value > 0.3:
                        boolio = True
        if not 'ReScaled_Paraboloid_Sigma_SplineMPEMod' in frame:
                print "Uh oh..."
        return boolio

tray.AddModule(CorrectedParaboloidCut,'ResolutionSlice')

tray.AddModule("I3Writer", "i3writer",
                  CompressionLevel = 0,
                  DropOrphanStreams = [],
                  Streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                  Filename = outfile,
                  )




tray.AddModule("TrashCan", "the can");


tray.Execute()
tray.Finish()
	

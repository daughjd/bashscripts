#!/usr/bin/env python

from I3Tray import *
from icecube import dataclasses, dataio, icetray
import os
import sys

runfiles = sys.argv[1:]
outfile = '/net/user/daughjd/data/DCStream_Level2_IC86.2012_data_Run00120811_AllSubs.i3'

tray = I3Tray()

tray.AddModule("I3Reader","i3reader")(
    ("FilenameList", runfiles),
    ("SkipKeys",['MCPMTResponseMap','CleanInIceRawData'])
    )


tray.AddModule("I3Writer", "i3writer",
                  CompressionLevel = 0,
                  DropOrphanStreams = [],
                  Streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                  Filename = outfile,
                  )




tray.AddModule("TrashCan", "the can");


tray.Execute()
tray.Finish()
	

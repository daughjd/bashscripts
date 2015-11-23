#!/usr/bin/env python
from I3Tray import *

from os.path import expandvars
from icecube import dataclasses, dataio, icetray, filter_tools
import os
import sys

load("libicepick")
load("libDomTools")

runfile = sys.argv[1]

writefile = runfile.replace('.i3','.LowEnScraps.i3')
writefile = writefile.replace('naoko/PnF24hRunResplit','daughjd/data')

print writefile

tray = I3Tray()

tray.AddModule("I3Reader","i3reader")(
    ("Filename", runfile),
    ("SkipKeys",[])
    )
    

def DCGrab(frame):
  if frame.Has('FilterMask'):
    if (frame['FilterMask'].get('DeepCoreFilter_12').condition_passed or frame['FilterMask'].get('DeepCoreFilter_TwoLayerExp_12').condition_passed):
      frame.Put('DCEvent',icetray.I3Bool(True))
    else:
      frame.Put('DCEvent',icetray.I3Bool(False))  
    frame.Put("DCV",icetray.I3Bool(frame['FilterMask'].get('DeepCoreFilter_12').condition_passed))
    frame.Put("TLDCV",icetray.I3Bool(frame['FilterMask'].get('DeepCoreFilter_TwoLayerExp_12').condition_passed))
tray.AddModule(DCGrab,'grabble')

tray.AddModule("I3IcePickModule<I3SimpleFilter>",'DCEvent_Filter',
    DecisionName = "Grub", 
    DiscardEvents = True,
    InputDecisionName = "DCEvent"                
    )
tray.AddModule("I3IcePickModule<I3SimpleFilter>",'DCEvent_DC',
    DecisionName = "Grubb",
    DiscardEvents = False,
    InputDecisionName = "DCV"
    )
tray.AddModule("I3IcePickModule<I3SimpleFilter>",'DCEvent_TL',
    DecisionName = "Grubbb",
    DiscardEvents = False,
    InputDecisionName = "TLDCV"
    )
tray.AddModule("I3Writer","writer")(
    ("filename", writefile),
    ("DropOrphanStreams", [icetray.I3Frame.DAQ]),
    ("SkipKeys", ['Grub','Grubb','Grubbb']),
    ("Streams", [])
    )

tray.AddModule("TrashCan", "the can");


tray.Execute()
tray.Finish()
	

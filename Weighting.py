#!/usr/bin/env python

from icecube import icetray,dataclasses,dataio,icepick,tableio,neutrinoflux
from I3Tray import *
from icecube import hdfwriter
from optparse import OptionParser
from ChkGRBWeighting import ChokedGRBFlux, ChokedGRBWeighting


parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-i", "--input", action="store", type="string", default="",
                  dest="INPUT", help="Input  3 file to process")
parser.add_option("-o", "--output", action="store", type="string", default="",
                  dest="OUTPUT", help="Output i3 file")
parser.add_option("-g", "--gcd", action="store", type="string", default="",
                  dest="GCD", help="GCD file for input i3 file")
parser.add_option("-p", "--prettyprint", action="store_true",
                  dest="PRETTY", help="Do nothing other than big tray dump")
parser.add_option("--chkgrbweight", action="store_true", default=False,
                  dest="CHKGRB", help="Use Chk GRB event weighting for Nu Simulation")
parser.add_option("-f","--flavor",action="store",type="string",default='numu',dest="FLAVOR",help="Neutrino Flavor (numu,nue,nutau)")

parser.add_option("--numsimfiles",type="int",action="store",default=1,dest="SIMFILES",help="Number of Simfiles in input file")

parser.add_option("-n", "--numevents",
                  type      = "int",
                  action    = "store",
                  default   = -1,
                  help      = "Number of events to process (default: all)",
                  )


# get parsed args
(options,args) = parser.parse_args()



tray = I3Tray()

icetray.load('libDomTools')
icetray.load('libmcsummary')
icetray.load('libneutrinoflux')
icetray.load('libatmo-weights')
icetray.load('libtrigger-splitter')
icetray.load('libjeb-filter-2012')

if len(options.GCD) > 1:
  runfiles = [options.GCD,options.INPUT]
else:
  runfiles = [options.INPUT]

writefile = options.OUTPUT
#outfile = runfiles[1].replace('.i3','.AtmoWeightCorrect.hdf5')
#tabler = hdfwriter.I3HDFTableService(outfile,1)


tray.AddModule("I3Reader","i3reader")(
    ("FilenameList", runfiles),
    ("SkipKeys",[])
    )

def PassedEitherDCFilter(frame):
        if frame.Has("FilterMask"):
            if frame["FilterMask"].get("DeepCoreFilter_12").condition_passed or frame["FilterMask"].get("DeepCoreFilter_TwoLayerExp_12").condition_passed:
                return True
            # end if()
        # end if()
        return False
    # end PassedDCFilter()

tray.AddModule(PassedEitherDCFilter,'DropOtherSplits')

tray.AddModule("I3MCTagger",'mctagger',
               MCTreeName = 'I3MCTree'
              )

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

if options.CHKGRB:
  tray.AddModule(ChokedGRBWeighting,'ChkGRBWeight',
       gamma = 3.0,
       nfiles = options.SIMFILES,
       weightdict = 'I3MCWeightDict',
       specweightname = 'ChkGRBSpectrumWeight',
       eventweightname = 'ChkGRBEventWeight',
       flavor = options.FLAVOR
       )


tray.AddModule("I3Writer", "i3writer",
                  CompressionLevel = 0,
                  DropOrphanStreams = [],
                  Streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                  Filename = writefile,
                  SkipKeys = [],
                  )
tray.AddModule('TrashCan','TRASHY')

if options.numevents < 0:
    tray.Execute()
else:
    tray.Execute(options.numevents)
tray.Finish ()



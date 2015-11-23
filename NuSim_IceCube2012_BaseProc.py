#!/usr/bin/env python

from I3Tray import *
from icecube import icetray, dataclasses, dataio, jeb_filter_2012, filter_tools, trigger_sim
from icecube import phys_services
from icecube.jeb_filter_2012 import filter_globals
from icecube.jeb_filter_2012.filter_globals import which_split
import os, sys, time
from math import log10, cos, radians
from optparse import OptionParser

## Import DOMS -> DeepCore string/DOM definition class
#from icecube.DeepCore_Filter import DOMS

# For the Stopwatch
#load("libpfdebugging")


start_time = time.asctime()

print 'Started:', start_time
 
# handling of command line arguments  
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)
parser.add_option("-i", "--input", action="store", type="string", default="", 
		  dest="INPUT", help="Input  3 file to process")
parser.add_option("-o", "--output", action="store", type="string", default="",
		  dest="OUTPUT", help="Output i3 file")
parser.add_option("-g", "--gcd", action="store", type="string", default="",
		  dest="GCD", help="GCD file for input i3 file")
parser.add_option("-d", "--decode", action="store_true", 
		  dest="DECODE", help="Add the decoder to the processing")
parser.add_option("--simdata", action="store_true", 
		  dest="SIMDATA", help="This is IC86 sim data, decode and triger")
parser.add_option("--qify", action="store_true", default=False,
		  dest="QIFY", help="Apply QConverter, use if file is P frame only")
parser.add_option("--no-retrigger", action="store_false", default=True,
		  dest="RETRIGGER", help="Apply trigger simulation (only applies to simulated data)")
parser.add_option("-p", "--prettyprint", action="store_true", 
		  dest="PRETTY", help="Do nothing other than big tray dump")
parser.add_option("--chkgrbweight", action="store_true", default=False, 
		  dest="CHKGRB", help="Use Chk GRB event weighting for Nu Simulation")

# get parsed args
(options,args) = parser.parse_args()

GCD = options.GCD
inputfile = options.INPUT
outputfile = options.OUTPUT
dodecode = options.DECODE
simdata = options.SIMDATA
prettyprint = options.PRETTY
grbevent = options.CHKGRB

outputfile = outputfile.replace('/data/sim/IceCube/2011/generated/GENIE-in-ice/numu/32/','/net/user/daughjd/data/GENIE_Filtered/')

print 'Opening file %s' % inputfile
 
print 'Preparing to write i3 file  %s' % outputfile 
 

#dbhost = "dbs2.icecube.wisc.edu"
#dbuser = "www"


########################################

class Skipper(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("num", "Num to skip at start of i3file", 0)
        self.AddOutBox("OutBox")
    

    def Configure(self):
        self.num = self.GetParameter("num")
        print 'Will skip the first events:', self.num
        self.count = 0

    def DAQ(self, frame):
        self.count += 1
        if(self.count < self.num):
#            print "skip", self.count, self.num
            pass
        else:
            self.PushFrame(frame)


tray = I3Tray() 
# Magic stopwatch                                                               
# sw_rootfile = 'StopWatch_' + outputfile + '.root'
# print 'Stopwatch file:', sw_rootfile
# tray.AddService("I3StopwatchServiceFactory","I3StopwatchService")(
#         ("InstallServiceAs","I3StopwatchService"),
#         ("OutputFileName", sw_rootfile),
#         ("UseHighTimeResolution",False),
#         )


tray.AddModule("I3Reader", "reader", filenamelist=[GCD,inputfile],
               SkipKeys = ['I3DST11',
                           'I3SuperDST',
                           'I3VEMCalData',
                           'PoleMuonLlhFit',
                           'PoleMuonLlhFitCutsFirstPulseCuts',
                           'PoleMuonLlhFitFitParams',
                           'CramerRaoPoleL2IpdfGConvolute_2itParams',
                           'CramerRaoPoleL2MPEFitParams',
                           'PoleL2IpdfGConvolute_2it',
                           'PoleL2IpdfGConvolute_2itFitParams',
                           'PoleL2MPEFit',
                           'PoleL2MPEFitCuts',
                           'PoleL2MPEFitFitParams',
                           'PoleL2MPEFitMuE',
                           ])
# move that old filterMask out of the way
tray.AddModule("Rename","filtermaskmover",Keys=["FilterMask","OrigFilterMask"])

# kick the stopwatch                                                            

#tray.AddModule("I3StopwatchStart","stopwatchstart")(
#        ("ID","I3ClientTime"),
#        ("StopwatchServiceInstalledAs","I3StopwatchService"),
#        )



if options.QIFY:
    tray.AddModule("QConverter", "qify", WritePFrame=False)
#tray.AddModule(Skipper, "skipper", num =2)

if simdata:
    #My IC86 simulation is P-framed, and has no triggers                                          

    if options.RETRIGGER:
        # some cleanup first
        tray.AddModule("Delete", "delete_triggerHierarchy", Keys=["I3TriggerHierarchy", "TimeShift"])
        time_shift_args = { "I3MCTreeNames": ['I3MCTree'],
                            "I3MCPMTResponseMapNames" : [],
                            "I3MCHitSeriesMapNames" : [] }
        gcd_file = dataio.I3File(GCD)
        tray.AddSegment(trigger_sim.TriggerSim, "trig", gcd_file = gcd_file ,
                        time_shift_args = time_shift_args)
    #Cleanout full waveforms from SLC-only readouts
    tray.AddModule("Rename","RenameslcCleanup",
                   keys=['InIceRawData','DirtyInIceRawData']
                   )
    tray.AddModule("I3SLCcleanup","CleanSLCForSim",
                   InIceInput = "DirtyInIceRawData",
                   InIceOutput ="InIceRawData"
                   )
    

## base processing requires:  GCD and frames being fed in by reader or Inlet
## base processing include:  decoding, TriggerCheck, Bad Dom cleaning, calibrators, Feature extraction,
##     pulse cleaning (seeded RT, and Time Window), PoleMuonLineit, PoleMuonLlh,
##     Cuts module and Mue on PoleMuonLlh
tray.AddSegment(jeb_filter_2012.BaseProcessing, "BaseProc",
		pulses=filter_globals.CleanedMuonPulses,
		decode=dodecode,simulation=simdata,
		)

#tray.AddSegment(jeb_filter_2012.IceTopVEMCal, "VEMCALStuff")
## Filters that use the InIce Trigger splitting

tray.AddSegment(jeb_filter_2012.MuonFilter, "MuonFilter",
                pulses = filter_globals.CleanedMuonPulses,
                If = which_split(split_name=filter_globals.InIceSplitter))

tray.AddSegment(jeb_filter_2012.CascadeFilter, "CascadeFilter",
                pulses = filter_globals.CleanedMuonPulses,
                muon_llhfit_name = filter_globals.muon_llhfit,
                If = which_split(split_name=filter_globals.InIceSplitter))	

tray.AddSegment(jeb_filter_2012.FSSFilter, "FSSFilter",
                pulses = filter_globals.SplitUncleanedInIcePulses,
                If = which_split(split_name=filter_globals.InIceSplitter)
                )

#tray.AddSegment(jeb_filter_2012.LowUpFilter, "LowUpFilter",
#                If = which_split(split_name=filter_globals.InIceSplitter)
#                )

#tray.AddSegment(jeb_filter_2012.ShadowFilter, "ShawdowFilters",
#                If = which_split(split_name=filter_globals.InIceSplitter)
#                )

tray.AddSegment(jeb_filter_2012.GCFilter, "GCFilter",
                If = which_split(split_name=filter_globals.InIceSplitter)
                )

tray.AddSegment(jeb_filter_2012.VEFFilter, "VEFFilter", 
                pulses = filter_globals.CleanedMuonPulses,
                If = which_split(split_name=filter_globals.InIceSplitter)
                )

#tray.AddSegment(jeb_filter_2012.OnlineL2Filter, "OnlineL2",
#                pulses = filter_globals.CleanedMuonPulses,
#                llhfit_name = filter_globals.muon_llhfit,
#                improved_linefit = True,
#                paraboloid = False,
#                If = which_split(split_name=filter_globals.InIceSplitter)
#                )

tray.AddSegment(jeb_filter_2012.DeepCoreFilter,"DeepCoreFilter",
                pulses = filter_globals.SplitUncleanedInIcePulses,
                If = which_split(split_name=filter_globals.InIceSplitter)
                )



## Filters on the Null split

#tray.AddSegment(jeb_filter_2012.EHEFilter, "EHEFilter",
#                If = which_split(split_name=filter_globals.NullSplitter)
#                )

#tray.AddSegment(jeb_filter_2012.MinBiasFilters, "MinBias",
#                If = which_split(split_name=filter_globals.NullSplitter)
#                )

#tray.AddSegment(jeb_filter_2012.SlopFilters, "SLOP",
#                If = which_split(split_name=filter_globals.NullSplitter)
#                )

#tray.AddSegment(jeb_filter_2012.FixedRateTrigFilter, "FixedRate",
#                If = which_split(split_name=filter_globals.NullSplitter)
#                )


## Filters on the CR split
#tray.AddSegment(jeb_filter_2012.CosmicRayFilter, "CosmicRayFilter", 
#                pulseMask = filter_globals.SplitUncleanedITPulses,
#                If = which_split(split_name=filter_globals.IceTopSplitter)
#                )


## not sure which split this gets.  (ALL?  Recos used are only on the InIce..)
tray.AddSegment(jeb_filter_2012.DST, "DSTFilter", 
                dstname  = "I3DST12",
                pulses   = filter_globals.CleanedMuonPulses,
                If = which_split(split_name=filter_globals.InIceSplitter))

# Generate filter Masks for all P frames
filter_mask_randoms = phys_services.I3GSLRandomService(9999)
print filter_globals.filter_pairs + filter_globals.sdst_pairs

tray.AddModule(filter_tools.FilterMaskMaker, "MakeFilterMasks",
               OutputMaskName = filter_globals.filter_mask,
               FilterConfigs = filter_globals.filter_pairs+ filter_globals.sdst_pairs,
               RandomService = filter_mask_randoms)

# Merge the FilterMasks:
tray.AddModule("OrPframeFilterMasks", "make_q_filtermask",
               InputName = filter_globals.filter_mask,
               OutputName = filter_globals.qfilter_mask)


#Q+P frame specific keep module needs to go first, as KeepFromSubstram
#will rename things, let's rename post keep.  

prekeeps = filter_globals.q_frame_keeps + [filter_globals.rawdaqdata,'JEBEventInfo'] + filter_globals.null_split_keeps + filter_globals.inice_split_keeps + filter_globals.icetop_split_keeps + filter_globals.onlinel2filter_keeps

#tray.AddModule("Keep","keep_before_merge",
#               keys = prekeeps)

## second set of prekeeps, conditional on filter content, based on newly created Qfiltermask
#Determine if we should apply harsh keep for events that failed to pass any filter
##  Note: excluding the sdst_streams entries

tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckAll",
               FilterNameList = filter_globals.filter_streams,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedAnyFilter",
               DiscardEvents = False, 
               Streams = [icetray.I3Frame.DAQ]
               )
def do_save_just_superdst(frame):
    if frame.Has("PassedAnyFilter"):
        if not frame["PassedAnyFilter"].value:
            return True    #  <- Event failed to pass any filter.  
        else:
            return False # <- Event passed some filter

    else:
        print "Failed to find key frame Bool!!"
        return False

#tray.AddModule("Keep", "KeepOnlySuperDSTs",
#               keys = filter_globals.keep_nofilterpass+["PassedAnyFilter"],
#               If = do_save_just_superdst
#               )

## Now clean up the events that not even the SuperDST filters passed on.
tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckSDST",
               FilterNameList = filter_globals.sdst_streams,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedKeepSuperDSTOnly",
               DiscardEvents = False,
               Streams = [icetray.I3Frame.DAQ]
               )

def dont_save_superdst(frame):
    if frame.Has("PassedKeepSuperDSTOnly") and frame.Has("PassedAnyFilter"):
        if frame["PassedAnyFilter"].value:
            return False  #  <- these passed a regular filter, keeper
        elif not frame["PassedKeepSuperDSTOnly"].value:
            return True    #  <- Event failed to pass SDST filter.  
        else:
            return False # <- Event passed some  SDST filter

    else:
        print "Failed to find key frame Bool!!"
        return False

#tray.AddModule("Keep", "KeepOnlyDSTs",
#               keys = filter_globals.keep_dst_only + ["PassedAnyFilter","PassedKeepSuperDSTOnly"],
#               If = dont_save_superdst
#               )


## Frames should now contain only what is needed.  now flatten, write/send to server
# Squish P frames back to single Q frame, one for each split:
#tray.AddModule("KeepFromSubstream","null_stream",
#               StreamName = filter_globals.NullSplitter,
#               KeepKeys = filter_globals.null_split_keeps,
#               )

#tray.AddModule("KeepFromSubstream","inice_split_stream",
#               StreamName = filter_globals.InIceSplitter,
#               KeepKeys = filter_globals.inice_split_keeps + filter_globals.onlinel2filter_keeps,
#               )

#tray.AddModule("KeepFromSubstream","icetop_split_stream",
#               StreamName = filter_globals.IceTopSplitter,
#               KeepKeys = filter_globals.icetop_split_keeps,
#               )

# Apply small keep list (SuperDST/SmallTrig/DST/FilterMask for non-filter passers


print 'Want I3DAQData: ',filter_globals.filters_keeping_allraw

# Remove I3DAQData object for events not passing one of the 'filters_keeping_allraw'
tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheck",
               FilterNameList = filter_globals.filters_keeping_allraw,
               FilterResultName = filter_globals.qfilter_mask,
               DecisionName = "PassedConventional",
               DiscardEvents = False,
               Streams = [icetray.I3Frame.DAQ]
               )

## Clean out the SuperDST errata when you have full I3DAQData
def I3DAQDataCleaner(frame):
    if frame.Has("PassedConventional"):
        if frame['PassedConventional'].value == False:
            frame.Delete(filter_globals.rawdaqdata)

## TODO: This needs to be done in the server.
tray.AddModule(I3DAQDataCleaner,"CleanErrataForConventional",
               Streams=[icetray.I3Frame.DAQ])

# Stop the stopwatch
#tray.AddModule("I3StopwatchStop","stopwatchstop")(
#        ("ID","I3ClientTime"),
#        ("StopwatchServiceInstalledAs","I3StopwatchService"),
#        ("ProcTimeFrameName","ExecutionTime")
#        )
### Add User Specified Event Picks ###
def TransientAnalysisEventGrab(frame):
	preciousfood = False
	if frame.Has('FilterMask'):
		if frame['FilterMask'].get('DeepCoreFilter_12').condition_passed or frame['FilterMask'].get('DeepCoreFilter_TwoLayerExp_12').condition_passed:
			preciousfood = True
	frame.Put('DCFilterStreamsEvent',icetray.I3Bool(preciousfood))

tray.AddModule(TransientAnalysisEventGrab,'grabule')

tray.AddModule("I3IcePickModule<I3SimpleFilter>",'DeepCoreStreamEvents',
		DecisionName = 'Fubar',
		InputDecisionName = 'DCFilterStreamsEvent',
		DiscardEvents = True
		)

tray.AddModule("Delete",'delibird',keys=['Fubar','CalibratedIceTopFADC_HLC','CalibratedIceTopATWD_SLC','CalibratedIceTopATWD_HLC',
	                                 'DirtyInIceRawData','CalibratedWaveforms','IceTopHLCVEMPulses','IceTopHLCPulseInfo',
					 'IceTopRawData','IceTopPulses_HLC','IceTopPulses_SLC','IceTopSLCVEMPulses','BaseProc_imprv_LF_linefit_final_rusage',
					 'CascadeFilter_CascadeLinefit_rusage','FSSFilter_pulses_hlc','FSSFilter_pulses_slc','InIceRawData',
                                         'IceTopCalibratedWaveforms','IceTopPulses','IceTopDSTPulses','CleanIceTopRawData','CleanInIceRawData','MCPMTResponseMap'])
### Add User Specified Further Processing ###

### Choked GRB Event Weighting ###

if grbevent:
  def MCWeightDictGrab(frame):
    if (frame.Has('I3MCWeightDict') or not frame.Has('PrimaryEnergy')):
      frame.Put('PrimaryEnergy',dataclasses.I3Double(frame['I3MCWeightDict'].get('PrimaryNeutrinoEnergy')))
      frame.Put('InteractionProb',dataclasses.I3Double(frame['I3MCWeightDict'].get('TotalInteractionProbabilityWeight'))) 
    else:
      pass
    
  tray.AddModule(MCWeightDictGrab,'MCGrab')  
'''
  def ChokedGRBFlux(energy,gamma,flavor,nevents,Agen,whichweight,IntProb):
    emaxp = (gamma*1.4e4/4)
    emaxk = (gamma*1.4e4/2)
    FRatio = {'numu':0.4,'nue':0.2,'nutau':0.4}
    pionic_nus = 45*(1/10.0-1/30.0) + 1350*(1/(2.0*900)-1/(2.0*10000)) + 135000*(1/(3.0*1000000)-1/(3.0*(emaxp)**3))
    kaonic_nus = 2*(1/10.0-1/200.0) + 400*(1/(2.0*(200)**2)-1/(2.0*(20000)**2)) + 400*20000*(1/(3.0*(20000)**3)-1/(3.0*(emaxk)**3))
    TotalNus = (kaonic_nus+pionic_nus)*Agen*FRatio[flavor]*10000
    WeightedMCNus = 0.83479014829153153*nevents
    ### Pionic ###
    if energy < 30:
      pionweight = 45.0
    elif energy < 100:
      pionweight = 45.0*(30*energy**-1)
    elif energy < emaxp:
      pionweight = 45.0*(30*100*energy**-2)
    else:
      pionweight = 0.0  
    ### Kaonic ###
    if energy < 200:
      kaonweight = (2.0)
    elif energy < 20000:
      kaonweight = (2.0)*(200*energy**-1)                        
    elif energy < emaxp:
      kaonweight = (2.0)*(200*20000*energy**-2) 
    else:
      kaonweight = 0.0
    spectralweight = (pionweight + kaonweight)/47.0
    grbevents = (spectralweight*(WeightedMCNus/TotalNus)**-1)*IntProb
  
    if whichweight == 'spectral':
      return spectralweight
    if whichweight == 'grbevents':  
      return grbevents

  def ChokedGRBWeighting(frame,gamma,flavor,nfiles,events_per_file,Agen):
    if frame.Has('PrimaryEnergy'):
      energy = frame['PrimaryEnergy'].value
      spectralweight = ChokedGRBFlux(energy,gamma,flavor,nfiles*events_per_file,Agen,'spectral',frame['InteractionProb'].value)
      grbeventweight = ChokedGRBFlux(energy,gamma,flavor,nfiles*events_per_file,Agen,'grbevents',frame['InteractionProb'].value)
      frame.Put("ChkGRBSpectrumWeight",dataclasses.I3Double(spectralweight))
      frame.Put("ChkGRBEventWeight",dataclasses.I3Double(grbeventweight))

  tray.AddModule(ChokedGRBWeighting,'grbweight',
    gamma=3,
    flavor='nue',
    nfiles = 1,
    Agen = 3.14159*1200**2,
    events_per_file = 100000
    )
'''

outputfile = inputfile.replace('/data/sim/IceCube/2011/generated/CORSIKA-in-ice/7437/03000-03999/','/net/user/daughjd/data/2012_L2/')
outputfile = outputfile.replace('2011_corsika','2012_Pro_Corsika')



# Write the physics and DAQ frames
tray.AddModule( "I3Writer", "EventWriter", filename=outputfile,
		Streams=[icetray.I3Frame.Physics,icetray.I3Frame.DAQ],
		DropOrphanStreams=[icetray.I3Frame.DAQ]
		)

tray.AddModule("TrashCan","byebye")

if prettyprint:
    print tray
    exit(0)

#tray.Execute()
tray.Execute()
#tray.Execute(3270)
#tray.Execute(200)

tray.Finish()

stop_time = time.asctime()

print 'Started:', start_time
print 'Ended:', stop_time

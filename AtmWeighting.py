import os
from icecube import icetray, dataio, dataclasses
from I3Tray import *

@icetray.traysegment
def NuOsc(tray, name):
	# """
	# Calibration and feature extraction.   
	# Parameters:
	# 	GCD
	# 	Description : Path to the GCD file that WaveCalibrator will use
	# 	Default     : 'noGCDdefined'
	# 	
	# """
	print "Got the nugen flag, going to add the flux calculations and tag 5 interesting particles"

	# # HACK HACK HACK shouldn't have to rerun this, but here we are. 
	# tray.AddModule("Delete","nugencorrection",keys=["AtmoWeight","AtmoWeightNoOsc"])

	#load("libmcsummary")
	#tray.AddModule( "I3MCTagger", "tagger",
	#	Mode = 2,
	#	MCTreeName = "I3MCTree",
	#	)
		
	#def isMySignal(frame):
	#	zenBool = (math.cos(frame["PrimaryNu"].dir.zenith) < 0)
	#
	#	position = frame["InteractionVertex"].pos
	#	x = position.x
	#	y = position.y
	#	z = position.z
	#	dist = math.sqrt(x*x + y*y)
	#	posBool = False
	#	if (dist<250 and  z<-100): posBool = True
	#	
	#	# because we don't have accurate NuFlux tables below 10Gev, Aug2012
	#	# and NeutrinoFlux will log_fatal if it sees an event with too low E
	#	fluxBool = frame["PrimaryNu"].energy>10
	#	
	#	myBool = zenBool and posBool and fluxBool
	#	frame["isMySignal"] = icetray.I3Bool(myBool)
	#
	#tray.AddModule(isMySignal,"superNu",Streams=[icetray.I3Frame.DAQ])

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
		If = lambda frame: frame["PrimaryNu"].energy>10
		)
	
	tray.AddService("NeutrinoFluxFactory", "numu-flux-service-bartol",
		InstallServiceAs="bartol_numu",
		)
	tray.AddService("NeutrinoFluxFactory", "nue-flux-service-bartol",
		ConventionalModelName = "bartol_nue",
		InstallServiceAs = "bartol_nue",
		PromptModelName = "naumov_rqpm_nue",
		)
	load("libatmo-weights")
	tray.AddModule("AtmoWeights", "bartolWeights",
		OutputName = "AtmoWeights_bartol",
		If = lambda frame: frame["PrimaryNu"].energy>10,
		DumpToFrame = False
	)


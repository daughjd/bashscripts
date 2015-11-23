from icecube import icetray, dataclasses
from icecube import dipolefit, clast, tensor_of_inertia, linefit, improvedLinefit, lilliput, cramer_rao
from icecube.icetray import I3Units
from icecube import finiteReco
from icecube import photonics_service
from icecube import coordinate_service
from icecube import gulliver_modules, MuonVariables
from icecube import paraboloid
from icecube import millipede
from icecube.icetray import I3Units
from icecube.common_variables import direct_hits

allstrings = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86]

@icetray.traysegment
def OfflineDeepCoreReco(tray, name, If = lambda f: True, suffix = '',
                        Pulses = '',
                                                ):

    tray.AddModule('I3DipoleFit', name + '_DipoleFit' + suffix,
                   Name = 'DipoleFit' + suffix,
                   InputRecoPulses = Pulses,   
                   AmpWeightPower = 0.0,
                   DipoleStep = 0,
                   MinHits = 5,   
                   If = If,
                   )

    tray.AddModule('I3CLastModule', name+'_CascadeLast' + suffix,
                   Name = 'CascadeLast' + suffix,
                   InputReadout = Pulses,
                   MinHits = 3, #Default 
                   AmplitudeOption = 1, #Default
                   AmplitudeWeight = 1.0, #Default
                   AmandaMode = False, #Default   
                   DirectHitRadius = 300.0*I3Units.m, #Default
                   DirectHitWindow = 200.0*I3Units.ns, #Default
                   DirectHitThreshold = 3, #Default
                   AmEnergyParam0 = 1.448, #Default
                   AmEnergyParam1 = -1.312, #Default
                   AmEnergyParam2 = 0.9294, #Default
                   AmEnergyParam3 = -0.1696, #Default
                   AmEnergyParam4 = 0.001123, #Default
                   AmEnergyParam5 = 0.0, #Default
                   EnergyParam0 = -0.431, #Default
                   EnergyParam1 = 1.842, #Default  
                   EnergyParam2 = -0.49199, #Default
                   EnergyParam3 = 0.07499, #Default 
                   EnergyParam4 = 0.0, #Default
                   EnergyParam5 = 0.0, #Default
                   InputSelection = '',
                   If = If,
                   )

    tray.AddModule('I3TensorOfInertia', name + '_ToI' + suffix,
                   InputReadout = Pulses,
                   Name = 'ToI'+suffix,  
                   MinHits = 3,
                   AmplitudeOption = 1,
                   AmplitudeWeight = 1,
                   InputSelection = '',
                   If = If,
                   )

    tray.AddSegment( improvedLinefit.simple,'LineFit'+suffix, inputResponse = Pulses, fitName = 'LineFit'+suffix, If = If)

    tray.AddSegment( lilliput.I3SinglePandelFitter, 'SPEFitSingle'+suffix, pulses = Pulses, seeds = ['LineFit'+suffix], If = If)
    
    tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit2'+suffix, pulses = Pulses, n_iterations = 2, seeds = [ 'SPEFitSingle'+suffix ], If = If)

    tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit6'+suffix, pulses = Pulses, n_iterations = 6, seeds = [ 'SPEFitSingle'+suffix ], If = If)

    #tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit16'+suffix, pulses = Pulses, n_iterations = 16, seeds = [ 'SPEFitSingle'+suffix ], If = If)

    #tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit32'+suffix, pulses = Pulses, n_iterations = 32, seeds = [ 'SPEFitSingle'+suffix ], If = If)

    #tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit64'+suffix, pulses = Pulses, n_iterations = 64, seeds = [ 'SPEFitSingle'+suffix ], If = If)

    #use only first hits.  Makes sense for an SPE likelihood
    tray.AddModule('CramerRao', 'SPE6'+suffix + '_SPEFitCramerRao' + suffix,
                   InputResponse = Pulses,
                   InputTrack = 'SPEFit6'+suffix,
                   OutputResult = 'SPEFit6CramerRao'+suffix,
                   AllHits = False, # ! doesn't make sense to use all hit for SPE pdf
                   DoubleOutput = False, # Default
                   z_dependent_scatter = True, # Default
                   If = If,
                   )

    #likelihood = lilliput.__add_pandel__(tray=tray, trayname = lilliput.default_trayname, pulses = Pulses, domllh = 'SPE1st')
    #reusing the minimizer from the SPEFit done earlier
    #minimizer = lilliput.__add_simplexminimizer__(tray = tray, trayname =  lilliput.default_trayname)

    #tray.AddService( "I3BasicSeedServiceFactory", 'SPEFit6' + "Seed"+suffix,
    #    FirstGuesses = ['SPEFit6'+suffix],
    #    )

    #tray.AddService( "I3BasicSeedServiceFactory", 'SPEFit6'+"SeedwithVertex"+suffix,
#	        FirstGuesses              =  ['SPEFit6'+suffix],
#	        InputReadout              = Pulses,              # ! Use pulses for vertex correction
#	        TimeShiftType             = "TFirst",            # ! Use TFirst for vertex correction
	        #ChargeFraction           = 0.9,                 # Default
	        #FixedEnergy              = float( "nan" ),      # Default
	        #MaxMeanTimeResidual      = 1000.0 * I3Units.ns, # Default
	        #NChEnergyGuessPolynomial = [],                  # Default
	        #SpeedPolice              = True,                # Default
	        #AddAlternatives          = "None",              # Default
	        #AltTimeShiftType         = "TFirst",            # Default
	        #OnlyAlternatives         = False,               # Default
	    )

#    tray.AddModule( "I3ParaboloidFitter", 'SPEFit6' + "Paraboloid" + suffix, 
#	        SeedService               = 'SPEFit6' + "Seed" + suffix,          # seedservice from above
#	        LogLikelihood             = likelihood,                      #likelihood service has been installed earlier - I guess you know how to do that!
#	        MaxMissingGridPoints      = 1,                               # ! Allow 1 missfit grid point
#	        VertexStepSize            = 5.0 * I3Units.m,                 # ! Use smallvertex size
#	        ZenithReach               = 2.0 * I3Units.degree,            # !
#	        AzimuthReach              = 2.0 * I3Units.degree,            # !
#	        GridpointVertexCorrection = 'SPEFit6'+"SeedwithVertex"+suffix,    # !Name of vertex correction service
#	        Minimizer                 = minimizer,                       # minimizerservice has to be installed
	        #NumberOfSamplingPoints   = 8,                               # Default
	        #NumberOfSteps            = 3,                               # Default
	        #MCTruthName              = "",                              # Default
#	        If                        = If
#	    )
    '''
    tray.AddModule("I3FirstPulsifier", "second-first-pulsify",
                    InputPulseSeriesMapName = Pulses,
                    OutputPulseSeriesMapName = Pulses+'_FIRSTONLY',
                    KeepOnlyFirstCharge = False,   # default
                    UseMask = False,               # default
		    If = If
                  )
    tray.AddModule("I3FirstPulsifier", "second-first-pulsify_all",
                    InputPulseSeriesMapName = 'SplitInIcePulses',
                    OutputPulseSeriesMapName = 'SplitInIcePulses_FIRSTONLY',
                    KeepOnlyFirstCharge = False,   # default
                    UseMask = False,               # default
                    If = If
                  )
    
    dh_defs = direct_hits.default_definitions


    tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'dha',
		    DirectHitsDefinitionSeries = dh_defs,
		    ParticleName = 'SPEFit6'+suffix,
		    OutputI3DirectHitsValuesBaseName = 'SPEFit6'+suffix+'_DirectHits',
		    PulseSeriesMapName =  Pulses+'_FIRSTONLY',
		    If = lambda f: f.Has('SPEFit6'+suffix)
                    )
    tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'dha_additional',
                    DirectHitsDefinitionSeries = dh_defs,
                    ParticleName = 'SPEFit6'+suffix,
                    OutputI3DirectHitsValuesBaseName = 'SPEFit6'+suffix+'_DirectHits_AllPulses',
                    PulseSeriesMapName =  'SplitInIcePulses_FIRSTONLY',
                    If = lambda f: f.Has('SPEFit6'+suffix)
                    )
    
    ### Add FiniteReco ###
    
    PhotonicsFiniteRecoTabledir   = "/nv/hp11/jdaughhetee3/data2/software/photonics-tables/SPICEMie_i3coords"
    PhotonicsFiniteRecoDriverdir  =  PhotonicsFiniteRecoTabledir+"/driverfiles"
    PhotonicsFiniteRecoDriverfile = "SPICEMie_i3coords_level2_muon_resampled.list"
   
    tray.AddService('I3GulliverMinuit2Factory', 'Minuit',
        Algorithm     = 'SIMPLEX',
        MaxIterations = 2000, 
        Tolerance     = 0.01
    )
    
   
    tray.AddService("I3PhotonicsServiceFactory","photonics_fr",
        PhotonicsTopLevelDirectory = PhotonicsFiniteRecoTabledir,
        DriverFileDirectory        = PhotonicsFiniteRecoDriverdir,
        PhotonicsLevel2DriverFile  = PhotonicsFiniteRecoDriverfile,
        PhotonicsTableSelection    = 2, 
        ServiceName                = "PhotonicsServiceFiniteReco"
    )
    
	   
	   
    tray.AddService("I3GulliverFinitePhPnhFactory", "PhPnh"+suffix,
        InputReadout   = Pulses,                       #put the HLC+SLC pulsemap here
        PhotorecName   = "PhotonicsServiceFiniteReco",
        ProbName       = "PhPnhPhotorec",
        #RCylinder     = 300,                          #the radius around the track, should be the same as below
        RCylinder      = 200,                          #should it be the same as what was used before?
        SelectStrings  = allstrings,     # Do we need this?
        StringLLH      = "true",
    )

    tray.AddModule("I3StartStopPoint", "qVertexReco"+suffix,
        Name            =  'SPEFit4'+suffix,             #put the SPE4Fit here
        InputRecoPulses =  Pulses,                    #be sure to take a HLC+SLC pulse series
        ExpectedShape   =   70,                       #contained track, this way the start AND stop point are reconstructed
        CylinderRadius  =  200,                       #only for DOMs within this radius around the track finiteReco calculates hit probabilities
        If = If
    )
	   
   
    ## #### keep the MC direction and only variate the vertex
    ## ## stop point reco : take start point and variate length
    tray.AddService("I3SimpleParametrizationFactory","simparStartVertex"+suffix,
        StepLinL     = 25.0*I3Units.m,
        BoundsLinL      = [0,2*I3Units.km],
        VertexMode   = "Default",
    )
	
    tray.AddService("I3BasicSeedServiceFactory","seedserveStartVertex"+suffix,
        FirstGuess    = 'SPEFit4'+suffix+"_Finite",
        InputReadout  =  Pulses,
        TimeShiftType = "TFirst" 
    )
	   
    tray.AddModule("I3SimpleFitter",'SPEFit4'+suffix+"_Finite_FitLengthTrack",
        SeedService     = "seedserveStartVertex"+suffix,
        Parametrization = "simparStartVertex"+suffix,
        LogLikelihood   = "PhPnh"+suffix,
        Minimizer       = "Minuit",
        StoragePolicy   = "OnlyBestFit",
        If = If
    )
    tray.AddService("I3SimpleParametrizationFactory","simparStopVertex"+suffix,
        StepLinL   = 25.0*I3Units.m,
        BoundsLinL = [0,2*I3Units.km],
        VertexMode = "Stop",
    )

   
    ## # seed service using the stop position -> set option + USE RESULTS OF STOP RECO HERE
    tray.AddService("I3BasicSeedServiceFactory","seedserveStopVertex"+suffix,
        FirstGuess    = 'SPEFit4'+suffix+"_Finite_FitLengthTrack",
        InputReadout  = Pulses,
        TimeShiftType =   "TFirst",
    )
	   
    tray.AddModule("I3SimpleFitter",'SPEFit4'+suffix+"_Finite_Start",
        SeedService     = "seedserveStopVertex"+suffix,
        Parametrization = "simparStopVertex"+suffix,
        LogLikelihood   = "PhPnh"+suffix,
        Minimizer       = "Minuit",
        StoragePolicy   = "OnlyBestFit",
        If = If
    )
	   
    tray.AddModule("I3StartStopLProb","qLLHR"+suffix, # gets + "_StartStopParams"
	        Name        = 'SPEFit4'+suffix+"_Finite_Start", 
	        ServiceName = "PhPnh"+suffix,
	        If = If
    )
    
    tray.AddSegment(millipede.MillipedeLowEnergyMuon,'millipede',
		    If=If,
#		    MuonPhotonicsService='PhotonicsServiceFiniteReco',
		    Pulses=Pulses,
		    ReadoutWindow='SplitUncleanedInIcePulsesTimeRange',
		    SeedTrack='SPEFit6'+suffix)
    '''	    

    tray.AddModule("Delete",name+"DC_Cleanup",keys=['ToI_DCEval2','ToI_DCEval3','LineFit_DC_debiasedPulses','LineFit_DC_linefit_final_rusage','DeepCoreL2Reco_DipoleFit_DC_rusage'])

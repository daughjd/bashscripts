## extract data from HDF5 format

import tables
import cPickle as pickle
import os

## Loading and saving
def load (filename):
    """Load `filename` using the pickle module."""
    with open (filename) as f:
        out = pickle.load (f)
        try:
            out.__cache_source_filename = filename
        except:
            pass
        return out

def save (obj, filename):
    """Dump `obj` to `filename` using the pickle module."""
    with open (filename, 'wb') as f:
        pickle.dump (obj, f, -1)


## Extracting useful parameters
def columns (fit, colnames):
    return [fit + '.' + colname for colname in colnames]

def direction (fit):
    return columns (fit, ['azimuth', 'zenith', 'fit_status'])

def energy (fit):
    return columns (fit, ['energy'])

def value (fit):
    return columns (fit, ['exists', 'value'])

def params (fit):
    # TODO: probably actually need .fit_status, not FitParams.exists
    ## relationship between .fit_status and FitParams.exists?
    return columns (fit, ['exists', 'logl', 'ndof', 'rlogl'])

def cscdllhparams (fit):
    ## In CscdLlh, seems to be logL convention instead of logl
    return columns (fit, ['exists', 'logL', 'rlogL'])

def fitstatus (fit):
    # good == 0
    return columns (fit, ['fit_status'])

def fitcuts (fit):
    return columns (fit, ['exists', 'l_dir', 'n_chan', 'n_dir'])

def toiparams (fit):
    return columns (fit, ['exists', 'mineval', 'eval2', 'eval3', 'evalratio'])

def flatten (x):
    out = []
    for xi in x:
        if isinstance (xi, list):
            for xij in xi:
                out.append (xij)
        else:
            out.append (xi)
    return out

## example of data to retrieve
column_specs_example = [
    # general stuff
    columns ('FilterMask',
        ['CascadeFilter_10', 'EHEFilter_10', 'MuonFilter_10', 'ICOnlineL2Filter_10']),
    'I3EventHeader',

    # cut variable keys
    rescaparams ('RescaVals'),  ## cascade cramer-rao
    toiparams ('ToIParams'),
    direction ('ToI'),
    direction ('TWSRTLineFit'),
    vertex ('TWSRTLineFit'),
    'TWSRTLineFit.speed',
    lfparams ('TWSRTLineFitParams'),
    direction ('CascadeLlhVertexFit'),
    vertex ('CascadeLlhVertexFit'),
    cscdllhparams ('CascadeLlhVertexFitParams'),
    direction ('SPEFit_SPEFitSingleseed_4iter'),  
    vertex ('SPEFit_SPEFitSingleseed_4iter'),
    params ('SPEFit_SPEFitSingleseed_4iterFitParams'),
    direction ('MpodFit_CredoFitSeed_5Iter'),
    energy ('MpodFit_CredoFitSeed_5Iter'),
    vertex ('MpodFit_CredoFitSeed_5Iter'),
    params ('MpodFit_CredoFitSeed_5IterFitParams'),
    'MpodFit_CredoFitSeed_5IterFitParams.qtotal',
    'MpodFit_CredoFitSeed_5IterFitParams.predicted_qtotal',
    ]

## I booked and named some of the values below from the I3MCTree
nugen_specs = sorted (flatten (column_specs + ['I3MCWeightDict', 'MCFirstInIceCascade', 
                                                     'MCPrimary',
                                                     'MCMaxCascade', 'MCMaxInIce']))

## See attached __init__.py file
### Put it in some 'hdfdataset' folder and put that in your python path (you can use the version at pub.icecube.wisc.edu:/home/hellauer/hdfdataset_example with the CMakeLists.txt and put it in your icerec build)
### This is code written by Mike Richman with additions by me and other UMD people to stack HDF5 data (hdf5.root.var.cols.stuff) into a python file.data.var.stuff format
### It also includes weighting for simulation
### For an example, check out pub.icecube.wisc.edu:/home/hellauer/diffsim_example
### Load the file in ipython with pickle (see load function above) and use <tab> to see what's in there

## If in icecube I3_BUILD dir:                          
from icecube import hdfdataset
from icecube.hdfdataset import I3HDFNugenDataSet

def trim_nugen (self):
    """Trim nugen data."""
    h5_files = sorted (glob ('{0}/Level3*nugen*.h5'.format (input_dir)))
    filename = 'IC86I.NuE.CscdL3.hdfdataset'
    mctruth = 'MCPrimary'
    h = I3HDFNugenDataSet (h5_files, load=nugen_specs, verbose=-2, truth=mctruth)
    save (h, filename)


## Extract only the necessary data and make variable calculations, if necessary
### Use arrays.py for dictionary object['key'] and attribute object.key functionality, also used in grbllh_roadmap.py (file.stuff format)
### You don't need to use this format for grbllh (it takes numpy arrays or lists as input), but it makes the datasets easy to deal with all at once
### For an example, check out pub.icecube.wisc.edu:/home/hellauer/diffsim_example
import arrays

def get_l3_cscd_arrays (dataset, season="", flavor=""):

    data = dataset.data

    out = arrays.Arrays ()

    linefit = data.TWSRTLineFit
    linefitparams = data.TWSRTLineFitParams

    if season == "86I":
        cscdllh = data.CascadeLlhVertexFit
        cscdllhparams = data.CascadeLlhVertexFitParams
        spe4 = data.SPEFit_SPEFitSingleseed_4iter  ## Calculated in L3 processing
        spe4params = data.SPEFit_SPEFitSingleseed_4iterFitParams
        out.spefit_zenith = spe4.zenith

        monopod = data.MpodFit_CredoFitSeed_5Iter
        monopodparams = data.MpodFit_CredoFitSeed_5IterFitParams

    rescaparams = data.RescaVals
    out.cr_theta = np.sqrt (rescaparams.CovMat_theta)
    out.cr_phi = np.sqrt (rescaparams.CovMat_phi)
    out.cramer_rao = 1 / np.sqrt (2) * np.sqrt (
        out.cr_theta**2 + out.cr_phi**2 * np.sin (monopod.zenith)**2)

    hitmultiplicity = data.HitMultiplicityValues

    out.charge_per_string = np.log10(monopodparams.qtotal) / hitmultiplicity.n_hit_strings

    if isinstance (dataset, hdfdataset.I3HDFNugenDataSet):
        ## True neutrino (direction, energy)
        if ("86" in season) and (flavor == "nue"):
            truth = data.MCPrimary
        else:
            truth = data.MCTruth
        out.true_zenith = truth.zenith
        out.true_azimuth = truth.azimuth
        out.true_energy = truth.energy
        out.type = truth.type  #for nu and nu_bar id
        out.oneweight = dataset.oneweight
        out.weight_E2 = dataset.get_weights ('power', -2)  ## see L921 in hdfdataset/python/__init__.py file for this function
        out.weight_E23 = dataset.get_weights ('power', -2.3)
        assert (flavor == "nue" or flavor == "numu" or flavor == "nutau")
        if flavor != 'nutau' and flavor != "":
            out.weight_atm = dataset.get_weights ('honda2006_{0}'.format(flavor))
        out.meta.n_gen = dataset.n_gen

    elif isinstance (dataset, hdfdataset.I3HDFExpDataSet):
        out.time = dataset.start_datetime
        out.meta.livetime = dataset.livetime
        out.weights = np.ones_like (out.lfv) / out.meta.livetime
        out.run = data.I3EventHeader.Run

    else:
        raise ValueError ('not prepared for {0}'.format (type (dataset)))

    return out

## make 'diffsim' object used in grbllh_roadmap.py
nugen_nue = cuts.get_l3_cscd_arrays (nugen_hdfdataset, season, "nue")
save (nugen_nue, 'nugen_nue.arrays')

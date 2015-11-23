## extract data from HDF5 format

import tables
import cPickle as pickle
import os

from scipy.interpolate import UnivariateSpline

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

def variable_map (fit):
    return columns(fit, ['NCh_Clean'])

def vertex (fit):
    return columns (fit, ['exists','x','y','z','length','speed'])

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
column_specs = [
    # general stuff
    columns ('FilterMask',
        ['CascadeFilter_13', 'MuonFilter_13', 'DeepCoreFilter_13']),
    'I3EventHeader',

    # cut variable keys
    direction ('SplineMPEMod'),
    vertex ('SplineMPEMod_Contained'),
    value ('ReScaled_Paraboloid_Sigma_SplineMPEMod'),
    direction ('LineFit_DC'),
    params ('SplineMPEModFitParams'),
    variable_map ('L4_Variables_LES'),
    variable_map ('L4_Variables_HES'),
    ]

## I booked and named some of the values below from the I3MCTree
nugen_specs = sorted (flatten (column_specs + ['I3MCWeightDict']))

## See attached __init__.py file
### Put it in some 'hdfdataset' folder and put that in your python path (you can use the version at pub.icecube.wisc.edu:/home/hellauer/hdfdataset_example with the CMakeLists.txt and put it in your icerec build)
### This is code written by Mike Richman with additions by me and other UMD people to stack HDF5 data (hdf5.root.var.cols.stuff) into a python file.data.var.stuff format
### It also includes weighting for simulation
### For an example, check out pub.icecube.wisc.edu:/home/hellauer/diffsim_example
### Load the file in ipython with pickle (see load function above) and use <tab> to see what's in there

## If in icecube I3_BUILD dir:                          
from icecube import hdfdataset
from icecube.hdfdataset import I3HDFNugenDataSet

def trim_nugen ():
    """Trim nugen data."""
    h5_files = "/data/user/daughjd/grbsearch/hdftabs/FinalSample_nugen_numu_IC86.2013.010090.Table.hdf5"
    filename = 'IC86II.NugenNuMu10090.FinalLevel.hdfdataset'
    mctruth = 'PrimaryNu'
    h = I3HDFNugenDataSet (h5_files, load=nugen_specs, verbose=-2, truth=mctruth)
    save (h, filename)
    return h

### Define New Spectral Weighting Options ###
def get_subphotospheric_weights(dataset):
	truth = dataset.truth
	fluxspline = pickle.load(open("data_objects/subphotofluxspline.pkl","r"))  ### Valid from 5-1000 GeV
	@np.vectorize
	def subphotoflu(energy):
		if energy < 1000.0:
			return fluxspline(energy)
		else:
			#return 0.0008144194258991092 * (energy/800.0)**-2
			return 0.0
	fluence = subphotoflu(truth.energy)
	weights = dataset.oneweight * fluence / dataset.n_gen
	return weights


def get_chkgrb_weights (dataset, gamma=3.0, Ej=10.0**51.5):
        """Get Choked GRB Weights."""
        truth = dataset.truth
        emaxp = gamma*1.4*10000/4
        hadcoolbreak_pionic = pow(gamma,5)*(3.7037037037037037*pow(10.0,50))/Ej
        radcoolbreak_pionic = 33.333333*gamma

        emaxk = gamma*1.4*10000/2
        hadcoolbreak_kaonic = pow(gamma,5)*(2.4691358024691357*pow(10.0,51))/Ej
        radcoolbreak_kaonic = 6666.66666*gamma

        ebreak1_pionic = 0.0
        ebreak2_pionic = 0.0
        ebreak1_kaonic = 0.0
        ebreak2_kaonic = 0.0

        if hadcoolbreak_pionic <= radcoolbreak_pionic:
                ebreak1_pionic = hadcoolbreak_pionic
                ebreak2_pionic = radcoolbreak_pionic
        if hadcoolbreak_pionic > radcoolbreak_pionic:
                ebreak1_pionic = radcoolbreak_pionic
                ebreak2_pionic = hadcoolbreak_pionic

        if hadcoolbreak_kaonic <= radcoolbreak_kaonic:
                ebreak1_kaonic = hadcoolbreak_kaonic
                ebreak2_kaonic = radcoolbreak_kaonic
        if hadcoolbreak_kaonic > radcoolbreak_kaonic:
                ebreak1_kaonic = radcoolbreak_kaonic
                ebreak2_kaonic = hadcoolbreak_kaonic
	@np.vectorize
        def individual_weight(ebreak1_pionic,ebreak2_pionic,ebreak1_kaonic,ebreak2_kaonic,emaxp,emaxk,energy):
                pionic = 0.0
                kaonic = 0.0
                if (energy > emaxk):
                    return 0

                if (energy > emaxp):
                    pionic = 0.0
                elif (energy < ebreak1_pionic):
                    pionic = np.power(energy,-2.0)
                elif (energy < ebreak2_pionic):
                    pionic = np.power(energy,-3.0)*ebreak1_pionic
                else:
                    pionic = np.power(energy,-4.0)*ebreak1_pionic*ebreak2_pionic

                if (energy > emaxk):
                    kaonic = 0.0
                elif (energy < ebreak1_kaonic):
                    kaonic = np.power(energy,-2.0)
                elif (energy < ebreak2_kaonic):
                    kaonic = np.power(energy,-3.0)*ebreak1_kaonic
                else:
                    kaonic = np.power(energy,-4.0)*ebreak1_kaonic*ebreak2_kaonic

                return 0.94*pionic + 0.06*kaonic

        fluence = individual_weight(ebreak1_pionic,ebreak2_pionic,ebreak1_kaonic,ebreak2_kaonic,emaxp,emaxk,truth.energy)
        weights = dataset.oneweight * fluence / dataset.n_gen
        return weights

## Extract only the necessary data and make variable calculations, if necessary
### Use arrays.py for dictionary object['key'] and attribute object.key functionality, also used in grbllh_roadmap.py (file.stuff format)
### You don't need to use this format for grbllh (it takes numpy arrays or lists as input), but it makes the datasets easy to deal with all at once
### For an example, check out pub.icecube.wisc.edu:/home/hellauer/diffsim_example

nug_numu = trim_nugen()

import arrays
import numpy as np

def get_simdata_arrays (dataset,simtype="nugen"):

    data = dataset.data

    out = arrays.Arrays ()

    linefit = data.LineFit_DC
    recofit = data.SplineMPEMod
    les_l4_vars = data.L4_Variables_LES
    hes_l4_vars = data.L4_Variables_LES
    out.nch = les_l4_vars.NCh_Clean + hes_l4_vars.NCh_Clean
    ## Reconstruction Info ##
    out.reco_zenith = recofit.zenith
    out.reco_azimuth = recofit.azimuth
    out.reco_err = data.ReScaled_Paraboloid_Sigma_SplineMPEMod.value
    ## True neutrino (direction, energy)
    truth = data.PrimaryNu
    out.true_zenith = truth.zenith
    out.true_azimuth = truth.azimuth
    out.true_energy = truth.energy
    #out.eproxy = 100.*(truth.energy * truth.energy**-1)
    out.eproxy = 5.0 * out.nch
    out.type = truth.type  #for nu and nu_bar id
    out.oneweight = dataset.oneweight
    out.weight_E2 = dataset.get_weights ('power', -2)  ## see L921 in hdfdataset/python/__init__.py file for this function
    out.weight_E3 = dataset.get_weights ('power', -3)
    out.weight_E25 = dataset.get_weights ('power', -2.5)
    #out.weight_E23 = dataset.get_weights ('power', -2.3)
    #assert (flavor == "nue" or flavor == "numu" or flavor == "nutau")
    #if flavor != 'nutau' and flavor != "":
    #    out.weight_atm = dataset.get_weights ('honda2006_{0}'.format(flavor))
    out.chkgrb_weight = get_chkgrb_weights(dataset,gamma=3.0, Ej=10.0**51.5)
    out.subphoto_weight = get_subphotospheric_weights(dataset)
    out.meta.n_gen = dataset.n_gen
    return out

## make 'diffsim' object used in grbllh_roadmap.py
nugen_numu = get_simdata_arrays (nug_numu,"nugen")
save (nugen_numu, 'data_objects/nugen_numu.arrays')


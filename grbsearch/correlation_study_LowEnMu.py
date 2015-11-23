import grbllh, grblist
import datetime,astrodate
import arrays
import pickle
import tables
from vars_class import Vars
from arrays import Arrays
from icecube import hdfdataset
import fakeps,time,os
import numpy as np

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


grbs = pickle.load(open("data_objects/grbarray.pkl",'r'))              ### GRB List between May 15th 2012 -> April 30th 2013
data_vars = pickle.load(open("data_objects/fulldatavars.pkl",'r'))     ### Full Data Set
exp = pickle.load(open("data_objects/bckgrdvars.pkl",'r'))             ### Data outside of 2 hr windows about GRBs
diffsim = pickle.load(open("data_objects/nugen_numu.arrays",'r'))     ### GENIE 1460 NuMu Simulation Events @ Final Level

### Define Analysis Type ###
llh_type = grbllh.LlhType.per_confchan

execfile("grbllh_roadmap.py")
llh_psarrays()
pssim=load("ps_numu.arrays")

llh_sources_events()
bg_sources = load('bg.sources')
bg_events = load('bg.events')
ps_sources= load('ps.sources')
ps_events = load('ps.events')
sig_events_diff = load('diff.events')

llh_pdfs()
pdf_space_bg = load('pdf_space_bg.pickle')
pdf_ratio_energy = load('pdf_ratio_energy.pickle')

llh_throwers()
#sig_thrower_numu = load('ps_numu_E2.thrower')
#sig_thrower_numu = load('ps_numu_ChkGRB.thrower')
sig_thrower_numu = load('ps_numu_SubPhoto.thrower')

bg_thrower = load('bg_thrower.pickle')

n_trials = 1e8


#llh_null_tsd()   ## Calculate null Test Statistic Distribution
null_tsd = load('tsd_1e8.pickle')

n_trials= 1e3

#llh_sensitivity_ps(beta=0.9,thresh=5,median=True)  ## Calculate Sensitivity

#sens = load('disc.pickle')


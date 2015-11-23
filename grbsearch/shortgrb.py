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
import pylab
import subprocess

from scipy.interpolate import UnivariateSpline

Norm = 1.0

def ShortGRBDistro(z):
        z_0=0.32
        return Norm * z**2 * np.e**(-z/z_0)


zs = np.linspace(0,8,10000)
dz = np.diff(zs)[0]

Norm = ((ShortGRBDistro(zs)*dz).sum())**-1

@np.vectorize
def CDFShortGRB(z):
        cumulative = np.linspace(0,z,10000)
        dizzy = np.diff(cumulative)[0]
        return (ShortGRBDistro(cumulative).sum() * dizzy) / ((ShortGRBDistro(zs).sum()*dz))


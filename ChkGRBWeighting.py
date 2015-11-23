import os
from icecube import icetray, dataio, dataclasses
from I3Tray import *

def ChokedGRBFlux(energy,gamma,flavor,nevents,Agen,whichweight,IntProb,inspectral,simemin,simemax,Ej=3*10**51):
  emaxp = (gamma*1.4e4/4)
  emaxk = (gamma*1.4e4/2)
  FRatio = {'numu':0.4,'nue':0.2,'nutau':0.4}
  hadcoolbreak_pionic = (gamma**5)*(3.7037037037037037e+50)/Ej
  radcoolbreak_pionic = 33.333333*gamma
  hadcoolbreak_kaonic = (gamma**5)*(2.4691358024691357e+51)/Ej
  radcoolbreak_kaonic = 6666.66666*gamma

  if simemax > emaxk:
    pionic_nus = 45*(1/simemin-1/30.0) + 1350*(1/(2.0*900)-1/(2.0*10000)) + 135000*(1/(3.0*1000000)-1/(3.0*(emaxp)**3))
    kaonic_nus = 2*(1/simemin-1/200.0) + 400*(1/(2.0*(200)**2)-1/(2.0*(20000)**2)) + 400*20000*(1/(3.0*(20000)**3)-1/(3.0*(emaxk)**3))
    Normalization_p=pionic_nus/(1/(simemin)**((inspectral-2)+1)-1/(emaxp)**((inspectral-2)+1))
    Normalization_k=kaonic_nus/(1/(simemin)**((inspectral-2)+1)-1/(emaxk)**((inspectral-2)+1))
  else:
    if simemax < 30.0:
      pionic_nus = 45*(1/simemin-1/simemax)
      kaonic_nus = 2*(1/simemin-1/simemax)
      Normalization_p=45*(1/(simemin)-1/(simemax))/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
      Normalization_k=2*(1/(simemin)-1/(simemax))/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
    elif simemax < 100:
      pionic_nus = 45*(1/simemin-1/30.0) + 1350*(1/(2.0*900)-1/(2.0*simemax**2))
      kaonic_nus = 2*(1/simemin-1/simemax)
      Normalization_p=pionic_nus/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
      Normalization_k=2*(1/(simemin)-1/(simemax))/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
    elif simemax < 200 and simemax < emaxp:
      pionic_nus = 45*(1/simemin-1/30.0) + 1350*(1/(2.0*900)-1/(2.0*10000)) + 135000*(1/(3.0*1000000)-1/(3.0*(simemax)**3))
      kaonic_nus = 2*(1/simemin-1/simemax)
      Normalization_p=pionic_nus/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
      Normalization_k=2*(1/(simemin)-1/(simemax))/(1/(simemin)**((inspectral-2)+1)-1/(simemax)**((inspectral-2)+1))
  ### Normalization ###
  Normalization=(Normalization_p+Normalization_k)/47.0
  

  TotalNus = (kaonic_nus+pionic_nus)*Agen*FRatio[flavor]*10000
  WeightedMCNus = Normalization*nevents
  ### Pionic ###
  pionweight=0.0
  kaonweight=0.0
  if energy < (gamma**5)*(3.7037037037037037e+50)/Ej:
    pionweight = 45.0*(energy**-(2-inspectral))
  elif energy < 33.333333*gamma:
    pionweight = 45.0*(30*energy**-(3-inspectral))
  elif energy < emaxp:
    pionweight = 45.0*(30*100*energy**-(4-inspectral))
  elif energy >= emaxp:
    pionweight = 0.0  
  ### Kaonic ###
  if energy < (gamma**5)*(2.4691358024691357e+51)/Ej:
    kaonweight = (2.0)*(energy**-(2-inspectral))
  elif energy < 6666.66666*gamma:
    kaonweight = (2.0)*(200*energy**-(3-inspectral))                        
  elif energy < emaxk:
    kaonweight = (2.0)*(200*20000*energy**-(4-inspectral)) 
  elif energy >= emaxk:
    kaonweight = 0.0
  spectralweight = (pionweight + kaonweight)/47.0
  grbevents = (spectralweight*(WeightedMCNus/TotalNus)**-1)*IntProb
  
  if whichweight == 'spectral':
    return spectralweight
  if whichweight == 'grbevents':  
    return grbevents

  


def ChokedGRBWeighting(frame,gamma,flavor,nfiles,weightdict,specweightname,eventweightname):
  if frame.Has('PrimaryNu'):
    energy = frame['PrimaryNu'].energy
    Agen=3.14159*(frame[weightdict]['InjectionSurfaceR'])**2
    simspectrumindex = frame[weightdict]['PowerLawIndex']
    simemin=10**(frame[weightdict]['MinEnergyLog'])
    simemax=10**(frame[weightdict]['MaxEnergyLog'])
    events_per_file = frame[weightdict]['NEvents']
    intprob = frame[weightdict]['TotalInteractionProbabilityWeight']
    spectralweight = ChokedGRBFlux(energy,gamma,flavor,nfiles*events_per_file,Agen,'spectral',intprob,simspectrumindex,simemin,simemax)
    grbeventweight = ChokedGRBFlux(energy,gamma,flavor,nfiles*events_per_file,Agen,'grbevents',intprob,simspectrumindex,simemin,simemax)
    frame.Put(specweightname,dataclasses.I3Double(spectralweight))
    frame.Put(eventweightname,dataclasses.I3Double(grbeventweight))


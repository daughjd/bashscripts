import sys,tables, numpy,pickle,pylab

from subprocess import call

### Data ###
datums  = tables.openFile('LowNchan_ResSlice_CompleteDataSet.330Days_FINAL.hdf5')

finalbool_data = datums.root.Final_Level_Bool.col('value')
splinemod_zen_data = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_data = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_data = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_data = datums.root.SplineMPEModParaboloid.col('azimuth')
finite_x_data = datums.root.SplineMPEMod_Contained.col('x')
finite_y_data = datums.root.SplineMPEMod_Contained.col('y')
finite_z_data = datums.root.SplineMPEMod_Contained.col('z')
pulse_time_data = datums.root.LET_RecoPulses.col('time')
pulse_event_data = datums.root.LET_RecoPulses.col('Event')
pulse_run_data = datums.root.LET_RecoPulses.col('Run')
pulse_charge_data = datums.root.LET_RecoPulses.col('charge')
rescaled_paraboloid_data = datums.root.ReScaled_Paraboloid_Sigma_SplineMPEMod.col('value')

datums.close()

### GENIE NuMu ###

datums  = tables.openFile('LowNchan_ResSlice_genie_ic.1460.FINAL.hdf5')

finalbool_numu = datums.root.Final_Level_Bool.col('value')
splinemod_zen_numu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_numu = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_numu = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_numu = datums.root.SplineMPEModParaboloid.col('azimuth')
finite_x_numu = datums.root.SplineMPEMod_Contained.col('x')
finite_y_numu = datums.root.SplineMPEMod_Contained.col('y')
finite_z_numu = datums.root.SplineMPEMod_Contained.col('z')
pulse_time_numu = datums.root.LET_RecoPulses.col('time')
pulse_event_numu = datums.root.LET_RecoPulses.col('Event')
pulse_run_numu = datums.root.LET_RecoPulses.col('Run')
pulse_charge_numu = datums.root.LET_RecoPulses.col('charge')
rescaled_paraboloid_numu = datums.root.ReScaled_Paraboloid_Sigma_SplineMPEMod.col('value')
primary_energy_numu = datums.root.PrimaryNu.col('energy')
primary_zenith_numu = datums.root.PrimaryNu.col('zenith')
primary_azi_numu = datums.root.PrimaryNu.col('azimuth')
atmo_numu = datums.root.AtmoWeight.col('value')

datums.close()

### Nugen NuMu ###

datums  = tables.openFile('LowNchan_ResSlice_nugen_numu_IC86.2013.010090.FINAL.hdf5')

finalbool_nugmu = datums.root.Final_Level_Bool.col('value')
splinemod_zen_nugmu = datums.root.SplineMPEMod.col('zenith')
splinemod_azi_nugmu = datums.root.SplineMPEMod.col('azimuth')
spline_para_zen_nugmu = datums.root.SplineMPEModParaboloid.col('zenith')
spline_para_azi_nugmu = datums.root.SplineMPEModParaboloid.col('azimuth')
finite_x_nugmu = datums.root.SplineMPEMod_Contained.col('x')
finite_y_nugmu = datums.root.SplineMPEMod_Contained.col('y')
finite_z_nugmu = datums.root.SplineMPEMod_Contained.col('z')
pulse_time_nugmu = datums.root.LET_RecoPulses.col('time')
pulse_event_nugmu = datums.root.LET_RecoPulses.col('Event')
pulse_run_nugmu = datums.root.LET_RecoPulses.col('Run')
pulse_charge_nugmu = datums.root.LET_RecoPulses.col('charge')
rescaled_paraboloid_nugmu = datums.root.ReScaled_Paraboloid_Sigma_SplineMPEMod.col('value')
primary_energy_nugmu = datums.root.PrimaryNu.col('energy')
primary_zenith_nugmu = datums.root.PrimaryNu.col('zenith')
primary_azi_nugmu = datums.root.PrimaryNu.col('azimuth')
atmo_nugmu = datums.root.AtmoWeight.col('value')

datums.close()


dataweight=1.0/(28475585.0)
dataweight=dataweight*numpy.ones(len(finalbool_data))

atmo_numu = atmo_numu*10**9/4000.
atmo_nugmu = atmo_nugmu*10**9/1000.

def AngleBetweenAngles(theta1,phi1,theta2,phi2):
   x1 = numpy.sin(theta1)*numpy.cos(phi1)
   y1 = numpy.sin(theta1)*numpy.sin(phi1)
   z1 = numpy.cos(theta1)
   x2 = numpy.sin(theta2)*numpy.cos(phi2)
   y2 = numpy.sin(theta2)*numpy.sin(phi2)
   z2 = numpy.cos(theta2)
   return numpy.arccos(x1*x2+y1*y2+z1*z2)


spline_primary_diff_numu = AngleBetweenAngles(primary_zenith_numu,primary_azi_numu,splinemod_zen_numu,splinemod_azi_numu)
spline_primary_diff_nugmu = AngleBetweenAngles(primary_zenith_nugmu,primary_azi_nugmu,splinemod_zen_nugmu,splinemod_azi_nugmu)


### MC Investigations ###

low_samplepulsearrays_data=numpy.array([])
low_samplepulsearrays_numu=numpy.array([])
low_samplepulsearrays_nugmu=numpy.array([])

#high_samplepulsearrays_data=numpy.array([])
#high_samplepulsearrays_numu=numpy.array([])
#high_samplepulsearrays_nugmu=numpy.array([])

pulse_id_data = pulse_event_data * pulse_run_data
pulse_id_nugmu = pulse_event_nugmu * pulse_run_nugmu
pulse_id_numu = pulse_event_numu * pulse_run_numu

q_nch_data = []
q_nch_nugmu = []
q_nch_numu = []

for i in numpy.unique(pulse_id_data)[:3000]:
	q_nch_data.append(pulse_charge_data[pulse_id_data==i].sum() / (1.0 * len(pulse_time_data[pulse_id_data==i])))
	low_samplepulsearrays_data = numpy.hstack([low_samplepulsearrays_data,pulse_time_data[pulse_id_data==i]-numpy.median(pulse_time_data[pulse_id_data==i])])

for i in numpy.unique(pulse_id_nugmu)[:3000]:
	q_nch_nugmu.append(pulse_charge_nugmu[pulse_id_nugmu==i].sum() / len(pulse_time_nugmu[pulse_id_nugmu==i]))
        low_samplepulsearrays_nugmu = numpy.hstack([low_samplepulsearrays_nugmu,pulse_time_nugmu[pulse_id_nugmu==i]-numpy.median(pulse_time_nugmu[pulse_id_nugmu==i])])

for i in numpy.unique(pulse_id_numu)[:3000]:
        q_nch_numu.append(pulse_charge_numu[pulse_id_numu==i].sum() / len(pulse_time_numu[pulse_id_numu==i]))
        low_samplepulsearrays_numu = numpy.hstack([low_samplepulsearrays_numu,pulse_time_numu[pulse_id_numu==i]-numpy.median(pulse_time_numu[pulse_id_numu==i])])

q_nch_data = numpy.array(q_nch_data)
q_nch_numu = numpy.array(q_nch_numu)
q_nch_nugmu = numpy.array(q_nch_nugmu)


'''
for i in numpy.unique(pulse_id_data[(res_slice_data)*(nolow_nch_real_final_level_data==1)])[:3000]:
        high_samplepulsearrays_data = numpy.hstack([high_samplepulsearrays_data,pulse_time_data[pulse_id_data==i]-numpy.median(pulse_time_data[pulse_id_data==i])])

for i in numpy.unique(pulse_id_nugmu[(res_slice_data)*(nolow_nch_real_final_level_nugmu==1)])[:3000]:
        high_samplepulsearrays_nugmu = numpy.hstack([high_samplepulsearrays_nugmu,pulse_time_nugmu[pulse_id_nugmu==i]-numpy.median(pulse_time_nugmu[pulse_id_nugmu==i])])

for i in numpy.unique(pulse_id_numu[(res_slice_data)*(nolow_nch_real_final_level_data==1)])[:3000]:
        high_samplepulsearrays_numu = numpy.hstack([high_samplepulsearrays_numu,pulse_time_numu[pulse_id_numu==i]-numpy.median(pulse_time_numu[pulse_id_numu==i])])

pylab.figure()
pylab.hist(q_nch_data,bins=numpy.linspace(0,4,20),histtype='step',lw=2,ec='k',label="Data",normed=True)
pylab.hist(q_nch_nugmu,bins=numpy.linspace(0,4,20),histtype='step',lw=2,ec='b',label="Nugen",normed=True)
pylab.hist(q_nch_numu,bins=numpy.linspace(0,4,20),histtype='step',lw=2,ec='g',label="GENIE",normed=True)
pylab.xlabel("Charge / Nch")
pylab.ylabel("Normed Counts")
pylab.grid()
pylab.legend()
pylab.savefig("ChargeNchRatioDistros")




	
'''



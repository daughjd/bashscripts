import grbllh, numpy, grblist
import datetime,astrodate
import arrays
import pickle
import tables
from vars_class import Vars
from arrays import Arrays
from icecube import hdfdataset
import jdutil

### Transient Analysis Time Frame ###
search_start = datetime.datetime(2012,05,12)
search_stop = datetime.datetime(2013,04,30)



grbs = grblist.get_grbs_from_grbweb_file(grb_file_name='lowen_trans_grblist.txt',t_min=search_start,t_max=search_stop)

def DataGrab(hdf):
	dat=tables.openFile(hdf)
	zen = dat.root.SplineMPEMod.col('zenith')
	azi = dat.root.SplineMPEMod.col('azimuth')
	err = dat.root.ReScaled_Paraboloid_Sigma_SplineMPEMod.col('value')
	nch = dat.root.FinalLevelNch.col('value')
	run = dat.root.I3EventHeader.col('Run')
	time_mjd = dat.root.timeMJD.col('value')
	time_jd = time_mjd + 2400000.5
	time=[]
	for temp in time_jd:
		time.append(jdutil.jd_to_datetime(temp))
	time = numpy.array(time)
	data_vars = Vars()
	data_vars.zenith = zen
	data_vars.azimuth = azi
	data_vars.error = err
	data_vars.nchan = nch
	data_vars.run = run
	data_vars.timeMJD = time_mjd
	data_vars.time = time
	data_vars.eproxy = 5.0*data_vars.nchan
	return data_vars


def BckGrdDataGrab(hdf):
	dat=tables.openFile(hdf)
	zen = dat.root.SplineMPEMod.col('zenith')
	azi = dat.root.SplineMPEMod.col('azimuth')
	err = dat.root.ReScaled_Paraboloid_Sigma_SplineMPEMod.col('value')
	nch = dat.root.FinalLevelNch.col('value')
	run = dat.root.I3EventHeader.col('Run')
	time_mjd = dat.root.timeMJD.col('value')
	time_jd = time_mjd + 2400000.5
	time=[]
	for temp in time_jd:
		time.append(jdutil.jd_to_datetime(temp))
	time = numpy.array(time)
	ontime_bool = numpy.ones(len(nch)) == 1.
	for indy,event_time in enumerate(time_mjd):
		flip_bool=False
		jdevent=event_time+astrodate.MJD_0
		for burst in grbs.grbs:
			jdburst_start = astrodate.JulianDate(burst.t_start)
			jdburst_end = astrodate.JulianDate(burst.t_end)
			cleared = False
			if jdevent - jdburst_start.jd < -0.0833:
				cleared=True
			if jdevent - jdburst_end.jd > 0.0833:
				cleared=True
			if not cleared:
				flip_bool = True
		if flip_bool:
			ontime_bool[indy] = False

	data_vars = Vars()
	data_vars.zenith = zen[ontime_bool]
	data_vars.azimuth = azi[ontime_bool]
	data_vars.error = err[ontime_bool]
	data_vars.nchan = nch[ontime_bool]
	data_vars.timeMJD = time_mjd[ontime_bool]
	data_vars.run = run[ontime_bool]
	data_vars.time = time[ontime_bool]
        data_vars.eproxy = 5.0*data_vars.nchan
	return data_vars

full_data = DataGrab("hdftabs/FinalLevel_IC86.2012.CompleteDataSet.L6Out.330Days.hdf5")
bg_data = BckGrdDataGrab("hdftabs/FinalLevel_IC86.2012.CompleteDataSet.L6Out.330Days.hdf5")

pickle.dump(grbs,open("data_objects/grbarray.pkl",'w'))
pickle.dump(full_data,open("data_objects/fulldatavars.pkl",'w'))
pickle.dump(bg_data,open("data_objects/bckgrdvars.pkl",'w'))



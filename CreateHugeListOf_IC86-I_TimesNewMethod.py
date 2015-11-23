#!/usr/bin/env python
"""
this script creates the list of times used for timescrambling
it is based on the previous script by Juanan (/net/user/aguilar/work/IceCube/psLab_RHEL_6.0_amd64/macro_llh/ic79/PrintOutAllEvTimesIC79.C)
it samples 1E6 times from a histogram of the times of the original (real data) events
I take the same original (real data) events as Jake in (/net/user/jfeintzeig/2012/PointSources/scripts/BDTs/CalcTotalLiveTime.py)
"""

import glob, numpy, tables
from ROOT import TH1D

def main():

    nbins = 316.*(24.*15.) #4min bins?
    tmin = 55694.4164795 #tmin and tmax were obtained running this script and looping over all events
    tmax = 56062.4181092

    minhisto = TH1D("minhisto","minhisto",int(nbins),tmin,tmax)
    minhisto1= TH1D("minhisto1","minhisto1",int(nbins),tmin,tmax)

    dates=[]
    runlist=file('Prelim_IC86-I_v1.5a.txt')
    for line in runlist.readlines()[77:]:
        dates.append([int(line.split()[0]),line.split()[1]])
   
    runs={}
    for item in dates:
        month=item[1].split('-')[1]
        day=item[1].split('-')[2]
        year=item[1].split('-')[0]
        runs[str(item[0])]=(year,month,day)

    for run in runs:
        print "filling times for run %s" % (run)
        # get date
        year,month,day=runs[run]
        l3files=glob.glob('/data/ana/Muon/level3/exp/%s/%s%s/*%s*.hd5' % (year,month,day,run))
        for item in l3files:
            try:
                f=tables.openFile(item)
                time_start_mjd_day=f.root.I3EventHeader.col('time_start_mjd_day')
                time_start_mjd_sec=f.root.I3EventHeader.col('time_start_mjd_sec')
                time_start_mjd_ns =f.root.I3EventHeader.col('time_start_mjd_ns')
                for i in range(0,len(time_start_mjd_day)):
                    time = time_start_mjd_day[i]+(time_start_mjd_sec[i]+time_start_mjd_ns[i]*1e-9)/86400.;
                    bin=minhisto1.Fill(time)
                    if minhisto.GetBinContent( bin )==0:
                        minhisto.Fill(time)
                f.close()
            except Exception, e:
                    print "problem with file ", item, " pr: ",e
                    
    outfile = open('HugeListOf_IC86-I_TimesNewMethod.txt', 'w')
    for i in range(0,1000000):
        rndtime=minhisto.GetRandom()
        outfile.write(str("%.11f" % rndtime)+"\n")
    outfile.close()


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    main()

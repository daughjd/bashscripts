#/bin/bash
for i in {1..100}

do
  echo "Processing injection with 3.0 Index"
  echo "#PBS -N DecHotSpotTestRA_ChkGRBSensitivity$i" > secondsubmit.pbs 
  echo "#PBS -l nodes=1:ppn=1" >> secondsubmit.pbs
  echo "#PBS -l walltime=40:00:00" >> secondsubmit.pbs
  echo "#PBS -q cygnusforce-6" >> secondsubmit.pbs
  echo "#PBS -l mem=4000mb" >> secondsubmit.pbs
  echo "#PBS -e /nv/hp11/jdaughhetee3/data2/logs/" >> secondsubmit.pbs
  echo "#PBS -o /nv/hp11/jdaughhetee3/data2/logs/" >> secondsubmit.pbs
  echo "/nv/hp11/jdaughhetee3/data2/software/psLab-LowE/start.sh<<EOF" >> secondsubmit.pbs
  echo "cd /nv/hp11/jdaughhetee3/data2/software/psLab-LowE/macro_llh/low-energy_ic86-2/qScripts/" >> secondsubmit.pbs
  echo 'root -l -q -b SensitivityGen_IC86II_LowEn.C\\("'$i'05","3.0","10.0**51.5","true"\\)' >> secondsubmit.pbs
  echo "EOF" >> secondsubmit.pbs
  msub -h secondsubmit.pbs
  sleep 0.5
done

  
# Ejet coefficients 0.1 , 0.1274275 , 0.16237767 , 0.20691381 , 0.26366509 ,   0.33598183 ,   0.42813324 ,   0.54555948 , 0.6951928 , 0.88586679 , 1.12883789 , 1.43844989 ,
#	            1.83298071 , 2.33572147 , 2.97635144 , 3.79269019 , 4.83293024 , 6.15848211 , 7.8475997  , 10.        


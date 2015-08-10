#!/bin/bash
while true; do
  RUNNING=`showq -u jdaughhetee3 | grep " Running " | awk '{print $1}'`

  NRUN=$(expr ${#RUNNING} + 1)
  let NRUN="NRUN / 12"

  HELD=`showq -u jdaughhetee3 | grep " UserHold " | awk '{print $1}'`
  NHELD=$(expr ${#HELD} + 1)
  let NHELD="NHELD / 12"

  echo $NRUN Jobs Running
  echo $NHELD Jobs Held

  if [ $NHELD -lt 1 ]
    then
      exit
  fi

  NDIF=0

  if [ $NRUN -lt 1200 ]
    then
      NDIF=$(expr 1200 - $NRUN)
      if [ $NDIF -gt $NHELD ]
        then
          let NDIF=$NHELD
      fi
  fi

  echo Submitting $NDIF more jobs!

  COUNT=0

  while [ $COUNT -lt $NDIF ]
  do
    HELD=`showq -u jdaughhetee3 | grep " UserHold " | awk '{print $1}'`
    jobby=$(echo $HELD | awk 'BEGIN {FS=" "}{print $1}')
    echo $jobby
    COUNT=$(($COUNT+1))
    mjobctl -u all $jobby
    sleep 2
  done

  if [ "$NHELD" -lt 1 ]
    then
      exit
  fi
  sleep 300
done



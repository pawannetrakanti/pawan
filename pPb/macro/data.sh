#!/bin/sh                                                                                                                                                                                                           
kspecies="ppb"
for (( j=0 ;  j < 6; j++ ))
  do
  leadb=$j
  for (( k=0 ;  k < 6; k++ ))
    do
    sleadb=$k
#    echo "Submitting process for $kspecies : Leading bin : $leadb and Sub-leading bin : $sleadb "
    ./rundata.sh $leadb $sleadb $kspecies 
#   runs for pT1 > 120 GeV/c  and pT2 > 30 GeV/c
#      ./rundata.sh 2 0 $kspecies  
    echo "----------------------"
  done
done


#!/bin/bash

POSNU_RAD=(3000 4000 5500 7000 9000 12000)
RAYSOL_RAN=13000
ENERGY_EXP=(16 17 18 19 20 21)
for j in $(seq 0 1 5)
do
  echo ${ENERGY_EXP[$j]}
  cp setupRNO_comparison.txt setupRNO_comparison_E"${ENERGY_EXP[$j]}".txt
  sed -i "s/EXPONENT=18/EXPONENT=${ENERGY_EXP[$j]}/g" setupRNO_comparison_E"${ENERGY_EXP[$j]}".txt
  sed -i "s/POSNU_RADIUS=4000/POSNU_RADIUS=${POSNU_RAD[$j]}/g" setupRNO_comparison_E"${ENERGY_EXP[$j]}".txt
  sed -i "s/DATA_SAVE_MODE=0/DATA_SAVE_MODE=1/g" setupRNO_comparison_E"${ENERGY_EXP[$j]}".txt

done

DATA_SAVE_MODE=0

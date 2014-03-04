#!/bin/bash

#this takes the file copied from the user guide references and makes it into a setup.txt file for AraSim. 

#argument is the textfile with the linebroken references (each reference should have an # in front of it)

#(cat $1 | tr -d '\n' | sed 's/\#/\n/g'  ; echo) > setup.txt

#(cat $1 | tr -s [:space:] | sed 's/\n/\t/g'  ; echo) > setup.txt


(cat $1 | sed ':begin;$!N;s/\n/ /;tbegin' | sed 's/\#/\n/g' ) > setup.txt

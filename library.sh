#!/bin/bash

#this creates the library libSim.so needed for any dynamically loaded Detector / Event etc... clases


# comm="gcc -shared -Wl,--no-as-needed -o libSim.so `ls | grep -E ^[^A]*\.o$ | tr '\n' ' '` `ls AraRoot | grep -E \.o$ | awk '{print "AraRoot/" $1}' | tr '\n' ' '` `ls AraRootFormat | grep -E \.o$ | awk '{print "AraRootFormat/" $1}' | tr '\n' ' '`"

comm="gcc -shared -Wl,--no-as-needed -L$ARA_UTIL_INSTALL_DIR/lib -lAraEvent -o libSim.so `ls | grep -E ^[^A]*[^s]o$ | tr '\n' ' '`"

echo "$comm"

eval "$comm"

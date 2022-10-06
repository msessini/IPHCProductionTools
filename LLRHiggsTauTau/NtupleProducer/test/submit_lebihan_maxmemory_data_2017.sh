#!/bin/bash

source /cvmfs/cms.cern.ch/crab3/crab.sh
eval 'scramv1 runtime -sh'


##MC
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py


crab submit crab3_DataB_2017.py

crab submit crab3_DataC_2017.py

crab submit crab3_DataD_2017.py

crab submit crab3_DataE_2017.py

crab submit crab3_DataF_2017.py


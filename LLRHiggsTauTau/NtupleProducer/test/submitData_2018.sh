#!/bin/bash

source /cvmfs/cms.cern.ch/crab3/crab.sh
eval 'scramv1 runtime -sh'


##MC
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py


crab submit config2018/data/crab3_DataA_2017.py

crab submit config2018/data/crab3_DataB_2017.py

crab submit config2018/data/crab3_DataC_2017.py

#crab submit crab3_DataE_2017.py



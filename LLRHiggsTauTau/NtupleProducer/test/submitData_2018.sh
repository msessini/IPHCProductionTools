#!/bin/bash

##Data
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py
sed -i 's/.*IsEmbed=.*/IsEmbed=False/g' analyzer.py
sed -i 's/.*PERIOD =.*/PERIOD = '\'' '\''/g' analyzer.py

crab submit config2018/data/crab3_SingleMuonA_2018.py

crab submit config2018/data/crab3_SingleMuonB_2018.py

crab submit config2018/data/crab3_SingleMuonC_2018.py

sed -i 's/.*PERIOD =.*/PERIOD = '\''D'\''/g' analyzer.py

crab submit config2018/data/crab3_SingleMuonD_2018.py



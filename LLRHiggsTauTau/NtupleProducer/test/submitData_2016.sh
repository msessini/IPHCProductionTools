#!/usr/bin/bash

##Data
sed -i 's/.*YEAR   =.*/YEAR   = 2016/g' analyzer.py
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py
sed -i 's/.*IsEmbed=.*/IsEmbed=False/g' analyzer.py
sed -i 's/.*PERIOD =.*/PERIOD = '\'' '\''/g' analyzer.py
sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string(" ")/g' ../python/HiggsTauTauProducer.py

crab submit config2016/data/crab3_SingleMuonB_2016.py

crab submit config2016/data/crab3_SingleMuonC_2016.py

crab submit config2016/data/crab3_SingleMuonD_2016.py

crab submit config2016/data/crab3_SingleMuonE_2016.py

crab submit config2016/data/crab3_SingleMuonF_2016.py

crab submit config2016/data/crab3_SingleMuonG_2016.py

crab submit config2016/data/crab3_SingleMuonH_2016.py




#!/bin/bash

##Embed
sed -i 's/.*YEAR   =.*/YEAR   = 2017/g' analyzer.py
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py
sed -i 's/.*IsEmbed=.*/IsEmbed=True/g' analyzer.py
sed -i 's/.*PERIOD =.*/PERIOD = '\'' '\''/g' analyzer.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2017/embedded/crab3_EmbedB_2017.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2017/embedded/crab3_EmbedC_2017.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2017/embedded/crab3_EmbedD_2017.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2017/embedded/crab3_EmbedE_2017.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2017/embedded/crab3_EmbedF_2017.py

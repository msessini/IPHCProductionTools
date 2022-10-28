#!/bin/bash

##Embed
sed -i 's/.*IsMC=.*/IsMC=False/g' analyzer.py
sed -i 's/.*IsEmbed=.*/IsEmbed=True/g' analyzer.py
sed -i 's/.*PERIOD =.*/PERIOD = '\'' '\''/g' analyzer.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/embedded/crab3_EmbedA_2018.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/embedded/crab3_EmbedB_2018.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/embedded/crab3_EmbedC_2018.py

sed -i 's/.*PERIOD =.*/PERIOD = '\''D'\''/g' analyzer.py
sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("DY_mutau_embedded")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/embedded/crab3_EmbedD_2018.py


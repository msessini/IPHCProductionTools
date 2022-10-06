sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_1qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/DYJets/crab3_DY1Jets_ll_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_2qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/DYJets/crab3_DY2Jets_ll_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_3qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/DYJets/crab3_DY3Jets_ll_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_4qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/DYJets/crab3_DY4Jets_ll_2018.py

##########

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_lnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/WJets/crab3_WJets_lnu_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_1qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/WJets/crab3_W1Jets_lnu_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_2qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/WJets/crab3_W2Jets_lnu_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_3qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/WJets/crab3_W3Jets_lnu_2018.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_4qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2018/WJets/crab3_W4Jets_lnu_2018.py


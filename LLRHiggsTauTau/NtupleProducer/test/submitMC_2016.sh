#!/usr/bin/bash

##MC
sed -i 's/.*YEAR   =.*/YEAR   = 2016/g' analyzer.py
sed -i 's/.*IsMC=.*/IsMC=True/g' analyzer.py
sed -i 's/.*IsEmbed=.*/IsEmbed=False/g' analyzer.py
sed -i 's/.*PERIOD =.*/PERIOD = '\'' '\''/g' analyzer.py

##TTbar
sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("ttbar_dilep")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/TTbar/crab3_TT_dilep_2016.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("ttbar_hadr")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/TTbar/crab3_TT_hadr_2016.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("ttbar_semilep")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/TTbar/crab3_TT_semilep_2016.py

########
##EWK

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("EWKWplus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWplus2Jets_Wlnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKWplus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWplus2Jets_Wlnu_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKWplus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWplus2Jets_Wlnu_v3_2016.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("EWKWminus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWminus2Jets_Wlnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKWminus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWminus2Jets_Wlnu_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKWminus2Jets_Wlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKWminus2Jets_Wlnu_v3_2016.py

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("EWKZ2Jets_Zll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKZ2Jets_Zll_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKZ2Jets_Zll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKZ2Jets_Zll_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("EWKZ2Jets_Zll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_EWKZ2Jets_Zll_v3_2016.py

#########
##WZ

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WZ_2l2q")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WZ/crab3_WZ_2l2q_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WZ_3l1nu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WZ/crab3_WZ_3l1nu_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WZ_1l3nu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WZ/crab3_WZ_1l3nu_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WZ_1l1nu2q")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WZ/crab3_WZ_1l1nu2q_2016.py  

#########
##ZZ

sed -i 's/.*DataMCType.*/				       DataMCType    = cms.untracked.string("ZZ_4l")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/ZZ/crab3_ZZ_4l_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("ZZ_2l2q")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/ZZ/crab3_ZZ_2l2q_2016.py

#########
##VV

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("VV_2l2nu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_VV_2l2nu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("VV_2l2nu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/EWK/crab3_VV_2l2nu_v2_2016.py

#########
##WW

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WW_1l1nu2q")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WW/crab3_WW_1l1nu2q_2016.py

#########
##singleTop

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("tw")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/singleTop/crab3_ST_tW_top_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("tbarw")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/singleTop/crab3_ST_tW_antitop_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("ST_tchannel_top")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/singleTop/crab3_ST_tchannel_top_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("ST_tchannel_antitop")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/singleTop/crab3_ST_tchannel_antitop_2016.py

#########
##signal

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WminusH_tautau")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/signal/crab3_WminusH_tautau_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("WplusH_tautau")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/signal/crab3_WplusH_tautau_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("ZH_tautau")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/signal/crab3_ZH_tautau_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("H_tautau_ggF")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/signal/crab3_H_tautau_ggF_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("H_tautau_VBF")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/signal/crab3_H_tautau_VBF_2016.py

##########
##DYJets

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_ll_10to50")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DYJets_ll_10to50_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_ll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DYJets_ll_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_ll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DYJets_ll_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_1qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DY1Jets_ll_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_2qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DY2Jets_ll_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_3qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DY3Jets_ll_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("DY_4qll")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/DYJets/crab3_DY4Jets_ll_2016.py

##########
##WJets

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_lnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_WJets_lnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_lnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_WJets_lnu_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_lnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_WJets_lnu_v3_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_lnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_WJets_lnu_v4_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_1qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W1Jets_lnu_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_2qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W2Jets_lnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_2qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W2Jets_lnu_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_3qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W3Jets_lnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_3qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W3Jets_lnu_v2_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_4qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W4Jets_lnu_v1_2016.py

sed -i 's/.*DataMCType.*/                                      DataMCType    = cms.untracked.string("W_4qlnu")/g' ../python/HiggsTauTauProducer.py
crab submit config2016/WJets/crab3_W4Jets_lnu_v2_2016.py
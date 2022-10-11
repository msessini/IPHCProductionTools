# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "H_tautau_VBF_2018"
config.General.workArea = "crab3_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.JobType.maxMemoryMB=4500
config.JobType.inputFiles = (['/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/data','/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/MVADM'])

config.section_("Data")
config.Data.inputDataset = "/VBFHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 1000
config.Data.totalUnits = 1000 #number of event
config.Data.outLFNDirBase= '/store/user/msessini/Prod_Octobre22'
config.Data.outputDatasetTag = "H_tautau_VBF_2018"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "TT_semilep_2016"
config.General.workArea = "crab3_production2016"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.maxMemoryMB=4500
config.JobType.inputFiles = (['/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/data','/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/MVADM'])

config.section_("Data")
config.Data.inputDataset = "/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 50000
config.Data.totalUnits = 44000000#-1
config.Data.outLFNDirBase= '/store/user/msessini/Prod_Octobre22_v2'
config.Data.outputDatasetTag = "TT_semilep_2016"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
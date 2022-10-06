# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "ST_tW_antitop_2018"
config.General.workArea = "crab3_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.JobType.maxMemoryMB=2500

config.section_("Data")
config.Data.inputDataset = "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 20000
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/msessini/Prod_2018_v3'
config.Data.outputDatasetTag = "ST_tW_antitop_2018"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

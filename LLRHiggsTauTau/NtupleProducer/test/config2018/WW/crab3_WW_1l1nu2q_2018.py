# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "crab3_WW_1l1nu2q_2018"
config.General.workArea = "crab3_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.JobType.maxMemoryMB=2500

config.section_("Data")
config.Data.inputDataset = "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 60000
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/msessini/Prod_2018_v3'
config.Data.outputDatasetTag = "WW_1l1nu2q_2018"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

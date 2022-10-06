# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "ST_tchannel_top_2018"
config.General.workArea = "crab3_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.maxMemoryMB=4500

config.section_("Data")
config.Data.inputDataset = "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 3
config.Data.totalUnits =700#-1 #number of event #number of files here...
config.Data.outLFNDirBase = '/store/user/msessini/Prod_2018_v3'
config.Data.outputDatasetTag = "ST_tchannel_top_2018"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

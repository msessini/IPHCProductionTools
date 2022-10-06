# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "EmbedC_2018"
config.General.workArea = "crab3_Embed_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.JobType.maxMemoryMB = 3500

config.section_("Data")
config.Data.inputDataset = "/EmbeddingRun2018C/MuTauFinalState-inputDoubleMu_102X_miniAOD-v1/USER"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 50000
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/msessini/EmbedProd_2018_v3'
config.Data.outputDatasetTag = "EmbedC_2018"
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
config.Site.whitelist = ['T2_DE_RWTH','T2_CH_CERN','T2_FR_IPHC']

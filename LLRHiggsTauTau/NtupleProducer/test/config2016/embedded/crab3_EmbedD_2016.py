# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "EmbedD_2016"
config.General.workArea = "crab3_Embed_production2016"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = (['/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/data','/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/MVADM'])

config.section_("Data")
config.Data.inputDataset = "/EmbeddingRun2016D/MuTauFinalState-inputDoubleMu_94X_Legacy_miniAOD-v5/USER"
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 40
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase= '/store/user/msessini/Prod_Embed2016'
config.Data.outputDatasetTag = "EmbedD_2016"
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
config.Site.whitelist = ['T2_DE_RWTH','T2_CH_CERN','T2_FR_IPHC']

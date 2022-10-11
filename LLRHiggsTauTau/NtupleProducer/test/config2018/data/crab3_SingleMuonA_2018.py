# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "SingleMuon_A_2018"
config.General.workArea = "crab3_Data_production2018"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py'

config.JobType.inputFiles = (['/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/data','/opt/sbg/cms/safe1/cms/msessini/ProductionTools/CMSSW_10_2_23/MVADM'])
config.section_("Data")
config.Data.inputDataset = "/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD"
config.Data.inputDBS = 'global'
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 2
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase= '/store/user/msessini/Prod_Octobre22'
config.Data.outputDatasetTag = "SingleMuon_A_2018"

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

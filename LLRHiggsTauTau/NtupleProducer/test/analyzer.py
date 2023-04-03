#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os, re
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

# Year/Period
YEAR   = 2018
PERIOD = ' '

#samples list (it could be moved to a cfg file for better reading
#samples = [
#]
#apply corrections?
APPLYMUCORR=False
APPLYELECORR=True
APPLYFSR=False #this is by far the slowest module (not counting SVFit so far)
#Cuts on the Objects (add more cuts with &&)
#MUCUT="(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && abs(eta)<2.4 && pt>8"
#ELECUT="abs(eta)<2.5 && gsfTrack.trackerExpectedHitsInner.numberOfHits<=1 && pt>10"
#TAUCUT="pt>15"
#JETCUT="pt>15"

USEPAIRMET=False # input to SVfit: true: MVA pair MET; false: PFmet (HF inclusion set using USE_NOHFMET)
APPLYMETCORR=False # flag to enable (True) and disable (False) Z-recoil corrections for MVA MET response and resolution
USE_NOHFMET = False # True to exclude HF and run on silver json


SVFITBYPASS=True # use SVFitBypass module, no SVfit computation, adds dummy userfloats for MET and SVfit mass
#USECLASSICSVFIT=True # if True use the ClassicSVfit package, if False use the SVFitStandAlone package

BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)
APPLYTESCORRECTION=False # shift the central value of the tau energy scale before computing up/down variations
COMPUTEUPDOWNSVFIT=False # compute SVfit for up/down TES variation
COMPUTEMETUPDOWNSVFIT=False # compute SVfit for up/down MET JES variation
doCPVariables=False # compute CP variables and PV refit
COMPUTEQGVAR = False # compute QG Tagger for jets
IsMC=False
IsEmbed=True
Is25ns=True
HLTProcessName='HLT' #Different names possible, check e.g. at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.
if not IsMC:
    HLTProcessName='HLT' #It always 'HLT' for real data
print "HLTProcessName: ",HLTProcessName

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byDeepTau2017v2p1VSjetraw"#"byIsolationMVA3oldDMwoLTraw"
PVERTEXCUT=""#"!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
MUCUT="isLooseMuon && pt>10"
ELECUT="pt>10"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="pt>20" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>0"
LLCUT="mass>-99999"
BCUT="pt>5"

# ------------------------
DO_ENRICHED=False # do True by default, both ntuples and enriched outputs are saved!
STORE_ENRICHEMENT_ONLY=True # When True and DO_ENRICHED=True only collection additional to MiniAOD standard are stored. They can be used to reproduce ntuples when used together with oryginal MiniAOD with two-file-solution
# ------------------------

execfile(PyFilePath+"python/HiggsTauTauProducer.py")


### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    # - Sync files -
    # Legacy 2016
    #'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/0026BF69-A58A-E811-BBA8-1866DA87AB31.root'
    #'/store/data/Run2016F/SingleMuon/MINIAOD/17Jul2018-v1/00000/002F631B-E98D-E811-A8FF-1CB72C1B64E6.root'
    #'/store/data/Run2016C/Tau/MINIAOD/17Jul2018-v1/40000/FC60B37F-468A-E811-8C7B-0CC47A7C3424.root'
    #'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/270000/EC774BCF-16E9-E811-B699-AC1F6B0DE2F4.root'
    #'/store/mc/RunIISummer16MiniAODv3/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/E8B2E181-FF3E-E911-B72D-0025904C7FC2.root'

    # Legacy 2017
    #'/store/data/Run2017D/SingleElectron/MINIAOD/31Mar2018-v1/90000/FEFBCFEA-5939-E811-99DA-0025905B85D6.root'
    #'/store/data/Run2017F/SingleMuon/MINIAOD/31Mar2018-v1/80004/EEC206A7-5737-E811-8E0A-0CC47A2B0388.root'
    #'/store/data/Run2017C/Tau/MINIAOD/31Mar2018-v1/90000/FE83BD44-CD37-E811-8731-842B2B180922.root'
    #'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/0056F0F4-4F44-E811-B415-FA163E5B5253.root'
    #'/store/mc/RunIIFall17MiniAODv2/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/A6C340CF-9588-E811-9E52-782BCB46E733.root',
    #'/store/mc/RunIIFall17MiniAODv2/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/941B9AB4-A189-E811-BA1C-0242AC130002.root'

    # Legacy 2018
    #'/store/data/Run2018D/EGamma/MINIAOD/22Jan2019-v2/110000/2FD9167C-A197-E749-908D-3809F7B83A23.root'
    #'/store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/0011A0F3-23A3-D444-AD4D-2A2FFAE23796.root'
    #'/store/data/Run2018D/Tau/MINIAOD/PromptReco-v2/000/320/497/00000/584D46DF-5F95-E811-9646-FA163E293146.root'
    #'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/B3F93EA2-04C6-E04E-96AF-CB8FAF67E6BA.root'
    #'/store/mc/RunIIAutumn18MiniAOD/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/250000/F22484D3-E820-F040-8E63-0864695D025B.root'

    #'/store/mc/RunIIAutumn18MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/062A981D-4A57-664A-A583-E803A658594B.root'
    #'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/B3F93EA2-04C6-E04E-96AF-CB8FAF67E6BA.root'
    #'/store/mc/RunIIAutumn18MiniAOD/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/1F27B94F-20C0-A341-A7AC-8D240C836C21.root'
    #'/store/mc/RunIIAutumn18MiniAOD/VBFHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/240000/F55690CF-BE27-D845-ADD0-4C794EDB202D.root'
    #'/store/mc/RunIIAutumn18MiniAOD/GluGluHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/9A8D0B74-11F2-A147-AE75-056DE50F66A2.root'
    #'/store/mc/RunIIAutumn18MiniAOD/ZHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/40000/AF84FB68-A30A-4D4A-9965-82CDCD7AC2E1.root'
    #'/store/mc/RunIIAutumn18MiniAOD/GluGluHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/11345215-5D59-3A48-8895-86A56D94D62F.root' #ths file for singular matrix 1:86577324
    #'/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/1/merged_100.root',
    #'/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/1/merged_1000.root'
    #'/store/user/aakhmets/gc_storage/MuTau_data_2017_CMSSW944_gridka/TauEmbedding_MuTau_data_2017_CMSSW944_Run2017F/1/merged_100.root'
    #'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/0034E0CE-45EA-2249-A559-62746A83F4F1.root'
    #'/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/59/merged_8958.root'
    #'/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/270000/161FC6D9-A2CE-F340-A45C-1BE5C1D68A01.root'
    #'/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/120000/40707AA5-7FB4-CB46-8763-7E92612DFB1C.root'
    #'/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/325/022/00000/75DF89C5-1538-204F-8E3F-CA967EA985EE.root'
    #'/store/mc/RunIIAutumn18MiniAOD/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/389B502B-B27C-B447-A3ED-F0C6062C8C03.root'
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_6547.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_6647.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_6747.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_6847.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_6947.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7047.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7147.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7247.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7347.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7447.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_747.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7547.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7647.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7747.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7847.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_7947.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8047.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8147.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8247.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8347.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8447.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_847.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8547.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8647.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8747.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8847.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_8947.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9047.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9147.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9247.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9347.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9447.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_947.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9547.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/48/merged_9647.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1048.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1148.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1248.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1348.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1448.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_148.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1548.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1648.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1748.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1848.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_1948.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2048.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2148.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2248.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2348.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2448.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_248.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2548.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2748.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2848.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_2948.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3048.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3148.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3248.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3348.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3448.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_348.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3548.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3648.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3748.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3848.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_3948.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4048.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4148.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4248.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4348.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4448.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_448.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4548.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4648.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4748.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_48.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4848.root', 
    '/store/user/aakhmets/gc_storage/MuTau_data_2018ABC_CMSSW1020/TauEmbedding_MuTau_data_2018ABC_CMSSW1020_Run2018B/49/merged_4948.root'
    ),
)

# process.source.skipEvents = cms.untracked.uint32(968)
# process.source.eventsToProcess = cms.untracked.VEventRange("1:4006-1:4006") # run only on event=2347130 (syntax= from run:evt - to run:evt)
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#    "325022:1152-325022:1152", "325022:1218-325022:1218", "325022:1345-325022:1345", "325022:1038-325022:1038", "325022:1150-325022:1150", 
#    "325022:1156-325022:1156", "325022:1168-325022:1168", "325022:1200-325022:1200", "325022:1226-325022:1226", "325022:1318-325022:1318", 
#    "325022:1092-325022:1092", "325022:1141-325022:1141", "325022:1205-325022:1205", "325022:1114-325022:1114", "325022:1144-325022:1144", 
#    "325022:1149-325022:1149", "325022:1154-325022:1154"
#)
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("319991:214-319991:214", "319991:181-319991:181")

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("319579:647-319579:647", "319579:514-319579:514")
#process.source.lumisToSkip = cms.untracked.VLuminosityBlockRange(
#	"317212:46-317212:46", "317212:167-317212:167", "317212:168-317212:168", "317527:331-317527:331", "317527:332-317527:332", "317527:336-317527:336",
#	"317182:481-317182:481", "317182:480-317182:480", "317182:482-317182:482", "317182:488-317182:488", "317320:446-317320:446", "317320:447-317320:447",
#	"317320:984-317320:984", "317696:207-317696:207", "317696:208-317696:208", "317696:210-317696:210", 
#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = -1

# JSON mask for data --> defined in the lumiMask file
# from JSON file
if not IsMC:
  if YEAR == 2016:
    execfile(PyFilePath+"python/lumiMask_2016.py")
  if YEAR == 2017:
    execfile(PyFilePath+"python/lumiMask_2017.py")
  if YEAR == 2018:
    execfile(PyFilePath+"python/lumiMask_2018.py")
  process.source.lumisToProcess = LUMIMASK

##
## Output file
##

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   TauSpinnerReco = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
    )
						  )

process.p = cms.Path(process.Candidates)

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *
from RecoTauTag.RecoTau.PATTauDiscriminationAgainstElectronMVA6_cfi import *


taus = "slimmedTaus"


# Refitted Vertices collection
process.load("HiggsCPinTauDecays.TauRefit.AdvancedRefitVertexProducer_cfi")

process.AdvancedRefitVertexBSProducer.srcLeptons = cms.VInputTag("softLeptons")
process.p *= (process.AdvancedRefitVertexBS)


process.AdvancedRefitVertexNoBSProducer.srcLeptons = cms.VInputTag("softLeptons")
process.p *= (process.AdvancedRefitVertexNoBS)

SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) ) 


process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))
#process.TFileService=cms.Service('TFileService',fileName=cms.string('refFiles/Mu16_sync.root'))

if DO_ENRICHED:
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('Enriched_miniAOD.root'),
        outputCommands = cms.untracked.vstring('keep *'),
        fastCloning     = cms.untracked.bool(False),
        #Compression settings from MiniAOD allowing to save about 10% of disc space compared to defults ->
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        overrideInputFileSplitLevels = cms.untracked.bool(True)

        # <-
    )
    if STORE_ENRICHEMENT_ONLY:
        # Store only additional collections compared to MiniAOD necessary to reproduce ntuples (basically MVAMET, lepton pairs with SVFit and corrected jets)
        # Size of about 10% of full EnrichedMiniAOD
        process.out.outputCommands.append('drop *')
        process.out.outputCommands.append('keep *_SVllCand_*_*')
        process.out.outputCommands.append('keep *_SVbypass_*_*')
        process.out.outputCommands.append('keep *_barellCand_*_*')
        process.out.outputCommands.append('keep *_corrMVAMET_*_*')
        process.out.outputCommands.append('keep *_MVAMET_*_*')
        process.out.outputCommands.append('keep *_jets_*_*')
        process.out.outputCommands.append('keep *_patJetsReapplyJEC_*_*')
        process.out.outputCommands.append('keep *_softLeptons_*_*')
        process.out.outputCommands.append('keep *_genInfo_*_*')
        #process.out.fileName = 'EnrichementForMiniAOD.root' #FIXME: change name of output file?
    process.end = cms.EndPath(process.out)


#process.p = cms.Path(process.Candidates)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

#processDumpFile = open('process.dump' , 'w')
#print >> processDumpFile, process.dumpPython()

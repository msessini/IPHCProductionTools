import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

#set this cut in the cfg file
try: APPLYELECORR
except NameError:
    APPLYELECORR="None"
ELECORRTYPE=APPLYELECORR
try: IsMC
except NameError:
    IsMC=True
try: IsEmbed
except NameError:
    IsEmbed=True
try: YEAR
except NameError:
    YEAR = 2018
try: PERIOD
except:
    PERIOD ="A"
print 'Year+Period:', str(YEAR)+PERIOD
try: doCPVariables
except NameError:
    doCPVariables=True       
try: LEPTON_SETUP
except NameError:
    LEPTON_SETUP=2012
try: APPLYFSR
except NameError:
    APPLYFSR=False
try: BUILDONLYOS
except NameError:
    BUILDONLYOS=False
try: Is25ns
except NameError:
    Is25ns=True

try: USE_NOHFMET
except NameError:
    USE_NOHFMET=False

#PFMetName = "slimmedMETs"
PFMetName = "slimmedMETsModifiedPuppiMET"
if USE_NOHFMET: PFMetName = "slimmedMETsNoHF"

try: APPLYMETCORR
except NameError:
    APPLYMETCORR=True

try: HLTProcessName
except NameError:
    HLTProcessName='HLT'


### ----------------------------------------------------------------------
### Trigger list
### ----------------------------------------------------------------------
if YEAR == 2016:
  print 'Using HLT trigger 2016'
  execfile(PyFilePath+"python/triggers_80X.py")  # 2016 triggers and filters
if YEAR == 2017:
  print 'Using HLT trigger 2017'
  execfile(PyFilePath+"python/triggers_92X.py")  # 2017 triggers and filters
if YEAR == 2018:
  print 'Using HLT trigger 2018'
  execfile(PyFilePath+"python/triggers_102X.py") # 2018 triggers and filters

if IsEmbed:
    if YEAR == 2016:
        print 'Using HLT trigger 2016'
        execfile(PyFilePath+"python/triggers_80X_Embed.py")  # 2016 triggers and filters
    if YEAR == 2017:
        print 'Using HLT trigger 2017'
        execfile(PyFilePath+"python/triggers_92X_Embed.py")  # 2017 triggers and filters
    if YEAR == 2018:
        print 'Using HLT trigger 2018'
        execfile(PyFilePath+"python/triggers_102X_Embed.py") # 2018 triggers and filters

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")    
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

if IsMC:
  if YEAR == 2016:
    process.GlobalTag.globaltag = '102X_mcRun2_asymptotic_v7'        # 2016 MC
  if YEAR == 2017:
    process.GlobalTag.globaltag = '102X_mc2017_realistic_v7'        # 2017 MC
  if YEAR == 2018:
    process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'  # 2018 MC
else :
  if YEAR == 2016:
    process.GlobalTag.globaltag = '102X_dataRun2_v12'                # 2016 Data
  if YEAR == 2017:
    process.GlobalTag.globaltag = '102X_dataRun2_v12'                # 2017 Data
  if YEAR == 2018:
    if PERIOD=="D":
        process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v15'    # 2018D Data
    else:
        process.GlobalTag.globaltag = '102X_dataRun2_v12'           # 2018ABC Data

print "GT: ",process.GlobalTag.globaltag

nanosec="25"
if not Is25ns: nanosec="50"

METfiltersProcess = "PAT" if IsMC else "RECO" # NB! this is not guaranteed to be true! the following is valid on 2015 Run C + Run D data. Check:
# NB: for MET filters, use PAT or RECO depending if the miniAOD was generated simultaneously with RECO or in a separated step
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

### ----------------------------------------------------------------------
### Counters 
### ----------------------------------------------------------------------
process.nEventsTotal = cms.EDProducer("EventCountProducer")       # don't change producer name
process.nEventsPassTrigger = cms.EDProducer("EventCountProducer") # these names are then "hard-coded" inside the ntuplizer plugin

### ----------------------------------------------------------------------
### Trigger bit Requests - filter 
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

if IsEmbed:
    process.hltFilter = hlt.hltHighLevel.clone(
        TriggerResultsTag = cms.InputTag("TriggerResults","","SIMembedding"),
        HLTPaths = TRIGGERLIST,
        andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
        throw = cms.bool(False) #if True: throws exception if a trigger path is invalid  
        )
else:
    process.hltFilter = hlt.hltHighLevel.clone(
        TriggerResultsTag = cms.InputTag("TriggerResults","",HLTProcessName),
        HLTPaths = TRIGGERLIST,
        andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
        throw = cms.bool(False) #if True: throws exception if a trigger path is invalid  
        )

### ----------------------------------------------------------------------
### L1ECALPrefiringWeightRecipe (for 2016 and 2017 MC only)
### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
### ----------------------------------------------------------------------
prefireEra = "2016BtoH"
if YEAR==2017: prefireEra = "2017BtoF"

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
  DataEra = cms.string(prefireEra),
  UseJetEMPt = cms.bool(False),
  PrefiringRateSystematicUncty = cms.double(0.2),
  SkipWarnings = False
)

# Trigger Unpacker Module
#process.patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
#                                            patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
#                                            triggerResults = cms.InputTag("TriggerResults","",HLTProcessName),
#                                            unpackFilterLabels = cms.bool(True)
#)

#MC stuff

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )


process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string(PVERTEXCUT),
  filter = cms.bool(False), # if True, rejects events . if False, produce emtpy vtx collection
)


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string(MUCUT),
)


process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
    vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    miniRelIsoChgCollection = cms.InputTag("isoForMu:miniIsoChg"), #nanoAOD
    miniRelIsoAllCollection = cms.InputTag("isoForMu:miniIsoAll"), #nanoAOD
    ptRatioCollection = cms.InputTag("ptRatioRelForMu:ptRatio"), #nanoAOD
    ptRelCollection = cms.InputTag("ptRatioRelForMu:ptRel"), #nanoAOD               
    jetNDauChargedMVASelCollection = cms.InputTag("ptRatioRelForMu:jetNDauChargedMVASel"), #nanoAOD
    sampleType = cms.int32(LEPTON_SETUP),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    cut = cms.string(""),
    flags = cms.PSet(
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isGood = cms.string(MUCUT)
    )
)

if IsEmbed:process.softMuons.genCollection=cms.InputTag("prunedGenParticles","","MERGE")
else: process.softMuons.genCollection=cms.InputTag("prunedGenParticles")

process.muons =  cms.Sequence(process.bareSoftMuons+ process.softMuons)


###
###
###


process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    
from RecoMET.METProducers.PFMET_cfi import pfMet
    
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(
      #  process,
     #   isData= (not IsMC),
    #    isEmbeddedSample=IsEmbed,
   #     fixEE2017 = bool(YEAR=='2017'),
  #      fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
 #       postfix="ModifiedMET",
#        )
    
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True );
runMetCorAndUncFromMiniAOD(process,
                           isData= (not IsMC),
                           #ACisEmbeddedSample=IsEmbed,
                           metType="Puppi",
                           fixEE2017 = bool(YEAR=='2017'),
                           fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
                           postfix="ModifiedPuppiMET",
                           jetFlavor="AK4PFPuppi",
                           )




###
### Electrons
###


# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")


#**********************
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#**********************

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if YEAR==2016:
  EgammaPostRecoSeq_ERA = '2016-Legacy'      # 2016 data
  rerunIDs = True
  rerunEnergyCorrections = False
if YEAR==2017:
  EgammaPostRecoSeq_ERA = '2017-Nov17ReReco' # 2017 data
  rerunIDs = True
  rerunEnergyCorrections = True
if YEAR == 2018:
  EgammaPostRecoSeq_ERA = '2018-Prompt'      # 2018 data
  rerunIDs = True
  rerunEnergyCorrections = True
setupEgammaPostRecoSeq(process,
                       runVID=rerunIDs,
                       runEnergyCorrections=rerunEnergyCorrections,
                       era=EgammaPostRecoSeq_ERA)
		       
process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("slimmedElectrons"),
   rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
   vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   #genCollection = cms.InputTag("prunedGenParticles"),
   sampleType = cms.int32(LEPTON_SETUP),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string(ELECUT)
   )
   
if IsEmbed:process.softElectrons.genCollection=cms.InputTag("prunedGenParticles","","MERGE")
else:process.softElectrons.genCollection=cms.InputTag("prunedGenParticles")

process.electrons = cms.Sequence(process.softElectrons+process.egammaPostRecoSeq)


### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("(isGlobalMuon || userFloat('isPFMuon'))"), #
           deltaR              = cms.double(0.05),  
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)



##
## Taus
##

# Davide first update for 2018 May 2019
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
                    updatedTauName = updatedTauName,
                    toKeep = ["MVADM_2016_v1","MVADM_2017_v1","deepTau2017v2p1","2017v2", 
                               ])

tauIdEmbedder.runTauID()

# old sequence starts here
process.bareTaus = cms.EDFilter("PATTauRefSelector",
   src = cms.InputTag("slimmedTausNewID"), 
   cut = cms.string(TAUCUT),
   )

##NOT USED FOR NOW, TBD Later
process.cleanTaus = cms.EDProducer("PATTauCleaner",
    src = cms.InputTag("bareTaus"),
    # preselection (any string-based cut on pat::Tau)
    preselection = cms.string(
            'tauID("decayModeFinding") > 0.5 &'
            ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
            ' tauID("againstMuonTight") > 0.5 &'
            ' tauID("againstElectronMedium") > 0.5'
        ),
    
   # overlap checking configurables
   checkOverlaps = cms.PSet(
      muons = cms.PSet(
          src       = cms.InputTag("cleanPatMuons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      electrons = cms.PSet(
          src       = cms.InputTag("cleanPatElectrons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      ),
        # finalCut (any string-based cut on pat::Tau)
        finalCut = cms.string(' '),
)

# TES corrections: https://indico.cern.ch/event/887196/contributions/3743090/attachments/1984772/3306737/TauPOG_TES_20200210.pdf

# EES corrections: https://indico.cern.ch/event/868279/contributions/3665970/attachments/1959265/3267731/FES_9Dec_explained.pdf

# NominalTESCorrection=-1#in percent\
APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data

# 2016 data - MVAoldDM2017v2
#NomTESUncDM0   = cms.double(1.0)  # in percent, up/down uncertainty of TES
#NomTESUncDM1   = cms.double(0.9)  # in percent, up/down uncertainty of TES
#NomTESUncDM10  = cms.double(1.1)  # in percent, up/down uncertainty of TES
#NomTESUncDM11  = --> Missing <--  # in percent, up/down uncertainty of TES
#NomTESCorDM0   = cms.double(-0.6) # DecayMode==0
#NomTESCorDM1   = cms.double(-0.5) # DecayMode==1
#NomTESCorDM10  = cms.double(0.0)  # DecayMode==10
#NomTESCorDM11  = --> Missing <--  # DecayMode==11

# 2017 data - MVAoldDM2017v2
#if YEAR == 2017:
#    NomTESUncDM0   = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM1   = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM10  = cms.double(0.9)  # in percent, up/down uncertainty of TES
#    NomTESUncDM11  = cms.double(1.0)  # in percent, up/down uncertainty of TES
#    NomTESCorDM0   = cms.double(0.7)  # DecayMode==0
#    NomTESCorDM1   = cms.double(-0.2) # DecayMode==1
#    NomTESCorDM10  = cms.double(0.1)  # DecayMode==10
#    NomTESCorDM11  = cms.double(-0.1) # DecayMode==11

# 2018 data - MVAoldDM2017v2
#if YEAR == 2018:
#    NomTESUncDM0   = cms.double(1.1)  # in percent, up/down uncertainty of TES
#    NomTESUncDM1   = cms.double(0.9)  # in percent, up/down uncertainty of TES
#    NomTESUncDM10  = cms.double(0.8)  # in percent, up/down uncertainty of TES
#    NomTESUncDM11  = --> Missing <--  # in percent, up/down uncertainty of TES
#    NomTESCorDM0   = cms.double(-1.3) # DecayMode==0
#    NomTESCorDM1   = cms.double(-0.5) # DecayMode==1
#    NomTESCorDM10  = cms.double(-1.2) # DecayMode==10
#    NomTESCorDM11  = --> Missing <--  # DecayMode==11

# 2016 data - DeepTau2017v2p1
NomTESUncDM0      = cms.double(0.8)  # in percent, up/down uncertainty of TES
NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
NomTESUncDM10     = cms.double(0.8)  # in percent, up/down uncertainty of TES
NomTESUncDM11     = cms.double(1.1)  # in percent, up/down uncertainty of TES
NomTESCorDM0      = cms.double(-0.9) # DecayMode==0
NomTESCorDM1      = cms.double(-0.1) # DecayMode==1
NomTESCorDM10     = cms.double(0.3)  # DecayMode==10
NomTESCorDM11     = cms.double(-0.2) # DecayMode==11

#EES BARREL
NomEFakeESCorDM0B     = cms.double(0.679) #DecayMode==0
NomEFakeESUncDM0BUp    = cms.double(0.806) #DecayMode==0
NomEFakeESUncDM0BDown  = cms.double(0.982) #DecayMode==0
NomEFakeESCorDM1B      = cms.double(3.389) #DecayMode==1
NomEFakeESUncDM1BUp    = cms.double(1.168) #DecayMode==1
NomEFakeESUncDM1BDown  = cms.double(2.475) #DecayMode==1
#EES ENDCAP
NomEFakeESCorDM0E      = cms.double(-3.5)   #DecayMode==0
NomEFakeESUncDM0EUp    = cms.double(1.808)  #DecayMode==0
NomEFakeESUncDM0EDown  = cms.double(1.102)  #DecayMode==0
NomEFakeESCorDM1E      = cms.double(5.)      #DecayMode==1
NomEFakeESUncDM1EUp    = cms.double(6.57)   #DecayMode==1
NomEFakeESUncDM1EDown  = cms.double(5.694)  #DecayMode==1

TESyear = "2016Legacy"

# 2017 data - DeepTau2017v2p1
if YEAR == 2017:
    NomTESUncDM0      = cms.double(1.0)  # in percent, up/down uncertainty of TES
    NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
    NomTESUncDM10     = cms.double(0.7)  # in percent, up/down uncertainty of TES
    NomTESUncDM11     = cms.double(1.4)  # in percent, up/down uncertainty of TES
    NomTESCorDM0      = cms.double(0.4)  # DecayMode==0
    NomTESCorDM1      = cms.double(0.2)  # DecayMode==1
    NomTESCorDM10     = cms.double(0.1)  # DecayMode==10
    NomTESCorDM11     = cms.double(-1.3) # DecayMode==1

    #EES BARREL
    NomEFakeESCorDM0B      = cms.double(0.911) #DecayMode==0
    NomEFakeESUncDM0BUp    = cms.double(1.343) #DecayMode==0
    NomEFakeESUncDM0BDown  = cms.double(0.882) #DecayMode==0
    NomEFakeESCorDM1B      = cms.double(1.154) #DecayMode==1
    NomEFakeESUncDM1BUp    = cms.double(2.162) #DecayMode==1
    NomEFakeESUncDM1BDown  = cms.double(0.973) #DecayMode==1
    #EES ENDCAP
    NomEFakeESCorDM0E      = cms.double(-2.604)   #DecayMode==0
    NomEFakeESUncDM0EUp    = cms.double(2.249)    #DecayMode==0
    NomEFakeESUncDM0EDown  = cms.double(1.43)     #DecayMode==0
    NomEFakeESCorDM1E      = cms.double(1.5)    #DecayMode==1
    NomEFakeESUncDM1EUp    = cms.double(6.461)      #DecayMode==1
    NomEFakeESUncDM1EDown  = cms.double(4.969)    #DecayMode==1

    TESyear = "2017ReReco"

# 2018 data - DeepTau2017v2p1
if YEAR == 2018:
    NomTESUncDM0          = cms.double(0.9)  # in percent, up/down uncertainty of TES
    NomTESUncDM1          = cms.double(0.5)  # in percent, up/down uncertainty of TES
    NomTESUncDM10         = cms.double(0.7)  # in percent, up/down uncertainty of TES
    NomTESUncDM11         = cms.double(1.2)  # in percent, up/down uncertainty of TES
    NomTESCorDM0          = cms.double(-1.6) # DecayMode==0
    NomTESCorDM1          = cms.double(-0.5) # DecayMode==1
    NomTESCorDM10         = cms.double(-1.2) # DecayMode==10
    NomTESCorDM11         = cms.double(-0.4) # DecayMode==11

    #EES BARREL
    NomEFakeESCorDM0B      = cms.double(1.362)    #DecayMode==0
    NomEFakeESUncDM0BUp    = cms.double(0.904)    #DecayMode==0
    NomEFakeESUncDM0BDown  = cms.double(0.474)    #DecayMode==0
    NomEFakeESCorDM1B      = cms.double(1.954)    #DecayMode==1
    NomEFakeESUncDM1BUp    = cms.double(1.226)    #DecayMode==1
    NomEFakeESUncDM1BDown  = cms.double(1.598)    #DecayMode==1
    #EES ENDCAP
    NomEFakeESCorDM0E      = cms.double(-3.097)   #DecayMode==0
    NomEFakeESUncDM0EUp    = cms.double(3.404)    #DecayMode==0
    NomEFakeESUncDM0EDown  = cms.double(1.25)     #DecayMode==0
    NomEFakeESCorDM1E      = cms.double(-1.5)     #DecayMode==1
    NomEFakeESUncDM1EUp    = cms.double(5.499)    #DecayMode==1
    NomEFakeESUncDM1EDown  = cms.double(4.309)    #DecayMode==1

    TESyear = "2018ReReco"

process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   PFCollection = cms.InputTag("packedPFCandidates"),
   offlinebeamSpot = cms.InputTag("offlineBeamSpot"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string(TAUDISCRIMINATOR),

   NominalTESUncertaintyDM0         = NomTESUncDM0,
   NominalTESUncertaintyDM1         = NomTESUncDM1,
   NominalTESUncertaintyDM10        = NomTESUncDM10,
   NominalTESUncertaintyDM11        = NomTESUncDM11,
   NominalTESCorrectionDM0          = NomTESCorDM0,
   NominalTESCorrectionDM1          = NomTESCorDM1,
   NominalTESCorrectionDM10         = NomTESCorDM10,
   NominalTESCorrectionDM11         = NomTESCorDM11,

   NominalEFakeESCorrectionDM0B      = NomEFakeESCorDM0B,
   NominalEFakeESUncertaintyDM0BUp   = NomEFakeESUncDM0BUp, 
   NominalEFakeESUncertaintyDM0BDown = NomEFakeESUncDM0BDown, 
   NominalEFakeESCorrectionDM1B      = NomEFakeESCorDM1B,
   NominalEFakeESUncertaintyDM1BUp   = NomEFakeESUncDM1BUp, 
   NominalEFakeESUncertaintyDM1BDown = NomEFakeESUncDM1BDown, 
   NominalEFakeESCorrectionDM0E      = NomEFakeESCorDM0E,
   NominalEFakeESUncertaintyDM0EUp   = NomEFakeESUncDM0EUp, 
   NominalEFakeESUncertaintyDM0EDown = NomEFakeESUncDM0EDown, 
   NominalEFakeESCorrectionDM1E      = NomEFakeESCorDM1E,
   NominalEFakeESUncertaintyDM1EUp   = NomEFakeESUncDM1EUp, 
   NominalEFakeESUncertaintyDM1EDown = NomEFakeESUncDM1EDown, 

   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   # ApplyTESUpDown = cms.bool(True if IsMC else False), # no shift computation when data
   flags = cms.PSet(
        isGood = cms.string("")
        ),

   year = cms.string(TESyear)
   )
   
if IsEmbed:process.softTaus.genCollection = cms.InputTag("prunedGenParticles","","MERGE")
else:process.softTaus.genCollection = cms.InputTag("prunedGenParticles")

process.taus=cms.Sequence(process.rerunMvaIsolationSequence + process.slimmedTausNewID + process.bareTaus + process.softTaus)

### ----------------------------------------------------------------------
### gen info, only from MC
### ----------------------------------------------------------------------
process.genInfo = cms.EDProducer("GenFiller",
         src = cms.InputTag("prunedGenParticles"),
         storeLightFlavAndGlu = cms.bool(True) # if True, store also udcs and gluons (first copy)
 )

if IsEmbed:process.genInfo.src = cms.InputTag("prunedGenParticles","","MERGE")
else:process.genInfo.src = cms.InputTag("prunedGenParticles")
                  
if IsMC : process.geninfo = cms.Sequence(process.genInfo)
else : process.geninfo = cms.Sequence()


### ----------------------------------------------------------------------
### Search for FSR candidates
### ----------------------------------------------------------------------
process.load("UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff")
process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),#cms.InputTag("cmgPhotonSel"),
    matchFSR = cms.bool(True)
    )

process.fsrSequence = cms.Sequence(process.fsrPhotonSequence + process.appendPhotons)
muString = "appendPhotons:muons"
eleString = "appendPhotons:electrons"
if not APPLYFSR : 
    process.fsrSequence = cms.Sequence()
    muString = "softMuons"
    eleString = "softElectrons"
    tauString = "softTaus"
#Leptons
process.softLeptons = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
)


#
#Jets
#

# apply new jet energy corrections and recompute btaggers
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

# JEC corrections
jecLevels = None

#if IsMC:
#    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
#else:
#    jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]

jecLevels = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]

# bTaggers
#btagVector = []
#
#if YEAR == 2018:
#    btagVector.append('None')
#
#if YEAR == 2017:
#    btagVector2017 = [
#        'pfDeepFlavourJetTags:probb',
#        'pfDeepFlavourJetTags:probbb',
#        'pfDeepFlavourJetTags:problepb',
#        'pfDeepFlavourJetTags:probc',
#        'pfDeepFlavourJetTags:probuds',
#        'pfDeepFlavourJetTags:probg'
#    ]
#    btagVector.extend(btagVector2017)
#
#if YEAR == 2016:
#    btagVector2016 = [
#        'pfDeepFlavourJetTags:probb',
#        'pfDeepFlavourJetTags:probbb',
#        'pfDeepFlavourJetTags:problepb',
#        'pfDeepFlavourJetTags:probc',
#        'pfDeepFlavourJetTags:probuds',
#        'pfDeepFlavourJetTags:probg',
#        'pfDeepCSVJetTags:probudsg',
#        'pfDeepCSVJetTags:probb',
#        'pfDeepCSVJetTags:probc',
#        'pfDeepCSVJetTags:probbb'
#    ]
#    btagVector.extend(btagVector2016)

# Update jet collection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),
   #btagDiscriminators = btagVector,
   labelName = 'UpdatedJEC'
)

# Update the jet sequences
#if YEAR == 2016:
#    process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC *
#                                       process.updatedPatJetsUpdatedJEC *
#                                       process.patJetCorrFactorsTransientCorrectedUpdatedJEC *
#                                       process.pfImpactParameterTagInfosUpdatedJEC *
#                                      process.pfInclusiveSecondaryVertexFinderTagInfosUpdatedJEC *
#                                      process.pfDeepCSVTagInfosUpdatedJEC *
#                                      process.pfDeepCSVJetTagsUpdatedJEC *
#                                       process.pfDeepFlavourTagInfosUpdatedJEC *
#                                       process.pfDeepFlavourJetTagsUpdatedJEC *
#                                       process.updatedPatJetsTransientCorrectedUpdatedJEC *
#                                      process.selectedUpdatedPatJetsUpdatedJEC)

#if YEAR == 2017:
#    process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC *
#                                       process.updatedPatJetsUpdatedJEC *
#                                       process.patJetCorrFactorsTransientCorrectedUpdatedJEC *
#                                       process.pfImpactParameterTagInfosUpdatedJEC *
#                                       process.pfInclusiveSecondaryVertexFinderTagInfosUpdatedJEC *
#                                       process.pfDeepCSVTagInfosUpdatedJEC *
#                                       process.pfDeepFlavourTagInfosUpdatedJEC *
#                                       process.pfDeepFlavourJetTagsUpdatedJEC *
#                                       process.updatedPatJetsTransientCorrectedUpdatedJEC *
#                                       process.selectedUpdatedPatJetsUpdatedJEC)

#if YEAR == 2018:
#    process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC *
#                                       process.updatedPatJetsUpdatedJEC *
#                                       process.selectedUpdatedPatJetsUpdatedJEC)

# Jet Selector after JEC and bTagging
#process.jets = cms.EDFilter("PATJetRefSelector",
#                            src = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
#                            cut = cms.string(JETCUT),
#)


 # ----------------------
 # Jet energy corrections
 # ----------------------
process.ak4PFL1FastjetCHS = cms.EDProducer("L1FastjetCorrectorProducer",
#    srcRho = cms.InputTag("kt6PFJets", "rho"),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet')
)
process.ak4PFL2RelativeCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)
process.ak4PFL3AbsoluteCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)
process.ak4PFResidualCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)

#Corrections applied to miniaod slimmedJets
pfJECS = cms.PSet(
  L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
  L2Relative = cms.string("ak4PFL2RelativeCHS"),
  L3Absolute = cms.string("ak4PFL3AbsoluteCHS")
)
if not IsMC: pfJECS = cms.PSet(
  L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
  L2Relative = cms.string("ak4PFL2RelativeCHS"),
  L3Absolute = cms.string("ak4PFL3AbsoluteCHS"),
  L2L3Residual = cms.string("ak4PFResidualCHS")
)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

# JER
# --------
process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")

process.slimmedJetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag("updatedPatJetsUpdatedJEC"),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK4PFchs'),
        algopt = cms.string('AK4PFchs_pt'),
        #resolutionFile = cms.FileInPath('Autumn18_V7_MC_PtResolution_AK4PFchs.txt'),
        #scaleFactorFile = cms.FileInPath('combined_SFs_uncertSources.txt'),

        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),

        debug = cms.untracked.bool(False),
        # Systematic variation
        # 0: Nominal
        # -1: -1 sigma (down variation)
        # 1: +1 sigma (up variation)
        variation = cms.int32(0),  # If not specified, default to 0
)

process.slimmedJetsSmearedDown = process.slimmedJetsSmeared.clone(variation=cms.int32(-1))
process.slimmedJetsSmearedUp = process.slimmedJetsSmeared.clone(variation=cms.int32(1))

# Jet Selector after JEC, JER and bTagging

if IsMC:
    process.jets = cms.EDFilter("PATJetRefSelector",
                                #src = cms.InputTag("slimmedJets"),
                                src = cms.InputTag("slimmedJetsSmeared"),
                                cut = cms.string(JETCUT)
                                )
else:
    process.jets = cms.EDFilter("PATJetRefSelector",
                                #src = cms.InputTag("slimmedJets"),
                                src = cms.InputTag("updatedPatJetsUpdatedJEC"),
                                cut = cms.string(JETCUT)
                                )




##
## QG tagging for jets
##

if COMPUTEQGVAR:

    QGlikelihood_tag = 'QGLikelihoodObject_v1_AK4PFchs'
    if YEAR == 2017 or YEAR == 2018:
      QGlikelihood_tag = 'QGLikelihoodObject_v1_AK4PFchs_2017'

    from CondCore.CondDB.CondDB_cfi import CondDB
     
    process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDB.clone(
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
      ),
      toGet = cms.VPSet(
        cms.PSet(
          record = cms.string('QGLikelihoodRcd'),
          tag    = cms.string(QGlikelihood_tag),
          label  = cms.untracked.string('QGL_AK4PFchs'),
        ),
      ),
    )
    process.es_prefer_qg = cms.ESPrefer("PoolDBESSource", "QGPoolDBESSource")

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag("slimmedJetsSmeared")    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
    process.jetSequence = cms.Sequence(process.slimmedJetsSmeared *  process.slimmedJetsSmearedDown * process.slimmedJetsSmearedUp * process.jets*process.QGTagger)

else:
    process.jetSequence = cms.Sequence(process.slimmedJetsSmeared *  process.slimmedJetsSmearedDown * process.slimmedJetsSmearedUp * process.jets)


# Add latest pileup jet ID
process.load("RecoJets.JetProducers.PileupJetID_cfi")
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_81x, _chsalgos_94x, _chsalgos_102x

if YEAR == 2016:
    PUalgo = _chsalgos_81x
if YEAR == 2017:
    PUalgo = _chsalgos_94x
if YEAR == 2018:
    PUalgo = _chsalgos_102x

process.pileupJetIdUpdated = process.pileupJetId.clone(
   jets = cms.InputTag("jets"),
   inputIsCorrected = True,
   applyJec = False,
   vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
   algos = PUalgo
)
process.jetSequence += cms.Sequence(process.pileupJetIdUpdated)


# il primo legge la collezione dei leptoni e stampa quali sono
#process.beforeLLcombiner = cms.EDFilter("beforeCombiner",
#    src = cms.InputTag("softLeptons")
#)

##
## Build ll candidates (here OS)
##
decayString="softLeptons softLeptons"
checkcharge=False
if BUILDONLYOS:
    decayString="softLeptons@+ softLeptons@-"
    checkcharge=True
process.barellCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string(decayString),
                                    cut = cms.string(LLCUT),
                                    checkCharge = cms.bool(checkcharge)
)

#il seconod legge le pairs e stampa quali sono e da chi sono composti
#process.afterLLcombiner = cms.EDFilter("afterCombiner",
#    srcPairs = cms.InputTag("barellCand")
#)

## ----------------------------------------------------------------------
## MVA MET
## ----------------------------------------------------------------------

process.METSequence = cms.Sequence()
if USEPAIRMET:
    print "Using pair MET (MVA MET)"
    from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
    runMVAMET(process, jetCollectionPF = "patJetsReapplyJEC")
    process.MVAMET.srcLeptons = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
    process.MVAMET.requireOS = cms.bool(False)
    process.MVAMET.permuteLeptonsWithinPlugin = cms.bool(False)
    process.MVAMET.leptonPermutations = cms.InputTag("barellCand")

    process.MVAMETInputs = cms.Sequence(
        process.slimmedElectronsTight + process.slimmedMuonsTight + process.slimmedTausLoose + process.slimmedTausLooseCleaned + process.patJetsReapplyJECCleaned +
        process.pfCHS + process.pfChargedPV + process.pfChargedPU + process.pfNeutrals + process.neutralInJets +
        process.pfMETCands + process.pfTrackMETCands + process.pfNoPUMETCands + process.pfPUCorrectedMETCands + process.pfPUMETCands +
        process.pfChargedPUMETCands + process.pfNeutralPUMETCands + process.pfNeutralPVMETCands + process.pfNeutralUnclusteredMETCands +
        process.pfChs +
        process.ak4PFCHSL1FastjetCorrector + process.ak4PFCHSL2RelativeCorrector + process.ak4PFCHSL3AbsoluteCorrector + process.ak4PFCHSResidualCorrector +
        process.ak4PFCHSL1FastL2L3Corrector + process.ak4PFCHSL1FastL2L3ResidualCorrector +
        process.tauDecayProducts + process.tauPFMET + process.tauMET + process.tausSignificance
    )
    for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET", "pfChargedPUMET", "pfNeutralPUMET", "pfNeutralPVMET", "pfNeutralUnclusteredMET"]:
        process.MVAMETInputs += getattr(process, met)
        process.MVAMETInputs += getattr(process, "ak4JetsFor"+met)
        process.MVAMETInputs += getattr(process, "corr"+met)
        process.MVAMETInputs += getattr(process, met+"T1")
        process.MVAMETInputs += getattr(process, "pat"+met)
        process.MVAMETInputs += getattr(process, "pat"+met+"T1")        

    process.METSequence += cms.Sequence(process.MVAMETInputs + process.MVAMET)


else:
    print "Using event pfMET (same MET for all pairs)"

    # patch to get a standalone MET significance collection
    
 #   process.METSignificance = cms.EDProducer ("ExtractMETSignificance",
#                                              srcMET=cms.InputTag("slimmedMETsModifiedMET","","TEST")
  #                                            )
    process.PUPPIMETSignificance = cms.EDProducer ("ExtractPUPPIMETSignificance",
                                              srcPUPPIMET=cms.InputTag("slimmedMETsModifiedPuppiMET")
                                              )
    # add variables with MET shifted for TES corrections
    process.ShiftMETforTES = cms.EDProducer ("ShiftMETforTES",
                                             srcMET  = cms.InputTag("slimmedMETsModifiedPuppiMET"),
                                             tauCollection = cms.InputTag("softTaus")
                                             )
    #process.ShiftMETforTES = cms.EDProducer ("ShiftMETforTES",
     #                                        srcMET  = cms.InputTag("slimmedMETsPuppi"),
      #                                      tauCollection = cms.InputTag("softTaus")
       #                                      )
    # add variables with MET shifted for EES corrections (E->tau ES)
    #process.ShiftMETforEES = cms.EDProducer ("ShiftMETforEES",
     #                                        srcMET  = cms.InputTag("slimmedMETsModifiedMET","","TEST"),
      #                                       tauCollection = cms.InputTag("softTaus")
       #                                      )
    process.ShiftMETforEES = cms.EDProducer ("ShiftMETforEES",
                                             srcMET  = cms.InputTag("slimmedMETsModifiedPuppiMET"),
                                             tauCollection = cms.InputTag("softTaus")
                                             )
    process.METSequence += process.egmPhotonIDSequence
    process.METSequence += process.puppiMETSequence
    process.METSequence += process.fullPatMetSequenceModifiedPuppiMET
   # process.METSequence += process.fullPatMetSequenceModifiedMET
    #process.METSequence += process.METSignificance
    process.METSequence += process.PUPPIMETSignificance
    process.METSequence += process.ShiftMETforTES
    process.METSequence += process.ShiftMETforEES
    
####

    # 2017 and 2018 ECAL bad calibration filter to be rerun, fix from:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

    # In 2017 and 2018 some problematic crystals --> pass list of crystals
    if YEAR == 2017 or YEAR == 2018:
        baddetEcallist = cms.vuint32(
            [872439604,872422825,872420274,872423218,
            872423215,872416066,872435036,872439336,
            872420273,872436907,872420147,872439731,
            872436657,872420397,872439732,872439339,
            872439603,872422436,872439861,872437051,
            872437052,872420649,872422436,872421950,
            872437185,872422564,872421566,872421695,
            872421955,872421567,872437184,872421951,
            872421694,872437056,872437057,872437313])

    # In 2016 no problem --> pass empty list
    if YEAR == 2016:
        baddetEcallist = cms.vuint32([])

    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
        "EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt        = cms.double(50.),
        baddetEcal    = baddetEcallist,
        taggingMode = cms.bool(True),
        debug = cms.bool(False)
        )


## ----------------------------------------------------------------------
## Z-recoil correction
## ----------------------------------------------------------------------

# corrMVAPairMET = []
if IsMC and APPLYMETCORR:
    if USEPAIRMET:
        process.selJetsForZrecoilCorrection = cms.EDFilter("PATJetSelector",
            src = cms.InputTag("jets"),                                      
            cut = cms.string("pt > 30. & abs(eta) < 4.7"), 
            filter = cms.bool(False)
        )
        process.corrMVAMET = cms.EDProducer("ZrecoilCorrectionProducer",                                                   
            srcPairs = cms.InputTag("barellCand"),
            srcMEt = cms.InputTag("MVAMET", "MVAMET"),
            srcGenParticles = cms.InputTag("prunedGenParticles"),
            srcJets = cms.InputTag("selJetsForZrecoilCorrection"),
            correction = cms.string("HTT-utilities/RecoilCorrections/data/MvaMET_MG_2016BCD.root")
        )
        process.METSequence += process.selJetsForZrecoilCorrection        
        process.METSequence += process.corrMVAMET

    else:
        raise ValueError("Z-recoil corrections for PFMET not implemented yet !!")


srcMETTag = None
if USEPAIRMET:
  srcMETTag = cms.InputTag("corrMVAMET") if (IsMC and APPLYMETCORR) else cms.InputTag("MVAMET", "MVAMET")
else:
  # MET corrected for central TES and EES shifts of the taus
  #srcMETTag = cms.InputTag("ShiftMETcentral")
  srcMETTag = cms.InputTag(PFMetName,"","TEST")

## ----------------------------------------------------------------------
## SV fit
## ----------------------------------------------------------------------
#if USECLASSICSVFIT:
#    print "Using CLASSIC_SV_FIT"
process.SVllCand = cms.EDProducer("ClassicSVfitInterface",
                                  srcPairs   = cms.InputTag("barellCand"),
                                  srcSig     = cms.InputTag("PUPPIMETSignificance", "PUPPIMETSignificance"),
                                  srcCov     = cms.InputTag("PUPPIMETSignificance", "PUPPIMETCovariance"),
                                  usePairMET = cms.bool(USEPAIRMET),
                                  srcMET     = srcMETTag,
                                  computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False),
                                  computeForUpDownMET = cms.bool(COMPUTEMETUPDOWNSVFIT if IsMC else False),
                                  METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                  METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                  METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                  METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
                                  METdxUP_EES   = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
                                  METdyUP_EES   = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
                                  METdxDOWN_EES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
                                  METdyDOWN_EES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES")
)
#else:
#    print "Using STANDALONE_SV_FIT"
#    process.SVllCand = cms.EDProducer("SVfitInterface",
#                                      srcPairs   = cms.InputTag("barellCand"),
#                                      srcSig     = cms.InputTag("METSignificance", "METSignificance"),
#                                      srcCov     = cms.InputTag("METSignificance", "METCovariance"),
#                                      usePairMET = cms.bool(USEPAIRMET),
#                                      srcMET     = srcMETTag,
#                                      computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT if IsMC else False)
#    )

## ----------------------------------------------------------------------
## SV fit BYPASS (skip SVfit, don't compute SVfit pair mass)
## ----------------------------------------------------------------------
process.SVbypass = cms.EDProducer ("SVfitBypass",
                                    srcPairs   = cms.InputTag("barellCand"),
                                    usePairMET = cms.bool(USEPAIRMET),
                                    srcMET     = srcMETTag,
                                    srcSig     = cms.InputTag("PUPPIMETSignificance", "PUPPIMETSignificance"),
                                    srcCov     = cms.InputTag("PUPPIMETSignificance", "PUPPIMETCovariance"),
                                    METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                    METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                    METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                    METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
                                    METdxUP_EES   = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
                                    METdyUP_EES   = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
                                    METdxDOWN_EES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
                                    METdyDOWN_EES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES")
)


## ----------------------------------------------------------------------
## Ntuplizer
## ----------------------------------------------------------------------
process.HTauTauTree = cms.EDAnalyzer("HTauTauNtuplizer",
                      fileName = cms.untracked.string ("CosaACaso"),
                      applyFSR = cms.bool(APPLYFSR),
                      IsMC = cms.bool(IsMC),
                      IsEmbed = cms.bool(IsEmbed),
                      do_MCSummary = cms.bool(IsMC),
                      do_MCComplete = cms.bool(IsMC),
                      year = cms.int32(YEAR),        
                      doCPVariables = cms.bool(doCPVariables),               
                      vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      secVtxCollection = cms.InputTag("slimmedSecondaryVertices"), # FRA
                      puCollection = cms.InputTag("slimmedAddPileupInfo"),
                      rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoMiniRelIsoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
                      rhoForJER = cms.InputTag("fixedGridRhoAll"), # FRA
                      PFCandCollection = cms.InputTag("packedPFCandidates"),
                      jetCollection = cms.InputTag("jets"),
                      #JECset = cms.untracked.string("patJetCorrFactors"),
                      JECset = cms.untracked.string("patJetCorrFactorsUpdatedJEC"),
                      SmearedJets = cms.InputTag("jets"),
                      SmearedJetsDown = cms.InputTag("slimmedJetsSmearedDown"),
                      SmearedJetsUp = cms.InputTag("slimmedJetsSmearedUp"),     
                      computeQGVar = cms.bool(COMPUTEQGVAR),
                      QGTagger = cms.InputTag("QGTagger", "qgLikelihood"),
                     # stage2TauCollection = cms.InputTag("caloStage2Digis","Tau","SIMembedding"),
                      #stage2JetCollection = cms.InputTag("caloStage2Digis","Jet","SIMembedding"),
                      ak8jetCollection = cms.InputTag("slimmedJetsAK8"),
                      lepCollection = cms.InputTag("softLeptons"),
                      lheCollection = cms.InputTag("LHEEventProduct"),
                      genCollection = cms.InputTag("generator"),
                      genericCollection = cms.InputTag("genInfo"),
                      genjetCollection = cms.InputTag("slimmedGenJets"),
                      totCollection = cms.InputTag("nEventsTotal"),
                      passCollection = cms.InputTag("nEventsPassTrigger"),
                      lhepCollection = cms.InputTag("externalLHEProducer"),
                     # triggerResultsLabel = cms.InputTag("TriggerResults", "", "SIMembedding"), #Different names for MiniAODv2 at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.                      
                      triggerSet = cms.InputTag("slimmedPatTrigger"),    # FRA
                      #triggerSet = cms.InputTag("selectedPatTrigger"),
                      triggerList = HLTLIST,
                      #metFilters = cms.InputTag ("TriggerResults","","MERGE"),
                      PUPPImetCollection = cms.InputTag("slimmedMETsModifiedPuppiMET"),
#                      srcPFMETCov = cms.InputTag("METSignificance", "METCovariance"),
 #                     srcPFMETSignificance = cms.InputTag("METSignificance", "METSignificance"),
                      srcPUPPIMETCov = cms.InputTag("PUPPIMETSignificance", "PUPPIMETCovariance"),
                      srcPUPPIMETSignificance = cms.InputTag("PUPPIMETSignificance", "PUPPIMETSignificance"),
                      HT = cms.InputTag("externalLHEProducer"),
                      beamSpot = cms.InputTag("offlineBeamSpot"),
                      #PrunedGenCollection = cms.InputTag("prunedGenParticles","","MERGE"),
                      #nBadMu = cms.InputTag("removeBadAndCloneGlobalMuons"),
                      genLumiHeaderTag = cms.InputTag("generator"),
                      #metERCollection = cms.InputTag("slimmedMETsModifiedMET"),
                     # metERCollection = cms.InputTag("slimmedMETsModifiedMET","","TEST"),        
                      ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),
                      L1prefireProb     = cms.InputTag("prefiringweight:nonPrefiringProb"),
                      L1prefireProbUp   = cms.InputTag("prefiringweight:nonPrefiringProbUp"),
                      L1prefireProbDown = cms.InputTag("prefiringweight:nonPrefiringProbDown"),       
                      TausNewIDCollection = cms.InputTag("slimmedTausNewID"),
                      RefitVtxBSCollection = cms.InputTag("AdvancedRefitVertexBSProducer"),
                      RefitVtxNoBSCollection = cms.InputTag("AdvancedRefitVertexNoBSProducer")
 
)


#if USE_NOHFMET:
#    process.HTauTauTree.metCollection = cms.InputTag("slimmedMETsNoHF")
#else:
    # MET corrected for central TES and EES shifts of the taus
#    process.HTauTauTree.metCollection = srcMETTag

#if YEAR == 2016 or YEAR == 2017:
#    process.HTauTauTree.JECset = cms.untracked.string("patJetCorrFactorsTransientCorrectedUpdatedJEC")
#if YEAR == 2018:
#    process.HTauTauTree.JECset = cms.untracked.string("patJetCorrFactorsUpdatedJEC")

if IsEmbed:
    process.HTauTauTree.stage2TauCollection = cms.InputTag("caloStage2Digis","Tau","SIMembedding")
    process.HTauTauTree.stage2JetCollection = cms.InputTag("caloStage2Digis","Jet","SIMembedding")
    process.HTauTauTree.triggerResultsLabel = cms.InputTag("TriggerResults", "", "SIMembedding")
    process.HTauTauTree.metFilters = cms.InputTag ("TriggerResults","","MERGE")
    process.HTauTauTree.PrunedGenCollection = cms.InputTag("prunedGenParticles","","MERGE")
    #L1prefireProb     = cms.InputTag("prefiringweight:nonPrefiringProb"),
    #L1prefireProbUp   = cms.InputTag("prefiringweight:nonPrefiringProbUp"),
    #L1prefireProbDown = cms.InputTag("prefiringweight:nonPrefiringProbDown"),     
else:
    process.HTauTauTree.stage2TauCollection = cms.InputTag("caloStage2Digis","Tau")
    process.HTauTauTree.stage2JetCollection = cms.InputTag("caloStage2Digis","Jet")
    process.HTauTauTree.triggerResultsLabel = cms.InputTag("TriggerResults", "", HLTProcessName)
    process.HTauTauTree.metFilters = cms.InputTag ("TriggerResults","",METfiltersProcess)
    process.HTauTauTree.PrunedGenCollection = cms.InputTag("prunedGenParticles")
    #L1prefireProb     = cms.InputTag("prefiringweight:nonPrefiringProb"),
    #L1prefireProbUp   = cms.InputTag("prefiringweight:nonPrefiringProbUp"),
    #L1prefireProbDown = cms.InputTag("prefiringweight:nonPrefiringProbDown"),     


if SVFITBYPASS:
    process.HTauTauTree.candCollection = cms.InputTag("SVbypass")
    process.SVFit = cms.Sequence (process.SVbypass)

else:
    process.HTauTauTree.candCollection = cms.InputTag("SVllCand")
    process.SVFit = cms.Sequence (process.SVllCand)


#print particles gen level - DEBUG purposes
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(10),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("prunedGenParticles")
)

##
## Paths
##
process.PVfilter = cms.Path(process.goodPrimaryVertices)
process.ecalBadCalib = cms.Path(process.ecalBadCalibReducedMINIAODFilter)
process.l1ECALPref = cms.Path(process.prefiringweight)

# Prepare lepton collections
process.Candidates = cms.Sequence(
    # process.printTree         + # just for debug, print MC particles
    process.nEventsTotal       +
    #process.hltFilter         + 
    process.nEventsPassTrigger +
    process.egammaPostRecoSeq  +
    process.muons              +
    process.electrons          + process.cleanSoftElectrons +
    process.taus               +
    process.fsrSequence        +
    process.softLeptons        + process.barellCand +
    process.jecSequence        + process.jetSequence +
    process.METSequence        +
    process.geninfo            +
    process.SVFit
    )
# always run ntuplizer
process.trees = cms.EndPath(process.HTauTauTree)


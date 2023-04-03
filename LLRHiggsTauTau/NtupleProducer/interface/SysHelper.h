#ifndef SysHelper_h
#define SysHelper_h

#include <string>
#include <TTree.h>
#include <TVector3.h>
#include <TVector2.h>
#include <RooWorkspace.h>
#include <boost/functional/hash.hpp>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <TauPolSoftware/SimpleFits/interface/PTObject.h>
#include <CondTools/BTau/interface/BTagCalibrationReader.h>
#include <CondFormats/BTauObjects/interface/BTagCalibration.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/FastMTT.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CPTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DataMCType.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/SelectionTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CorrectionTools.h>
#include <HiggsCPinTauDecays/ImpactParameter/interface/IpCorrection.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>


typedef std::map<std::string, std::unique_ptr<JetCorrectionUncertainty>> myJECMap;

class SysHelper{
 public:
  SysHelper(Int_t theYear, std::string dataMCstring);
  ~SysHelper();
  void ResetVariables();
  void MakeBranches(TTree *tree, bool isNominal);
  void GetCollections(const edm::View<pat::CompositeCandidate>* cands_, const edm::View<reco::Candidate>* daus_, const edm::View<pat::Jet> *jets_, const edm::View<pat::Jet> *jetsUp_, const edm::View<pat::Jet> *jetsDown_);
  void GetGenInfo(edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag_, std::map<std::string, double> theomap, Int_t dataMCtype);
  void FillGenTaus(std::vector<std::vector<unsigned int>> signal_Tauidx, std::vector<std::vector<std::vector<double>>> tauandprod_p4, std::vector<std::vector<int>> tauandprod_charge, std::vector<std::vector<std::vector<double>>> tauandprod_vtx, std::vector<std::vector<int>> taudandprod_pdgid, std::vector<unsigned int> tau_JAK);
  void GetEventInfo(bool isEmbed, bool isMC, ULong64_t runNumber, Float_t nPU, ULong64_t evtidx, Int_t lumi);
  void GetDecayProducts(std::vector<LorentzVectorParticle> a1lvp, std::vector<TrackParticle> muontrack, std::vector<TrackParticle> tautrack, std::vector<std::vector<std::vector<double>>> pi_P4, std::vector<std::vector<double>> pi_charges);
  void GetPV(std::vector<float> _RefitPVBS_x, std::vector<float> _RefitPVBS_y, std::vector<float> _RefitPVBS_z, std::vector<std::vector<std::vector<double>>> _RefitPVBS_Cov, std::vector<size_t> _VertexHashBS1, std::vector<size_t> _VertexHashBS2, std::vector<size_t> _LeptonHash);
  bool FillPV(int tauIndex, int muIndex);
  void GetMETCov(float cov00, float cov10, float cov11);
  void GetJECUnc(std::map<std::string, std::vector<Float_t>> JECmapUp, std::map<std::string, std::vector<Float_t>> JECmapDown, myJECMap* JECmap, JetCorrectionUncertainty* JECunc);
  void GetPrefiringWeights(double pref, double prefUp, double prefDown);
  void GetTauSpinnerWeights(const double wEven, const double wOdd, const double wMM);
  void FillTree(TTree *tree, std::string sysType, std::string var, const edm::Event& event, bool trig, std::vector<Long64_t> _daughters_trgMatched, std::vector<math::XYZTLorentzVector> LeptonP4, const pat::MET &PUPPImet_, int Npartons);

 private:
  //Tau
  int _tauIndex;
  int _tauGenMatch;
  float _tauDM;
  float _tauIPx;
  float _tauIPy;
  float _tauIPz;
  float _tauIPsignificance;
  float _tauPt;
  float _tauEta;
  float _tauPhi;
  float _tauE;
  float _tauSVx;
  float _tauSVy;
  float _tauSVz;
  float _GEFtauE;
  float _GEFtauPt;
  float _GEFtauPhi;
  float _GEFtauEta;
  //Muon
  int _muIndex;
  int _muGenMatch;
  float _muPt;
  float _muEta;
  float _muPhi;
  float _muE;
  double _muIso;
  float _muIPx;
  float _muIPy;
  float _muIPz;
  float _muIPsignificance;
  //Gen
  math::XYZTLorentzVector genMuonP4;
  math::XYZTLorentzVector genTauP4;
  float _genTaupx;
  float _genTaupy;
  float _genTaupz;
  float _genTauE;
  float _genTauSVx;
  float _genTauSVy;
  float _genTauSVz;
  float _genMuonpx;
  float _genMuonpy;
  float _genMuonpz;
  float _genMuonE;
  float _genMuonIPx;
  float _genMuonIPy;
  float _genMuonIPz;
  float _genPVx;
  float _genPVy;
  float _genPVz;
  double _gendpPhiCP;
  double _genpvPhiCP;
  //Triggers
  bool _trgIsoMu;
  bool _extraIsoMu;
  bool _trgXMuTau;
  //Jets, Pair and MET
  bool _isOSpair;
  bool _isIso;
  bool _isMediumID;
  bool _isTightJetID;
  double _pairvisMass;
  int _Njets;
  int _Nbtags;
  double _leadingjetPt;
  double _leadingjetEta;
  double _leadingjetPhi;
  double _leadingjetE;
  double _subleadingjetPt;
  double _subleadingjetEta;
  double _subleadingjetPhi;
  double _subleadingjetE;
  double _subsubleadingjetPt;
  double _subsubleadingjetEta;
  double _subsubleadingjetPhi;
  double _subsubleadingjetE;
  double _dijetMass;
  double _dijetPt;
  double _dijetdeltaEta;
  double _ditauPt;
  double _fastMTTmass;
  float _muMETmt;
  float _PUPPImet;
  float _PUPPImetphi;
  float _PUPPIMETCov00;
  float _PUPPIMETCov10;
  float _PUPPIMETCov11;
  double _pvPhiCP;
  double _dpPhiCP;
  //Weights
  double _wEven;
  double _wOdd;
  double _wMM;
  double _wPrefiring;
  double _wPrefiringUp;
  double _wPrefiringDown;
  double _wIDvsJet;
  double _wIDvsJetUp;
  double _wIDvsJetDown;
  double _wIDvsEle;
  double _wIDvsEleUp;
  double _wIDvsEleDown;
  double _wIDvsMu;
  double _wIDvsMuUp;
  double _wIDvsMuDown;
  double _wTrg;
  double _wTrgUp;
  double _wTrgDown;
  double _wIDMu;
  double _wTrkMu;
  double _wPU;
  double _wZpT;
  double _wZpTUp;
  double _wZpTDown;
  double _wToppT;
  double _wToppTUp;
  double _wToppTDown;
  double _wBtag;
  double _wBtagUp;
  double _wBtagDown;
  double _wPSISRUp;
  double _wPSISRDown;
  double _wPSFSRUp;
  double _wPSFSRDown;
  double _wScaleUp;
  double _wScaleDown;
  double _wMC;
  double _wSignal;
  double _wEmb;
  double _wEmbUp;
  double _wEmbDown;
  double _wTot;
  //Maps
  std::map<std::string, double> _TauIDSFmap;
  std::map<std::string, double> _MuonIDTrkSFmap;
  std::map<std::string, double> _TriggerSFmap;
  std::map<std::string, double> _ZpTreweightingmap;
  std::map<std::string, double> _ToppTreweightingmap;
  std::map<std::string, double> _TheoreticalUncmap;
  std::map<std::string, double> _EmbeddingSFmap;
  //A1 LVP & MUON TRACK & TAU TRACK
  std::vector<LorentzVectorParticle> A1LVP;
  std::vector<TrackParticle> MuonTrack;
  std::vector<TrackParticle> TauTrack;
  //Refit Pions
  std::vector<std::vector<std::vector<double>>> RefitPionsP4;
  std::vector<std::vector<double>> RefitPionsCharge;
  //PV Refit BS
  float _pvx;
  float _pvy;
  float _pvz;
  float _pvCov00;
  float _pvCov11;
  float _pvCov22;
  float _pvCov01;
  float _pvCov02;
  float _pvCov12;
  std::vector<std::vector<double>> pvcov;
  //Event
  std::string _Idstr;
  Int_t _Id;
  int _Npartons;
  bool _isEmbed;
  bool _isData;
  bool _isMC;
  Int_t _theYear;
  ULong64_t _runNumber;
  ULong64_t _evt;
  Int_t _lumi;
  float _nPU;
  //MC type
  bool _isZ;
  bool _isW;
  bool _isSignal;
  bool _isQCD;
  bool _isVV;
  bool _isTTbar;
  bool _isSingleTop;
  //Data
  std::shared_ptr<RooWorkspace> _w;
  std::string _Label;
  TH1* _histTES;
  TGraph* _histFES;
  RecoilCorrector* _recoilPuppiMetCorrector;
  MEtSys* _recoilPuppiMetShifter;
  IpCorrection* _IPcorrector;
  TH2D* _histbtagEfficiency;
  const BTagCalibration* _btagCalib;
  BTagCalibrationReader* _btagReader;
  TH1D* _PU_data;
  TH1D* _PU_mc;
  //
  std::vector<float> _RefitPVBS_x;
  std::vector<float> _RefitPVBS_y;
  std::vector<float> _RefitPVBS_z;
  std::vector<std::vector<std::vector<double>>> _RefitPVBS_Cov;
  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_x;
  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_y;
  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_z;
  std::vector<size_t> _VertexHashBS1;
  std::vector<size_t> _VertexHashBS2;
  std::vector<size_t> _LeptonHash;
  const edm::View<pat::CompositeCandidate>* cands;
  const edm::View<reco::Candidate>* daus;
  const edm::View<pat::Jet> *jets;
  const edm::View<pat::Jet> *jetsUp;
  const edm::View<pat::Jet> *jetsDown;
  std::map<std::string, std::vector<Float_t>> JECmapUp;
  std::map<std::string, std::vector<Float_t>> JECmapDown;
  myJECMap* JECmap;
  JetCorrectionUncertainty* JECunc;
  edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag;
  pat::CompositeCandidate SelectedPair;
  //
  bool SelectPair(std::string sysType, std::string var, const edm::Event& event, std::vector<math::XYZTLorentzVector> LeptonP4, bool trig, std::vector<Long64_t> _daughters_trgMatched);
  std::vector<math::XYZTLorentzVector> SelectJets(const reco::Candidate* cand1, const reco::Candidate* cand2, const edm::View<pat::Jet>* jets, int pTcut, std::string sysType);
  std::pair<std::vector<math::XYZTLorentzVector>, TVector2> CorrectedJetsandMET(std::string sysType, std::string var, const pat::MET& PUPPImet, const reco::Candidate* cand1, const reco::Candidate* cand2, myJECMap* JECmap);
  void GetSampleType();
};
#endif

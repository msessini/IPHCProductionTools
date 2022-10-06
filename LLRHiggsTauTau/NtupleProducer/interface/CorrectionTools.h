#ifndef CorrectionTools_h
#define CorrectionTools_h

#include <cmath>
#include <TMVA/Reader.h>
#include <RooFunctor.h>
#include <RooWorkspace.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <TauPOG/TauIDSFs/interface/TauIDSFTool.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/PileUp.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DataMCType.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/SelectionTools.h>
#include <CondTools/BTau/interface/BTagCalibrationReader.h>
#include <CondFormats/BTauObjects/interface/BTagCalibration.h>
#include <HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h>
#include <HTT-utilities/RecoilCorrections/interface/MEtSys.h>
#include <HiggsCPinTauDecays/ImpactParameter/interface/IpCorrection.h>

namespace corrector{
  void IPCorrection(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, TVector3 IP, double Eta, IpCorrection ipCorrector);
  void METRecoilCorrection(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, const pat::MET& PUPPImet, int Njets, float shiftedMETx, float shiftedMETy, std::string sysType, std::string var, RecoilCorrector* recoilPuppiMetCorrector, MEtSys* recoilPuppiMetShifter);
  math::XYZTLorentzVector P4Corrected(math::XYZTLorentzVector p4, int genmatch, int DM, std::string Unc, TH1* histTES, TGraph* histFES);
  math::XYZTLorentzVector MuP4Corrected(math::XYZTLorentzVector p4, std::string Unc);
};

namespace weight{
  std::map<std::string, double> ZpTreweighting(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, std::string sysType, std::shared_ptr<RooWorkspace> w);
  std::map<std::string, double> ToppTreweighting(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, std::string sysType);
  std::map<std::string, double> PrefiringWeight(double pref, double prefUp, double prefDown, Int_t theYear, std::string sysType);
  std::map<std::string, double> TauIDSF(int genmatch, float DM, math::XYZTLorentzVector p4, bool isEmbed, std::string sysType, std::shared_ptr<RooWorkspace> w, std::string Label);
  std::map<std::string, double> MuonIDTrkSF(math::XYZTLorentzVector p4, double m_iso, bool isEmbed, std::shared_ptr<RooWorkspace> w);
  std::map<std::string, double> TriggerSF(int genmatch, float DM, math::XYZTLorentzVector taup4, double m_iso, math::XYZTLorentzVector mup4, Int_t theYear, bool isEmbed, std::string sysType, std::shared_ptr<RooWorkspace> w);
  std::map<std::string, double> BTaggingSF(std::vector<const pat::Jet*> selectedJets, Int_t theYear, std::string sysType, TH2D* btagEfficiency, const BTagCalibration *calib, BTagCalibrationReader *reader);
  //
  double SignalReweighting(Int_t theYear, TString id);
  double PileUpreweighting(float nPU, TH1D* PU_data, TH1D* PU_mc);
  double Stitching(Int_t theYear, int id, int taugenmatch, int mugenmatch, double mvis, bool isZ, bool isW, int Npartons);
};

#endif

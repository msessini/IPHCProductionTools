#ifndef CPTools_h
#define CPTools_h

#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/TrackParticle.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LorentzVectorParticle.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CorrectionTools.h>
#include <TauPolSoftware/SimpleFits/interface/GlobalEventFit.h>
#include <HiggsCPinTauDecays/TauDecaysInterface/interface/SCalculator.h>
#include <HiggsCPinTauDecays/ImpactParameter/interface/IpCorrection.h>
#include <HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h>

class CPTools : public SCalculator{
 public:
  CPTools(std::vector<std::vector<double>> pionsP4, std::vector<double> pionsCharge, math::XYZTLorentzVector muonp4, TVector3 muonref);
  ~CPTools(){};
  //
  PTObject setMET(float MET_px, float MET_py, TMatrixD MET_cov);
  void initGEF(LorentzVectorParticle a1, TrackParticle muon, TVector3 PV, std::vector<std::vector<double>> pvcov, PTObject MET);
  float getTau(std::string var);
  void setImpactParameter();
  void correctIP(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, double muEta, IpCorrection ipcorrector);
  TVector3 getImpactParameter();
  float getIPsignificance(std::vector<std::vector<double>> pvcov);
  double getPhiCPwithDP();
  double getPhiCPwithPV();
  double getPhiCPwithPV(TLorentzVector hadTauP4);

 private:
  std::string _hadDecay;
  std::vector<TLorentzVector> _pionsP4;
  std::vector<double> _pionsCharge;
  TLorentzVector _hadronicTauP4;
  TLorentzVector _muonP4;
  TVector3 _PV;
  TVector3 _muonRef;
  TVector3 _muonIP;

};

#endif

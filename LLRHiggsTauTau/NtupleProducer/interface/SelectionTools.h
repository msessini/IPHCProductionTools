#ifndef SelectionTools_h
#define SelectionTools_h

#include <cmath>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/TauDecay_CMSSW.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/PDGInfo.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>

namespace seltools{

  bool EleVeto(const reco::Candidate* cand);
  bool MuVeto(const reco::Candidate* cand);
  bool DiEle(const reco::Candidate* cand1, const reco::Candidate* cand2);
  bool DiMuon(const reco::Candidate* cand1, const reco::Candidate* cand2);
  int GenMatch(const reco::Candidate* lep, const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag);
  bool CompareLegs(const reco::Candidate *i, const reco::Candidate *j);
  bool CHECK_BIT(unsigned long long var, int pos);
  bool ComparePairsbyIso(const pat::CompositeCandidate i, const pat::CompositeCandidate j);
  float ComputeMT(math::XYZTLorentzVector VisP4, float METx, float METy);

};
#endif
  

#include <LLRHiggsTauTau/NtupleProducer/interface/SelectionTools.h>

bool seltools::EleVeto(const reco::Candidate* cand)
{
  if(cand->isElectron()){
    if(cand->pt()>10 && std::abs(cand->eta())<2.5){
      if(std::abs(userdatahelpers::getUserFloat(cand,"dxy")) < 0.045 && std::abs(userdatahelpers::getUserFloat(cand,"dz")) < 0.2){
        if(userdatahelpers::getUserFloat(cand,"isEleNoIsoID90") == 1){
          if(userdatahelpers::getUserInt(cand,"isConversionVeto") == 1){
            if(userdatahelpers::getUserInt(cand,"missingHit")<=1){
              if(userdatahelpers::getUserFloat(cand,"combRelIsoPF")<(0.3*cand->pt())){
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}

bool seltools::MuVeto(const reco::Candidate* cand)
{
  if(cand->isMuon()){
    if(cand->pt()>10 && std::abs(cand->eta())<2.4){
      if(std::abs(userdatahelpers::getUserFloat(cand,"dxy")) < 0.045 && std::abs(userdatahelpers::getUserFloat(cand,"dz")) < 0.2){
        if(seltools::CHECK_BIT(userdatahelpers::getUserInt(cand,"muonID"),2)){
          if(userdatahelpers::getUserFloat(cand,"combRelIsoPF")<(0.3*cand->pt())){
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool seltools::DiEle(const reco::Candidate* cand1, const reco::Candidate* cand2)
{
  bool kin=false, vertex=false, isele=false, iso=false;

  if(cand1->isElectron() && cand2->isElectron() &&  ((cand1->charge()/abs(cand1->charge())) != (cand2->charge()/abs(cand2->charge())))){
    {
      if(deltaR(cand1->p4(),cand2->p4())>0.15)
        {
          kin = (cand1->pt()>15 && cand2->pt()>15 && fabs(cand1->eta())<2.5 && fabs(cand2->eta())<2.5);
          vertex = (fabs(userdatahelpers::getUserFloat(cand1,"dxy"))<0.045 && fabs(userdatahelpers::getUserFloat(cand2,"dxy"))< 0.045 && fabs(userdatahelpers::getUserFloat(cand1,"dz")) < 0.2 && fabs(userdatahelpers::getUserFloat(cand2,"dz")) < 0.2);
          isele=(userdatahelpers::getUserFloat(cand1,"isEleNoIsoID90") && userdatahelpers::getUserFloat(cand2,"isEleNoIsoID90"));
          iso= (userdatahelpers::getUserFloat(cand1,"combRelIsoPF") < (0.3*cand1->pt()) && userdatahelpers::getUserFloat(cand2,"combRelIsoPF")<(0.3*cand2->pt()));
          if((kin && vertex && isele && iso)==1)return true;
        }
    }
  }
  return false;
}

bool seltools::DiMuon(const reco::Candidate* cand1, const reco::Candidate* cand2)
{
  bool kin=false, vertex=false, isele=false, iso=false;

  if(cand1->isMuon() && cand2->isMuon() &&  ((cand1->charge()/abs(cand1->charge())) != (cand2->charge()/abs(cand2->charge())))){
    {
      if(deltaR(cand1->p4(),cand2->p4())>0.15)
        {
          kin = (cand1->pt()>15 && cand2->pt()>15 && fabs(cand1->eta())<2.4 && fabs(cand2->eta())<2.4);
          vertex = (fabs(userdatahelpers::getUserFloat(cand1,"dxy"))<0.045 && fabs(userdatahelpers::getUserFloat(cand2,"dxy"))< 0.045 && fabs(userdatahelpers::getUserFloat(cand1,"dz")) < 0.2 && fabs(userdatahelpers::getUserFloat(cand2,"dz")) < 0.2);
          isele=((userdatahelpers::getUserFloat(cand1,"isPFMuon") && userdatahelpers::getUserFloat(cand2,"isPFMuon")) && (userdatahelpers::getUserFloat(cand1,"isGlobalMuon") && userdatahelpers::getUserFloat(cand2,"isGlobalMuon")) && (userdatahelpers::getUserFloat(cand1,"isTrackerMuon") && userdatahelpers::getUserFloat(cand2,"isTrackerMuon")));
          iso= (userdatahelpers::getUserFloat(cand1,"combRelIsoPF") < (0.3*cand1->pt()) && userdatahelpers::getUserFloat(cand2,"combRelIsoPF")<(0.3*cand2->pt()));
          if((kin && vertex && isele && iso)==1)return true;
        }
    }
  }
  return false;
}

bool seltools::CHECK_BIT(unsigned long long var, int pos){
  unsigned long long CHECK1=(unsigned long long)(var & ((unsigned long long)(1) << pos));
  unsigned long long CHECK2=(unsigned long long)((unsigned long long)(1) << pos);
  return ( CHECK1==CHECK2 );
}

int seltools::GenMatch(const reco::Candidate* lep, const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag)
{

  int GenMatch = 6;
  float minDR = 0.2;

  edm::Handle<edm::View<pat::GenericParticle> > candHandle;
  event.getByToken(generictag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();
  //
  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen) {
    unsigned pdgid = abs(igen->pdgId());
    unsigned flag = igen->userInt ("generalGenFlags");
    //
    if (deltaR(igen->p4(),lep->p4())<minDR) {
      minDR = deltaR(igen->p4(),lep->p4());

      bool type1 = pdgid==11 && seltools::CHECK_BIT(flag,0) && igen->pt()>8;
      bool type2 = pdgid==13 && seltools::CHECK_BIT(flag,0) && igen->pt()>8;
      bool type3 = pdgid==11 && seltools::CHECK_BIT(flag,5) && igen->pt()>8;
      bool type4 = pdgid==13 && seltools::CHECK_BIT(flag,5) && igen->pt()>8;
      bool type5 = pdgid==15 && seltools::CHECK_BIT(flag,0) && igen->pt()>15;

      if (type1) GenMatch = 1;
      else if (type2) GenMatch = 2;
      else if (type3) GenMatch = 3;
      else if (type4) GenMatch = 4;
      else if (type5) GenMatch = 5;
    }
  }
  return GenMatch;
}

bool seltools::CompareLegs(const reco::Candidate *i, const reco::Candidate *j){
  int iType=2,jType=2;

  if(i->isElectron())iType=1;
  else if(i->isMuon())iType=0;

  if(j->isElectron())jType=1;
  else if(j->isMuon())jType=0;

  if(iType>jType) return false;
  else if(iType==jType && i->pt()<j->pt()) return false;

  return true;
}


bool seltools::ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j){

  //Second criteria: ISO
  float isoi=999,isoj=999;
  int cand1j=-1,cand1i=-1;

  if(seltools::CompareLegs(i.daughter(0),i.daughter(1)))cand1i=0;
  else cand1i=1;
  if(seltools::CompareLegs(j.daughter(0),j.daughter(1)))cand1j=0;
  else cand1j=1;

  //step 1, leg 1 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(cand1j),"combRelIsoPF");
  if (!i.daughter(cand1i)->isMuon() && !i.daughter(cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(cand1i),"byDeepTau2017v2p1VSjetraw");
  if (!j.daughter(cand1j)->isMuon() && !j.daughter(cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(cand1j),"byDeepTau2017v2p1VSjetraw");
  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 2, leg 1 Pt
  if(i.daughter(cand1i)->pt()>j.daughter(cand1j)->pt()) return true;
  else if(i.daughter(cand1i)->pt()<j.daughter(cand1j)->pt()) return false;

  //step 3, leg 2 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"combRelIsoPF");
  if (!i.daughter(1-cand1i)->isMuon() && !i.daughter(1-cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(1-cand1i),"byDeepTau2017v2p1VSjetraw");
  if (!j.daughter(1-cand1j)->isMuon() && !j.daughter(1-cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(1-cand1j),"byDeepTau2017v2p1VSjetraw");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 4, leg 2 Pt
  if(i.daughter(1-cand1i)->pt()>j.daughter(1-cand1j)->pt()) return true;

  return false;
}

float seltools::ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy)
{
  math::XYZTLorentzVector METP4 (METx, METy, 0, 0); // I only care about transverse plane
  double dphi = deltaPhi(visP4.phi(), METP4.phi());
  return sqrt(2.*visP4.Pt()*METP4.Pt()*(1.-cos(dphi)));
}




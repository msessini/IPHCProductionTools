#include <LLRHiggsTauTau/NtupleProducer/interface/CPTools.h>

CPTools::CPTools(std::vector<std::vector<double>> pionsp4, std::vector<double> pionscharge, math::XYZTLorentzVector muonp4, TVector3 muonref) {

  int piSize = pionsp4.size();
  //
  if(piSize == 3) _hadDecay = "a1";
  else if(piSize == 2) _hadDecay = "rho";
  else if(piSize == 1) _hadDecay = "pion";
  else _hadDecay = "none";
  //
  for(std::vector<std::vector<double>>::iterator iPion=pionsp4.begin(); iPion!=pionsp4.end(); ++iPion) {
    _pionsP4.push_back(TLorentzVector(iPion->at(1), iPion->at(2), iPion->at(3), iPion->at(0)));
  }
  //
  _pionsCharge = pionscharge;
  _muonP4 = TLorentzVector(muonp4.px(), muonp4.py(), muonp4.pz(), muonp4.energy());
  _muonRef = muonref;
  //
}

void CPTools::correctIP(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, double muEta, IpCorrection ipcorrector) {
  corrector::IPCorrection(event, generictag, _muonIP, muEta, ipcorrector);
}

PTObject CPTools::setMET(float MET_px, float MET_py, TMatrixD MET_cov) {

  TMatrixTSym<double> covMET(2);
  covMET[0][0] = MET_cov[0][0];
  covMET[0][1] = MET_cov[0][1];
  covMET[1][0] = MET_cov[1][0];
  covMET[1][1] = MET_cov[0][0];

  TMatrixT<double> MET(2,1);
  MET[0][0] = MET_px;
  MET[1][0] = MET_py;

  PTObject METinput(MET,covMET);
  return METinput;
}

void CPTools::initGEF(LorentzVectorParticle a1, TrackParticle muon, TVector3 pv, std::vector<std::vector<double>> pvcov, PTObject METinput) {

  TMatrixTSym<double> _PVCov(3);
  for(int j=0; j<3; j++){
    for(int k=j; k<3; k++){
      _PVCov(j,k) = pvcov.at(j).at(k);
      _PVCov(k,j) = pvcov.at(j).at(k);
    }
  }
  _PV = pv;
  //
  GlobalEventFit GEF(muon,a1,METinput,_PV,_PVCov);
  GEF.setMassConstraint(125.3);
  //GEF.setMinimizer("Minuit2"); // only in more recent version of GEF
  GEF.SetCorrectPt(false);
  GEFObject FitTaus = GEF.Fit();
  if(FitTaus.Fitconverged()) {
    _hadronicTauP4 = FitTaus.getTauH().LV();
  }
  else _hadronicTauP4 = TLorentzVector(0., 0., 0., 0.);
}

float CPTools::getTau(std::string var) {
  if(var=="E") return _hadronicTauP4.E();
  else if(var=="Pt") return _hadronicTauP4.Pt();
  else if(var=="Phi") return _hadronicTauP4.Phi();
  else if(var=="Eta") return _hadronicTauP4.Eta();
  else return 0;
}

void CPTools::setImpactParameter() {

  TVector3 k, p, IP;
  k.SetXYZ(_muonRef.x() - _PV.x(), _muonRef.y() - _PV.y(), _muonRef.z() - _PV.z());
  p.SetXYZ(_muonP4.Px(), _muonP4.Py(), _muonP4.Pz());
  if (p.Mag() != 0) IP = k - (p.Dot(k) / p.Mag2()) * p;
  else IP.SetXYZ(-999, -999, -999);
  _muonIP = IP; 
}

float CPTools::getIPsignificance(std::vector<std::vector<double>> pvcov) {

  TVector3 n = _muonIP.Unit();
  ROOT::Math::SVector<double, 3> Sn(n.x(),n.y(),n.z());
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> cov_PV;
  for(int j=0; j<3; j++){
    for(int k=j; k<3; k++){
      cov_PV(j,k) = pvcov.at(j).at(k);
      cov_PV(k,j) = pvcov.at(j).at(k);
    }
  }
  double alpha = TMath::Sqrt((ROOT::Math::Dot(Sn,cov_PV*Sn)));
  double sig = _muonIP.Mag()/alpha;
  return sig;
}

TVector3 CPTools::getImpactParameter() {
  return _muonIP;
}

double CPTools::getPhiCPwithDP() {
  return AcopAngle_DPIP(_hadDecay, "muon", _pionsP4, _muonP4, _muonIP);
}

double CPTools::getPhiCPwithPV() {
  return AcopAngle_PVIP(_hadDecay, "muon", _hadronicTauP4, _pionsP4, _pionsCharge, _muonP4, _muonIP);
}

double CPTools::getPhiCPwithPV(TLorentzVector hadTauP4) {
  return AcopAngle_PVIP(_hadDecay, "muon", hadTauP4, _pionsP4, _pionsCharge, _muonP4, _muonIP);
}

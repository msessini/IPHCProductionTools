#include <LLRHiggsTauTau/NtupleProducer/interface/SysHelper.h>

SysHelper::SysHelper(Int_t theYear, std::string dataMCstring)
{
  _theYear = theYear;
  _Idstr = dataMCstring;
  //Tau
  _tauIndex = -99;
  _tauGenMatch = -99;
  _tauPt = -99;
  _tauEta = -99;
  _tauPhi = -99;
  _tauE = -99;
  _tauDM = -99;
  _tauIPx = -99;
  _tauIPy = -99;
  _tauIPz = -99;
  _tauIPsignificance = -99;
  _tauSVx = -99;
  _tauSVy = -99;
  _tauSVz = -99;
  _GEFtauE = -99;
  _GEFtauPt = -99;
  _GEFtauPhi = -99;
  _GEFtauEta = -99;
  //Muon
  _muIndex = -99;
  _muGenMatch = -99;
  _muPt = -99;
  _muEta = -99;
  _muPhi = -99;
  _muE = -99;
  _muIso = -99;
  _muIPx = -99;
  _muIPy = -99;
  _muIPz = -99;
  _muIPsignificance = -99;
  //Gen
  _genTaupx = -99;
  _genTaupy = -99;
  _genTaupz = -99;
  _genTauE = -99;
  _genTauSVx = -99;
  _genTauSVy = -99;
  _genTauSVz = -99;
  _genMuonpx = -99;
  _genMuonpy = -99;
  _genMuonpz = -99;
  _genMuonE = -99;
  _genMuonIPx = -99;
  _genMuonIPy = -99;
  _genMuonIPz = -99;
  _genPVx = -99;
  _genPVy = -99;
  _genPVz = -99;
  _gendpPhiCP = -99;
  _genpvPhiCP = -99;
  //Jets, Pair and MET
  _isOSpair = false;
  _isIso = false;
  _isMediumID = false;
  _pairvisMass = -99;
  _Njets = -99;
  _Nbjets = -99;
  _leadingjetPt = -99;
  _leadingjetEta = -99;
  _leadingjetPhi = -99;
  _leadingjetE = -99;
  _subleadingjetPt = -99;
  _subleadingjetEta = -99;
  _subleadingjetPhi = -99;
  _subleadingjetE = -99;
  _subsubleadingjetPt = -99;
  _subsubleadingjetEta = -99;
  _subsubleadingjetPhi = -99;
  _subsubleadingjetE = -99;
  _dijetMass = -99;
  _dijetPt = -99;
  _dijetdeltaEta = -99;
  _ditauPt = -99;
  _fastMTTmass = -99;
  _muMETmt = -99;
  _PUPPImet = -99;
  _PUPPImetphi = -99;
  _PUPPIMETCov00 = -99;
  _PUPPIMETCov10 = -99;
  _PUPPIMETCov11 = -99;
  _pvPhiCP = -99;
  _dpPhiCP = -99;
  //Weights
  _wEven = 1.;
  _wOdd = 1.;
  _wMM = 1.;
  _wPrefiring = 1.;
  _wPrefiringUp = 1.;
  _wPrefiringDown = 1.;
  _wIDvsJet = 1.;
  _wIDvsJetUp = 1.;
  _wIDvsJetDown = 1.;
  _wIDvsEle = 1.;
  _wIDvsEleUp = 1.;
  _wIDvsEleDown = 1.;
  _wIDvsMu = 1.;
  _wIDvsMuUp = 1.;
  _wIDvsMuDown = 1.;
  _wTrg = 1.;
  _wTrgUp = 1.;
  _wTrgDown = 1.;
  _wIDMu = 1.;
  _wTrkMu = 1.;
  _wPU = 1.;
  _wZpT = 1.;
  _wToppT = 1.;
  _wBtag = 1.;
  _wBtagUp = 1.;
  _wBtagDown = 1.;
  _wZpTUp = 1.;
  _wZpTDown = 1.;
  _wToppTUp = 1.;
  _wToppTDown = 1.;
  _wPSISRUp = 1.;
  _wPSISRDown = 1.;
  _wPSFSRUp = 1.;
  _wPSFSRDown = 1.;
  _wScaleUp = 1.;
  _wScaleDown = 1.;
  _wMC = 1.;
  _wSignal = 1.;
  //Maps
  _TauIDSFmap.clear();
  _MuonIDTrkSFmap.clear();
  _TriggerSFmap.clear();
  _ZpTreweightingmap.clear();
  _ToppTreweightingmap.clear();
  //A1 LVP & MUON TRACK
  A1LVP.clear();
  MuonTrack.clear();
  //Refit Pions
  RefitPionsP4.clear();
  RefitPionsCharge.clear();
  //PV 
  _pvx = -99;
  _pvy = -99;
  _pvz = -99;
  _pvCov00 = -99;
  _pvCov11 = -99;
  _pvCov22 = -99;
  _pvCov01 = -99;
  _pvCov02 = -99;
  _pvCov12 = -99;
  pvcov.clear();
  //Event
  _Npartons = -99;
  _isEmbed = false;
  _isData = false;
  _isMC = false;
  _nPU = -99;
  //Sample type
  _Id = -99;
  _isZ = false;
  _isW = false;
  _isSignal = false;
  _isQCD = false;
  _isVV = false;
  _isTTbar = false;
  _isSingleTop = false;
  //Data
  auto year = std::to_string(theYear);
  if(theYear == 2016) _Label = "2016Legacy";
  else _Label = year+"ReReco";
  //
  TFile *WorkSpace=TFile::Open(("$CMSSW_BASE/data/energy_scale/htt_scalefactors_legacy_"+year+".root").c_str(), "READ");
  _w = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  WorkSpace->Close();
  delete WorkSpace;
  //
  TFile *TES = new TFile(("$CMSSW_BASE/data/energy_scale/TauES_dm_DeepTau2017v2p1VSjet_"+_Label+".root").c_str(),"READ");
  _histTES = dynamic_cast<TH1*>((const_cast<TFile*>(TES))->Get("tes"));
  TES->Close();
  delete TES;
  //
  TFile *FES= new TFile(("$CMSSW_BASE/data/energy_scale/TauFES_eta-dm_DeepTau2017v2p1VSe_"+_Label+".root").c_str(),"READ");
  _histFES = dynamic_cast<TGraph*>((const_cast<TFile*>(FES))->Get("fes"));
  FES->Close();
  delete FES;
  //
  _recoilPuppiMetCorrector = new RecoilCorrector(("$CMSSW_BASE/data/recoil/Type1_PuppiMET_"+year+".root").c_str());
  _recoilPuppiMetShifter = new MEtSys(("$CMSSW_BASE/data/recoil/PuppiMETSys_"+year+".root").c_str());
  //
  TFile *btagFile = new TFile(("$CMSSW_BASE/data/b_tag/eta_pt_btagEff_"+year+".root").c_str(),"READ");
  _histbtagEfficiency = dynamic_cast<TH2D*>((const_cast<TFile*>(btagFile))->Get("eta_pt_eff"));
  btagFile->Close();
  delete btagFile;
  //
  if(theYear == 2016) {
    //_btagCalib = new const BTagCalibration("deepcsv", (std::string)getenv("CMSSW_BASE") + "/data/b_tag/DeepCSV_2016LegacySF_V1.csv"); 
    //PRODUIRE LES EFF MAP
  }
  else if(theYear == 2017) {
    //_btagCalib = new const BTagCalibration("deepcsv", (std::string)getenv("CMSSW_BASE") + "/data/b_tag/DeepCSV_94XSF_V5_B_F.csv"); 
    //PRODUIRE LES EFF MAP
  }
  else if(theYear == 2018) {
    _btagCalib = new const BTagCalibration("deepcsv", (std::string)getenv("CMSSW_BASE") + "/data/b_tag/DeepCSV_102XSF_V1.csv");
  }
  _btagReader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
  _btagReader->load(*_btagCalib, BTagEntry::FLAV_B, "comb");
  //
  _IPcorrector = new IpCorrection("$CMSSW_BASE/data/IP/ip_2018.root");
  //
  if(theYear == 2016) {
    TFile* filePUdistribution_data=TFile::Open("$CMSSW_BASE/data/pileup/Data_Pileup_2016_271036-284044_80bins.root", "READ");
    TFile* filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/MC_Moriond17_PU25ns_V1.root", "READ");
    //
    _PU_data=(TH1D*)filePUdistribution_data->Get("pileup");
    _PU_mc=(TH1D*)filePUdistribution_MC->Get("pileup");
    //
    filePUdistribution_data->Close();
    filePUdistribution_MC->Close();
  }
  if(theYear == 2017) {
    TFile* filePUdistribution_data=TFile::Open("$CMSSW_BASE/data/pileup/pu_distributions_data_2017.root", "READ");
    TFile* filePUdistribution_MC = nullptr;
    if(_Idstr=="dy_ll") filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_DYJetsToLL-LO.root", "READ");
    else if(_Idstr=="w_3qlnu") filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_W3JetsToLNu-LO.root", "READ");
    else if(_Idstr=="ww_1l1nu2q") filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_WWTo1L1Nu2Q.root", "READ");
    else if(_Idstr=="wminush_tautau") filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_WminusHToTauTau_M-125.root", "READ");
    else if(_Idstr=="wplush_tautau") filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_WplusHToTauTau_M-125.root", "READ");
    else filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileup_2017_DYJetsToLL-ext.root", "READ");
    //
    _PU_data=(TH1D*)filePUdistribution_data->Get("pileup");
    _PU_mc=(TH1D*)filePUdistribution_MC->Get("pileup");
    //
    filePUdistribution_data->Close();
    filePUdistribution_MC->Close();
  }
  if(theYear == 2018) {
    TFile* filePUdistribution_data=TFile::Open("$CMSSW_BASE/data/pileup/pileUp_data_Autumn18.root", "READ");
    TFile* filePUdistribution_MC=TFile::Open("$CMSSW_BASE/data/pileup/pileUp_MC_Autumn18.root", "READ");
    //
    _PU_data=(TH1D *)filePUdistribution_data->Get("pileup");
    _PU_mc=(TH1D *)filePUdistribution_MC->Get("pileup");
    //
    filePUdistribution_data->Close();
    filePUdistribution_MC->Close();
  }
  //
}

SysHelper::~SysHelper() 
{
  delete _recoilPuppiMetCorrector;
  delete _recoilPuppiMetShifter;
  delete _histTES;
  delete _histFES;
  delete _histbtagEfficiency;
  delete _btagCalib;
  delete _btagReader;
  delete _IPcorrector;
}


void SysHelper::ResetVariables()
{
  //Tau
  _tauIndex = -99;
  _tauGenMatch = -99;
  _tauPt = -99;
  _tauEta = -99;
  _tauPhi = -99;
  _tauE = -99;
  _tauDM = -99;
  _tauIPx = -99;
  _tauIPy = -99;
  _tauIPz = -99;
  _tauIPsignificance = -99;
  _tauSVx = -99;
  _tauSVy = -99;
  _tauSVz = -99;
  _GEFtauE = -99;
  _GEFtauPt = -99;
  _GEFtauPhi = -99;
  _GEFtauEta = -99;
  //Muon
  _muIndex = -99;
  _muGenMatch = -99;
  _muPt = -99;
  _muEta = -99;
  _muPhi = -99;
  _muE = -99;
  _muIso = -99;
  _muIPx = -99;
  _muIPy = -99;
  _muIPz = -99;
  _muIPsignificance = -99;
  //Jets, Pair and MET
  _isOSpair = false;
  _isIso = false;
  _isMediumID = false;
  _pairvisMass = -99;
  _Njets = -99;
  _Nbjets = -99;
  _leadingjetPt = -99;
  _leadingjetEta = -99;
  _leadingjetPhi = -99;
  _leadingjetE = -99;
  _subleadingjetPt = -99;
  _subleadingjetEta = -99;
  _subleadingjetPhi = -99;
  _subleadingjetE = -99;
  _subsubleadingjetPt = -99;
  _subsubleadingjetEta = -99;
  _subsubleadingjetPhi = -99;
  _subsubleadingjetE = -99;
  _dijetMass = -99;
  _dijetPt = -99;
  _dijetdeltaEta = -99;
  _ditauPt = -99;
  _fastMTTmass = -99;
  _muMETmt = -99;
  _PUPPImet = -99;
  _PUPPImetphi = -99;
  _pvPhiCP = -99;
  _dpPhiCP = -99;
  //Weights
  _wPrefiring = 1.;
  _wPrefiringUp = 1.;
  _wPrefiringDown = 1.;
  _wIDvsJet = 1.;
  _wIDvsJetUp = 1.;
  _wIDvsJetDown = 1.;
  _wIDvsEle = 1.;
  _wIDvsEleUp = 1.;
  _wIDvsEleDown = 1.;
  _wIDvsMu = 1.;
  _wIDvsMuUp = 1.;
  _wIDvsMuDown = 1.;
  _wTrg = 1.;
  _wTrgUp = 1.;
  _wTrgDown = 1.;
  _wIDMu = 1.;
  _wTrkMu = 1.;
  _wPU = 1.;
  _wZpT = 1.;
  _wToppT = 1.;
  _wZpTUp = 1.;
  _wZpTDown = 1.;
  _wToppTUp = 1.;
  _wToppTDown = 1.;
  _wBtag = 1.;
  _wBtagUp = 1.;
  _wBtagDown = 1.;
  _wPSISRUp = 1.;
  _wPSISRDown = 1.;
  _wPSFSRUp = 1.;
  _wPSFSRDown = 1.;
  _wScaleUp = 1.;
  _wScaleDown = 1.;
  _wMC = 1.;
  _wSignal = 1.;
  //Maps
  _TauIDSFmap.clear();
  _MuonIDTrkSFmap.clear();
  _TriggerSFmap.clear();
  _ZpTreweightingmap.clear();
  _ToppTreweightingmap.clear();
}

void SysHelper::GetCollections(const edm::View<pat::CompositeCandidate>* cands_, const edm::View<reco::Candidate>* daus_, const edm::View<pat::Jet>* jets_, const edm::View<pat::Jet>* jetsUp_, const edm::View<pat::Jet>* jetsDown_)
{

  cands = cands_;
  daus = daus_;
  jets = jets_;
  jetsUp = jetsUp_;
  jetsDown = jetsDown_;
}

void SysHelper::GetGenInfo(edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag_, std::map<std::string, double> theomap, Int_t dataMCtype) {

  generictag = generictag_;
  _TheoreticalUncmap = theomap;
  _Id = dataMCtype;
}

void SysHelper::FillGenTaus(std::vector<std::vector<unsigned int>> signal_Tauidx, std::vector<std::vector<std::vector<double>>> tauandprod_p4, std::vector<std::vector<int>> tauandprod_charge, std::vector<std::vector<std::vector<double>>> tauandprod_vtx, std::vector<std::vector<int>> taudandprod_pdgid, std::vector<unsigned int> tau_JAK) {

  if(signal_Tauidx.size()!=0) {
    for(unsigned int i = 0; i<signal_Tauidx.size(); i++) {
      if(signal_Tauidx.at(i).size()!=2) continue;
      int tauIdx, muIdx;
      if(tau_JAK.at(0) == 2) {
	tauIdx = 1;
	muIdx = 0;
      }
      else if(tau_JAK.at(1) == 2) {
        tauIdx = 0;
        muIdx = 1;
      }
      else return;
      //Fill gen PV
      _genPVx = tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(0).at(0);
      _genPVy = tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(0).at(1);
      _genPVz = tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(0).at(2);
      //Fill gen tau vis
      _genTaupx = tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(0).at(1);
      _genTaupy = tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(0).at(2);
      _genTaupz = tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(0).at(3);
      _genTauE = tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(0).at(0);
      _genTauSVx = tauandprod_vtx.at(signal_Tauidx.at(i).at(tauIdx)).at(1).at(0);
      _genTauSVy = tauandprod_vtx.at(signal_Tauidx.at(i).at(tauIdx)).at(1).at(1);
      _genTauSVz = tauandprod_vtx.at(signal_Tauidx.at(i).at(tauIdx)).at(1).at(2);
      //Fill gen mu
      _genMuonpx = tauandprod_p4.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(1);
      _genMuonpy = tauandprod_p4.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(2);
      _genMuonpz = tauandprod_p4.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(3);
      _genMuonE = tauandprod_p4.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(0);	  
      TVector3 genMuRef(tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(0) - _genPVx,
		        tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(1) - _genPVy,
		        tauandprod_vtx.at(signal_Tauidx.at(i).at(muIdx)).at(2).at(2) - _genPVz);
      //Fill gen CP
      if(tau_JAK.at(tauIdx) == 5) {
	math::XYZTLorentzVector genMuonP4(_genMuonpx, _genMuonpy, _genMuonpz, _genMuonE);
	TLorentzVector genHadTauP4(_genTaupx, _genTaupy, _genTaupz, _genTauE);
	std::vector<std::vector<double>> genPionsP4;
	std::vector<double> genPionsCharge;
	for(unsigned int j = 0; j<tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).size(); j++) {
	  if(abs(taudandprod_pdgid.at(signal_Tauidx.at(i).at(tauIdx)).at(j)) == 211 || taudandprod_pdgid.at(signal_Tauidx.at(i).at(tauIdx)).at(j) == 111) {
	    std::vector<double> genPionP4;
	    genPionP4.push_back(tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(j).at(0));
	    genPionP4.push_back(tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(j).at(1));
	    genPionP4.push_back(tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(j).at(2));
	    genPionP4.push_back(tauandprod_p4.at(signal_Tauidx.at(i).at(tauIdx)).at(j).at(3));
	    genPionsP4.push_back(genPionP4);
	    genPionsCharge.push_back(tauandprod_charge.at(signal_Tauidx.at(i).at(tauIdx)).at(j));
	  }
	}
	CPTools genCP(genPionsP4, genPionsCharge, genMuonP4, genMuRef);
	genCP.setImpactParameter();
	TVector3 genMuIP = genCP.getImpactParameter();
	_genMuonIPx = genMuIP.X();
	_genMuonIPy = genMuIP.Y();
	_genMuonIPz = genMuIP.Z();
	_gendpPhiCP = genCP.getPhiCPwithDP();
	_genpvPhiCP = genCP.getPhiCPwithPV(genHadTauP4);
      }
      else return;
    }
  }
  else return;
}

void SysHelper::GetEventInfo(bool isEmbed, bool isMC, ULong64_t runNumber, Float_t nPU, ULong64_t evtidx, Int_t lumi) {

  _isEmbed = isEmbed;
  _isMC = isMC;
  _isData = !(isEmbed || isMC);
  _runNumber = runNumber;
  _nPU = nPU;
  _evt = evtidx;
  _lumi = lumi;
}

void SysHelper::GetDecayProducts(std::vector<LorentzVectorParticle> a1lvp, std::vector<TrackParticle> muontrack, std::vector<std::vector<std::vector<double>>> pi_P4, std::vector<std::vector<double>> pi_charges) {

  A1LVP = a1lvp;
  MuonTrack = muontrack;
  RefitPionsP4 = pi_P4;
  RefitPionsCharge = pi_charges;
}

void SysHelper::GetPV(std::vector<float> x, std::vector<float> y, std::vector<float> z, std::vector<std::vector<std::vector<double>>> cov, std::vector<size_t> BS1, std::vector<size_t> BS2, std::vector<size_t> Hash) {

  _RefitPVBS_x = x;
  _RefitPVBS_y = y;
  _RefitPVBS_z = z;
  _RefitPVBS_Cov = cov;
  _VertexHashBS1 = BS1;
  _VertexHashBS2 = BS2;
  _LeptonHash = Hash;
}

bool SysHelper::FillPV(int tauIndex, int muonIndex) {

  if(_RefitPVBS_x.size()>0 && muonIndex!=-99  && tauIndex!=-99)
    {
      std::vector<size_t> hashes;
      size_t hash = 0;
      boost::hash_combine(hash, _LeptonHash.at(tauIndex));
      boost::hash_combine(hash, _LeptonHash.at(muonIndex));

      hashes.push_back(hash);
      hash = 0;

      boost::hash_combine(hash, _LeptonHash.at(muonIndex));
      boost::hash_combine(hash, _LeptonHash.at(tauIndex));
      hashes.push_back(hash);

      for (unsigned int ivertex =0; ivertex<_VertexHashBS1.size(); ivertex++){
	size_t selectionHash = 0;
	boost::hash_combine(selectionHash, _VertexHashBS1.at(ivertex));
	boost::hash_combine(selectionHash, _VertexHashBS2.at(ivertex));
	if ( std::find(hashes.begin(), hashes.end(), selectionHash) != hashes.end() ){
	  _pvx = _RefitPVBS_x.at(ivertex);
	  _pvy = _RefitPVBS_y.at(ivertex);
	  _pvz = _RefitPVBS_z.at(ivertex);
	  pvcov = _RefitPVBS_Cov.at(ivertex);
	}
      } // loop over refitted vertices collection
      return true;
    }
  else return false;
}

void SysHelper::GetMETCov(float cov00, float cov10, float cov11) {

  _PUPPIMETCov00 = cov00;
  _PUPPIMETCov10 = cov10;
  _PUPPIMETCov11 = cov11;
}

void SysHelper::GetJECUnc(std::map<std::string, std::vector<Float_t>> JECmapUp_, std::map<std::string, std::vector<Float_t>> JECmapDown_, myJECMap* JECmap_, JetCorrectionUncertainty* JECunc_) {

  JECmapUp = JECmapUp_;
  JECmapDown = JECmapDown_;
  JECmap = JECmap_;
  JECunc = JECunc_;
}

void SysHelper::GetPrefiringWeights(double pref, double prefUp, double prefDown) {

  _wPrefiring = pref;
  _wPrefiringUp = prefUp;
  _wPrefiringDown = prefDown;
}

void SysHelper::GetTauSpinnerWeights(const double wEven, const double wOdd, const double wMM) {

  _wEven = wEven;
  _wOdd = wOdd;
  _wMM = wMM;
}

void SysHelper::MakeBranches(TTree *tree, bool isNominal) {
  
  tree->Branch("run", &_runNumber);
  tree->Branch("evt", &_evt);
  tree->Branch("lumi", &_lumi);
  tree->Branch("year", &_theYear);
  tree->Branch("tauIndex", &_tauIndex);
  tree->Branch("tauGenMatch", &_tauGenMatch);
  tree->Branch("tauPt", &_tauPt);
  tree->Branch("tauEta", &_tauEta);
  tree->Branch("tauPhi", &_tauPhi);
  tree->Branch("tauE", &_tauE);
  tree->Branch("GEFtauPt", &_GEFtauPt);
  tree->Branch("GEFtauEta", &_GEFtauEta);
  tree->Branch("GEFtauPhi", &_GEFtauPhi);
  tree->Branch("GEFtauE", &_GEFtauE);
  tree->Branch("tauDM", &_tauDM);
  tree->Branch("tauIPx", &_tauIPx);
  tree->Branch("tauIPy", &_tauIPy);
  tree->Branch("tauIPz", &_tauIPz);
  tree->Branch("tauIPsignificance", &_tauIPsignificance);
  tree->Branch("tauSVx", &_tauSVx);
  tree->Branch("tauSVy", &_tauSVy);
  tree->Branch("tauSVz", &_tauSVz);
  tree->Branch("muIndex", &_muIndex);
  tree->Branch("muGenMatch", &_muGenMatch);
  tree->Branch("muPt", &_muPt);
  tree->Branch("muEta", &_muEta);
  tree->Branch("muPhi", &_muPhi);
  tree->Branch("muE", &_muE);
  tree->Branch("muIso", &_muIso);
  tree->Branch("muIPx", &_muIPx);
  tree->Branch("muIPy", &_muIPy);
  tree->Branch("muIPz", &_muIPz);
  tree->Branch("muIPsignificance", &_muIPsignificance);
  tree->Branch("pvx", &_pvx);
  tree->Branch("pvy", &_pvy);
  tree->Branch("pvz", &_pvz);
  tree->Branch("pvCov00", &_pvCov00);
  tree->Branch("pvCov11", &_pvCov11);
  tree->Branch("pvCov22", &_pvCov22);
  tree->Branch("pvCov01", &_pvCov01);
  tree->Branch("pvCov02", &_pvCov02);
  tree->Branch("pvCov12", &_pvCov12);
  tree->Branch("isOSpair", &_isOSpair);
  tree->Branch("isIso", &_isIso);
  tree->Branch("isMediumID", &_isMediumID);
  tree->Branch("pairvisMass", &_pairvisMass);
  tree->Branch("Njets", &_Njets);
  tree->Branch("leadingjetPt", &_leadingjetPt);
  tree->Branch("leadingjetEta", &_leadingjetEta);
  tree->Branch("leadingjetPhi", &_leadingjetPhi);
  tree->Branch("leadingjetE", &_leadingjetE);
  tree->Branch("subleadingjetPt", &_subleadingjetPt);
  tree->Branch("subleadingjetEta", &_subleadingjetEta);
  tree->Branch("subleadingjetPhi", &_subleadingjetPhi);
  tree->Branch("subleadingjetE", &_subleadingjetE);
  tree->Branch("subsubleadingjetPt", &_subsubleadingjetPt);
  tree->Branch("subsubleadingjetEta", &_subsubleadingjetEta);
  tree->Branch("subsubleadingjetPhi", &_subsubleadingjetPhi);
  tree->Branch("subsubleadingjetE", &_subsubleadingjetE);
  tree->Branch("dijetPt", &_dijetPt);
  tree->Branch("dijetMass", &_dijetMass);
  tree->Branch("dijetdeltaEta", &_dijetdeltaEta);
  tree->Branch("ditauPt", &_ditauPt);
  tree->Branch("fastMTTmass", &_fastMTTmass);
  tree->Branch("muMETmt", &_muMETmt);
  tree->Branch("PUPPImet", &_PUPPImet);
  tree->Branch("PUPPImetphi", &_PUPPImetphi);
  tree->Branch("PUPPIMETCov00", &_PUPPIMETCov00);
  tree->Branch("PUPPIMETCov10", &_PUPPIMETCov10);
  tree->Branch("PUPPIMETCov11", &_PUPPIMETCov11);
  tree->Branch("pvPhiCP", &_pvPhiCP);
  tree->Branch("dpPhiCP", &_dpPhiCP);
  tree->Branch("MCId", &_Id);
  tree->Branch("Npartons", &_Npartons);
  tree->Branch("isData", &_isData);
  tree->Branch("isZ", &_isZ);
  tree->Branch("isW", &_isW);
  tree->Branch("isSignal", &_isSignal);
  tree->Branch("isQCD", &_isQCD);
  tree->Branch("isVV", &_isVV);
  tree->Branch("isTTbar", &_isTTbar);
  tree->Branch("isSingleTop", &_isSingleTop);
  tree->Branch("wEven", &_wEven);
  tree->Branch("wOdd", &_wOdd);
  tree->Branch("wMM", &_wMM);
  tree->Branch("wPrefiring", &_wPrefiring);
  tree->Branch("wIDvsJet", &_wIDvsJet);
  tree->Branch("wIDvsEle", &_wIDvsEle);
  tree->Branch("wIDvsMu", &_wIDvsMu);
  tree->Branch("wTrg", &_wTrg);
  tree->Branch("wIDMu", &_wIDMu);
  tree->Branch("wTrkMu", &_wTrkMu);
  tree->Branch("wPU", &_wPU);
  tree->Branch("wZpT", &_wZpT);
  tree->Branch("wToppT", &_wToppT);
  tree->Branch("wBtag", &_wBtag);
  tree->Branch("wMC", &_wMC);
  tree->Branch("wSignal", &_wSignal);
  if(isNominal) {
    tree->Branch("wPrefiringUp", &_wPrefiringUp);
    tree->Branch("wPrefiringDown", &_wPrefiringDown);
    tree->Branch("wIDvsJetUp", &_wIDvsJetUp);
    tree->Branch("wIDvsJetDown", &_wIDvsJetDown);
    tree->Branch("wIDvsEleUp", &_wIDvsEleUp);
    tree->Branch("wIDvsEleDown", &_wIDvsEleDown);
    tree->Branch("wIDvsMuUp", &_wIDvsMuUp);
    tree->Branch("wIDvsMuDown", &_wIDvsMuDown);
    tree->Branch("wTrgUp", &_wTrgUp);
    tree->Branch("wTrgDown", &_wTrgDown);
    tree->Branch("wZpTUp", &_wZpTUp);
    tree->Branch("wZpTDown", &_wZpTDown);
    tree->Branch("wToppTUp", &_wToppTUp);
    tree->Branch("wToppTDown", &_wToppTDown);
    tree->Branch("wBtagUp", &_wBtagUp);
    tree->Branch("wBtagDown", &_wBtagDown);
    tree->Branch("wPSISRUp", &_wPSISRUp);
    tree->Branch("wPSISRDown", &_wPSISRDown);
    tree->Branch("wPSFSRUp", &_wPSFSRUp);
    tree->Branch("wPSFSRDown", &_wPSFSRDown);
    tree->Branch("wScaleUp", &_wScaleUp);
    tree->Branch("wScaleDown", &_wScaleDown);
    tree->Branch("genPVx", &_genPVx);
    tree->Branch("genPVy", &_genPVy);
    tree->Branch("genPVz", &_genPVz);
    tree->Branch("genTaupx", &_genTaupx);
    tree->Branch("genTaupy", &_genTaupy);
    tree->Branch("genTaupz", &_genTaupz);
    tree->Branch("genTauE", &_genTauE);
    tree->Branch("genTauSVx", &_genTauSVx);
    tree->Branch("genTauSVy", &_genTauSVy);
    tree->Branch("genTauSVz", &_genTauSVz);
    tree->Branch("genMuonpx", &_genMuonpx);
    tree->Branch("genMuonpy", &_genMuonpy);
    tree->Branch("genMuonpz", &_genMuonpz);
    tree->Branch("genMuonE", &_genMuonE);
    tree->Branch("genMuonIPx", &_genMuonIPx);
    tree->Branch("genMuonIPy", &_genMuonIPy);
    tree->Branch("genMuonIPz", &_genMuonIPz);
    tree->Branch("gendpPhiCP", &_gendpPhiCP);
    tree->Branch("genpvPhiCP", &_genpvPhiCP);
  }
}

void SysHelper::FillTree(TTree *tree, std::string sysType, std::string var, const edm::Event& event, bool trig, std::vector<Long64_t> _daughters_trgMatched, std::vector<math::XYZTLorentzVector> LeptonP4, const pat::MET &PUPPImet, int Npartons) {

  // Set sample type boolean
  GetSampleType();
  //
  bool PairSelection = SelectPair(sysType, var, event, LeptonP4, trig, _daughters_trgMatched);
  if(PairSelection == false) return;
  //
  double dRmin1=0.00001;
  double dRmin2=0.00001;
  for(unsigned int i=0;i<LeptonP4.size();i++) {
    if(abs(deltaR(LeptonP4.at(i),SelectedPair.daughter(0)->p4()))<dRmin1){_muIndex=i;dRmin1=deltaR(LeptonP4.at(i),SelectedPair.daughter(0)->p4());}
    if(abs(deltaR(LeptonP4.at(i),SelectedPair.daughter(1)->p4()))<dRmin2){_tauIndex=i;dRmin2=deltaR(LeptonP4.at(i),SelectedPair.daughter(1)->p4());}
  }
  //
  bool PVFilling = FillPV(_tauIndex, _muIndex);
  if(PVFilling == false) return;
  TVector3 PV(_pvx, _pvy, _pvz);
  //
  const reco::Candidate *tau = SelectedPair.daughter(1);
  const reco::Candidate *mu = SelectedPair.daughter(0);

  // Tau and Mu ES are corrected from here
  _tauDM = userdatahelpers::getUserFloat(tau,"MVADM2017v1");
  if(_tauDM<0) return;
  _tauPt = tau->pt();
  _tauEta = tau->eta();
  _tauPhi = tau->phi();
  _tauE = tau->energy();
  if(!_isData) {
    _tauGenMatch = seltools::GenMatch(SelectedPair.daughter(1), event, generictag);
    _muGenMatch = seltools::GenMatch(SelectedPair.daughter(0), event, generictag);
  }
  _muPt = mu->pt();
  _muEta = mu->eta();
  _muPhi = mu->phi();
  _muE = mu->energy();
  _muIso = userdatahelpers::getUserFloat(mu,"combRelIsoPF");
  _pairvisMass = (tau->p4() + mu->p4()).M();
  _Npartons = Npartons;
  //
  std::pair<std::vector<math::XYZTLorentzVector>, TVector2> JetsandMET = CorrectedJetsandMET(sysType, var, PUPPImet, mu, tau, JECmap);
  // Jets are corrected from here
  std::vector<math::XYZTLorentzVector> SelectedJets = JetsandMET.first;
  _Njets = SelectedJets.size();
  if(_Njets>=1){
    _leadingjetPt = SelectedJets.at(0).pt();
    _leadingjetEta = SelectedJets.at(0).eta();
    _leadingjetPhi = SelectedJets.at(0).phi();
    _leadingjetE = SelectedJets.at(0).energy();
    if(_Njets>=2){
      _subleadingjetPt = SelectedJets.at(1).pt();
      _subleadingjetEta = SelectedJets.at(1).eta();
      _subleadingjetPhi = SelectedJets.at(1).phi();
      _subleadingjetE = SelectedJets.at(1).energy();
      _dijetPt = (SelectedJets.at(0) + SelectedJets.at(1)).pt();
      _dijetMass = (SelectedJets.at(0) + SelectedJets.at(1)).M();
      _dijetdeltaEta = std::abs(SelectedJets.at(0).eta() - SelectedJets.at(1).eta());
      if(_Njets>=3){
        _subsubleadingjetPt = SelectedJets.at(2).pt();
        _subsubleadingjetEta = SelectedJets.at(2).eta();
        _subsubleadingjetPhi = SelectedJets.at(2).phi();
        _subsubleadingjetE = SelectedJets.at(2).energy();
      }
    }
  }
  //
  float ShiftedPUPPImet_px;
  float ShiftedPUPPImet_py;
  TVector2 thePUPPImet = JetsandMET.second;
  if(sysType == "TES" || sysType == "MES") {
    ShiftedPUPPImet_px = thePUPPImet.Px() + LeptonP4.at(_tauIndex).px() + LeptonP4.at(_muIndex).px() - tau->px() - mu->px();
    ShiftedPUPPImet_py = thePUPPImet.Py() + LeptonP4.at(_tauIndex).py() + LeptonP4.at(_muIndex).py() - tau->py() - mu->py();
  }
  else {
    ShiftedPUPPImet_px = thePUPPImet.Px();
    ShiftedPUPPImet_py = thePUPPImet.Py();
  }
  if((_isSignal || _isW || _isZ) && !_isEmbed){
    corrector::METRecoilCorrection(event, generictag, PUPPImet, _Njets, ShiftedPUPPImet_px, ShiftedPUPPImet_py, sysType, var, _recoilPuppiMetCorrector, _recoilPuppiMetShifter);
  }
  // MET is corrected from here
  _ditauPt = std::sqrt(std::pow((tau->px() + mu->px() + ShiftedPUPPImet_px),2) + std::pow((tau->py() + mu->py() + ShiftedPUPPImet_py),2));
  _muMETmt = seltools::ComputeMT(mu->p4(),ShiftedPUPPImet_px,ShiftedPUPPImet_py);
  if(_muMETmt>50) return;
  _PUPPImet = std::sqrt(std::pow(ShiftedPUPPImet_px,2) + std::pow(ShiftedPUPPImet_py,2));
  _PUPPImetphi = (TVector3(ShiftedPUPPImet_px,ShiftedPUPPImet_py,0.)).Phi();

  // Compute FastMTT ditau mass here
  TMatrixD PUPPImetCovMatrix(2,2);
  PUPPImetCovMatrix[0][0] = _PUPPIMETCov00;
  PUPPImetCovMatrix[1][0] = _PUPPIMETCov10;
  PUPPImetCovMatrix[0][1] = _PUPPIMETCov10;
  PUPPImetCovMatrix[1][1] = _PUPPIMETCov11;
  //
  FastMTT FastMTTAlgo;
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
  //
  classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToMuDecay, _muPt, _muEta, _muPhi, mu->mass(), 0);
  classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, _tauPt, _tauEta, _tauPhi, tau->mass(), _tauDM);
  //
  measuredTauLeptons.push_back(lep1);
  measuredTauLeptons.push_back(lep2);
  FastMTTAlgo.run(measuredTauLeptons, ShiftedPUPPImet_px, ShiftedPUPPImet_py, PUPPImetCovMatrix);
  //
  classic_svFit::LorentzVector tau1P4mtt = FastMTTAlgo.getTau1P4();
  classic_svFit::LorentzVector tau2P4mtt = FastMTTAlgo.getTau2P4();
  //
  _fastMTTmass = (tau1P4mtt + tau2P4mtt).mass();
   
  // Compute PhiCP here for a1mu channel
  if(_tauDM == 10 && RefitPionsP4.size()==3 && A1LVP.at(_tauIndex).LV().P()!=0 && pvcov.size()!=0) {
    _pvCov00 = pvcov[0][0];
    _pvCov11 = pvcov[1][1];
    _pvCov22 = pvcov[2][2];
    _pvCov01 = pvcov[0][1];
    _pvCov02 = pvcov[0][2];
    _pvCov12 = pvcov[1][2];

    TVector3 muRef(mu->vx(), mu->vy(), mu->vz());
    CPTools CP(RefitPionsP4.at(_tauIndex), RefitPionsCharge.at(_tauIndex), mu->p4(), muRef);
    PTObject METinput = CP.setMET(ShiftedPUPPImet_px, ShiftedPUPPImet_py, PUPPImetCovMatrix);
    CP.initGEF(A1LVP.at(_tauIndex), MuonTrack.at(_muIndex), PV, pvcov, METinput);
    _GEFtauE = CP.getTau("E");
    _GEFtauPt = CP.getTau("Pt");
    _GEFtauPhi = CP.getTau("Phi");
    _GEFtauEta = CP.getTau("Eta");
    CP.setImpactParameter();
    if(_isMC) {
      IpCorrection ipcorrector = (*_IPcorrector);
      CP.correctIP(event, generictag, mu->eta(), ipcorrector);
    }
    TVector3 muIP = CP.getImpactParameter();
    //
    _muIPx = muIP.X();
    _muIPy = muIP.Y();
    _muIPz = muIP.Z();
    _muIPsignificance = CP.getIPsignificance(pvcov);
    _pvPhiCP = CP.getPhiCPwithPV();
    _dpPhiCP = CP.getPhiCPwithDP();

    //Calculate tau IP and its significance, needed for Fake Factors
    LorentzVectorParticle TauLVP = A1LVP.at(_tauIndex);
    TrackParticle TauTrack(TauLVP.getParMatrix(), TauLVP.getCovMatrix(), TauLVP.PDGID(), TauLVP.Mass(), TauLVP.Charge(), TauLVP.BField());

    std::vector<float> h_param = {float(TauTrack.Parameter(TrackParticle::kappa)),
                                  float(TauTrack.Parameter(TrackParticle::lambda)),
                                  float(TauTrack.Parameter(TrackParticle::phi)),
                                  float(TauTrack.Parameter(TrackParticle::dxy)),
                                  float(TauTrack.Parameter(TrackParticle::dz))};

    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(TauLVP.Parameter(0),TauLVP.Parameter(1),TauLVP.Parameter(2));
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pvtx(PV.X(),PV.Y(),PV.Z());

    ImpactParameter tauIP;
    TVector3 ip = tauIP.CalculatePCA(TauTrack.BField(), h_param, ref, pvtx);
    _tauIPx = ip.x();
    _tauIPy = ip.y();
    _tauIPz = ip.z();
    _tauSVx = TauLVP.Vertex().X();
    _tauSVy = TauLVP.Vertex().Y();
    _tauSVz = TauLVP.Vertex().Z();

    ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov;
    TMatrixTSym<double> cov = TauLVP.getCovMatrix();
    SMatrixSym3D SigmaPrV;
    for(int i=0; i<5; i++) {
      for(int j=0; j<5; j++) {
        helixCov(i,j) = cov(i,j);
        if(i<3 && j<3) SigmaPrV(i,j) = pvcov[i][j];
      }
    }

    ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = tauIP.CalculatePCACovariance(helixCov, SigmaPrV);

    double mag = ip.Mag();
    ROOT::Math::SVector<double, 3> ip_svec;
    ip_svec(0) = ip.X();
    ip_svec(1) = ip.Y();
    ip_svec(2) = ip.Z();
    ip_svec = ip_svec.Unit();

    double uncert = sqrt(ROOT::Math::Dot( ip_svec, ip_cov * ip_svec));
    _tauIPsignificance = mag/uncert;
  }
  // Various Weights
  if((_isMC && !_isQCD) || _isEmbed) {
    _TauIDSFmap = weight::TauIDSF(_tauGenMatch, _tauDM, tau->p4(), _isEmbed, sysType, _w, _Label);
    _MuonIDTrkSFmap = weight:: MuonIDTrkSF(mu->p4(), _muIso, _isEmbed, _w);
    _TriggerSFmap = weight::TriggerSF(_tauGenMatch, _tauDM, tau->p4(), _muIso, mu->p4(), _theYear, _isEmbed, sysType, _w);
    _wIDvsJet = _TauIDSFmap["wIDvsJet"];
    _wIDvsEle = _TauIDSFmap["wIDvsEle"];
    _wIDvsMu = _TauIDSFmap["wIDvsMu"];
    _wIDMu = _MuonIDTrkSFmap["wID"];
    _wTrkMu = _MuonIDTrkSFmap["wTrk"];
    _wTrg = _TriggerSFmap["wTrg"];
    _wMC = _TheoreticalUncmap["wMC"];
    if(_isEmbed && _wMC > 10000.) _wMC *= 0.000000001;
    //
    if(sysType == "Nominal") {
      _wIDvsJetUp = _TauIDSFmap["wIDvsJetUp"];
      _wIDvsJetDown = _TauIDSFmap["wIDvsJetDown"];
      _wIDvsEleUp = _TauIDSFmap["wIDvsEleUp"];
      _wIDvsEleDown = _TauIDSFmap["wIDvsEleDown"];
      _wIDvsMuUp = _TauIDSFmap["wIDvsMuUp"];
      _wIDvsMuDown = _TauIDSFmap["wIDvsMuDown"];
      _wTrgUp = _TriggerSFmap["wTrgUp"];
      _wTrgDown = _TriggerSFmap["wTrgDown"];
    }
    //
    if(!_isEmbed) {
      _wPU = weight::PileUpreweighting(_nPU, _PU_data, _PU_mc);
      if(_isZ) {
	_ZpTreweightingmap = weight::ZpTreweighting(event, generictag, sysType, _w);
	_wZpT = _ZpTreweightingmap["wZpT"];
	if(sysType == "Nominal") {
	  _wZpTUp = _ZpTreweightingmap["wZpTUp"];
	  _wZpTDown = _ZpTreweightingmap["wZpTDown"];
	}
      }
      if(_isTTbar) {
	_ToppTreweightingmap = weight::ToppTreweighting(event, generictag, sysType);
	_wToppT = _ToppTreweightingmap["wToppT"];
	if(sysType == "Nominal") {
	  _wToppTUp = _ToppTreweightingmap["wToppTUp"];
	  _wToppTDown = _ToppTreweightingmap["wToppTDown"];
	}
      }
      if(_isSignal) {
	_wSignal = weight::SignalReweighting(_theYear, _Idstr); 
	if(sysType == "Nominal") {
	  _wPSISRUp = _TheoreticalUncmap["wPSISRUp"];
	  _wPSISRDown = _TheoreticalUncmap["wPSISRDown"];
	  _wPSFSRUp = _TheoreticalUncmap["wPSFSRUp"];
	  _wPSFSRDown = _TheoreticalUncmap["wPSFSRDown"];
	  _wScaleUp = _TheoreticalUncmap["wScaleUp"];
	  _wScaleDown = _TheoreticalUncmap["wScaleDown"];
	}
      }
    }
  }
  tree->Fill();
}


/////////////////////////////
//End of public definitions//
/////////////////////////////

std::pair<std::vector<math::XYZTLorentzVector>, TVector2> SysHelper::CorrectedJetsandMET(std::string sysType, std::string var, const pat::MET& PUPPImet, const reco::Candidate* cand1, const reco::Candidate* cand2, myJECMap* JECmap) {

  std::vector<std::string> JESsys = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_YEAR",
    "EC2_YEAR",
    "Absolute_YEAR",
    "HF_YEAR",
    "RelativeSample_YEAR"
  };
  bool JES = false;
  if(std::find(JESsys.begin(), JESsys.end(), sysType) != JESsys.end()) JES = true;

  auto year = std::to_string(_theYear);
  auto foundYear = sysType.find("YEAR");
  if(foundYear != std::string::npos) sysType.replace(sysType.end()-4, sysType.end(), year);

  std::pair<std::vector<math::XYZTLorentzVector>, TVector2> JetsandMET;

  int pTcut = 20; //Looser cut for Up and Down JES variations
  if(!JES) pTcut = 30;

  //Jets
  std::vector<math::XYZTLorentzVector> JetsP4Vect = SelectJets(cand1, cand2, jets, pTcut, sysType);
  std::vector<math::XYZTLorentzVector> JetsP4VectJES;
  math::XYZTLorentzVector JetsP4Sum;
  math::XYZTLorentzVector JetsP4SumJES;
  for(unsigned int i=0; i<JetsP4Vect.size(); i++) {
    JetsP4Sum+=JetsP4Vect.at(i);
    if(JES) {
      (*JECmap)[sysType]->setJetEta(JetsP4Vect.at(i).eta());
      (*JECmap)[sysType]->setJetPt(JetsP4Vect.at(i).pt());
      float uncertainty = 0;
      if(var == "Up") uncertainty = (*JECmap)[sysType]->getUncertainty(true);
      if(var == "Down") uncertainty = -(*JECmap)[sysType]->getUncertainty(false);
      //Now cut at pT>30 as it should be
      if((JetsP4Vect.at(i)*(1+uncertainty)).pt() > 30) {
        JetsP4SumJES += JetsP4Vect.at(i)*(1+uncertainty);
	JetsP4VectJES.push_back(JetsP4Vect.at(i)*(1+uncertainty));
      }
    }
  }
  if(JES) {
    JetsP4Vect.clear();
    JetsP4Vect = JetsP4VectJES;
  }

  math::XYZTLorentzVector JetsP4SumJER;
  if(sysType == "JER") {
    JetsP4Vect.clear();
    if(var == "Up") JetsP4Vect = SelectJets(cand1, cand2, jetsUp, pTcut, sysType);
    if(var == "Down") JetsP4Vect = SelectJets(cand1, cand2, jetsDown, pTcut, sysType);
    for(unsigned int i=0; i<JetsP4Vect.size(); i++) {
      JetsP4SumJER+=JetsP4Vect.at(i);
    }
  }
  if(JetsP4Vect.size()>0){
    struct {
      bool operator()(math::XYZTLorentzVector a, math::XYZTLorentzVector b) const {return a.pt() > b.pt();}
    } SortbyPt;
    std::sort(JetsP4Vect.begin(), JetsP4Vect.end(), SortbyPt);
  }
  JetsandMET.first = JetsP4Vect;

  // Shift MET accordingly
  double PUPPImet_px;
  double PUPPImet_py;
  if(JES) {
    PUPPImet_px = PUPPImet.px() + JetsP4Sum.Px() - JetsP4SumJES.Px();
    PUPPImet_py = PUPPImet.py() + JetsP4Sum.Py() - JetsP4SumJES.Py();
  }
  else if (sysType == "JER") {
    PUPPImet_px = PUPPImet.px() + JetsP4Sum.Px() - JetsP4SumJER.Px();
    PUPPImet_py = PUPPImet.py() + JetsP4Sum.Py() - JetsP4SumJER.Py();
  }
  else {
    PUPPImet_px = PUPPImet.px();
    PUPPImet_py = PUPPImet.py();
  }

  JetsandMET.second = TVector2(PUPPImet_px, PUPPImet_py);

  return JetsandMET;
}

std::vector<math::XYZTLorentzVector> SysHelper::SelectJets(const reco::Candidate* cand1, const reco::Candidate* cand2, const edm::View< pat::Jet>* jets, int pTcut, std::string sysType) {

  std::vector<const pat::Jet*> SelectedPATJets;
  std::vector<math::XYZTLorentzVector> JetsP4Vect;
  //
  for(edm::View<pat::Jet>::const_iterator ijet=jets->begin(); ijet!=jets->end(); ++ijet) {

    float jecFactor = ijet->jecFactor("Uncorrected");
    float jetRawPt = jecFactor * ijet->pt();
    if(_theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;

    // PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    int NumNeutralParticles =ijet->neutralMultiplicity();
    int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());
    bool tightJetID = false;

    // 2016 data
    if(_theYear == 2016) {
      if(absjeta <= 2.7) tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
      else if(absjeta <= 3.0) tightJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
      else tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    }
    // 2017 data
    else if(_theYear == 2017) {
      if(absjeta <= 2.7) tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0) || absjeta>2.4) );
      else if(absjeta <= 3.0) tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
      else tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 );
    }
    // 2018 data
    else if(_theYear == 2018) {
      if(absjeta <= 2.6) tightJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 );
      else if(absjeta>2.6 && absjeta <= 2.7) tightJetID = ( NHF<0.90 && NEMF<0.99 && CHM>0 );
      else if(absjeta>2.7 && absjeta <= 3.0) tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
      else tightJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );
    }

    if(ijet->pt()<pTcut || abs(ijet->eta())>4.7 || deltaR(ijet->p4(),cand1->p4())<0.5 || deltaR(ijet->p4(),cand2->p4())<0.5 || !tightJetID) continue;
    JetsP4Vect.push_back(ijet->p4());
    SelectedPATJets.push_back(&(*ijet));
  }

  // BTagging
  if((_isMC && !_isQCD) || _isEmbed) {
    std::map<std::string, double> BTaggingSFmap = weight::BTaggingSF(SelectedPATJets, _theYear, sysType, _histbtagEfficiency, _btagCalib, _btagReader);
    _wBtag = BTaggingSFmap["wBtag"];
    if(sysType == "Nominal") {
      _wBtagUp = BTaggingSFmap["wBtagUp"];
      _wBtagDown = BTaggingSFmap["wBtagDown"];
    }
  }
  //
  return JetsP4Vect;
}

bool SysHelper::SelectPair(std::string sysType, std::string var, const edm::Event& event, std::vector<math::XYZTLorentzVector> LeptonP4, bool trig, std::vector<Long64_t> _daughters_trgMatched)
{
  int Npairs = 0;
  math::XYZTLorentzVector TauP4Corrected, MuonP4;
  std::vector<pat::CompositeCandidate> candVector;

  for(edm::View<pat::CompositeCandidate>::const_iterator candi = cands->begin(); candi!=cands->end();++candi) {
    pat::CompositeCandidate cand = (*candi);

    float decayMode=-1.;
    int decayModeFindingNewDMs=-1;
    int genMatch=-1;
    int byVVLooseDeepTau2017v2p1VSe=-1,byTightDeepTau2017v2p1VSmu=-1, byMediumDeepTau2017v2p1VSjet=-1, byVVVLooseDeepTau2017v2p1VSjet=-1;

    //MuTau
    if(cand.daughter(0)->isMuon() && !cand.daughter(0)->isElectron() && !cand.daughter(1)->isMuon() && !cand.daughter(1)->isElectron()){
      decayMode=userdatahelpers::getUserFloat (cand.daughter(1), "decayMode");
      decayModeFindingNewDMs=userdatahelpers::getUserInt (cand.daughter(1), "decayModeFindingNewDMs");
      if(!_isData){
        genMatch = seltools::GenMatch(cand.daughter(1), event, generictag);
      }
      else genMatch = 6;
      //!\/
      if(_isMC){
        if(sysType == "TES") TauP4Corrected = corrector::P4Corrected(cand.daughter(1)->p4(),genMatch,decayMode, var, _histTES, _histFES);
        else TauP4Corrected = corrector::P4Corrected(cand.daughter(1)->p4(),genMatch,decayMode, "Nom", _histTES, _histFES);
        cand.daughter(1)->setP4(TauP4Corrected);
      }
      else TauP4Corrected = cand.daughter(1)->p4();
      //
      byVVVLooseDeepTau2017v2p1VSjet=userdatahelpers::getUserInt (cand.daughter(1), "byVVVLooseDeepTau2017v2p1VSjet");
      byMediumDeepTau2017v2p1VSjet=userdatahelpers::getUserInt (cand.daughter(1), "byMediumDeepTau2017v2p1VSjet");
      byVVLooseDeepTau2017v2p1VSe=userdatahelpers::getUserInt (cand.daughter(1), "byVVLooseDeepTau2017v2p1VSe");
      byTightDeepTau2017v2p1VSmu=userdatahelpers::getUserInt (cand.daughter(1), "byTightDeepTau2017v2p1VSmu");
      //
      if(_isMC && sysType == "MES"){
        MuonP4 = corrector::MuP4Corrected(cand.daughter(0)->p4(), var);
        cand.daughter(0)->setP4(MuonP4);
      }
      else MuonP4 = cand.daughter(0)->p4();
    }
    else continue;
    //di-lepton && third lepton
    bool DiLepton = false, ThirdLepton = false;
    for(edm::View<reco::Candidate>::const_iterator daui = daus->begin(); daui!=daus->end();++daui){
      const reco::Candidate* candi = &(*daui);
      if(deltaR(candi->p4(),TauP4Corrected)>0.5 && deltaR(candi->p4(),MuonP4)>0.5){
	if(seltools::EleVeto(candi) || seltools::MuVeto(candi)) ThirdLepton = true;
      }
      for(edm::View<reco::Candidate>::const_iterator dauj = daus->begin(); dauj!=daus->end();++dauj){
	const reco::Candidate* candj = &(*dauj);
	if(candi == candj) continue;
	if(seltools::DiEle(candi,candj) || seltools::DiMuon(candi,candj)) DiLepton = true;
      }
    }
    if(DiLepton == true || ThirdLepton == true) return false; 
    //pT
    if(TauP4Corrected.Pt()<20 || MuonP4.Pt()<20) continue;
    //eta
    if(std::abs(MuonP4.Eta())>2.4 || std::abs(TauP4Corrected.Eta())>2.3) continue;
    //dz
    if(std::abs(userdatahelpers::getUserFloat(cand.daughter(0),"dz")) > 0.2 || std::abs(userdatahelpers::getUserFloat(cand.daughter(1),"dz")) > 0.2) continue;
    //dxy
    if(std::abs(userdatahelpers::getUserFloat(cand.daughter(0),"dxy")) > 0.045) continue;
    //decay mode
    if(decayModeFindingNewDMs<0.5 || decayModeFindingNewDMs == 5 || decayModeFindingNewDMs == 6) continue;
    //ID
    if(byVVVLooseDeepTau2017v2p1VSjet<0.5 || byVVLooseDeepTau2017v2p1VSe<0.5 || byTightDeepTau2017v2p1VSmu<0.5) continue;
    if(byVVVLooseDeepTau2017v2p1VSjet>0.5 && byMediumDeepTau2017v2p1VSjet>0.5) _isMediumID = true;
    if(!seltools::CHECK_BIT(userdatahelpers::getUserInt(cand.daughter(0),"muonID"),2)) continue;
    //muon isolation
    if(userdatahelpers::getUserFloat(cand.daughter(0),"combRelIsoPF") < 0.15) _isIso = true;
    //charge
    if(std::abs(cand.daughter(0)->charge())!=1 || std::abs(cand.daughter(1)->charge())!=1) continue;
    //We remove OS cut to keep SS region for DNN training and fake factors, save boolean instead for analysis
    if(cand.daughter(1)->charge()*cand.daughter(0)->charge()<0) _isOSpair = true; 
    //delta R
    if(deltaR(TauP4Corrected,MuonP4)<0.5) continue;
    int tauIndex=-99,muonIndex=-99;
    double dRmin1=0.00001;
    double dRmin2=0.00001;
    for(unsigned int i=0;i<LeptonP4.size();i++)
      {
	if(abs(deltaR(LeptonP4.at(i),cand.daughter(0)->p4()))<dRmin1){muonIndex=i;dRmin1=deltaR(LeptonP4.at(i),cand.daughter(0)->p4());}
	if(abs(deltaR(LeptonP4.at(i),cand.daughter(1)->p4()))<dRmin2){tauIndex=i;dRmin2=deltaR(LeptonP4.at(i),cand.daughter(1)->p4());}
      }

    if(tauIndex!=-99 && muonIndex!=-99)
      {
        if(trig && _theYear == 2016){
	  if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),2) && MuonP4.pt()>23 && std::abs(MuonP4.eta())<2.1);
	  else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),3) && MuonP4.pt()>23 && std::abs(MuonP4.eta())<2.1);
	  else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),4) && MuonP4.pt()>23 && std::abs(MuonP4.eta())<2.1);
	  else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),6) && MuonP4.pt()>23 && std::abs(MuonP4.eta())<2.1);
	  else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),19) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),19));
          else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),20) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),20));
	  else continue;
	}
        else if(trig && _theYear == 2017){
	  if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),0) && MuonP4.pt()>25);
          else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),1) && MuonP4.pt()>25);
          else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),4) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),4) && MuonP4.pt()>21 && TauP4Corrected.pt()>32 && std::abs(TauP4Corrected.eta())<2.1);
          else continue;
	}
        else if(trig && _theYear == 2018){
	  if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),0) && MuonP4.pt()>25);
          else if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),2) && MuonP4.pt()>25);
	  else if(_isMC || (_isData && _runNumber>317509)){
	    if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),8) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),8) && MuonP4.pt()>21 && std::abs(MuonP4.eta())<2.1 && TauP4Corrected.pt()>32 && std::abs(TauP4Corrected.eta())<2.1);
	    else continue;
          }
          else if(_isData && _runNumber<=317509){
	    if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),13) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),7) && MuonP4.pt()>21 && std::abs(MuonP4.eta())<2.1 && TauP4Corrected.pt()>32 && std::abs(TauP4Corrected.eta())<2.1);
	    else continue;
	  }
          else if(_isEmbed){
            if(seltools::CHECK_BIT(_daughters_trgMatched.at(muonIndex),7) && seltools::CHECK_BIT(_daughters_trgMatched.at(tauIndex),7) && MuonP4.pt()>21 && std::abs(MuonP4.eta())<2.1 && TauP4Corrected.pt()>32 && std::abs(TauP4Corrected.eta())<2.1);
            else continue;
          }
        }
	else continue;
      }
    candVector.push_back(cand);
    Npairs++;
  }
  if(Npairs!=0 && (TauP4Corrected+MuonP4).M()>40) {
    std::sort(candVector.begin(),candVector.end(),seltools::ComparePairsbyIso);
    SelectedPair = candVector.at(0);
    return true;
  }
  else return false;
}

void SysHelper::GetSampleType() {

  if(_Idstr=="dy_ll" || _Idstr=="dy_1qll" || _Idstr=="dy_2qll" || _Idstr=="dy_3qll" || _Idstr=="dy_4qll" || _Idstr=="dy_ll_10to50") _isZ = true;
  if(_Idstr=="w_lnu" || _Idstr=="w_1qlnu" || _Idstr=="w_2qlnu" || _Idstr=="w_3qlnu" || _Idstr=="w_4qlnu") _isW = true;
  if(_Idstr=="h_tautau_ggf" || _Idstr=="h_tautau_vbf" || _Idstr=="zh_tautau" || _Idstr=="wplush_tautau" || _Idstr=="wminush_tautau") _isSignal = true;
  if(_Idstr=="qcd") _isQCD = true;
  if(_Idstr=="wz_2l2q" || _Idstr=="wz_3l1nu" || _Idstr=="wz_1l3nu" || _Idstr=="wz_1l1nu2q" || _Idstr=="zz_4l" || _Idstr=="zz_2l2nu" || _Idstr=="zz_2l2q" || _Idstr=="ww_2l2nu" || _Idstr=="ww_1l1nu2q") _isVV = true;
  if(_Idstr=="ttbar_dilep" || _Idstr=="ttbar_hadr" || _Idstr=="ttbar_semilep") _isTTbar = true;
  if(_Idstr=="tw" || _Idstr=="tbarw" || _Idstr=="st_tchannel_top" || _Idstr=="st_tchannel_antitop") _isSingleTop = true;
}


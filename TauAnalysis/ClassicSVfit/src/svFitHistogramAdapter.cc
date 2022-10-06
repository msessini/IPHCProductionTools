#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TMath.h>
#include <TFile.h>
#include <TObject.h>
#include <TLorentzVector.h>

#include <numeric>

#include <boost/algorithm/string/replace.hpp>

using namespace classic_svFit;

TH1* HistogramTools::compHistogramDensity(TH1 const* histogram)
{
  TH1* histogram_density = static_cast<TH1*>(histogram->Clone((std::string(histogram->GetName())+"_density").c_str()));
  histogram_density->Scale(1.0, "width");
  return histogram_density;
}

void HistogramTools::extractHistogramProperties(
    TH1 const* histogram,
    double& xMaximum,
    double& xMaximum_interpol,
    double& xMean,
    double& xQuantile016,
    double& xQuantile050,
    double& xQuantile084
)
{
  // compute median, -1 sigma and +1 sigma limits on reconstructed mass

  if ( histogram->Integral() > 0. ) {
    Double_t q[3];
    Double_t probSum[3];
    probSum[0] = 0.16;
    probSum[1] = 0.50;
    probSum[2] = 0.84;
    (const_cast<TH1*>(histogram))->GetQuantiles(3, q, probSum);
    xQuantile016 = q[0];
    xQuantile050 = q[1];
    xQuantile084 = q[2];
  } else {
    xQuantile016 = 0.;
    xQuantile050 = 0.;
    xQuantile084 = 0.;
  }

  xMean = histogram->GetMean();

  TH1* histogram_density = HistogramTools::compHistogramDensity(histogram);
  if ( histogram_density->Integral() > 0. ) {
    int binMaximum = histogram_density->GetMaximumBin();
    xMaximum = histogram_density->GetBinCenter(binMaximum);
    double yMaximum = histogram_density->GetBinContent(binMaximum);
    if ( binMaximum > 1 && binMaximum < histogram_density->GetNbinsX() ) {
      int binLeft       = binMaximum - 1;
      double xLeft      = histogram_density->GetBinCenter(binLeft);
      double yLeft      = histogram_density->GetBinContent(binLeft);

      int binRight      = binMaximum + 1;
      double xRight     = histogram_density->GetBinCenter(binRight);
      double yRight     = histogram_density->GetBinContent(binRight);

      double xMinus     = xLeft - xMaximum;
      double yMinus     = yLeft - yMaximum;
      double xPlus      = xRight - xMaximum;
      double yPlus      = yRight - yMaximum;

      xMaximum_interpol = xMaximum + 0.5*(yPlus*square(xMinus) - yMinus*square(xPlus))/(yPlus*xMinus - yMinus*xPlus);
    } else {
      xMaximum_interpol = xMaximum;
    }
  } else {
    xMaximum = 0.;
    xMaximum_interpol = 0.;
  }
  delete histogram_density;
}

double HistogramTools::extractValue(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double value = maximum;
  return value;
}

double HistogramTools::extractUncertainty(TH1 const* histogram)
{
  double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084;
  HistogramTools::extractHistogramProperties(histogram, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084);
  double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum, 2.) + TMath::Power(maximum - quantile016, 2.)));
  return uncertainty;
}

double HistogramTools::extractLmax(TH1 const* histogram)
{
  TH1* histogram_density = HistogramTools::compHistogramDensity(histogram);
  double Lmax = histogram_density->GetMaximum();
  delete histogram_density;
  return Lmax;
}

TH1* HistogramTools::makeHistogram(const std::string& histogramName, double xMin, double xMax, double logBinWidth)
{
  if ( xMin <= 0. ) xMin = 0.1;
  int numBins = 1 + TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
  TArrayF binning(numBins + 1);
  binning[0] = 0.;
  double x = xMin;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    binning[idxBin] = x;
    x *= logBinWidth;
  }
  TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, binning.GetArray());  
  return histogram;
}

int SVfitQuantity::nInstances = 0;

SVfitQuantity::SVfitQuantity() :
  uniqueName("_SVfitQuantity_"+std::to_string(++SVfitQuantity::nInstances))
{
}

SVfitQuantity::~SVfitQuantity()
{
  if (histogram_ != nullptr) delete histogram_;
}

void SVfitQuantity::bookHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  if (histogram_ != nullptr) delete histogram_;
  histogram_ = createHistogram(vis1P4, vis2P4, met);
  histogram_->SetName((histogram_->GetName()+uniqueName).c_str());
}

const TH1* SVfitQuantity::getHistogram() const { return histogram_;}

void SVfitQuantity::writeHistogram() const
{
  if (histogram_ != nullptr)
  {
    std::string histogramName = histogram_->GetName();
    boost::replace_all(histogramName, uniqueName, "");
    histogram_->Write(histogramName.c_str(), TObject::kWriteDelete);
  }
}

void SVfitQuantity::fillHistogram(const double & value, const double & weight)
{
  histogram_->Fill(value, weight);
}

void SVfitQuantity::fillHistogram(
    const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
    const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met
)
{
  histogram_->Fill(fitFunction(tau1P4, tau2P4, tauSumP4, vis1P4, vis2P4, met));
}

double SVfitQuantity::extractValue() const
{
  return HistogramTools::extractValue(histogram_);
}

double SVfitQuantity::extractUncertainty() const
{
  return HistogramTools::extractUncertainty(histogram_);
}

double SVfitQuantity::extractLmax() const
{
  return HistogramTools::extractLmax(histogram_);
}

bool SVfitQuantity::isValidSolution() const
{
  return (extractLmax() > 0.0);
}

TH1* DiTauSystemPtSVfitQuantity::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return HistogramTools::makeHistogram("ClassicSVfitIntegrand_histogramPt", 1., 1.e+3, 1.025);
}

double DiTauSystemPtSVfitQuantity::fitFunction(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                               const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return tauSumP4.pt();
}

TH1* DiTauSystemEtaSVfitQuantity::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return new TH1D("ClassicSVfitIntegrand_histogramEta", "ClassicSVfitIntegrand_histogramEta", 198, -9.9, +9.9);
}

double DiTauSystemEtaSVfitQuantity::fitFunction(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                                const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return tauSumP4.eta();
}

TH1* DiTauSystemPhiSVfitQuantity::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return new TH1D("ClassicSVfitIntegrand_histogramPhi", "ClassicSVfitIntegrand_histogramPhi", 180, -TMath::Pi(), +TMath::Pi());
}

double DiTauSystemPhiSVfitQuantity::fitFunction(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                                const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return tauSumP4.phi();
}

TH1* DiTauSystemMassSVfitQuantity::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  double visMass = (vis1P4 + vis2P4).mass();
  double minMass = visMass/1.0125;
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  return HistogramTools::makeHistogram("ClassicSVfitIntegrand_histogramMass", minMass, maxMass, 1.025);
}

double DiTauSystemMassSVfitQuantity::fitFunction(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                                 const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  return tauSumP4.mass();
}

TH1* TransverseMassSVfitQuantity::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  classic_svFit::LorentzVector measuredDiTauSystem = vis1P4 + vis2P4;
  double visTransverseMass2 = square(vis1P4.Et() + vis2P4.Et()) - (square(measuredDiTauSystem.px()) + square(measuredDiTauSystem.py()));
  double visTransverseMass = TMath::Sqrt(TMath::Max(1., visTransverseMass2));
  double minTransverseMass = visTransverseMass/1.0125;
  double maxTransverseMass = TMath::Max(1.e+4, 1.e+1*minTransverseMass);
  return HistogramTools::makeHistogram("ClassicSVfitIntegrand_histogramTransverseMass", minTransverseMass, maxTransverseMass, 1.025);
}

double TransverseMassSVfitQuantity::fitFunction(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                                const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{

  double transverseMass2 = square(tau1P4.Et() + tau2P4.Et()) - (square(tauSumP4.px()) + square(tauSumP4.py()));
  return TMath::Sqrt(TMath::Max(1., transverseMass2));
}

TauSVfitQuantity::TauSVfitQuantity(size_t tauIndex) :
  classic_svFit::SVfitQuantity(),
  m_tauIndex(tauIndex),
  m_tauLabel("Tau"+std::to_string(tauIndex+1))
{
}

TauESVfitQuantity::TauESVfitQuantity(size_t tauIndex) : TauSVfitQuantity(tauIndex)
{
}
TH1* TauESVfitQuantity::createHistogram(const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  double visEnergy = (m_tauIndex == 0 ? vis1P4 : vis2P4).E();
  return classic_svFit::HistogramTools::makeHistogram(std::string("svfitAlgorithm_histogram"+m_tauLabel+"E").c_str(), visEnergy/1.025, TMath::Max(1.e+3, 1.e+1*visEnergy/1.025), 1.025);
}
double TauESVfitQuantity::fitFunction(const classic_svFit::LorentzVector& tau1P4, const classic_svFit::LorentzVector& tau2P4, const classic_svFit::LorentzVector& tauSumP4,
                                      const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return (m_tauIndex == 0 ? tau1P4 : tau2P4).E();
}

TauERatioSVfitQuantity::TauERatioSVfitQuantity(size_t tauIndex) : TauSVfitQuantity(tauIndex)
{
}
TH1* TauERatioSVfitQuantity::createHistogram(const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return new TH1D(std::string("svfitAlgorithm_histogram"+m_tauLabel+"ERatio").c_str(), std::string("svfitAlgorithm_histogram"+m_tauLabel+"ERatio").c_str(), 200, 0.0, 1.0);
}
double TauERatioSVfitQuantity::fitFunction(const classic_svFit::LorentzVector& tau1P4, const classic_svFit::LorentzVector& tau2P4, const classic_svFit::LorentzVector& tauSumP4,
                                           const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  double visEnergy = (m_tauIndex == 0 ? vis1P4 : vis2P4).E();
  double tauEnergy = (m_tauIndex == 0 ? tau1P4 : tau2P4).E();
  return (tauEnergy != 0.0 ? visEnergy / tauEnergy : 0.0);
}

TauPtSVfitQuantity::TauPtSVfitQuantity(size_t tauIndex) : TauSVfitQuantity(tauIndex)
{
}
TH1* TauPtSVfitQuantity::createHistogram(const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return classic_svFit::HistogramTools::makeHistogram(std::string("svfitAlgorithm_histogram"+m_tauLabel+"Pt").c_str(), 1., 1.e+3, 1.025);
}
double TauPtSVfitQuantity::fitFunction(const classic_svFit::LorentzVector& tau1P4, const classic_svFit::LorentzVector& tau2P4, const classic_svFit::LorentzVector& tauSumP4,
                                       const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return (m_tauIndex == 0 ? tau1P4 : tau2P4).pt();
}

TauEtaSVfitQuantity::TauEtaSVfitQuantity(size_t tauIndex) : TauSVfitQuantity(tauIndex)
{
}
TH1* TauEtaSVfitQuantity::createHistogram(const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return new TH1D(std::string("svfitAlgorithm_histogram"+m_tauLabel+"Eta").c_str(), std::string("svfitAlgorithm_histogram"+m_tauLabel+"Eta").c_str(), 198, -9.9, +9.9);
}
double TauEtaSVfitQuantity::fitFunction(const classic_svFit::LorentzVector& tau1P4, const classic_svFit::LorentzVector& tau2P4, const classic_svFit::LorentzVector& tauSumP4,
                                        const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return (m_tauIndex == 0 ? tau1P4 : tau2P4).eta();
}

TauPhiSVfitQuantity::TauPhiSVfitQuantity(size_t tauIndex) : TauSVfitQuantity(tauIndex)
{
}
TH1* TauPhiSVfitQuantity::createHistogram(const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return new TH1D(std::string("svfitAlgorithm_histogram"+m_tauLabel+"Phi").c_str(), std::string("svfitAlgorithm_histogram"+m_tauLabel+"Phi").c_str(), 180, -TMath::Pi(), +TMath::Pi());
}
double TauPhiSVfitQuantity::fitFunction(const classic_svFit::LorentzVector& tau1P4, const classic_svFit::LorentzVector& tau2P4, const classic_svFit::LorentzVector& tauSumP4,
                                        const classic_svFit::LorentzVector& vis1P4, const classic_svFit::LorentzVector& vis2P4, const classic_svFit::Vector& met) const
{
  return (m_tauIndex == 0 ? tau1P4 : tau2P4).phi();
}


HistogramAdapter::HistogramAdapter(std::vector<SVfitQuantity*> const& quantities) :
  quantities_(quantities)
{}

HistogramAdapter::~HistogramAdapter()
{
  /*
  for (std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); quantity != quantities_.end(); ++quantity)
  {
    delete *quantity;
  }
  */
}

void HistogramAdapter::setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  vis1P4_ = vis1P4;
  vis2P4_ = vis2P4;
  met_ = met;
}

void HistogramAdapter::setTau1And2P4(const LorentzVector& tau1P4, const LorentzVector& tau2P4) {
  tau1P4_ = tau1P4;
  tau2P4_ = tau2P4;
  tauSumP4_ = tau1P4_ + tau2P4_;
}

unsigned int HistogramAdapter::registerQuantity(SVfitQuantity* quantity)
{
  quantities_.push_back(quantity);
  return getNQuantities() - 1;
}

const SVfitQuantity* HistogramAdapter::getQuantity(unsigned int iQuantity) const
{
  if(iQuantity>=getNQuantities()) return 0;
  else return  quantities_[iQuantity];
}

void HistogramAdapter::bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met)
{
  for (std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); quantity != quantities_.end(); ++quantity)
  {
    (*quantity)->bookHistogram(vis1P4, vis2P4, met);
  }
}

void HistogramAdapter::fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tauSumP4,
                                      const LorentzVector& vis1P4, const LorentzVector& vis2P4, const Vector& met) const
{
  for (std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); quantity != quantities_.end(); ++quantity)
  {
    (*quantity)->fillHistogram(tau1P4, tau2P4, tauSumP4, vis1P4, vis2P4, met);
  }
}

void HistogramAdapter::writeHistograms(const std::string& likelihoodFileName) const
{
  TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
  likelihoodFile->cd();

  for (std::vector<SVfitQuantity*>::iterator quantity = quantities_.begin(); quantity != quantities_.end(); ++quantity)
  {
    (*quantity)->writeHistogram();
  }

  likelihoodFile->Write();
  likelihoodFile->Close();
  delete likelihoodFile;
}

double HistogramAdapter::DoEval(const double* x) const
{
  fillHistograms(tau1P4_, tau2P4_, tauSumP4_, vis1P4_, vis2P4_, met_);
  return 0.;
}

double HistogramAdapter::extractValue(size_t index) const
{
  return quantities_.at(index)->extractValue();
}

double HistogramAdapter::extractUncertainty(size_t index) const
{
  return quantities_.at(index)->extractUncertainty();
}

double HistogramAdapter::extractLmax(size_t index) const
{
  return quantities_.at(index)->extractLmax();
}

std::vector<double> HistogramAdapter::extractValues() const
{
  std::vector<double> results;
  std::transform(quantities_.begin(), quantities_.end(), results.begin(),
                 [](SVfitQuantity* quantity) { return quantity->extractValue(); });
  return results;
}

std::vector<double> HistogramAdapter::extractUncertainties() const
{
  std::vector<double> results;
  std::transform(quantities_.begin(), quantities_.end(), results.begin(),
                 [](SVfitQuantity* quantity) { return quantity->extractUncertainty(); });
  return results;
}

std::vector<double> HistogramAdapter::extractLmaxima() const
{
  std::vector<double> results;
  std::transform(quantities_.begin(), quantities_.end(), results.begin(),
                 [](SVfitQuantity* quantity) { return quantity->extractLmax(); });
  return results;
}

bool HistogramAdapter::isValidSolution() const
{
  return std::accumulate(quantities_.begin(), quantities_.end(), true,
                         [](bool result, SVfitQuantity* quantity) { return result && quantity->isValidSolution(); });
}

DiTauSystemHistogramAdapter::DiTauSystemHistogramAdapter(std::vector<SVfitQuantity*> const& quantities) :
  HistogramAdapter(quantities)
{
  indexPt_ = registerQuantity(new DiTauSystemPtSVfitQuantity());
  indexEta_ = registerQuantity(new DiTauSystemEtaSVfitQuantity());
  indexPhi_ = registerQuantity(new DiTauSystemPhiSVfitQuantity());
  indexMass_ = registerQuantity(new DiTauSystemMassSVfitQuantity());
  indexTransverseMass_ = registerQuantity(new TransverseMassSVfitQuantity());
}

double DiTauSystemHistogramAdapter::getPt() const
{
  return extractValue(indexPt_);
}

double DiTauSystemHistogramAdapter::getPtErr() const
{
  return extractUncertainty(indexPt_);
}

double DiTauSystemHistogramAdapter::getPtLmax() const
{
  return extractLmax(indexPt_);
}

double DiTauSystemHistogramAdapter::getEta() const
{
  return extractValue(indexEta_);
}

double DiTauSystemHistogramAdapter::getEtaErr() const
{
  return extractUncertainty(indexEta_);
}

double DiTauSystemHistogramAdapter::getEtaLmax() const
{
  return extractLmax(indexEta_);
}

double DiTauSystemHistogramAdapter::getPhi() const
{
  return extractValue(indexPhi_);
}

double DiTauSystemHistogramAdapter::getPhiErr() const
{
  return extractUncertainty(indexPhi_);
}

double DiTauSystemHistogramAdapter::getPhiLmax() const
{
  return extractLmax(indexPhi_);
}

double DiTauSystemHistogramAdapter::getMass() const
{
  return extractValue(indexMass_);
}

double DiTauSystemHistogramAdapter::getMassErr() const
{
  return extractUncertainty(indexMass_);
}

double DiTauSystemHistogramAdapter::getMassLmax() const
{
  return extractLmax(indexMass_);
}

double DiTauSystemHistogramAdapter::getTransverseMass() const
{
  return extractValue(indexTransverseMass_);
}

double DiTauSystemHistogramAdapter::getTransverseMassErr() const
{
  return extractUncertainty(indexTransverseMass_);
}

double DiTauSystemHistogramAdapter::getTransverseMassLmax() const
{
  return extractLmax(indexTransverseMass_);
}

TauTauHistogramAdapter::TauTauHistogramAdapter(std::vector<classic_svFit::SVfitQuantity*> const& quantities) :
  DiTauSystemHistogramAdapter(quantities)
{
  indexTau1Pt = registerQuantity(new TauPtSVfitQuantity(0));
  indexTau1Eta = registerQuantity(new TauEtaSVfitQuantity(0));
  indexTau1Phi = registerQuantity(new TauPhiSVfitQuantity(0));
  indexTau2Pt = registerQuantity(new TauPtSVfitQuantity(1));
  indexTau2Eta = registerQuantity(new TauEtaSVfitQuantity(1));
  indexTau2Phi = registerQuantity(new TauPhiSVfitQuantity(1));
}

classic_svFit::LorentzVector TauTauHistogramAdapter::GetFittedHiggsLV() const
{
  TLorentzVector momentum;
  momentum.SetPtEtaPhiM(getPt(), getEta(), getPhi(), getMass());
  return classic_svFit::LorentzVector(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
}

classic_svFit::LorentzVector TauTauHistogramAdapter::GetFittedTau1LV() const
{
  TLorentzVector momentum;
  momentum.SetPtEtaPhiM(extractValue(indexTau1Pt), extractValue(indexTau1Eta), extractValue(indexTau1Phi), classic_svFit::tauLeptonMass);
  return classic_svFit::LorentzVector(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
}

classic_svFit::LorentzVector TauTauHistogramAdapter::GetFittedTau2LV() const
{
  TLorentzVector momentum;
  momentum.SetPtEtaPhiM(extractValue(indexTau2Pt), extractValue(indexTau2Eta), extractValue(indexTau2Phi), classic_svFit::tauLeptonMass);
  return classic_svFit::LorentzVector(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
}

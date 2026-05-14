#include "CATS.h"
#include "CommonAnaFunctions.h"

#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"

#include "DLM_Integration.h"
#include "CATSconstants.h"
#include "DLM_SubPads.h"
#include "DLM_Bessel.h"
#include "DLM_Fitters.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFractionFitter.h"
#include "TGenPhaseSpace.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"

#include <math.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TGaxis.h"
#include <complex>

const double Mass_He3 = 2808.39142;
const int He3_PDG = 1000020030;
const int p_PDG = 2212;

double kstarToInvariantMass(const double kstar) {
  return std::sqrt(Mass_He3 * Mass_He3 + kstar * kstar) + std::sqrt(Mass_p * Mass_p + kstar * kstar);
}

//***************************************************************************
void he3p_square_well()
{
  TFile *OutputFile;

  // SET THE GAUSSIAN SOURCE in fm
  CATSparameters SOURCE_PARS(CATSparameters::tSource, 1, true);

  // SET THE BINNING
  double kmin = 4.;
  double kmax = 700;
  double binwidth = 0.5;
  unsigned NumMomBins = int((kmax - kmin) / binwidth); // 1 MeV binning
  const double mmin = 3747;
  const double mmax = 4247;

  OutputFile = new TFile("output/pHe3_square_well.root", "recreate");

  //std::vector<double> radii = {1.5, 2., 2.5, 3., 4., 5., 5.24, 5.59, 5.76, 6., 6.23, 6.39, 7.01};
  //std::vector<double> radii = {7.12, 7.06, 7.16};
  
  // order: nominal value, nominal - lower uncertainty, nominal + upper uncertainty
  std::vector<double> radii = {6.11, 6.11 - 0.57, 6.11 + 0.57,
                               5.15, 5.15 - 0.60, 5.15 + 0.58,
                               4.02, 4.02 - 0.52, 4.02 + 0.53};

  for (const auto& radius : radii) {
  
    // Check if the correlation function has already been computed for this radius value
    if (OutputFile->Get(Form("r=%.2f_fm/hhe3_Sigma_plus_CF", radius))) {
      std::cout << "Correlation function for R = " << radius << " fm already exists. Skipping computation.\n";
      continue;
    }
    
    std::cout << "Computing the correlation function for R = " << radius << " fm\n";

    SOURCE_PARS.SetParameter(0, radius);

    CATS cats;
    cats.SetMomBins(NumMomBins, kmin, kmax);
    cats.SetAnaSource(GaussSource, SOURCE_PARS);
    cats.SetAnaSource(0, SOURCE_PARS.GetParameter(0));
    cats.SetUseAnalyticSource(true);
    cats.SetMomentumDependentSource(false);
    cats.SetThetaDependentSource(false);
    cats.SetExcludeFailedBins(false);
    cats.SetQ1Q2(2);
    cats.SetPdgId(He3_PDG, p_PDG);
    cats.SetRedMass((Mass_He3 * Mass_p) / (Mass_He3 + Mass_p));
    cats.SetNumChannels(2);
    cats.SetChannelWeight(0, 0.25);
    cats.SetChannelWeight(1, 0.75);
    cats.SetQuantumStatistics(0);
    
    // set s wave on the singlet
    cats.SetNumPW(0, 1);
    // set s wave on the triplet
    cats.SetNumPW(1, 1);

    // Strong interaction
    // L=0, S=0
    double squareWellSingletPars[4] = {-2.3806, 1.796947, -114.5162, 2.049168};
    CATSparameters catsSquareWellSingletPars(CATSparameters::tPotential, 4, true);
    catsSquareWellSingletPars.SetParameters(squareWellSingletPars);
    cats.SetShortRangePotential(0, 0, DoubleSquareWell, catsSquareWellSingletPars);

    //double squareWllSingletPars[2] = {114.5162, 2.049168};
    //double squareWllSingletPars[2] = {114.5162, 1.049168};
    //CATSparameters catsSquareWellSingletPars(CATSparameters::tPotential, 2, true);
    //catsSquareWellSingletPars.SetParameters(squareWllSingletPars);
    //cats.SetShortRangePotential(0, 0, SquareWell, catsSquareWellSingletPars);

    //double doubleGausPars[4] = {-18.170, 1.689, -16.754, 1.820};
    //CATSparameters catsDoubleGausPars(CATSparameters::tPotential, 4, true);
    //catsDoubleGausPars.SetParameters(doubleGausPars);
    //cats.SetShortRangePotential(0, 0, DoubleGaussSum, catsDoubleGausPars);

    // L=0, S=1
    //double squareWelTripletPars[4] = {-6.3256, 0.922117, -10000., 10.0};
    double squareWelTripletPars[4] = {-6.3256, 0.922117, -155.1059, 1.175850};
    CATSparameters catsSquareWellTripletPars(CATSparameters::tPotential, 4, true);
    catsSquareWellTripletPars.SetParameters(squareWelTripletPars);
    cats.SetShortRangePotential(1, 0, DoubleSquareWell, catsSquareWellTripletPars);
    
    //double squareWelTripletPars[2] = {-155.1059, 1.175850};
    //CATSparameters catsSquareWellTripletPars(CATSparameters::tPotential, 2, true);
    //catsSquareWellTripletPars.SetParameters(squareWelTripletPars);
    //cats.SetShortRangePotential(1, 0, SquareWell, catsSquareWellTripletPars);
    ////////////////////////////////////

    cats.KillTheCat();

    TH1F cats_CF("hcats_CF", "cats_CF; #it{k}* (MeV/#it{c})", NumMomBins, kmin, kmax);
    TH1F cats_InvMass("hcats_InvMass", "cats_InvMass; #it{m}_{inv} (MeV/#it{c}^{2})", NumMomBins, mmin, mmax);
    for (unsigned uBin = 0; uBin < NumMomBins; uBin++)
    {
      double CF = cats.GetCorrFun(uBin);
      double CF_err = cats.GetCorrFunErr(uBin);
      cats_CF.SetBinContent(uBin, CF);
      cats_CF.SetBinError(uBin, CF_err);
      double kStar = cats_CF.GetBinCenter(uBin);
      double invMass = kstarToInvariantMass(kStar);
      cats_InvMass.SetBinContent(cats_InvMass.FindBin(invMass), CF);
      cats_InvMass.SetBinError(cats_InvMass.FindBin(invMass), CF_err);
    }

    TDirectory* OutputDir = OutputFile->mkdir(Form("r=%.2f_fm", radius));
    OutputDir->cd();
    cats_CF.Write();
    cats_InvMass.Write();

  }
  delete OutputFile;
}

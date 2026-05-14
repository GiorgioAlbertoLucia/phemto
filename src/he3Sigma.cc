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
const double Mass_Sigma = 1189.37;
const int He3_PDG = 1000020030;
const int p_PDG = 2212;
const int d_PDG = 1000010020;
const int Sigma_PDG = 3212;

//***************************************************************************
void he3Sigma()
{
  TFile *OutputFile;

  // SET THE GAUSSIAN SOURCE in fm
  CATSparameters SOURCE_PARS(CATSparameters::tSource, 1, true);

  double kmin = 4.;
  double kmax = 700;
  double binwidth = 0.5;
  unsigned NumMomBins = int((kmax - kmin) / binwidth); // 1 MeV binning

  OutputFile = new TFile("output/he3Sigma.root", "recreate");

  //std::vector<double> radii = {2., 3., 4., 5., 5.24, 5.76, 6., 6.39};
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

    CATS he3_Sigma_plus;
    he3_Sigma_plus.SetMomBins(NumMomBins, kmin, kmax);
    he3_Sigma_plus.SetAnaSource(GaussSource, SOURCE_PARS);
    he3_Sigma_plus.SetAnaSource(0, SOURCE_PARS.GetParameter(0));
    he3_Sigma_plus.SetUseAnalyticSource(true);
    he3_Sigma_plus.SetMomentumDependentSource(false);
    he3_Sigma_plus.SetThetaDependentSource(false);
    he3_Sigma_plus.SetExcludeFailedBins(false);
    he3_Sigma_plus.SetQ1Q2(2);
    he3_Sigma_plus.SetPdgId(Sigma_PDG, He3_PDG);
    he3_Sigma_plus.SetRedMass((Mass_Sigma * Mass_He3) / (Mass_Sigma + Mass_He3));
    he3_Sigma_plus.SetNumChannels(2);
    he3_Sigma_plus.SetChannelWeight(0, 0.25);
    he3_Sigma_plus.SetChannelWeight(1, 0.75);
    he3_Sigma_plus.SetQuantumStatistics(0);
    
  
    ////////////////////////////////////

    he3_Sigma_plus.KillTheCat();

    TH1F he3_Sigma_plus_CF("hhe3_Sigma_plus_CF", "he3_Sigma_plus_CF; #it{k}* (MeV/#it{c})", NumMomBins, kmin, kmax);
    for (unsigned uBin = 0; uBin < NumMomBins; uBin++)
    {
      double CF = he3_Sigma_plus.GetCorrFun(uBin);
      double CF_err = he3_Sigma_plus.GetCorrFunErr(uBin);
      he3_Sigma_plus_CF.SetBinContent(uBin, CF);
      he3_Sigma_plus_CF.SetBinError(uBin, CF_err);
      double kStar = he3_Sigma_plus_CF.GetBinCenter(uBin);
    }

    TDirectory* OutputDir = OutputFile->mkdir(Form("r=%.2f_fm", radius));
    OutputDir->cd();
    he3_Sigma_plus_CF.Write();

  }
  delete OutputFile;
}

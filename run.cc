// Intended only to load all the necessary headers for the library and scripts
// 
// Run only from the root terminal with the following command:
// .x run.cc
// Do NOT load this file in the ROOT interpreter

#include "TSystem.h"
#include "TROOT.h"

void run() {
    gInterpreter->AddIncludePath("/home/galucia/DLM/install/include/");
    //gSystem->Load("/home/galucia/DLM/install/CMake/libCATS.dylib");
    gSystem->Load("/home/galucia/gsl-2.8/.libs/libgsl.so");
    gSystem->Load("/home/galucia/DLM/install/lib/libCATS.so");
    
    //gROOT->ProcessLine(Form(".L %s/src/he3p.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("he3p()");

    gROOT->ProcessLine(Form(".L %s/src/he3Sigma.cc", gSystem->WorkingDirectory()));
    gROOT->ProcessLine("he3Sigma()");

    //gROOT->ProcessLine(Form(".L %s/src/he3p_square_well.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("he3p_square_well()");

    //gROOT->ProcessLine(Form(".L %s/src/he3p_coulomb.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("he3p_coulomb()");

    //gROOT->ProcessLine(Form(".L %s/src/he3p_square_well_pwave.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("he3p_square_well()");

    //gROOT->ProcessLine(Form(".L %s/src/pd_square_well.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("pd_square_well()");

    //gROOT->ProcessLine(Form(".L %s/src/draw_square_well_potential.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("draw_square_well_potential()");

    //gROOT->ProcessLine(Form(".L %s/src/correlationFunctionRadii.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("correlationFunctionRadii()");

    //gROOT->ProcessLine(Form(".L %s/src/coulombCorrelation_pd.cpp");
    //gROOT->ProcessLine("coulombCorrelation_pd()");

    //gROOT->ProcessLine(Form(".L %s/src/pp_studies.cpp");
    //gROOT->ProcessLine("pp_studies()");

    //gROOT->ProcessLine(Form(".L %s/src/he3p_strong_LL.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("he3p_strong_LL()");

    //gROOT->ProcessLine(Form(".L %s/src/pd_strong_LL.cc", gSystem->WorkingDirectory()));
    //gROOT->ProcessLine("pd_strong_LL()");

    gROOT->ProcessLine(".q");
}

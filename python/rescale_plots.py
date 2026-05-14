6#23fode to convert a TGraph to a TH1F

from ROOT import TGraph, TH1F, TFile

# utils 

def scale_x_axis(hist, scale_factor):
    hist.SetTitle(hist.GetTitle() + f" ({hist.GetXaxis().GetTitle()} in GeV)")
    hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
    hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() / scale_factor, hist.GetXaxis().GetXmax() / scale_factor)
    return hist


def convert_and_rebin_histogram(infile_path:str, hist_name:str, outfile: TFile, suffix:str):
    '''
        Convert the x scale from MeV to GeV and rebin the histogram
    '''

    print(f"Processing histogram: {hist_name} from file: {infile_path}")
    infile = TFile(infile_path, "READ")
    hist = infile.Get(hist_name)
    scale_x_axis(hist, scale_factor=1000)
    hist.SetName(hist_name)

    mmin, mmax, bin_width = 0.003, 0.8, 0.002
    if 'InvMass' in hist_name:
        mmin, mmax, bin_width = 3.747, 4.247, 0.0005

    nbins = int((mmax - mmin) / bin_width)
    tmp_hist = TH1F("tmp_hist", "tmp_hist", nbins, mmin, mmax)
    for i in range(1, nbins+1):
        tmp_hist.SetBinContent(i, hist.GetBinContent(i))
        tmp_hist.SetBinError(i, hist.GetBinError(i))

    if "/" in hist_name:
        hist_name = hist_name.split("/")[1]

    outfile.cd()
    tmp_hist.Write(hist_name+suffix)

    infile.Close()

def convert_and_rebin(infile_path:str, outfile: TFile, suffix:str, hist_names:list=None):
    '''
        Convert the x scale from MeV to GeV and rebin the histogram
    '''

    if hist_names is None:
        hist_names = ['hHe3_p_Coul_CF', 'hHe3_p_Coul_InvMass']
    for hist_name in hist_names:
        convert_and_rebin_histogram(infile_path, hist_name, outfile, suffix)
    
if __name__ == "__main__":
    
    radii = set([6.11, 6.11 - 0.57, 6.11 + 0.57, 5.15, 5.15 - 0.60, 5.15 + 0.58, 4.02, 4.02 - 0.52, 4.02 + 0.53])
    print(f"Radii to process: {radii}")
    
    outfile = TFile("../output/pHe3_square_well_rescaled.root", "RECREATE")
    for radius in radii:
        convert_and_rebin("../output/pHe3_square_well.root", outfile, f"_r={radius:.2f}_fm", [f"r={radius:.2f}_fm/hcats_CF"])
    outfile.Close()


    outfile = TFile("../output/he3Sigma_rescaled.root", "RECREATE")
    for radius in radii:
        convert_and_rebin("../output/he3Sigma.root", outfile, f"_r={radius:.2f}_fm", [f"r={radius:.2f}_fm/hhe3_Sigma_plus_CF"])
    outfile.Close()
import argparse
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain

ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

## include a common.h file
ROOT.gROOT.LoadMacro('inc/Common.h++')
from ROOT import nsigmaHe3

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
args = parser.parse_args()
if args.config_file == '':
    print('** No config file provided. Exiting. **')
    exit()

with open(args.config_file, 'r') as stream:
  confFile = yaml.safe_load(stream)
  conf = confFile["data_preparation"]

bin_width = conf["bin_width"]
pt_max = conf["max_pt"]
pt_min = conf["min_pt"]
pt_bins = np.arange(pt_min, pt_max + bin_width, bin_width)
max_abs_dca = conf["max_abs_dca"]
n_bins_dca = conf["n_bins_dca"]
input_data = conf["input_data"]
input_mc = conf["input_mc"]
output_file = conf["output_file"]
outFile = ROOT.TFile(output_file, "RECREATE")

base_sels = "fTPCnCls >= 110 && std::abs(fEta) < 0.8 && std::abs(fDCAxy) < 0.7 && pt > 1 && pt < 9.0 && clSize > 5 && matter==0 && nITScls == 7"

file_data_list = input_data if isinstance(input_data, list) else [input_data]
chainData = TChain("O2nucleitable")
for fileName in file_data_list:
  fileData = TFile(fileName)
  for key in fileData.GetListOfKeys() :
    keyName = key.GetName()
    if 'DF_' in keyName :
        chainData.Add(f'{fileName}/{keyName}/O2nucleitable')

mc_data_list = input_mc if isinstance(input_mc, list) else [input_mc]
chainMC = TChain("O2nucleitablemc")
for fileName in mc_data_list:
    fileData = TFile(fileName)
    for key in fileData.GetListOfKeys() :
        keyName = key.GetName()
        if 'DF_' in keyName :
            chainMC.Add(f'{fileName}/{keyName}/O2nucleitablemc')

mc_lambdab_data_list = input_mc if isinstance(input_mc, list) else [input_mc]
chainMCLambdab = TChain("O2nucleitablemc")
for fileName in mc_lambdab_data_list:
    fileData = TFile(fileName)
    for key in fileData.GetListOfKeys() :
        keyName = key.GetName()
        if 'DF_' in keyName :
            chainMCLambdab.Add(f'{fileName}/{keyName}/O2nucleitablemc')
############################################################################################################################################################################
rdf = ROOT.ROOT.RDataFrame(chainData) \
.Define("ptUncorr", "2 * std::abs(fPt)") \
.Define("pt", "ptUncorr + 0.0343554 + 0.96161 * std::exp(-1.51286 * ptUncorr)") \
.Define("p", "pt * cosh(fEta)") \
.Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)") \
.Define("matter", "fPt > 0") \
.Define("pidForTracking", "fFlags >> 12") \
.Define("nsigmaTPC", 'nsigmaHe3(fTPCInnerParam, fTPCsignal)') \
.Define("clSize", 'averageClusterSize(fITSclusterSizes)') \
.Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)") \
.Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)") \
.Filter(base_sels)

rdfMC = ROOT.ROOT.RDataFrame(chainMC) \
.Define("ptUncorr", "2 * std::abs(fPt)") \
.Define("pt", "ptUncorr + 0.0343554 + 0.96161 * std::exp(-1.51286 * ptUncorr)") \
.Define("p", "pt * cosh(fEta)") \
.Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)") \
.Define("matter", "fPt > 0") \
.Define("pidForTracking", "fFlags >> 12") \
.Define("nsigmaTPC", 'nsigmaHe3(fTPCInnerParam * 2, fTPCsignal)') \
.Define("clSize", 'averageClusterSize(fITSclusterSizes)') \
.Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)") \
.Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)") \
.Define("isPrimary", "fFlags & (1 << 9)") \
.Filter(base_sels)

rdfMCLambdab = ROOT.ROOT.RDataFrame(chainMCLambdab) \
.Define("ptUncorr", "2 * std::abs(fPt)") \
.Define("pt", "ptUncorr + 0.0343554 + 0.96161 * std::exp(-1.51286 * ptUncorr)") \
.Define("p", "pt * cosh(fEta)") \
.Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)") \
.Define("matter", "fPt > 0") \
.Define("pidForTracking", "fFlags >> 12") \
.Define("nsigmaTPC", 'nsigmaHe3(fTPCInnerParam * 2, fTPCsignal)') \
.Define("clSize", 'averageClusterSize(fITSclusterSizes)') \
.Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)") \
.Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)") \
.Define("isPrimary", "fFlags & (1 << 9)") \
.Filter(base_sels)

rdfMCSecondaries = rdfMC.Filter("fPDGcode == -1000020030 && std::abs(fMotherPDGcode) == 1010010030") 
rdfMC = rdfMC.Filter("fPDGcode == -1000020030 && isPrimary")
############################################################################################################################################################################

hPtNSigmaTPC = rdf.Histo2D(("hPtNSigmaTPC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "pt", "nsigmaTPC")
hPtNSigmaTPCMC = rdfMC.Histo2D(("hPtNSigmaTPCMC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -5, 5), "pt", "nsigmaTPC")
hPtNSigmaTPCMCLambdab = rdfMCLambdab.Histo2D(("hPtNSigmaTPCMCLambdab", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -5, 5), "pt", "nsigmaTPC")
hPtDCAxyMC = rdfMC.Histo2D(("hPtDCAxyMC", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hPtDCAxyMCSecondaries = rdfMCSecondaries.Histo2D(("hPtDCAxyMCSecondaries", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hPtDCAxyMCLambdab = rdfMCLambdab.Histo2D(("hPtDCAxyMCLambdab", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hDCAxyDCaz = rdf.Histo2D(("hDCAxyDCaz", ";DCA_{xy} (cm);DCA_{z} (cm)", n_bins_dca, -max_abs_dca, max_abs_dca, n_bins_dca, -max_abs_dca, max_abs_dca), "fDCAxy", "fDCAz")
hDCAxyDCazMC = rdfMC.Histo2D(("hDCAxyDCazMC", ";DCA_{xy} (cm);DCA_{z} (cm)", n_bins_dca, -max_abs_dca, max_abs_dca, n_bins_dca, -max_abs_dca, max_abs_dca), "fDCAxy", "fDCAz")

hNSigmaTPCClSize = rdf.Histo2D(("hNSigmaTPCClSize", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")
hNSigmaTPCClSizeMC = rdfMC.Histo2D(("hNSigmaTPCClSizeMC", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")
hNSigmaTPCClSizeMCLambdab = rdfMCLambdab.Histo2D(("hNSigmaTPCClSizeMCLambdab", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")

roo_nsigmaTPC = ROOT.RooRealVar("nsigmaTPC", "n#sigma_{TPC}", -4, 4)
mu_nsigmaTPC = ROOT.RooRealVar("mu_nsigmaTPC", "#mu_{n#sigma_{TPC}}", 0, -1, 1)
sigma_nsigmaTPC = ROOT.RooRealVar("sigma_nsigmaTPC", "#sigma_{n#sigma_{TPC}}", 1, 1.e-4, 3)
tau_nsigmaTPC = ROOT.RooRealVar("tau_nsigmaTPC", "#tau_{n#sigma_{TPC}}", 1, -4, 4)
roo_gaus = ROOT.RooGaussian("gaus_nsigmatpc", "gaus", roo_nsigmaTPC, mu_nsigmaTPC, sigma_nsigmaTPC)
roo_exp = ROOT.RooExponential("exp_nsigmatpc", "exp", roo_nsigmaTPC, tau_nsigmaTPC)
sgn_counts = ROOT.RooRealVar("sgn_counts", "sgn_counts", 1000, 0, 1.e6)
bkg_counts = ROOT.RooRealVar("bkg_counts", "bkg_counts", 1000, 0, 1.e6)
total_pdf = ROOT.RooAddPdf("total_pdf", "total_pdf", ROOT.RooArgList(roo_gaus, roo_exp), ROOT.RooArgList(sgn_counts, bkg_counts))
hPurity = ROOT.TH1F("hPurity", ";#it{p}_{T} (GeV/#it{c});Purity", len(pt_bins) - 1, pt_bins)
outFile.mkdir("nsigma_fits")
outFile.cd("nsigma_fits")

for iPtBin in range(1, hPtNSigmaTPC.GetXaxis().GetNbins() + 1):
  print(f"Pt bin {iPtBin}")
  hSlice = hPtNSigmaTPC.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_data_nsigma = ROOT.RooDataHist(f"roo_data_nsigma_{iPtBin}", f"roo_data_nsigma_{iPtBin}", ROOT.RooArgList(roo_nsigmaTPC), hSlice)
  total_pdf.fitTo(roo_data_nsigma)
  plot = roo_nsigmaTPC.frame()
  plot.SetName(f"frame_nsigma_{iPtBin}")
  roo_data_nsigma.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  total_pdf.plotOn(plot, ROOT.RooFit.Name("model"))
  total_pdf.paramOn(plot, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
  plot.Write()
  ## calc purity between -2 and 2
  roo_nsigmaTPC.setRange("sig_region", -2, 2)
  sgn_integral = roo_gaus.createIntegral(ROOT.RooArgSet(roo_nsigmaTPC), ROOT.RooArgSet(roo_nsigmaTPC), "sig_region")
  bkg_integral = roo_exp.createIntegral(ROOT.RooArgSet(roo_nsigmaTPC), ROOT.RooArgSet(roo_nsigmaTPC), "sig_region")
  sgn_counts_2s = sgn_integral.getVal() * sgn_counts.getVal()
  bkg_counts_2s = bkg_integral.getVal() * bkg_counts.getVal()
  purity = sgn_counts_2s / (sgn_counts_2s + bkg_counts_2s)
  purity_err = np.sqrt(purity * (1 - purity) / (sgn_counts_2s + bkg_counts_2s))
  hPurity.SetBinContent(iPtBin, purity)
  hPurity.SetBinError(iPtBin, purity_err)

rdf_bkg = rdf.Filter("std::abs(nsigmaTPC) > 3")
rdf = rdf.Filter("std::abs(nsigmaTPC) < 2")
hPtDCAxy = rdf.Histo2D(("hPtDCAxy", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hPtDCAxyBkg = rdf_bkg.Histo2D(("hPtDCAxyBkg", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")

outFile.cd()

hDCAxyDCaz.Write()
hDCAxyDCazMC.Write()
hPtDCAxy.Write()
hPtDCAxyBkg.Write()
hPurity.Write()
hPtNSigmaTPC.Write()
hNSigmaTPCClSize.Write()
hPtNSigmaTPCMC.Write()
hNSigmaTPCClSizeMC.Write()
hPtDCAxyMC.Write()
hPtDCAxyMCSecondaries.Write()
hPtNSigmaTPCMCLambdab.Write()
hNSigmaTPCClSizeMCLambdab.Write()
hPtDCAxyMCLambdab.Write()


## compute he3 <-- h3l and he3 efficiencies
rdfMCGen = ROOT.ROOT.RDataFrame(chainMC).Define("isPrimary", "fFlags & (1 << 9)")
rdfMCGenPrim = rdfMCGen.Filter("fPDGcode == -1000020030 && isPrimary && fMotherPDGcode == 0")
rdfMCGenSec = rdfMCGen.Filter("fPDGcode == -1000020030 && abs(fMotherPDGcode) == 1010010030")
hPtMCGenPrim = rdfMCGenPrim.Histo1D(("hPtMCGenPrim", ";#it{p}_{T} (GeV/#it{c})", len(pt_bins) - 1, pt_bins), "fgPt")
hPtMCGenPrim_Clone = hPtMCGenPrim.Clone("hPtMCGenPrim_Clone")
hPtMCGenSec = rdfMCGenSec.Histo1D(("hPtMCGenSec", ";#it{p}_{T} (GeV/#it{c})", len(pt_bins) - 1, pt_bins), "fgPt")
hMothRadVsPtGenSec = rdfMCGenSec.Histo2D(("hMothRadVsPtGenSec", ";#it{p}_{T} (GeV/#it{c});Mother radius (cm)", len(pt_bins) - 1, pt_bins, 40, 0, 20), "fgPt", "fMotherDecRad")
hPtMCGenSec_Clone = hPtMCGenSec.Clone("hPtMCGenSec_Clone")
hPtMCRecPrim = rdfMC.Histo1D(("hPtMCRecPrim", ";#it{p}_{T} (GeV/#it{c})", len(pt_bins) - 1, pt_bins), "fgPt")
hPtMCRecSec = rdfMCSecondaries.Histo1D(("hPtMCRecSec", ";#it{p}_{T} (GeV/#it{c})", len(pt_bins) - 1, pt_bins), "fgPt")
hEffPrim = hPtMCRecPrim.Clone("hEffPrim")
hEffPrim.Divide(hPtMCGenPrim_Clone)
hEffSec = hPtMCRecSec.Clone("hEffSec")
hEffSec.Divide(hPtMCGenSec_Clone)

h3l_spectrum_file = ROOT.TFile("utils/h3l_spectrum.root")
he3_spectrum_file = ROOT.TFile("utils/he3_spectrum.root")
## Get the total He3 spectrum and its systematic and statistical uncertainties
total_he3_spectrum_stat = he3_spectrum_file.Get("fStatTPCA")
total_he3_spectrum_stat.SetDirectory(0)
## fit he3 spectrum with mt_expo
mt_expo_he3 = ROOT.TF1('mtexpo_he3', '[2]*x*exp(-TMath::Sqrt(([0]*[0]+x*x))/[1])', 0., 10)
mt_expo_he3.FixParameter(0, 2.80839)  # He3 mass
mt_expo_he3.SetParLimits(1, 0.1, 1)
mt_expo_he3.SetParLimits(2, 1.e-08, 1)
cv_he3 = ROOT.TCanvas("cv_he3", "He3 spectrum fit", 800, 600)
total_he3_spectrum_stat.Fit(mt_expo_he3, 'R')
total_he3_spectrum_syst = he3_spectrum_file.Get("fSystTPCA")
total_he3_spectrum_stat.Draw("pe")
outFile.cd()
cv_he3.Write()
total_he3_spectrum_syst.SetDirectory(0)
## Get the total H3L spectrum and its systematic and statistical uncertainties
total_h3l_spectrum_stat = h3l_spectrum_file.Get("hStat")
total_h3l_spectrum_stat.SetDirectory(0)
## fit h3l spectrum with mt_expo
mt_expo_h3l = ROOT.TF1('mtexpo_h3l', '[2]*x*exp(-TMath::Sqrt(([0]*[0]+x*x))/[1])', 0., 10)
mt_expo_h3l.FixParameter(0, 2.99131)  # H3L mass
mt_expo_h3l.SetParLimits(1, 0.1, 1)
mt_expo_h3l.SetParLimits(2, 1.e-08, 1)
total_h3l_spectrum_stat.Fit(mt_expo_h3l, 'R')
total_h3l_spectrum_syst = h3l_spectrum_file.Get("hSyst")
total_h3l_spectrum_syst.SetDirectory(0)
mt_expo_h3l = h3l_spectrum_file.Get("pt_expo")

h3l_lv = ROOT.TLorentzVector()
he3_dau_lv = ROOT.TLorentzVector()
h3l_mass = 2.99131
he3_mass = 2.80839
pi_mass = 0.13957
phase_space = ROOT.TGenPhaseSpace()

n_trials = 1e5
h3l_th1 = ROOT.TH1F("h3l_th1", "H3L spectrum", len(pt_bins) - 1, pt_bins)
he3_from_hyp_th1 = ROOT.TH1F("he3_from_hyp_th1", "Non-prompt He spectrum", len(pt_bins) - 1, pt_bins) 
for i in range(int(n_trials)):
    h3l_lv.SetPtEtaPhiM(mt_expo_h3l.GetRandom(), ROOT.gRandom.Uniform(-1, 1), ROOT.gRandom.Uniform(0, ROOT.TMath.TwoPi()), h3l_mass)
    phase_space.SetDecay(h3l_lv, 2, np.array([he3_mass, pi_mass]))
    phase_space.Generate()
    he3_from_hyp_th1.Fill(phase_space.GetDecay(0).Pt())
    h3l_th1.Fill(h3l_lv.Pt())

hHe3ToHypConv = he3_from_hyp_th1.Clone("hHe3ToHypConv")
hHe3ToHypConv.Divide(h3l_th1)
hyp_br = 0.25  # Branching ratio for H3L -> He3 + pi
hyp_br_unc = 0.1  # Relative uncertainty on the branching ratio

hExpFraction = ROOT.TH1F("hExpFraction", ";#it{p}_{T} (GeV/#it{c});Expected fraction", len(pt_bins) - 1, pt_bins)
hHypToTHe3Ratio = hExpFraction.Clone("hHypToTHe3Ratio")
hHypToTHe3Ratio.SetTitle(";#it{p}_{T} (GeV/#it{c});H3L / He3 ratio")

for iPtBin in range(1, hExpFraction.GetNbinsX() + 1):
   yield_h3l = mt_expo_h3l.Integral(hExpFraction.GetXaxis().GetBinLowEdge(iPtBin), hExpFraction.GetXaxis().GetBinUpEdge(iPtBin))
   yield_h3l *= hHe3ToHypConv.GetBinContent(iPtBin)
   yield_he3_total = mt_expo_he3.Integral(hExpFraction.GetXaxis().GetBinLowEdge(iPtBin), hExpFraction.GetXaxis().GetBinUpEdge(iPtBin))
   hHypToTHe3Ratio.SetBinContent(iPtBin, yield_h3l / yield_he3_total)
   yield_ratio = yield_h3l * hyp_br / yield_he3_total

   ## get stat and syst uncertainties for he3 and h3l yields.
   bin_for_h3l = total_h3l_spectrum_stat.FindBin(hExpFraction.GetXaxis().GetBinCenter(iPtBin))
   unc_hyp = np.sqrt(total_h3l_spectrum_stat.GetBinError(bin_for_h3l)**2 + total_h3l_spectrum_syst.GetBinError(bin_for_h3l)**2) / total_h3l_spectrum_stat.GetBinContent(bin_for_h3l)
   bin_for_he3 = total_he3_spectrum_stat.FindBin(hExpFraction.GetXaxis().GetBinCenter(iPtBin))
   unc_he3 = np.sqrt(total_he3_spectrum_stat.GetBinError(bin_for_he3)**2 + total_he3_spectrum_syst.GetBinError(bin_for_he3)**2) / total_he3_spectrum_stat.GetBinContent(bin_for_he3)
   yield_ratio_unc = yield_ratio * np.sqrt(unc_hyp**2 + unc_he3**2 + hyp_br_unc**2)
   eff_ratio = hEffSec.GetBinContent(iPtBin) / hEffPrim.GetBinContent(iPtBin)
   hExpFraction.SetBinContent(iPtBin, yield_ratio * eff_ratio)
   hExpFraction.SetBinError(iPtBin, yield_ratio_unc * eff_ratio)

outFile.mkdir("h3l_he3_fraction")
outFile.cd("h3l_he3_fraction")
mt_expo_h3l.Write()
mt_expo_he3.Write()
hHe3ToHypConv.Write()
hPtMCGenPrim.Write()
hPtMCGenSec.Write()
hPtMCRecPrim.Write()
hPtMCRecSec.Write()
hEffPrim.Write()
hEffSec.Write()
hHypToTHe3Ratio.Write()
hExpFraction.Write()
hMothRadVsPtGenSec.Write()
outFile.Close() 
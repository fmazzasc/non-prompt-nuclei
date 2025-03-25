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

base_sels = "fTPCnCls >= 110 && std::abs(fEta) < 0.8 && std::abs(fDCAxy) < 0.7 && std::abs(fDCAxy) > 1e-7 && pt > 1 && pt < 9.0 && clSize > 5 && matter==0 && nITScls == 7 "

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

rdfMCSecondaries = rdfMC.Filter("std::abs(fPDGcode) == 1000020030 && !isPrimary")
rdfMC = rdfMC.Filter("std::abs(fPDGcode) == 1000020030 && isPrimary")
############################################################################################################################################################################

hPtNSigmaTPC = rdf.Histo2D(("hPtNSigmaTPC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "pt", "nsigmaTPC")
hPtNSigmaTPCMC = rdfMC.Histo2D(("hPtNSigmaTPCMC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -5, 5), "pt", "nsigmaTPC")
hPtDCAxyMC = rdfMC.Histo2D(("hPtDCAxyMC", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hPtDCAxyMCSecondaries = rdfMCSecondaries.Histo2D(("hPtDCAxyMCSecondaries", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, n_bins_dca, -max_abs_dca, max_abs_dca), "pt", "fDCAxy")
hNSigmaTPCClSize = rdf.Histo2D(("hNSigmaTPCClSize", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")
hNSigmaTPCClSizeMC = rdfMC.Histo2D(("hNSigmaTPCClSizeMC", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")

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
hPtDCAxy.Write()
hPtDCAxyBkg.Write()
hPurity.Write()
hPtNSigmaTPC.Write()
hNSigmaTPCClSize.Write()
hPtNSigmaTPCMC.Write()
hNSigmaTPCClSizeMC.Write()
hPtDCAxyMC.Write()
hPtDCAxyMCSecondaries.Write()

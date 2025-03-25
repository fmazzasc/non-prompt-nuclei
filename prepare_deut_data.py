import numpy as np
import ROOT
from ROOT import TFile, TChain
import argparse
import yaml
import os

ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

## include a common.h file
ROOT.gROOT.LoadMacro('inc/Common.h++')
from ROOT import nsigmaDeu


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
input_data = conf["input_data"]
input_mc = conf["input_mc"]
output_file = conf["output_file"]
outFile = ROOT.TFile(output_file, "RECREATE")

############################################################################################################################################################################
if conf['analyse_tree']:
  base_sels = "fTPCnCls >= 110 && std::abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && pt > 0.6 && pt < 9.0 && std::abs(nsigmaTPC) < 2"
  base_sels_mc = "fTPCnCls >= 110 && std::abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && pt > 0.6 && pt < 9.0 && matter==0"
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

  rdf = ROOT.ROOT.RDataFrame(chainData) \
  .Define("pt", "std::abs(fPt)") \
  .Define("p", "pt * cosh(fEta)") \
  .Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * sqrt(1.f / (fBeta * fBeta) - 1.f) - 1.8756129") \
  .Define("matter", "fPt > 0") \
  .Define("pidForTracking", "fFlags >> 12") \
  .Define("nsigmaTPC", 'nsigmaDeu(fTPCInnerParam, fTPCsignal)') \
  .Define("clSize", 'averageClusterSize(fITSclusterSizes)') \
  .Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)") \
  .Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)") \
  .Define("hasTRD", "fFlags & (1 << 6)") \
  .Filter(base_sels)


  rdfMC = ROOT.ROOT.RDataFrame(chainMC) \
  .Define("pt", "std::abs(fPt)") \
  .Define("p", "pt * cosh(fEta)") \
  .Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * sqrt(1.f / (fBeta * fBeta) - 1.f) - 1.8756129") \
  .Define("matter", "fPt > 0") \
  .Define("pidForTracking", "fFlags >> 12") \
  .Define("nsigmaTPC", 'nsigmaDeu(fTPCInnerParam, fTPCsignal)') \
  .Define("clSize", 'averageClusterSize(fITSclusterSizes)') \
  .Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)") \
  .Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)") \
  .Define("hasTRD", "fFlags & (1 << 6)") \
  .Define("isPrimary", "fFlags & (1 << 9)") \
  .Filter(base_sels_mc)

  rdfMCSecondaries = rdfMC.Filter("std::abs(fPDGcode) == 1000010020 && !isPrimary")
  rdfMC = rdfMC.Filter("std::abs(fPDGcode) == 1000010020 && isPrimary")
  hPtNSigmaTPC = rdf.Histo2D(("hPtNSigmaTPC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "pt", "nsigmaTPC")
  hPtTPCSignal = rdf.Histo2D(("hPtTPCSignal", ";#it{p}_{T} (GeV/#it{c});TPC signal", len(pt_bins) - 1, pt_bins, 100, 0, 1000), "pt", "fTPCsignal")
  hPtTofMass = rdf.Histo2D(("hPtTofMass", ";#it{p}_{T} (GeV/#it{c});TOF mass (GeV/#it{c}^{2})", len(pt_bins) - 1, pt_bins, 100, -0.5, 0.5), "pt", "tofMass")
  hPtNSigmaTPCMC = rdfMC.Histo2D(("hPtNSigmaTPCMC", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -10, 10), "pt", "nsigmaTPC")
  hPtTofMassMC = rdfMC.Histo2D(("hPtTofMassMC", ";#it{p}_{T} (GeV/#it{c});TOF mass (GeV/#it{c}^{2})", len(pt_bins) - 1, pt_bins, 100, -0.5, 0.5), "pt", "tofMass")
  hPtDCAxyMC = rdfMC.Histo2D(("hPtDCAxyMC", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.02, 0.02), "pt", "fDCAxy")
  hPtDCAxyMCSecondaries = rdfMCSecondaries.Histo2D(("hPtDCAxyMCSecondaries", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.02, 0.02), "pt", "fDCAxy")
  hNSigmaTPCClSize = rdf.Histo2D(("hNSigmaTPCClSize", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")
  hNSigmaTPCClSizeMC = rdfMC.Histo2D(("hNSigmaTPCClSizeMC", ";n#sigma_{TPC};ITS cluster size", 100, -3, 3, 30, 0, 15), "nsigmaTPC", "clSize")

############################################################################################################################################################################
else:
  kPt = 0
  kDCAxy = 1
  kDCAz = 2
  kNSigmaTPC = 3
  kTOFmass = 4
  kNITSClus = 5
  kNTPCClus = 6
  kCorrectPV = 7
  kIsSecondary = 8
  kFromWeakDecay = 9
  
  thn_file = ROOT.TFile.Open(input_data)
  thn_sparse_data = thn_file.Get("nuclei-spectra/spectra/hDCAHistsA_deuteron")
  thn_file_mc = ROOT.TFile.Open(input_mc)
  thn_sparse_mc = thn_file_mc.Get("nuclei-spectra/spectra/hDCAHistsA_deuteron")

  thn_sparse_data.GetAxis(kNITSClus).SetRange(thn_sparse_data.GetAxis(5).FindBin(7), thn_sparse_data.GetAxis(5).FindBin(7))
  thn_sparse_data.GetAxis(kNTPCClus).SetRange(thn_sparse_data.GetAxis(6).FindBin(110), thn_sparse_data.GetAxis(6).FindBin(150))
  thn_sparse_data.GetAxis(kPt).SetRange(thn_sparse_data.GetAxis(0).FindBin(pt_min), thn_sparse_data.GetAxis(0).FindBin(pt_max))
  thn_sparse_mc.GetAxis(kNITSClus).SetRange(thn_sparse_mc.GetAxis(5).FindBin(7), thn_sparse_mc.GetAxis(5).FindBin(7))
  thn_sparse_mc.GetAxis(kNTPCClus).SetRange(thn_sparse_mc.GetAxis(6).FindBin(110), thn_sparse_mc.GetAxis(6).FindBin(150))
  thn_sparse_mc.GetAxis(kPt).SetRange(thn_sparse_mc.GetAxis(0).FindBin(pt_min), thn_sparse_mc.GetAxis(0).FindBin(pt_max))
  ## rebin pt axis to match the binning of the config
  bin_width_orig = thn_sparse_data.GetAxis(0).GetBinWidth(1)
  rebin_factor = int(bin_width / bin_width_orig)
  rebin_factor_mc = int(bin_width / thn_sparse_mc.GetAxis(0).GetBinWidth(1))
  
  thn_sparse_data.GetAxis(kNSigmaTPC).SetRange(thn_sparse_data.GetAxis(3).FindBin(-2), thn_sparse_data.GetAxis(3).FindBin(3))

  hPtNSigmaTPC = thn_sparse_data.Projection(kNSigmaTPC, kPt)
  hPtNSigmaTPC.RebinX(rebin_factor)
  hPtNSigmaTPC.SetName("hPtNSigmaTPC")
  hPtTofMass = thn_sparse_data.Projection(kTOFmass, kPt)
  hPtTofMass.RebinX(rebin_factor)
  hPtTofMass.SetName("hPtTofMass")
  hPtNSigmaTPCMC = thn_sparse_mc.Projection(kNSigmaTPC, kPt)
  hPtNSigmaTPCMC.RebinX(rebin_factor)
  hPtNSigmaTPCMC.SetName("hPtNSigmaTPCMC")
  hPtTofMassMC = thn_sparse_mc.Projection(kTOFmass, kPt)
  hPtTofMassMC.RebinX(rebin_factor)
  hPtTofMassMC.SetName("hPtTofMassMC")

###########################################################################################################################################################################
roo_tofmass = ROOT.RooRealVar("tofmass", "M_{TOF}", -4, 4)
mu_tof = ROOT.RooRealVar("mu_tof", "#mu_{tof}", 0, -0.5, 0.5)
sigma_tof = ROOT.RooRealVar("sigma_tof", "#sigma_{tof}", 0.1, 0.01, 0.5)
roo_signal_pdf = ROOT.RooGaussian("roo_signal_pdf", "roo_signal_pdf", roo_tofmass, mu_tof, sigma_tof)


c1_bkg = ROOT.RooRealVar("c1_bkg", "c1_bkg", 0, -1, 1)
c2_bkg = ROOT.RooRealVar("c2_bkg", "c2_bkg", 0, -1, 1)
roo_bkg_pdf = ROOT.RooChebychev("roo_bkg_pdf", "roo_bkg_pdf", roo_tofmass, ROOT.RooArgList(c1_bkg, c2_bkg))

sgn_counts = ROOT.RooRealVar("sgn_counts", "sgn_counts", 1000, 0, 1.e6)
bkg_counts = ROOT.RooRealVar("bkg_counts", "bkg_counts", 1000, 0, 1.e6)
total_pdf = ROOT.RooAddPdf("total_pdf", "total_pdf", ROOT.RooArgList(roo_signal_pdf, roo_bkg_pdf), ROOT.RooArgList(sgn_counts, bkg_counts))
hPurity = ROOT.TH1F("hPurity", ";#it{p}_{T} (GeV/#it{c});Purity", len(pt_bins) - 1, pt_bins)

outFile.mkdir("tofmass_fits")
outFile.cd("tofmass_fits")

for iPtBin in range(1, hPtTofMass.GetXaxis().GetNbins() + 1):
  hSlice = hPtTofMass.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_data_tof = ROOT.RooDataHist(f"roo_data_tof_{iPtBin}", f"roo_data_tof_{iPtBin}", ROOT.RooArgList(roo_tofmass), hSlice)
  ## fit between -0.4 and 0.4
  total_pdf.fitTo(roo_data_tof, ROOT.RooFit.Range(-0.4, 0.4))
  plot = roo_tofmass.frame()
  plot.SetName(f"frame_tof_{iPtBin}")
  roo_data_tof.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  total_pdf.plotOn(plot, ROOT.RooFit.Name("total"))
  total_pdf.plotOn(plot, ROOT.RooFit.Components("roo_signal_pdf"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("signal"))
  total_pdf.plotOn(plot, ROOT.RooFit.Components("roo_bkg_pdf"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("background"))
  total_pdf.paramOn(plot, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
  plot.Write()
  ## calc purity between -0.15 and 0.15
  roo_tofmass.setRange("sig_region", -0.15, 0.15)
  sgn_integral = roo_signal_pdf.createIntegral(ROOT.RooArgSet(roo_tofmass), ROOT.RooArgSet(roo_tofmass), "sig_region")
  bkg_integral = roo_bkg_pdf.createIntegral(ROOT.RooArgSet(roo_tofmass), ROOT.RooArgSet(roo_tofmass), "sig_region")
  sgn_counts_2s = sgn_integral.getVal() * sgn_counts.getVal()
  bkg_counts_2s = bkg_integral.getVal() * bkg_counts.getVal()
  purity = sgn_counts_2s / (sgn_counts_2s + bkg_counts_2s)
  purity_err = np.sqrt(purity * (1 - purity) / (sgn_counts_2s + bkg_counts_2s))
  hPurity.SetBinContent(iPtBin, purity)
  hPurity.SetBinError(iPtBin, purity_err)


if conf['analyse_tree']:
  rdf_bkg = rdf.Filter("tofMass < -0.5")
  rdf = rdf.Filter("tofMass > -0.15 && tofMass < 0.15")
  hPtDCAxy = rdf.Histo2D(("hPtDCAxy", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.02, 0.02), "pt", "fDCAxy")
  hPtDCAxyBkg = rdf_bkg.Histo2D(("hPtDCAxyBkg", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.02, 0.02), "pt", "fDCAxy")

else:
  thn_sparse_data.GetAxis(kTOFmass).SetRange(thn_sparse_data.GetAxis(kTOFmass).FindBin(-0.15), thn_sparse_data.GetAxis(kTOFmass).FindBin(0.15))
  hPtDCAxy = thn_sparse_data.Projection(kDCAxy, kPt)
  hPtDCAxy.SetName("hPtDCAxy")
  hPtDCAxy.RebinX(rebin_factor)
  thn_sparse_data.GetAxis(kTOFmass).SetRange(1, thn_sparse_data.GetAxis(kTOFmass).FindBin(-0.5))
  hPtDCAxyBkg = thn_sparse_data.Projection(kDCAxy, kPt)
  hPtDCAxyBkg.SetName("hPtDCAxyBkg")
  hPtDCAxyBkg.RebinX(rebin_factor)

  thn_sparse_mc.GetAxis(kIsSecondary).SetRange(thn_sparse_mc.GetAxis(kIsSecondary).FindBin(0), thn_sparse_mc.GetAxis(kIsSecondary).FindBin(0))
  hPtDCAxyMC = thn_sparse_mc.Projection(kDCAxy, kPt)
  hPtDCAxyMC.SetName("hPtDCAxyMC")
  hPtDCAxyMC.RebinX(rebin_factor_mc)

  thn_sparse_mc.GetAxis(kIsSecondary).SetRange(thn_sparse_mc.GetAxis(kIsSecondary).FindBin(1), thn_sparse_mc.GetAxis(kIsSecondary).FindBin(1))
  thn_sparse_mc.GetAxis(kFromWeakDecay).SetRange(thn_sparse_mc.GetAxis(kFromWeakDecay).FindBin(1), thn_sparse_mc.GetAxis(kFromWeakDecay).FindBin(1))
  hPtDCAxyMCSecondaries = thn_sparse_mc.Projection(kDCAxy, kPt)
  hPtDCAxyMCSecondaries.SetName("hPtDCAxyMCSecondaries")
  hPtDCAxyMCSecondaries.RebinX(rebin_factor_mc)


outFile.cd()
hPtDCAxy.Write()
hPtTofMass.Write()
hPtDCAxyBkg.Write()
hPurity.Write()
hPtNSigmaTPC.Write()
hPtTofMassMC.Write()
hPtNSigmaTPCMC.Write()
hPtDCAxyMC.Write()
hPtDCAxyMCSecondaries.Write()

if conf['analyse_tree']:
  hNSigmaTPCClSize.Write()
  hNSigmaTPCClSizeMC.Write()

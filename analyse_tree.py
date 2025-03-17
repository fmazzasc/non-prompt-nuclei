import numpy as np
import ROOT
from ROOT import TFile, TChain

ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

## include a common.h file
ROOT.gROOT.LoadMacro('inc/Common.h+')

from ROOT import nsigmaHe3

pt_bins = np.array(np.arange(1.5, 5, 0.4) , dtype=float)
base_sels = "fTPCnCls >= 110 && std::abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && pt > 1 && pt < 9.0 && clSize > 5 && matter==0"


file_data_list = ['/data3/fmazzasc/he3/2024_data/al/AO2D.root', '/data3/fmazzasc/he3/2024_data/an/AO2D.root']
chainData = TChain("O2nucleitable")
for fileName in file_data_list:
  fileData = TFile(fileName)
  for key in fileData.GetListOfKeys() :
    keyName = key.GetName()
    if 'DF_' in keyName :
        chainData.Add(f'{fileName}/{keyName}/O2nucleitable')

mc_data_list = ['/data3/fmazzasc/he3/mc/25a3/AO2D.root']
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
hPtDCAxy = rdf.Histo2D(("hPtDCAxy", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.04, 0.04), "pt", "fDCAxy")
hPtDCAxyMC = rdfMC.Histo2D(("hPtDCAxyMC", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.04, 0.04), "pt", "fDCAxy")
hPtDCAxyMCSecondaries = rdfMCSecondaries.Histo2D(("hPtDCAxyMCSecondaries", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", len(pt_bins) - 1, pt_bins, 100, -0.04, 0.04), "pt", "fDCAxy")
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
outFile = ROOT.TFile("output.root", "recreate")
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

outFile.cd()
hPurity.Write()
hPtNSigmaTPC.Write()
hNSigmaTPCClSize.Write()
hPtDCAxy.Write()
hPtNSigmaTPCMC.Write()
hNSigmaTPCClSizeMC.Write()
hPtDCAxyMC.Write()


## template of primary He3
dcaxy = ROOT.RooRealVar("dcaxy", "DCA_{xy}", 0., -0.04, 0.04, "cm")
muDCAxy = ROOT.RooRealVar("muDCAxy", "#mu_{DCA_{xy}}", 0., -0.04, 0.04, "cm")
sigmaDCAxy = ROOT.RooRealVar("sigmaDCAxy", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-6, 1.e-2, "cm")
alpha = ROOT.RooRealVar("alpha", "#alpha", 2, 0.1, 5)
n = ROOT.RooRealVar("n", "n", 2, 0, 10)
cbShape = ROOT.RooCrystalBall("cbShape", "cbShape", dcaxy, muDCAxy, sigmaDCAxy, alpha, n, True)
hSigmaCB = ROOT.TH1F("hSigmaCB", ";#it{p}_{T} (GeV/#it{c});#sigma_{CB} (cm)", len(pt_bins) - 1, pt_bins)
hAlphaCB = ROOT.TH1F("hAlphaCB", ";#it{p}_{T} (GeV/#it{c});#alpha_{CB}", len(pt_bins) - 1, pt_bins)
hNCB = ROOT.TH1F("hNCB", ";#it{p}_{T} (GeV/#it{c});n_{CB}", len(pt_bins) - 1, pt_bins)

## template of secondaries He3 from H3L
muDCAxyH3L = ROOT.RooRealVar("muDCAxyH3L", "#mu_{DCA_{xy}}", 0., -0.04, 0.04, "cm")
sigmaDCAxyH3L = ROOT.RooRealVar("sigmaDCAxyH3L", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-3, 1.e-2, "cm")
alphaH3L = ROOT.RooRealVar("alphaH3L", "#alpha", 2, 0.5, 5)
nH3L = ROOT.RooRealVar("nH3L", "n", 2, 0, 10)
cbShapeH3L = ROOT.RooCrystalBall("cbShapeH3L", "cbShapeH3L", dcaxy, muDCAxy, sigmaDCAxyH3L, alphaH3L, nH3L, True)
hSigmaCBH3L = ROOT.TH1F("hSigmaCBH3L", ";#it{p}_{T} (GeV/#it{c});#sigma_{CB} (cm)", len(pt_bins) - 1, pt_bins)
hAlphaCBH3L = ROOT.TH1F("hAlphaCBH3L", ";#it{p}_{T} (GeV/#it{c});#alpha_{CB}", len(pt_bins) - 1, pt_bins)
hNCBH3L = ROOT.TH1F("hNCBH3L", ";#it{p}_{T} (GeV/#it{c});n_{CB}", len(pt_bins) - 1, pt_bins)

outFile.mkdir("dcaxy_fits_primaries")
outFile.mkdir("dcaxy_fits_secondaries")
for iPtBin in range(1, hPtDCAxyMC.GetXaxis().GetNbins() + 1):
  outFile.cd("dcaxy_fits_primaries")
  hSlice = hPtDCAxyMC.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_mc_dcaxy = ROOT.RooDataHist(f"roo_mc_dcaxy_{iPtBin}", f"roo_mc_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
  cbShape.fitTo(roo_mc_dcaxy)
  plot = dcaxy.frame()
  plot.SetName(f"frame_dcaxy_{iPtBin}")
  roo_mc_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  cbShape.plotOn(plot, ROOT.RooFit.Name("model"))
  cbShape.paramOn(plot, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
  plot.Write()
  hSigmaCB.SetBinContent(iPtBin, sigmaDCAxy.getVal())
  hSigmaCB.SetBinError(iPtBin, sigmaDCAxy.getError())
  hAlphaCB.SetBinContent(iPtBin, alpha.getVal())
  hAlphaCB.SetBinError(iPtBin, alpha.getError())
  hNCB.SetBinContent(iPtBin, n.getVal())
  hNCB.SetBinError(iPtBin, n.getError())
  outFile.cd("dcaxy_fits_secondaries")
  hSlice = hPtDCAxyMCSecondaries.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_mc_dcaxy = ROOT.RooDataHist(f"roo_mc_secondaries_dcaxy_{iPtBin}", f"roo_mc_secondaries_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
  cbShapeH3L.fitTo(roo_mc_dcaxy)
  plot = dcaxy.frame()
  plot.SetName(f"frame_secondaries_dcaxy_{iPtBin}")
  roo_mc_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  cbShapeH3L.plotOn(plot, ROOT.RooFit.Name("model"))
  cbShapeH3L.paramOn(plot, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
  plot.Write()
  hSigmaCBH3L.SetBinContent(iPtBin, sigmaDCAxyH3L.getVal())
  hSigmaCBH3L.SetBinError(iPtBin, sigmaDCAxyH3L.getError())
  hAlphaCBH3L.SetBinContent(iPtBin, alphaH3L.getVal())
  hAlphaCBH3L.SetBinError(iPtBin, alphaH3L.getError())
  hNCBH3L.SetBinContent(iPtBin, nH3L.getVal())
  hNCBH3L.SetBinError(iPtBin, nH3L.getError())


outFile.cd('dcaxy_fits_primaries')
hSigmaCB.Write()
hAlphaCB.Write()
hNCB.Write()
outFile.cd('dcaxy_fits_secondaries')
hSigmaCBH3L.Write()
hAlphaCBH3L.Write()
hNCBH3L.Write()


## now we fit the data with the MC templates, convolved with a Gaussian
outFile.cd()
outFile.mkdir("dcaxy_fits_data")
outFile.cd("dcaxy_fits_data")

gaus_sigma = ROOT.RooRealVar("gaus_sigma", "#sigma_{Gaus}", 1e-4, 1e-6, 1e-2)
gaus_reso = ROOT.RooGaussian("gaus_reso", "gaus_reso", dcaxy, muDCAxy, gaus_sigma)

hGausSigma = ROOT.TH1F("hGausSigma", ";#it{p}_{T} (GeV/#it{c});#sigma_{Gaus} (cm)", len(pt_bins) - 1, pt_bins)
hPrimaryFrac = ROOT.TH1F("hPrimaryFrac", ";#it{p}_{T} (GeV/#it{c});Primary fraction", len(pt_bins) - 1, pt_bins)
hSecondaryFrac = ROOT.TH1F("hSecondaryFrac", ";#it{p}_{T} (GeV/#it{c});Secondary fraction", len(pt_bins) - 1, pt_bins)

for iPtBin in range(1, hPtDCAxy.GetXaxis().GetNbins() + 1):
  hSlice = hPtDCAxy.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_data_dcaxy = ROOT.RooDataHist(f"roo_data_dcaxy_{iPtBin}", f"roo_data_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
  
  sigmaDCAxy.setVal(hSigmaCB.GetBinContent(iPtBin))
  sigmaDCAxy.setConstant(True)
  alpha.setVal(hAlphaCB.GetBinContent(iPtBin))
  alpha.setConstant(True)
  n.setVal(hNCB.GetBinContent(iPtBin))
  n.setConstant(True)
  sigmaDCAxyH3L.setVal(hSigmaCBH3L.GetBinContent(iPtBin))
  sigmaDCAxyH3L.setConstant(True)
  alphaH3L.setVal(hAlphaCBH3L.GetBinContent(iPtBin))
  alphaH3L.setConstant(True)
  nH3L.setVal(hNCBH3L.GetBinContent(iPtBin))
  nH3L.setConstant(True)

  convP = ROOT.RooFFTConvPdf(f"convP_{iPtBin}", f"convP_{iPtBin}", dcaxy, cbShape, gaus_reso)
  convS = ROOT.RooFFTConvPdf(f"convS_{iPtBin}", f"convS_{iPtBin}", dcaxy, cbShapeH3L, gaus_reso)

  fracPrim = ROOT.RooRealVar("fracPrim", "fracPrim", 0.5, 0, 1)
  total_pdf = ROOT.RooAddPdf(f"total_pdf_{iPtBin}", "total_pdf", ROOT.RooArgList(convP, convS), ROOT.RooArgList(fracPrim))

  total_pdf.fitTo(roo_data_dcaxy)
  plot = dcaxy.frame()
  plot.SetName(f"frame_dcaxy_data_{iPtBin}")
  roo_data_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  total_pdf.plotOn(plot, ROOT.RooFit.Name("model"))
  total_pdf.paramOn(plot, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
  plot.Write()

  hGausSigma.SetBinContent(iPtBin, gaus_sigma.getVal())
  hGausSigma.SetBinError(iPtBin, gaus_sigma.getError())
  hPrimaryFrac.SetBinContent(iPtBin, fracPrim.getVal())
  hPrimaryFrac.SetBinError(iPtBin, fracPrim.getError())
  hSecondaryFrac.SetBinContent(iPtBin, 1 - fracPrim.getVal())
  hSecondaryFrac.SetBinError(iPtBin, fracPrim.getError())





outFile.cd()
hGausSigma.Write()
hPrimaryFrac.Write()
hSecondaryFrac.Write()

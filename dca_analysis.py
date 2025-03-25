import ROOT
import yaml
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
args = parser.parse_args()
if args.config_file == '':
    print('** No config file provided. Exiting. **')
    exit()

with open(args.config_file, 'r') as stream:
  confFile = yaml.safe_load(stream)
  conf_prep = confFile["data_preparation"]
  conf_dca = confFile["dca_analysis"]

pt_max = conf_prep["max_pt"]
pt_min = conf_prep["min_pt"]
bin_width = conf_prep["bin_width"]
pt_bins = np.arange(pt_min, pt_max + bin_width, bin_width)

input_file = ROOT.TFile.Open(conf_prep["output_file"], "read")
hPurity = input_file.Get("hPurity")
hPurity.SetDirectory(0)
hPtDCAxy = input_file.Get("hPtDCAxy")
hPtDCAxy.SetDirectory(0)
hPtDCAxyMC = input_file.Get("hPtDCAxyMC")
hPtDCAxyMC.SetDirectory(0)
hPtDCAxyMCSecondaries = input_file.Get("hPtDCAxyMCSecondaries")
hPtDCAxyMCSecondaries.SetDirectory(0)
hPtDCAxyBkg = input_file.Get("hPtDCAxyBkg")
hPtDCAxyBkg.SetDirectory(0)

outFile = ROOT.TFile.Open(conf_dca["output_file"], "recreate")

dcaxy = ROOT.RooRealVar("dcaxy", "DCA_{xy}", 0., -0.04, 0.04, "cm")
muDCAxy = ROOT.RooRealVar("muDCAxy", "#mu_{DCA_{xy}}", 0., -0.04, 0.04, "cm")

## template of background from He3 sidebands
sigmaDCAxyBkg = ROOT.RooRealVar("sigmaDCAxyBkg", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
alphaBkg = ROOT.RooRealVar("alphaBkg", "#alpha", 2, 0.2, 5)
nBkg = ROOT.RooRealVar("nBkg", "n", 2, 1., 10)
cbShapeBkg = ROOT.RooCrystalBall("cbShapeBkg", "cbShapeBkg", dcaxy, muDCAxy, sigmaDCAxyBkg, alphaBkg, nBkg, True)
hSigmaCBBkg = ROOT.TH1F("hSigmaCBBkg", ";#it{p}_{T} (GeV/#it{c});#sigma_{CB} (cm)", len(pt_bins) - 1, pt_bins)
hAlphaCBBkg = ROOT.TH1F("hAlphaCBBkg", ";#it{p}_{T} (GeV/#it{c});#alpha_{CB}", len(pt_bins) - 1, pt_bins)
hNCBBkg = ROOT.TH1F("hNCBBkg", ";#it{p}_{T} (GeV/#it{c});n_{CB}", len(pt_bins) - 1, pt_bins)

outFile.mkdir("dcaxy_fits_bkg")
outFile.cd("dcaxy_fits_bkg")
for iPtBin in range(1, hPtDCAxyBkg.GetXaxis().GetNbins() + 1):
  hSlice = hPtDCAxyBkg.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_data_dcaxy = ROOT.RooDataHist(f"roo_data_bkg_dcaxy_{iPtBin}", f"roo_data_bkg_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
  cbShapeBkg.fitTo(roo_data_dcaxy)
  plot = dcaxy.frame()
  plot.SetName(f"frame_bkg_dcaxy_{iPtBin}")
  roo_data_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  cbShapeBkg.plotOn(plot, ROOT.RooFit.Name("model"))
  cbShapeBkg.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
  plot.Write()
  hSigmaCBBkg.SetBinContent(iPtBin, sigmaDCAxyBkg.getVal())
  hSigmaCBBkg.SetBinError(iPtBin, sigmaDCAxyBkg.getError())
  hAlphaCBBkg.SetBinContent(iPtBin, alphaBkg.getVal())
  hAlphaCBBkg.SetBinError(iPtBin, alphaBkg.getError())
  hNCBBkg.SetBinContent(iPtBin, nBkg.getVal())
  hNCBBkg.SetBinError(iPtBin, nBkg.getError())


## template of primary He3 from MC
sigmaDCAxy = ROOT.RooRealVar("sigmaDCAxy", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-6, 1.e-2, "cm")
alpha = ROOT.RooRealVar("alpha", "#alpha", 2, 0.5, 5)
n = ROOT.RooRealVar("n", "n", 2, 1., 10)
cbShape = ROOT.RooCrystalBall("cbShape", "cbShape", dcaxy, muDCAxy, sigmaDCAxy, alpha, n, True)
hSigmaCB = ROOT.TH1F("hSigmaCB", ";#it{p}_{T} (GeV/#it{c});#sigma_{CB} (cm)", len(pt_bins) - 1, pt_bins)
hAlphaCB = ROOT.TH1F("hAlphaCB", ";#it{p}_{T} (GeV/#it{c});#alpha_{CB}", len(pt_bins) - 1, pt_bins)
hNCB = ROOT.TH1F("hNCB", ";#it{p}_{T} (GeV/#it{c});n_{CB}", len(pt_bins) - 1, pt_bins)

## template of secondaries He3 from H3L MC
muDCAxyH3L = ROOT.RooRealVar("muDCAxyH3L", "#mu_{DCA_{xy}}", 0., -0.02, 0.02, "cm")
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
  cbShape.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
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
  cbShapeH3L.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
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

gaus_sigma = ROOT.RooRealVar("gaus_sigma", "#sigma_{Gaus}", 1e-4, 1e-4, 1e-2)
gaus_reso = ROOT.RooGaussian("gaus_reso", "gaus_reso", dcaxy, muDCAxy, gaus_sigma)

hGausSigma = ROOT.TH1F("hGausSigma", ";#it{p}_{T} (GeV/#it{c});#sigma_{Gaus} (cm)", len(pt_bins) - 1, pt_bins)
hPrimaryFrac = ROOT.TH1F("hPrimaryFrac", ";#it{p}_{T} (GeV/#it{c});Primary fraction", len(pt_bins) - 1, pt_bins)
hSecondaryFrac = ROOT.TH1F("hSecondaryFrac", ";#it{p}_{T} (GeV/#it{c});Secondary fraction", len(pt_bins) - 1, pt_bins)

for iPtBin in range(1, hPtDCAxy.GetXaxis().GetNbins() + 1):
  hSlice = hPtDCAxy.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
  roo_data_dcaxy = ROOT.RooDataHist(f"roo_data_dcaxy_{iPtBin}", f"roo_data_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
  
  ## fix the MC primary template
  sigmaDCAxy.setVal(hSigmaCB.GetBinContent(iPtBin))
  sigmaDCAxy.setConstant(True)
  alpha.setVal(hAlphaCB.GetBinContent(iPtBin))
  alpha.setConstant(True)
  n.setVal(hNCB.GetBinContent(iPtBin))
  n.setConstant(True)
  ## fix the MC secondary template
  sigmaDCAxyH3L.setVal(hSigmaCBH3L.GetBinContent(iPtBin))
  sigmaDCAxyH3L.setConstant(True)
  alphaH3L.setVal(hAlphaCBH3L.GetBinContent(iPtBin))
  alphaH3L.setConstant(True)
  nH3L.setVal(hNCBH3L.GetBinContent(iPtBin))
  nH3L.setConstant(True)
  ## fix the background template
  sigmaDCAxyBkg.setVal(hSigmaCBBkg.GetBinContent(iPtBin))
  sigmaDCAxyBkg.setConstant(True)
  alphaBkg.setVal(hAlphaCBBkg.GetBinContent(iPtBin))
  alphaBkg.setConstant(True)
  nBkg.setVal(hNCBBkg.GetBinContent(iPtBin))
  nBkg.setConstant(True)

  frac_sig = ROOT.RooRealVar("frac_sig", "frac_sig", 0.5, 0, 1)
  ## fix sig fraction to the purity
  frac_sig.setVal(hPurity.GetBinContent(iPtBin))
  frac_sig.setConstant(True)
  
  fracPrim = ROOT.RooRealVar("fracPrim", "fracPrim", 0.8, 0, 1)
  convP = ROOT.RooFFTConvPdf(f"convP_{iPtBin}", f"convP_{iPtBin}", dcaxy, cbShape, gaus_reso)
  convS = ROOT.RooFFTConvPdf(f"convS_{iPtBin}", f"convS_{iPtBin}", dcaxy, cbShapeH3L, gaus_reso)
  signal_pdf = ROOT.RooAddPdf(f"signal_pdf_{iPtBin}", "signal_pdf", ROOT.RooArgList(convP, convS), ROOT.RooArgList(fracPrim))
  
  full_pdf = ROOT.RooAddPdf(f"full_pdf_{iPtBin}", "full_pdf", ROOT.RooArgList(signal_pdf, cbShapeBkg), ROOT.RooArgList(frac_sig))
  full_pdf.fitTo(roo_data_dcaxy)

  plot = dcaxy.frame()
  plot.SetName(f"frame_dcaxy_data_{iPtBin}")
  roo_data_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  full_pdf.plotOn(plot, ROOT.RooFit.Name("model"))
  ## plot the components
  full_pdf.plotOn(plot, ROOT.RooFit.Components(f"convP_{iPtBin}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("prim"))
  full_pdf.plotOn(plot, ROOT.RooFit.Components(f"convS_{iPtBin}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("sec"))
  full_pdf.plotOn(plot, ROOT.RooFit.Components(f"cbShapeBkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("bkg"))
  full_pdf.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
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

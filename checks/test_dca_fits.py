import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("../inc/RooGausDExp.cxx+")


ffile = ROOT.TFile("../results/prepare_data_he3.root", "READ")
hPrim = ffile.Get("hPtDCAxyMCLambdab")
hPrim.SetDirectory(0)


dcaxy = ROOT.RooRealVar("dcaxy", "DCA_{xy}", 0., -0.04, 0.04, "cm")
muDCAxy = ROOT.RooRealVar("muDCAxy", "#mu_{DCA_{xy}}", 0., -0.04, 0.04, "cm")

## Double sided Crystal Ball parameters
sigmaDCAxy = ROOT.RooRealVar("sigmaDCAxy", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-6, 1.e-2, "cm")
alpha = ROOT.RooRealVar("alpha", "#alpha", 2, 0.5, 7)
n = ROOT.RooRealVar("n", "n", 2, .5, 20)
cbShape = ROOT.RooCrystalBall("cbShape", "cbShape", dcaxy, muDCAxy, sigmaDCAxy, alpha, n, True)

## double gaussian parameters
sigmaDCAxy1 = ROOT.RooRealVar("sigmaDCAxy1", "#sigma_{DCA_{xy}}^{1}", 1.e-3, 1.e-4, 1.e-1, "cm")
sigmaDCAxy2 = ROOT.RooRealVar("sigmaDCAxy2", "#sigma_{DCA_{xy}}^{2}", 3.e-3, 1.e-4, 1.e-1, "cm")
gausShape1 = ROOT.RooGaussian("gausShape1", "gausShape1", dcaxy, muDCAxy, sigmaDCAxy1)
gausShape2 = ROOT.RooGaussian("gausShape2", "gausShape2", dcaxy, muDCAxy, sigmaDCAxy2)
fracGaus = ROOT.RooRealVar("fracGaus", "f_{gaus}", 0.9, 0., 1.)
dGausShape = ROOT.RooAddPdf("gausShape", "gausShape", ROOT.RooArgList(gausShape1, gausShape2), ROOT.RooArgList(fracGaus))

## roogausdexp
alphaR = ROOT.RooRealVar("alphaR", "#alpha_{R}", 1.5, 0.01, 10.)
alphaL = ROOT.RooFormulaVar("alphaL", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaR))
gausExpMCxy = ROOT.RooGausDExp("gausExpMCxy", "GausDExp", dcaxy, muDCAxy, sigmaDCAxy, alphaL, alphaR)

## cbshape + gausShape
fracGausCB = ROOT.RooRealVar("fracGausCB", "f_{gausCB}", 0.5, 0., 1.)
cbShapeGaus = ROOT.RooAddPdf("cbShapeGaus", "cbShapeGaus", ROOT.RooArgList(cbShape, gausShape1), ROOT.RooArgList(fracGausCB))

outfile = ROOT.TFile("fit_results.root", "RECREATE")
outfile.mkdir("cb_fits")
outfile.mkdir("dGaus_fits")
outfile.mkdir("gausExp_fits")
outfile.mkdir("cbShapeGaus_fits")

for iPtBin in range(1, hPrim.GetXaxis().GetNbins() + 1):
    hSlice = hPrim.ProjectionY(f"slice_{iPtBin}", iPtBin, iPtBin)
    roo_mc_dcaxy = ROOT.RooDataHist(f"roo_mc_dcaxy_{iPtBin}", f"roo_mc_dcaxy_{iPtBin}", ROOT.RooArgList(dcaxy), hSlice)
    cbShape.fitTo(roo_mc_dcaxy)
    plot = dcaxy.frame()
    plot.SetName(f"frame_dcaxy_{iPtBin}")
    roo_mc_dcaxy.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
    cbShape.plotOn(plot, ROOT.RooFit.Name("model"))
    cbShape.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
    outfile.cd("cb_fits")
    plot.Write()
    dGausShape.fitTo(roo_mc_dcaxy)
    plotGaus = dcaxy.frame()
    plotGaus.SetName(f"frame_dcaxy_dGaus_{iPtBin}")
    roo_mc_dcaxy.plotOn(plotGaus, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
    dGausShape.plotOn(plotGaus, ROOT.RooFit.Name("model"))
    dGausShape.paramOn(plotGaus, ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plotGaus.chiSquare("model", "data")}'))
    outfile.cd("dGaus_fits")
    plotGaus.Write()
    gausExpMCxy.fitTo(roo_mc_dcaxy)
    plotExp = dcaxy.frame()
    plotExp.SetName(f"frame_dcaxy_gausExp_{iPtBin}")
    roo_mc_dcaxy.plotOn(plotExp, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
    gausExpMCxy.plotOn(plotExp, ROOT.RooFit.Name("model"))
    gausExpMCxy.paramOn(plotExp, ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plotExp.chiSquare("model", "data")}'))
    outfile.cd("gausExp_fits")
    plotExp.Write()
    cbShapeGaus.fitTo(roo_mc_dcaxy)
    plotCBGaus = dcaxy.frame()
    plotCBGaus.SetName(f"frame_dcaxy_cbShapeGaus_{iPtBin}")
    roo_mc_dcaxy.plotOn(plotCBGaus, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
    cbShapeGaus.plotOn(plotCBGaus, ROOT.RooFit.Name("model"))
    cbShapeGaus.paramOn(plotCBGaus, ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plotCBGaus.chiSquare("model", "data")}'))
    outfile.cd("cbShapeGaus_fits")
    plotCBGaus.Write()
    hSlice.Delete()



outfile.Close()



import ROOT
ROOT.gROOT.SetBatch(True)

file_lb = ROOT.TFile("../results/prepare_data_he3.root")
hPtDCAxy_lb = file_lb.Get("hPtDCAxyMCLambdab")
hPtDCAxy_lb.SetDirectory(0)

file_prompt = ROOT.TFile("../results/prepare_data_he3.root")
hPtDCAxy_prompt = file_prompt.Get("hPtDCAxyMC")
hPtDCAxy_prompt.SetDirectory(0)

## get projections over x-axis and plot them in different canvas
outfile = ROOT.TFile("out_compare.root", "RECREATE")
hPtDCAxy_lb.Write()
hPtDCAxy_prompt.Write()

for iBin in range(1, hPtDCAxy_lb.GetNbinsX() + 1):
    binContent_lb = hPtDCAxy_lb.ProjectionY(f"proj_lb_{iBin}", iBin, iBin)
    binContent_prompt = hPtDCAxy_prompt.ProjectionY(f"proj_prompt_{iBin}", iBin, iBin)
 
    c = ROOT.TCanvas(f"c_{iBin}", f"c_{iBin}", 800, 600)
    binContent_lb.SetLineColor(ROOT.kRed)
    binContent_prompt.SetLineColor(ROOT.kBlue)
    binContent_prompt.DrawNormalized()
    binContent_lb.DrawNormalized("same")
    leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    leg.AddEntry(binContent_lb, "From #Lambda_{b}", "l")
    leg.AddEntry(binContent_prompt, "Prompt", "l")
    leg.Draw()
    c.Write()

outfile.Close()
    

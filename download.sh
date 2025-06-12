############################################################################################################

mkdir -p he3_data
mkdir -p he3_data/data
mkdir -p he3_data/mc
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0041/413783/110250/AO2D.root file:he3_data/data/AO2D_24_skimmed.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0042/421683/114406/AO2D.root file:he3_data/mc/AO2D_25a3.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0041/413369/110166/AO2D.root file:he3_data/mc/AO2D_25e1_lambdab.root
############################################################################################################
# mkdir -p d_data
# mkdir -p d_data/data
# mkdir -p d_data/mc
# mkdir -p d_data/data/thnsparse
# mkdir -p d_data/mc/thnsparse
# alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/373337/91935/AO2D.root file:d_data/data/AO2D_24_minbias.root
# alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/370965/90732/AO2D.root file:d_data/mc/AO2D_25a3.root
# alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/374335/92500/AnalysisResults.root file:d_data/data/thnsparse/AnalysisResults_24_minbias.root  
# alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/378130/94721/AnalysisResults.root file:d_data/mc/thnsparse/AnalysisResults_25a3.root
############################################################################################################
mkdir -p utils
alien_cp alien:/alice/cern.ch/user/f/fmazzasc/he3_from_lambdab_analysis/h3l_spectrum.root file:utils/h3l_spectrum.root
alien_cp alien:/alice/cern.ch/user/f/fmazzasc/he3_from_lambdab_analysis/he3_spectrum.root file:utils/he3_spectrum.root
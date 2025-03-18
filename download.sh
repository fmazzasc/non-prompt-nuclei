mkdir -p he3_data
mkdir -p he3_data/data
mkdir -p he3_data/mc
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0036/367418/88263/AO2D.root file:he3_data/data/AO2D_an.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0029/290850/48265/AO2D.root file:he3_data/data/AO2D_al.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0033/337353/73325/AO2D.root file:he3_data/mc/AO2D_25a3.root

mkdir -p d_data
mkdir -p d_data/data
mkdir -p d_data/mc
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/372539/91519/AO2D.root file:d_data/data/AO2D_minbias.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0037/370965/90732/AO2D.root file:d_data/mc/AO2D_25a3.root

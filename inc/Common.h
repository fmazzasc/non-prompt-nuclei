#ifndef COMMON_H
#define COMMON_H

#include "ROOT/RDataFrame.hxx"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>

using std::map;
using std::string;
using std::vector;

double bb(double bg, double kp1, double kp2, double kp3, double kp4, double kp5)
{
  double beta = bg / std::sqrt(1. + bg * bg);
  double aa = std::pow(beta, kp4);
  double bb = std::pow(1. / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

float bbHe3(float mom) { return bb(mom / 2.80839, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
float nsigmaHe3(float mom, float sig) { return (sig / bbHe3(mom * 2) - 1. + 2.20376e-02) / 0.055; }

float bbH3(float mom) { return bb(mom / 2.80892f, -136.71, 0.441, 0.2269, 1.347, 0.8035); }
float nsigmaH3(float mom, float sig) { return (sig / bbH3(mom) - 1.) / 0.07; }

float bbDeut(float mom) { return bb(mom / 1.8756129f, 11.56475, 3.86889, 0.0138162, 2.29009, 0.912842); }
float nsigmaDeu(float mom, float sig) { return (sig / bbDeut(mom) - 1.) / 0.07; }

float bbHe4(float mom) { return bb(mom / 3.72738f, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
float nsigmaHe4(float mom, float sig) { return (sig / bbHe4(mom * 2) - 1.) / 0.07; }

float DCAxyCut(float pt, float nsigma)
{
  float invPt = 1.f / pt;
  return (7.62783e-04 + 4.59326e-03 * invPt + 6.89163e-03 * invPt * invPt) * nsigma;
}
float nSigmaDCAxy(double pt, float dcaxy) {
  return dcaxy / DCAxyCut(pt, 1);
}

float DCAzCut(float pt, float nsigma)
{
  float invPt = 1.f / pt;
  return (5.00000e-04 + 8.73690e-03 * invPt + 9.62329e-04 * invPt * invPt) * nsigma;
}

float averageClusterSize(uint32_t itsClusterSizes)
{
  float sum = 0;
  int nclusters = 0;
  int max = 0;
  for (int layer = 0; layer < 7; layer++) {
    int clsize = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (clsize > 0) {
      nclusters++;
      sum += clsize;
      if (clsize > max) {
        max = clsize;
      }
    }
  }
  if (nclusters == 0) {
    return 0;
  }
  // truncated mean
  return (sum - max) / (nclusters - 1);
};

float nSigmaDCAz(double pt, float dcaz) {
  return dcaz / DCAzCut(pt, 1);
}

#endif
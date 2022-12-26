#include "utils/defs.hpp"
#include "utils/ioroutines.hpp"
#include "graphmodel/graphmodel.hpp"
#include "simulators/hidden_cascade.hpp"
#include "simulators/industry_cascade.hpp"
#include "algorithms/vacc_minimizer.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace util;
using namespace graphmodel;
using namespace simulator;

#define DETGRAPHFILE "plot_output_det.csv"
#define VACGRAPHFILE "plot_output_vac.csv"
#define DIAGNOSTICSFILE "optimization_diagnostics.txt"
#define RAWGRAP "plot_output_raw.csv"
#define N_VARIANCEFILE "variances.csv"


void DiagnoseVariance(algorithm::InfluenceMinimizer& minc, int n, int k, std::ostream& out = std::cerr)
{
  for (double p = 0.01; p < 0.99; p = p + 0.98/k)
  {
    out << "NewData for p:" << p << "\n"; 
    minc.SetConstantProb(p);
    minc.DiagnosePerformance(n,out,true);
    out.flush();
  }
}

int main()
{
  // Initialize a metagraph.      
  /*
  IoR zed;
  zed.SetAnnouncementDirectory(ADIR);
  zed.SetCompanyDictionaryFile(CDFILE);
  zed.SetNodeDictionaryFile(NDFILE);
  zed.ReadCompanyDictionary();
  zed.ReadNodeDictionary();
  zed.ReadAnnouncements();
  */

  std::ofstream out;
  out.open(DIAGNOSTICSFILE);

  std::cerr << "Intialize IndustryCascade for " << VREFDATE << "\n";
  IndustryCascade foo(VREFDATE);

  //std::ofstream out2;
  //out2.open(COMPANYFILE);
  //foo.PrintCompanies(out2);
  //out2.close();


  std::cerr << " Running initial diagnostics \n";
  foo.GraphDiagnostics(out);
  algorithm::InfluenceMinimizer minim(foo);
  std::cerr << " Running variance diagnostics \n";
  //std::ofstream outv;
  //outv.open(N_VARIANCEFILE);
  auto start = std::chrono::high_resolution_clock::now(); 
  //DiagnoseVariance(minim, 120, 20, outv);
  //outv.close();
  auto stop = std::chrono::high_resolution_clock::now(); 
  double count = std::chrono::duration<double>(stop-start).count();
  std::cerr << "Diagnostics took " << count << "s \n";

  std::cerr << "Running singleton influence check";
  

  // Calclulate the optimimus with greefy algorithm for n = 1 2 3 4 5 6 7 8 9 10
  start = std::chrono::high_resolution_clock::now();
  foo.SetConstantProb(0.15);
  std::ofstream outd;
  outd.open(DETGRAPHFILE);
  algorithm::InfluenceMinimizer inf(foo);
  auto refz = inf.DefaultInfluence(20);
  outd << "0 " << refz.first << " " << refz.second <<  "\n";
  outd.flush();
  IndustryCascade fooz(foo);
  for (int n = 1; n < 21; ++n)
  {
    algorithm::InfluenceMinimizer inf(fooz);
    auto minz = inf.FindMinimalNodeComp(20);
    outd << n << " " << std::get<2>(minz) << " " << std::get<3>(minz) << "\n";
    outd.flush();
    int node = std::get<0>(minz);
    int comp = std::get<1>(minz);
    fooz.DeactivateFromInside(node,comp);
  }
  stop = std::chrono::high_resolution_clock::now(); 
  count = std::chrono::duration<double>(stop-start).count();
  std::cerr << "20 optimization points took " << count << "s \n";

  return 0;
}
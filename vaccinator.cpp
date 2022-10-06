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

#define RAWGRAPHFILE "plot_output_raw.csv"
#define DIAGNOSTICSFILE "optimization_diagnostics.txt"


void GenerateRawGraphs(IndustryCascade& ind, const double& p0, const double& pn, const int& n)
{
  algorithm::InfluenceMinimizer minim(ind);    
  
  const std::string filename = RAWGRAPHFILE;
  std::ofstream out;
  out.open(filename);
  const std::vector<double> weights = ind.GetSimulationWeights();
  std::cerr << "Starting the simulation cycles \n";
  for (double p = p0; p <= pn; p += (pn -p0)/n)
  {
    ind.SetConstantProb(p);
    auto ref_inf = ind.RunTotal(weights);
    std::cerr << "Ref: " << ref_inf.first << " std: " << sqrt(ref_inf.second) << "\n";
    auto res = minim.FindMinimalNodeComp();
    out << p << "," << ref_inf.first << "," <<  sqrt(ref_inf.second) << "," << std::get<2>(res) << "," << std::get<3>(res) << "\n";
    out.flush();
  }
  out.close();
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

  std::cerr << " Running initial diagnostics \n";
  foo.GraphDiagnostics(out);

  std::cerr << "Running singleton influence check";

  
  // Start timing.
  auto start = std::chrono::high_resolution_clock::now(); 
  foo.SetConstantProb(0.2);
  algorithm::InfluenceMinimizer minim(foo);

  minim.DiagnoseBetweenMinimal(out);

  auto stop = std::chrono::high_resolution_clock::now(); 
  double count = std::chrono::duration<double>(stop-start).count();
  std::cerr << "Singleton influence check took " << count << "s \n";
  return 0;
}
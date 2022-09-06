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
  std::cerr << "Intialize IndustryCascade for " << VREFDATE << "\n";
  IndustryCascade foo(VREFDATE);
  std::cerr << "all created.\nStarting minimization procedure.";

  // Start timing.
  auto start = std::chrono::high_resolution_clock::now(); 
  algorithm::InfluenceMinimizer minim(foo);    
  
  const std::string filename = "plot_output.csv";
  std::ofstream out;
  out.open(filename);
  const std::vector<double> weights = foo.GetSimulationWeights();
  for (double p = 0.05; p < 1; p += 0.05)
  {
    foo.SetConstantProb(p);
    auto ref_inf = foo.RunTotal(weights);
    std::cerr << "Ref: " << ref_inf.first << " std: " << sqrt(ref_inf.second) << "\n";
    auto res = minim.FindMinimalNodeComp();
    out << p << "," << ref_inf.first << "," <<  sqrt(ref_inf.second) << "," << std::get<2>(res) << "," << std::get<3>(res) << "\n";
    out.flush();
  }
  out.close();

  auto stop = std::chrono::high_resolution_clock::now(); 
  double count = std::chrono::duration<double>(stop-start).count();
  std::cerr << "Minimisation took " << count << "s \n";
  return 0;
}
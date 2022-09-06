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


double sum(const std::vector<double>& input)
{
  return std::accumulate(input.cbegin(), input.cend(), 0.0);
}

double avg(const std::vector<double>& input)
{
  return sum(input) / input.size();
}

double w_avg(const std::vector<double>& input, const std::vector<double>& weights)
{
  BOOST_ASSERT(input.size() == weights.size());
  double ws = 0;
  double div = 0;
  for (size_t i = 0; i < input.size(); ++i) 
  {
    ws += input[i]*weights[i];
    div += weights[i];
  }
  BOOST_ASSERT(div > 0);
  return ws/div;
}

double st_error(std::vector<double> input)
{
  const int& n = static_cast<int>(input.size());
  const double& avgs = avg(input);
  double err = 0.0;
  #pragma omp parallel for reduction(+:err)
  for(int i = 0; i <= n-1;++i)
  {
      err += (avgs - input[i])*(avgs - input[i]);
  }
  err /= (n-1);
  return std::sqrt(err); 
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
    double ref_inf = foo.RunTotal(weights);
    std::cerr << "Ref: " << ref_inf << "\n";
    auto res = minim.FindMinimalNodeComp();
    out << p << "," << ref_inf << "," << std::get<2>(res) << "\n";
    out.flush();
  }
  out.close();

  auto stop = std::chrono::high_resolution_clock::now(); 
  double count = std::chrono::duration<double>(stop-start).count();
  std::cerr << "Minimisation took " << count << "s \n";
  return 0;
}
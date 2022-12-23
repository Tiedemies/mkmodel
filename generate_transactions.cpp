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

#define ALTERN_GRAPH "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/Synthetics/generated_graph.csv"

int main()
{

  std::cerr << "Intialize IndustryCascade for " << GREFDATE << "\n";
  IndustryCascade foo(GREFDATE);

  std::cerr << "Read alternate graph";
  foo.ReadExternalGraph(ALTERN_GRAPH);
  for (auto tr: foo.RunGeneration())
  {
    std::cout << foo.to_string(tr) << "\n";
    std::cout.flush();
  }

  return 0;
}
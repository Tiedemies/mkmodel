#ifndef VACCMIN
#define VACCMIN
#include "../simulators/hidden_cascade.hpp"

#include <chrono>
#include <set>
#include <vector>
#include <map>
#include <random>

namespace algorithm
{
  class InfluenceMinimizer
  {
    public:
      // Constructor
      InfluenceMinimizer(simulator::IndustryCascade& ind);
      // Copy constructor
      InfluenceMinimizer(const InfluenceMinimizer& rhs);
      // Destructor 
      ~InfluenceMinimizer();
      // Find a minimal set via greedy algorithm of size n 
      std::set<int> GreedyMinimalSet(int n);
      // Find single minimal node. 
      std::tuple<int,int,double,double> FindMinimalNodeComp();
      std::set<std::tuple<int,int>> MinizeDiffBetween(int n);
    private:
      simulator::IndustryCascade& ind_; 
  };
}
// End
#endif
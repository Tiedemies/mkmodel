#ifndef VACCMIN
#define VACCMIN
#include "../simulators/hidden_cascade.hpp"
#include "../simulators/industry_cascade.hpp"
#include "../simulators/announcement_cascade.hpp"

#include <chrono>
#include <set>
#include <vector>
#include <map>
#include <random>
#include <iostream>

namespace algorithm
{
  class InfluenceMinimizer
  {
    public:
      // Constructor
      InfluenceMinimizer(simulator::IndustryCascade& ind);
      InfluenceMinimizer(simulator::AnnouncementCascade& anc);
      // Copy constructor
      InfluenceMinimizer(const InfluenceMinimizer& rhs);
      // Destructor 
      ~InfluenceMinimizer();
      // Find a minimal set via greedy algorithm of size n 
      std::set<int> GreedyMinimalSet(int n);
      // Find single minimal node. 
      std::tuple<int,int,double,double> FindMinimalNodeComp(int n = 4);
      std::tuple<int,double,double> FindMinimalNode(int n = 4);
      std::pair<double,double> DefaultInfluence(int n = 4);
      std::set<std::tuple<int,int>> MinizeDiffBetween(int n);

      void DiagnoseBetweenMinimal(std::ostream& = std::cerr);

      // Run a variance diagnostic up to n iterations      
      void DiagnosePerformance(int n, std::ostream& = std::cerr, bool anti = true);

      void SetConstantProb(double p);

    private:
      // simulator::IndustryCascade& ind_;
      simulator::AnnouncementCascade anc_; 
  };
}
// End
#endif
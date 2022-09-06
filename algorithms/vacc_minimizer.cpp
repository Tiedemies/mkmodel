// Influence Minimizer class implementation
#include "vacc_minimizer.hpp"


#include <chrono>
#include <set>
#include <vector>
#include <map>
#include <random>
#include <boost/assert.hpp>
#include "../simulators/industry_cascade.hpp"
#include "../simulators/hidden_cascade.hpp"

// Open namespace
namespace algorithm
{
  /* Constructor */
  InfluenceMinimizer::InfluenceMinimizer(simulator::IndustryCascade& ind) : ind_(ind)
  {
    BOOST_ASSERT(true); 
    // void
  }
  
  /* Destructor */
  InfluenceMinimizer::~InfluenceMinimizer()
  {
    BOOST_ASSERT(true);
    // void
  }

  /* Find a minimal node */
  std::tuple<int,int,double,double>
  InfluenceMinimizer::FindMinimalNodeComp()
  {
    // Precondition
    BOOST_ASSERT(ind_.hc_.GetN() > 0);
    
    int min_node = -1;
    int min_comp = -1;
    int max_node = -1;
    int max_comp = -1;
    double min_influence = std::numeric_limits<double>::max();
    double max_influence = std::numeric_limits<double>::min();
    double min_var = 0.0;
    /* Find the node */
    int simcount = 0;
    const auto weights = ind_.GetSimulationWeights();
    for (int comp = 0; comp < ind_.n_comp_;++comp)
    {
      for (auto node_it = ind_.insiders_.at(comp).begin(); node_it != ind_.insiders_.at(comp).end();++node_it)
      {
        int node = *node_it; 
        ind_.DeactivateFromInside(node,comp);
        auto res_pair = ind_.RunTotal(weights);
        const double& inf = res_pair.first;
        ind_.ReactivateInside();
        if (inf < min_influence)
        {
          min_influence = inf;
          min_node = node;
          min_comp = comp; 
          min_var = res_pair.second;
        }
        if (inf > max_influence)
        {
          max_influence = inf;
          max_node = node;
          max_comp = comp; 
        }
        std::cerr << ++simcount << " pairs tried \n" << "Current min: " << min_influence << "\n";
      }
     
    }

    // Postcondition
    BOOST_ASSERT(max_node > 0);
    BOOST_ASSERT(min_node > 0);
    BOOST_ASSERT(max_comp > 0);
    BOOST_ASSERT(min_comp > 0);
    std::cerr << "max influence: " << max_influence << " vs min " << min_influence <<"\n";
    return std::make_tuple(min_node, min_comp, min_influence, min_var);
  }



}
// DoneInfluenceMinimizer
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
        #pragma omp critical
        {
        if (inf < min_influence)
        {
          min_influence = inf;
          min_node = node;
          min_comp = comp; 
          min_var = res_pair.second;
        }
        }

        #pragma omp critical
        {
        if (inf > max_influence)
        {
          max_influence = inf;
          max_node = node;
          max_comp = comp; 
        }
        }
        #pragma omp atomic
        ++simcount; 

        std::cerr << simcount << " pairs tried \n" << "Current min: " << min_influence << "\n";
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

  /* Find a minimal node/comp using betweenness-heuristics  */
  std::set<std::tuple<int,int>> 
  InfluenceMinimizer::MinizeDiffBetween(int n)
  {
    // Precondition
    BOOST_ASSERT(ind_.hc_.GetN() > 0);
    const auto weights = ind_.GetSimulationWeights();
    std::set<std::tuple<int,int>> togo;

    // TODO ACTUAL IMPLEMENTATION 

    //Postcondition
    BOOST_ASSERT(togo.size() == n);
    return togo;
  }


  // DIagnostic minimization procedure
  void
  InfluenceMinimizer::DiagnoseBetweenMinimal(std::ostream& out)
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
    int simcount = 0;
    const auto weights = ind_.GetSimulationWeights();
    const auto& c_cent = ind_.GetCompCentrality();
    const auto& n_cent = ind_.GetNodeCentrality();
    /* First run the initial */
    ind_.hc_.gather_statistic_ = true;
    auto res_p_init = ind_.RunTotal();
    ind_.hc_.gather_statistic_ = false;

    const auto n_c_map = ind_.GetNeighbourActivation();
    
    for (int comp = 0; comp < ind_.n_comp_;++comp)
    {
      
      for (auto node_it = ind_.insiders_.at(comp).begin(); node_it != ind_.insiders_.at(comp).end();++node_it)
      {
        int node = *node_it; 
        ind_.DeactivateFromInside(node,comp);
        auto res_pair = ind_.RunTotal(weights);
        const double& inf = res_pair.first;
        const double& stdv = res_pair.second;
        out << comp << "," << node << "," << inf << "," << stdv << "," << c_cent.at(comp) << "," << n_cent.at(node) 
            << "," << ind_.comp_adj_.at(comp).size() << "," << ind_.hc_.adj_.at(node).size() << ","  
            << ind_.insiders_.at(comp).size() << ","<< ind_.foo_.GetGraph(ind_.date_)->GetInsiderOf(node).size() 
            << "," << n_c_map.at(ind_.hc_.Key(node,comp)) << "\n";
        ind_.ReactivateInside();
      }
    }
    // Postcondition
    BOOST_ASSERT(max_node > 0);
    BOOST_ASSERT(min_node > 0);
    BOOST_ASSERT(max_comp > 0);
    BOOST_ASSERT(min_comp > 0);
  }

}
// DoneInfluenceMinimizer
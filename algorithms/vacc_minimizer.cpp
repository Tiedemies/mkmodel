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
  InfluenceMinimizer::InfluenceMinimizer(simulator::IndustryCascade& ind) : anc_(simulator::AnnouncementCascade(ind))
  {
    BOOST_ASSERT(true); 
    // void
  }

  InfluenceMinimizer::InfluenceMinimizer(simulator::AnnouncementCascade& anc) : anc_(anc)
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
  InfluenceMinimizer::FindMinimalNodeComp(int n)
  {
    // Precondition
    BOOST_ASSERT(anc_.ic_.hc_.GetN() > 0);
    int min_node = -1;
    int min_comp = -1;
    double min_influence = std::numeric_limits<double>::max();
    double min_var = 0.0;
    /* Find the node */
    int simcount = 0;
    for (int comp = 0; comp < anc_.ic_.n_comp_;++comp)
    {
      for (auto node_it = anc_.ic_.insiders_.at(comp).begin(); node_it != anc_.ic_.insiders_.at(comp).end();++node_it)
      {
        int node = *node_it; 
        if (anc_.ic_.inside_of_.at(node).size() < 2)
        {
          continue;
        }
        std::vector<double> c_vec_anti(n,0.0);
        anc_.ic_.DeactivateFromInside(node,comp);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
          c_vec_anti[i] = anc_.RunSingleCascade(true);
        }
        anc_.ic_.ReactivateInside(node, comp);
        double inf = util::avg(c_vec_anti);
        if (inf < min_influence)
        {
          min_influence = inf;
          min_node = node;
          min_comp = comp; 
          min_var = util::st_error(c_vec_anti);
        }     
        ++simcount; 

        //std::cerr << simcount << " pairs tried \n" << "Current min: " << min_influence << "\n";
      }
     
    }

    // Postcondition
 
    BOOST_ASSERT(min_node > 0);
    BOOST_ASSERT(min_comp > 0);
    // std::cerr << "min " << min_influence <<" \n";
    return std::make_tuple(min_node, min_comp, min_influence, min_var);
  }

  /* Find a minimal node/comp using betweenness-heuristics  */
  std::set<std::tuple<int,int>> 
  InfluenceMinimizer::MinizeDiffBetween(int n)
  {
    // Precondition
    BOOST_ASSERT(anc_.ic_.hc_.GetN() > 0);
    const auto weights = anc_.ic_.GetSimulationWeights();
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
    BOOST_ASSERT(anc_.ic_.hc_.GetN() > 0);
    
    const auto weights = anc_.ic_.GetSimulationWeights();
    const auto& c_cent = anc_.ic_.GetCompCentrality();
    const auto& n_cent = anc_.ic_.GetNodeCentrality();
    /* First run the initial */
    anc_.ic_.hc_.gather_statistic_ = true;
    const auto& n_c_map = anc_.ic_.GetNeighbourOutActivation();
    const auto& n_c_map2 = anc_.ic_.GetNeighbourInActivation();
    anc_.ic_.hc_.gather_statistic_ = false;
    
    for (int comp = 0; comp < anc_.ic_.n_comp_;++comp)
    {
      
      for (auto node_it = anc_.ic_.insiders_.at(comp).begin(); node_it != anc_.ic_.insiders_.at(comp).end();++node_it)
      {
        int node = *node_it; 
        if (anc_.ic_.inside_of_.at(node).size() < 2)
        {
          continue;
        }
        anc_.ic_.DeactivateFromInside(node,comp);
        auto res_pair = anc_.ic_.RunTotal(weights);
        const double& inf = res_pair.first;
        const double& stdv = res_pair.second;
        out << comp << "," << node << "," << inf << "," << stdv << "," << c_cent.at(comp) << "," << n_cent.at(node) 
            << "," << anc_.ic_.comp_adj_.at(comp).size() << "," << anc_.ic_.hc_.adj_.at(node).size() << ","  
            << anc_.ic_.insiders_.at(comp).size() << ","<< anc_.ic_.foo_.GetGraph(anc_.ic_.date_)->GetInsiderOf(node).size() 
            << "," << n_c_map.at(anc_.ic_.hc_.Key(node,comp)) << "," << n_c_map2.at(anc_.ic_.hc_.Key(node,comp)) << "\n";
        anc_.ic_.ReactivateInside();
      }
    }
    // Postcondition
   
  }

  void 
  InfluenceMinimizer::DiagnosePerformance(int n, std::ostream& out, bool anti)
  {
    for (int nr = 1; nr < n; nr += (nr == 1)?9:10)
    {
      int nn = anti?nr:2*nr;
      std::vector<double> c_vec(24,0.0);
      #pragma omp parallel for
      for (int ik = 0; ik < nn*24; ++ik)
      {
        int k = ik%24;
        double temp = anc_.RunSingleCascade(anti)/(double)nn;
        #pragma omp critical
        {
          c_vec[k] += temp;
        } 
      }
      out << nn << "," << util::avg(c_vec) << "," << util::st_error(c_vec) << "\n";
    }
  }

  std::pair<double,double> 
  InfluenceMinimizer::DefaultInfluence(int n)
  {
    std::vector<double> c_vec_anti(n,0.0);
    //#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
      c_vec_anti[i] = anc_.RunSingleCascade(true);
    }
    return std::make_pair(util::avg(c_vec_anti), util::st_error(c_vec_anti));
  }

  void 
  InfluenceMinimizer::SetConstantProb(double p)
  {
    anc_.ic_.SetConstantProb(p);
  }  

  std::tuple<int,double,double>
  InfluenceMinimizer::FindMinimalNode(int n)
  {
    // Precondition
    BOOST_ASSERT(anc_.ic_.hc_.GetN() > 0);
    int min_node = -1;
    double min_influence = std::numeric_limits<double>::max();
    double min_var = 0.0;
    /* Find the node */
    int simcount = 0;

    for (int node = 0; node < anc_.ic_.inside_of_.size(); ++node)
    {  
      // Skip vaccinating those that are not inside in more than 1. 
      if (anc_.ic_.inside_of_.at(node).size() < 2 || anc_.ic_.hc_.is_disabled_.at(node))
      {
        continue;
      }
      anc_.ic_.DeactivateNode(node);
      std::vector<double> c_vec_anti(n,0.0);
      #pragma omp parallel for
      for (int i = 0; i < n; ++i)
      {
        c_vec_anti[i] = anc_.RunSingleCascade(true);
      }
      anc_.ic_.ReactivateInside(node);
      double inf = util::avg(c_vec_anti);
      if (inf < min_influence)
      {
        min_influence = inf;
        min_node = node;
        min_var = util::st_error(c_vec_anti);
      }     
      ++simcount; 
        //std::cerr << simcount << " pairs tried \n" << "Current min: " << min_influence << "\n";
    }
    // Postcondition
    BOOST_ASSERT(min_node > 0);
    // std::cerr << "min " << min_influence <<" \n";
    return std::make_tuple(min_node, min_influence, min_var);
  }

}
// DoneInfluenceMinimizer
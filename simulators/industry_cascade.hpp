#ifndef I_CASCADE
#define I_CASCADE
#include "hidden_cascade.hpp"
#include "../utils/defs.hpp"
#include "../utils/h_random.hpp"
#include "../utils/ioroutines.hpp"

#include "../graphmodel/graphmodel.hpp"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>
#include<list>

namespace algorithm
{
  class InfluenceMinimizer;
}

namespace simulator
{
  using namespace util;
  using namespace graphmodel;
/* 
 * The IndustryCascade class is a wrapper for a hidden cascade that also includes the information 
 * about announcement numbers and insider sets, and we can run the whole system through them.
 */
class IndustryCascade
{
  public:
    IndustryCascade(HiddenCascade h);
    IndustryCascade(int date, int days = 365);
    // Copy constructor:
    IndustryCascade(const IndustryCascade& rhs);
    ~IndustryCascade();
    std::pair<double,double> RunTotal();
    std::pair<double,double> RunTotal(const std::vector<double>& weights);
    std::vector<double> RunSingle(int i);
    void RandomizeCoefs(double mu);
    std::vector<double> GetSimulationWeights();
    void SetSimulationCoef(int n);
    void SetConstantProb(double p);
    void RunDiagnostics(); 
    void GraphDiagnostics(std::ostream& = std::cerr);
    void DeactivateFromInside(int node, int comp);
    void ReactivateInside(int node, int comp);
    void ReactivateInside(int node);
    void ReactivateInside();
    const double& GetTreeRank(int node);
    const double& GetTreeDiff(int node, int comp);
    const std::unordered_map<size_t,double>& GetNeighbourInActivation();
    const std::unordered_map<size_t,double>& GetNeighbourOutActivation();
    const std::unordered_map<size_t,double>& GetRecursiveNeighbourActivation();
    bool IsCached(int node, int comp);
    const std::vector<double>& GetNodeCentrality();
    const std::vector<double>& GetCompCentrality(); 

   
    struct Transaction
    {
      int node_;
      int comp_;
      int date_;
      double price_;
      double vol_;
      bool inside_;
      bool operator<(const Transaction& rhs) const;
    };
    // Generate transactions
    std::vector<Transaction> RunGeneration();

    // Check if the prices is rising (true) or dropping (false). 
    bool IsBear(int comp, int day);
    friend class algorithm::InfluenceMinimizer; 

  private:
    void InitializeFoo();
    void InitializeAnnouncements();
    void InitializeWindows();
    void InitializeTransactions();
    void InitializeProbabilities(); 
    void InitializePrices(); 

    void InitializeCompAdj();
    int simulation_coefficient_; 
    int window_size_;
    int date_; 
    int days_;
    int days_begin_;
    int days_end_;
    int n_node_;
    int n_comp_; 
    double def_p_;
    IoR foo_; 
    HiddenCascade hc_;
    Random rnd_;
    std::vector<std::vector<int>> insiders_;
    std::vector<std::vector<int>> inside_of_;
    std::vector<int> num_announcements_;
    std::vector<std::vector<int>> announcement_days_; 
    std::vector<std::vector<bool>> in_window_days_;
    std::vector<std::vector<bool>> trading_days_; 
    std::vector<std::vector<double>> prices_; 
    std::vector<Transaction> transacts_;
    std::vector<std::vector<double>> out_trading_prob_;
    std::vector<std::vector<double>> in_trading_prob_;
    std::vector<int> n_in_days_;
    std::vector<int> n_out_days_;
    std::vector<double> tree_rank_; 
    std::unordered_map<size_t, double> neighbour_out_activation_; 
    std::unordered_map<size_t, double> neighbour_in_activation_; 

    struct DisableCacheEntry
    {
      int node_; 
      int comp_;
      DisableCacheEntry(int node, int comp): node_(node), comp_(comp) {};
      bool operator==(const DisableCacheEntry& rhs) const { return rhs.node_ == node_ && rhs.comp_ == comp_; };  
    };

    std::list<DisableCacheEntry> disable_cache_; 
    std::vector<double> comp_centrality_; 
    std::vector<std::vector<int>> comp_adj_; 
};

// Helper functions
double InterPolate(const std::vector<double>& vec, int i);
}
#endif
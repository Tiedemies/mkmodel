
#ifndef H_CASCADE
#define H_CASCADE

#include "../utils/h_random.hpp"
#include "../utils/defs.hpp"
#include "../graphmodel/graphmodel.hpp"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>

// Forward declaration
namespace algorithm
{
  class InfluenceMinimizer;
}

namespace simulator
{
    using namespace graphmodel;
    using namespace util;
class HiddenCascade
{
    public:
    // 
    // Constructor and destructor
        HiddenCascade() = delete;
        HiddenCascade(const MonoGraph* mg, const double& p, const double& fp, const double& tp);
        ~HiddenCascade();
    // Simulate. 
        double Simulate(const std::vector<int>& inside);
    // Enhanced Simulate. 
        double SimulateEnc(const std::vector<int>& inside);
    // Change parameters: 
        void SetTradingProp(const std::vector<double>& props);
        void SetExpectedVolumes(const std::vector<double>& vols);
        void Disable(int node);
        void SetSimulationN(size_t n);
    // Generate a datapoint:
        std::vector<double> Generate(const std::vector<int>& inside, bool window_day);
    // Test for success of infomation transfer between u and v:
        bool IsSuccess(const int& u, const int& v) const;
    // Maximum  and  minimum number of activated nodes in a sim
        int GetMaxActivated() const;
        int GetMinActivated() const;
    // Deactivate or "vaccinate" a given node:
        void DeactivateNode(int u);
    // Re-activate node:
        void ReactivateNode(int u);
    // Deactivate a given connection from u to v
        void DeactivateConnection(int u, int v);
    // Re-acticate the connection, given a probability (default probability if not given)
        void ActivateConnection(int u, int v, double p = -1);
    // Randomize connection weights to expected value/deviation. 
        void RandomizeWeights(double mu, double sigma = -1);
        void SetConstantProb(double p);
    // Get number of nodes
        int GetN() const;
    // Query the latest variation of nodes 
        double LastVar() const;

        friend class IndustryCascade; 
        friend class algorithm::InfluenceMinimizer;
  
    private:
        inline size_t Key(const int& u, const int& v) const;
        size_t sim_n_; 
        std::vector<double> false_positive_prob_;
        std::vector<double> true_positive_prob_;
        std::unordered_map<size_t, double> p_map_;
        std::vector<std::vector<int>> adj_; 
        
        double uni_p_; 

        mutable Random rnd_; 
        std::vector<bool> is_disabled_;
        std::set<int> disabled_; 
        std::vector<double> simulated_activations_; 
        std::vector<double> volume_; 
        double num_activated_;
        int max_activated_;
        int min_activated_; 

        std::vector<int> infected_by_sim_; 

};
}

#endif
 
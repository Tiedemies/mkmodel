
#ifndef H_CASCADE
#define H_CASCADE
#include "defs.hpp"
#include "graphmodel.hpp"
#include "h_random.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>

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
    // Change parameters: 
        void SetTradingProp(const std::vector<double>& props);
        void SetExpectedVolumes(const std::vector<double>& vols);
        void Disable(int node);
        void SetSimulationN(size_t n);
    // Generate a datapoint:
        std::vector<double> Generate(const std::vector<int>& inside, bool window_day);
    // Test for success of infomation transfer between u and v:
        bool IsSuccess(const int& u, const int& v) const;
        int GetMaxActivated() const;
        int GetMinActivated() const;

        void DeactivateConnection(int u, int v);
        void ActivateConnection(int u, int v, double p = -1);
  
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
        std::vector<double> simulated_profits_;
        std::vector<double> simulated_activations_; 
        std::vector<double> trading_prob_;
        std::vector<double> volume_; 
        double num_activated_;
        int max_activated_;
        int min_activated_;  

};

#endif
 
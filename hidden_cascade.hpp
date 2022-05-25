
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
        HiddenCascade(const MonoGraph* mg, const double& p, const double& fp, const double& tp);
        ~HiddenCascade();
        void Simulate(const std::vector<int>& inside, bool simulate_profits);
        void Disable(int node);
        void SetSimulationN(size_t n);
        bool IsSuccess(const int& u, const int& v) const;
    private:
        inline size_t Key(const int& u, const int& v) const;
        size_t sim_n_; 
        std::vector<double> false_positive_prob_;
        std::vector<double> true_positive_prob_;
        std::unordered_map<size_t, double> p_map_;
        std::unordered_map<int,Alist2> adj_;
        mutable Random rnd_; 
        std::vector<bool> is_disabled_;
        std::set<int> disabled_; 
        std::vector<double> simulated_profits_;
        std::vector<double> simulated_activations_; 
};

#endif
 
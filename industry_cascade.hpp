#ifndef I_CASCADE
#define I_CASCADE
#include "hidden_cascade.hpp"
#include "defs.hpp"
#include "graphmodel.hpp"
#include "h_random.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>

/* 
 * The IndustryCascade class is a wrapper for a hidden cascade that also includes the information 
 * about announcement numbers and insider sets, and we can run the whole system through them.
 */
class IndustryCascade
{
    public:
        IndustryCascade(HiddenCascade h);
        ~IndustryCascade();
        std::vector<double> RunTotal();
        std::vector<double> RunSingle(int i);
        void SetSimulationCoef(int n);

    private:
        std::vector<std::vector<int>> insiders_;
        std::vector<int> num_announcements_; 
        HiddenCascade& hc_; 
        int simulation_coefficient_; 
};



#endif
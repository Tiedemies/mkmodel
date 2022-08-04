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

/* 
 * The IndustryCascade class is a wrapper for a hidden cascade that also includes the information 
 * about announcement numbers and insider sets, and we can run the whole system through them.
 */
class IndustryCascade
{
    public:
        IndustryCascade(HiddenCascade h);
        IndustryCascade(int year);
        ~IndustryCascade();
        std::vector<double> RunTotal();
        std::vector<double> RunSingle(int i);
        void SetSimulationCoef(int n);

    private:
        int simulation_coefficient_; 
        double def_p_;
        IoR foo_; 
        HiddenCascade hc_;
        std::vector<std::vector<int>> insiders_;
        std::vector<int> num_announcements_;
};



#endif
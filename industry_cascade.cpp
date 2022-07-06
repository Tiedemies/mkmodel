#include"industry_cascade.hpp"
#include"graphmodel.hpp"
#include"hidden_cascade.hpp"
#include"h_random.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<set>
#include<cmath>
#include<algorithm>

void cumulate(std::vector<double>&lhs, const std::vector<double>& rhs)
{
    if (lhs.empty())
    {
        lhs = rhs;
        return;
    }
    for (size_t i = 0; i < rhs.size(); ++i)
    {
        lhs[i] += rhs[i];
    }
} 

IndustryCascade::IndustryCascade(HiddenCascade h) : hc_(h), simulation_coefficient_(10)
{
    // void
}

IndustryCascade::~IndustryCascade()
{
    // void
}

std::vector<double> 
IndustryCascade::RunTotal()
{
    if (insiders_.empty())
    {
        throw std::logic_error("Insiders not set.");
    }
    double cumulative_min = 0.0;
    double cumulative_max = 0.0;
    std::vector<double> totals;
    for (int i = 0; i < static_cast<int>(insiders_.size()); ++i)
    {
        int n = simulation_coefficient_*num_announcements_[i];
        hc_.SetSimulationN(n);
        totals.push_back(hc_.Simulate(insiders_.at(i)));
    }
    return totals;
}
// std::vector<double> RunSingle(int i);


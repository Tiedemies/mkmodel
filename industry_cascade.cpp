#include"industry_cascade.hpp"
#include"graphmodel.hpp"
#include"hidden_cascade.hpp"
#include"h_random.hpp"
#include"ioroutines.hpp"
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


IndustryCascade::IndustryCascade(HiddenCascade h) : simulation_coefficient_(10), foo_(false), hc_(h)
{
   //void
}

IndustryCascade::IndustryCascade(int year) : simulation_coefficient_(10), def_p_(0.1), foo_(year),
    hc_(HiddenCascade(foo_.GetGraph(year),def_p_,def_p_,def_p_))
{
    foo_.ReadAnnounceTable();
    foo_.ReadCompanyDictionary();
    // We initialize the announcement numbers first
    num_announcements_.resize(foo_.cnames_.size()+1,0);
    insiders_.resize(foo_.cnames_.size());
    MonoGraph* graph = foo_.GetGraph(year);
    for (auto is_name: foo_.cnames_)
    {
        //DEBUG:
        //std::cerr << is_name.first << "," << is_name.second << "\n";
        int cnum = is_name.second;
        const auto& inside = graph->GetInsider(cnum);
        insiders_.at(cnum).insert(insiders_.at(cnum).end(), inside.cbegin(), inside.cend());
    }
    for (auto an_pair: foo_.an_table_)
    {
        std::string isin = an_pair.first;
        std::string cname = foo_.isin_company_[isin];
        const auto& anlist = an_pair.second;
        int cnum = foo_.cnames_[cname];
        for (auto an_it = std::lower_bound(an_pair.second.cbegin(),an_pair.second.cend(), year - 9999); *an_it <= year; ++an_it)
        {
            ++num_announcements_[cnum];
        }
        std::cerr << "cname: " << num_announcements_[cnum] << "\n"; 
    }
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
        // Skip the insiders that empty
        if (insiders_[i].empty())
        {
            continue;
        }
        int n = simulation_coefficient_*num_announcements_[i];
        hc_.SetSimulationN(n);
        totals.push_back(hc_.Simulate(insiders_.at(i)));
    }
    return totals;
}
// std::vector<double> RunSingle(int i);


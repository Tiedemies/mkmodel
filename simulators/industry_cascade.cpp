#include"industry_cascade.hpp"
#include"hidden_cascade.hpp"

#include"../graphmodel/graphmodel.hpp"

#include"../utils/h_random.hpp"
#include"../utils/ioroutines.hpp"
#include"../utils/defs.hpp"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<string>
#include<set>
#include<cmath>
#include<algorithm>
#include "boost/date_time/gregorian/gregorian.hpp"

IndustryCascade::IndustryCascade(HiddenCascade h) : simulation_coefficient_(10), foo_(false), hc_(h)
{
  //void
}

IndustryCascade::IndustryCascade(int year) : simulation_coefficient_(100), def_p_(0.1), foo_(year),
  hc_(HiddenCascade(foo_.GetGraph(year),def_p_,def_p_,def_p_))
{
  foo_.ReadAnnounceTable();
  foo_.ReadCompanyDictionary();
  using  namespace boost::gregorian; 
  date ref_time(from_simple_string(REFDAY));
  int yy = VREFDATE / 10000;
  int mm = (VREFDATE - yy*10000) / 100;
  int dd = VREFDATE%100; 
  date bgtime(yy,mm,dd);
  int days_begin = (bgtime - ref_time).days();
  int days_end = days_begin + 365; 

  // We initialize the announcement numbers first
  num_announcements_.resize(foo_.cnames_.size()+1,0);
  announcement_days_.resize(foo_.cnames_.size()+1);
  insiders_.resize(foo_.cnames_.size()+1);
  MonoGraph* graph = foo_.GetGraph(year);
  //std::cerr << "grap ok " << graph << "\n";
  for (auto is_name: foo_.cnames_)
  {
    //DEBUG:
    //std::cerr << is_name.first << "," << is_name.second << "\n";
    int cnum = is_name.second;
    const auto& inside = graph->GetInsider(cnum);
    insiders_.at(cnum).insert(insiders_.at(cnum).end(), inside.cbegin(), inside.cend());
  }
  int total = 0;
  int comps = 0;
  int max = 0;
  int min = std::numeric_limits<int>::max(); 
  //std::cerr << "insiders initialized\n";
  for (auto an_pair: foo_.an_table_)
  {
    std::string isin = an_pair.first;
    //std::cerr << "isin: " << isin << "\n";
    if(foo_.isin_company_.find(isin) == foo_.isin_company_.end())
    {
      continue;
    }
    std::string cname = foo_.isin_company_[isin];
    //std::cerr << "company: " << cname << "\n";
    if (foo_.cnames_.find(cname) == foo_.cnames_.end())
    {
      continue;
    }
    ++comps;
    int cnum = foo_.cnames_[cname];
    //std::cerr << "Number: " << cnum << "\n";
    for (auto an: an_pair.second)
    {
      if (days_begin < an && an <= days_end)
      {
        ++num_announcements_[cnum];
        announcement_days_[cnum].push_back(an-days_begin);
        ++total;
      }
    }
    min = std::min(num_announcements_[cnum], min);
    max = std::max(num_announcements_[cnum], max);
  }
  std::cerr << total << " announcements for " << comps << " companies in time window, min:" 
            << min << ", max:" << max << "\n"; 
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
  [[maybe_unused]] double cumulative_min = 0.0;
  [[maybe_unused]] double cumulative_max = 0.0;
  std::vector<double> totals;
  totals.resize(insiders_.size(),0.0);
  for (int i = 0; i < static_cast<int>(insiders_.size()); ++i)
  {
    // Skip the insiders that empty
    if (insiders_[i].empty())
    {
      continue;
    }
    int n = simulation_coefficient_*(num_announcements_[i] + 1);
    hc_.SetSimulationN(n);
    totals[i] = (hc_.Simulate(insiders_.at(i)));
    BOOST_ASSERT(!std::isnan(totals[i]));
  }
  return totals;
}
// std::vector<double> RunSingle(int i);

void 
IndustryCascade::RandomizeCoefs(double mu)
{
  hc_.RandomizeWeights(mu);
}


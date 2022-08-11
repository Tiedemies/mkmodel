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

IndustryCascade::IndustryCascade(int tdate, int days): simulation_coefficient_(100), date_(tdate), days_(days), n_node_(-1), n_comp_(-1), 
  def_p_(0.1), foo_(tdate), hc_(HiddenCascade(foo_.GetGraph(tdate),def_p_,def_p_,def_p_))
{
  InitializeFoo();
  InitializeAnnouncements();
  InitializeWindows();
  InitializeTransactions();
  InitializeProbabilities();
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
  //[[maybe_unused]] double cumulative_min = 0.0;
  //[[maybe_unused]] double cumulative_max = 0.0;
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

bool
IndustryCascade::Transaction::operator<(const IndustryCascade::Transaction& rhs) const
{
  // First priority is company number
  if (comp_ < rhs.comp_)
  {
    return true;
  }
  if (comp_ > rhs.comp_)
  {
    return false;
  }
  // Second priority is date
  if (date_ < rhs.date_)
  {
    return true;
  }
  if (date_ > rhs.date_)
  {
    return false;
  }
  // Third is node 
  if (node_ < rhs.node_)
  {
    return true;
  }
  if (node_ > rhs.node_)
  {
    return false;
  }  
  return false;

}

NodeTransactionTable 
IndustryCascade::RunGeneration()
{
  // Precondition
  BOOST_ASSERT(true);
  for(int comp = 0; comp < insiders_.size(); ++comp)
  {
    auto an_it = announcement_days_[comp].begin();
    bool inside = false;
    for(int d = 1; d < days_; ++d)
    {
      while (an_it != announcement_days_[comp].end() && *an_it <= d)
      {
        ++an_it;
      }
      if (an_it != announcement_days_[comp].end() && *an_it < d+5)
      {
        inside = true; 
      }
      else
      {
        inside = false; 
      }
      for(int node = 0; node < out_trading_prob_.size(); ++node)
      {

      }
    }
  }
}

void
IndustryCascade::InitializeFoo()
{
  // We will use announcements, companies, and transactions. 
  // We initialize the I/O and insiders list. 
  foo_.ReadAnnounceTable();
  foo_.ReadCompanyDictionary();
  foo_.ReadTransactionTable();

  // We initialize the announcement numbers first
  n_comp_ = foo_.cnames_.size() + 1;
  num_announcements_.resize(n_comp_,0);
  announcement_days_.resize(n_comp_);
  insiders_.resize(n_comp_);
  const MonoGraph* graph = foo_.GetGraph(date_);
  //std::cerr << "grap ok " << graph << "\n";
  for (auto is_name: foo_.cnames_)
  {
    //DEBUG:
    //std::cerr << is_name.first << "," << is_name.second << "\n";
    int cnum = is_name.second;
    const auto& inside = graph->GetInsider(cnum);
    insiders_.at(cnum).insert(insiders_.at(cnum).end(), inside.cbegin(), inside.cend());
  }
}

void
IndustryCascade::InitializeAnnouncements()
{
  using namespace boost::gregorian; 
  date ref_time(from_simple_string(REFDAY));
  int yy = VREFDATE / 10000;
  int mm = (VREFDATE - yy*10000) / 100;
  int dd = VREFDATE%100; 
  date bgtime(yy,mm,dd);
  int days_begin = (bgtime - ref_time).days();
  int days_end = days_begin + days_; 
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
      int cnum = foo_.cnames_[cname];
      //std::cerr << "Number: " << cnum << "\n";
      for (auto an: an_pair.second)
      {
        if (days_begin < an && an <= days_end)
        {
          ++num_announcements_[cnum];
          announcement_days_[cnum].push_back(an-days_begin);
        }
      }
    }
}

void 
IndustryCascade::InitializeTransactions()
{
  using namespace boost::gregorian; 
  date ref_time(from_simple_string(REFDAY));
  int yy = VREFDATE / 10000;
  int mm = (VREFDATE - yy*10000) / 100;
  int dd = VREFDATE%100; 
  date bgtime(yy,mm,dd);
  int days_begin = (bgtime - ref_time).days();
  int days_end = days_begin + days_; 
  for (auto tran: foo_.tr_table_)
  {
    int node = tran.first;
    n_node_ = std::max(n_node_, node+1);
    for(auto trans_v: tran.second.pt_)
    {
      const std::string& isin = trans_v.first;
      const std::string& cname = foo_.isin_company_[isin];
      const int cnum = foo_.cnames_[cname];
      auto announce_it = announcement_days_[cnum].begin();
      // Go through the transactions. 
      for (auto trans: trans_v.second)
      {
        int date = std::get<0>(trans);
        // Skip those that are not in the data.
        if (date <= days_begin || date > days_end)
        {
          continue;
        }
        int n_date = date-days_begin;
        Transaction tr;
        tr.node_ = node;
        tr.comp_ = cnum;
        tr.date_ = n_date;
        tr.price_ = std::get<1>(trans); 
        tr.vol_ = std::get<2>(trans);
        tr.inside_ = in_window_days_[cnum][n_date]; 
        transacts_.push_back(tr);
      }
    }
  }
  std::sort(transacts_.begin(), transacts_.end());
}

void
IndustryCascade::InitializeWindows()
{
  in_window_days_.resize(n_comp_, std::vector<bool>(days_+1,false));
  for (int comp = 0; comp < n_comp_; ++comp)
  {
    auto an_it = announcement_days_[comp].cbegin(); 
    if (an_it == announcement_days_[comp].cend())
    {
      continue;
    }
    for (int day = 1; day < days_; ++day)
    {
      while (an_it != announcement_days_[comp].end() && *an_it <= day)
      {
        ++an_it;
      }
      if (an_it  == announcement_days_[comp].end() || *an_it > day+5)
      {
        in_window_days_[comp][day] = false;
      }
      else
      {
        in_window_days_[comp][day] = true;
      }
    }
  }
}

void 
IndustryCascade::InitializeProbabilities()
{
  out_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));
  in_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));
  for (int comp = 0; comp < n_comp_; ++comp)
  {
    for(int day = 1; day < days_; ++day)
    {

    }
  }
}
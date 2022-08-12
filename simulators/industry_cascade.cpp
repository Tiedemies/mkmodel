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

IndustryCascade::IndustryCascade(int tdate, int days): simulation_coefficient_(100), window_size_(7), date_(tdate), days_(days), n_node_(-1), n_comp_(-1), 
  def_p_(0.1), foo_(tdate), hc_(HiddenCascade(foo_.GetGraph(tdate),def_p_,def_p_,def_p_))
{
  InitializeFoo();
  InitializePrices();
  InitializeAnnouncements();
  InitializeWindows();
  InitializeTransactions();
  InitializeProbabilities();
  //std::cerr << "Randomizer gives: " << rnd_.get() << "," << rnd_.get() << "," << rnd_.get() << "\n"; 
}


IndustryCascade::~IndustryCascade()
{
  // void
}

void
IndustryCascade::RunDiagnostics()
{
  std::cerr << "Running Diagnostics. \n";
  std::cerr << "Announcements sanity check:  \n";
  /*
  for (int comp = 0; comp < n_comp_; ++comp)
  {
    std::cerr << "*** Company " << comp << ": " << announcement_days_[comp].size() << "announcements";
    if (announcement_days_[comp].empty())
    {
      std::cerr << ".\n"; 
      continue;
    }
    std::cerr << ", first: " << announcement_days_[comp].front() << ", last: " << announcement_days_[comp].back()
      << ".\n"; 
  }
  */
  std::cerr << "*************\n";
  std::cerr << "Probabilities sanity check: \n";
  for (int node = 0; node < n_node_; ++node)
  {
    for(int comp = 0; comp < n_comp_; ++comp)
    {
      if (in_trading_prob_[node][comp] > 1 || in_trading_prob_[node][comp] < 0)
      {
        std::cerr << "bad probability at " << node << "," << comp << ": " << in_trading_prob_[node][comp] << "\n";
        throw std::logic_error("bad probability"); 
      }   
      if (out_trading_prob_[node][comp] > 1 || in_trading_prob_[node][comp] < 0)
      {
        std::cerr << "bad probability at " << node << "," << comp << ": " << out_trading_prob_[node][comp] << "\n";
        throw std::logic_error("bad probability"); 
      }   
    }
  }
  std::cerr << "passed \n"; 

}

std::vector<double> 
IndustryCascade::RunTotal()
{
  if (insiders_.empty())
  {
    throw std::logic_error("Insiders not set.");
  }
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

std::vector<IndustryCascade::Transaction>
IndustryCascade::RunGeneration()
{
  // Precondition
  std::vector<Transaction> result; 
  BOOST_ASSERT(true);
  bool was_in = false;
  hc_.SetSimulationN(1); 
  for(int comp = 0; comp < insiders_.size(); ++comp)
  {
    for(int d = 1; d < days_; ++d)
    {
      // No trades if day is not trading day
      if (std::isnan(prices_[comp][d]))
      {
        continue;
      }
      const bool& in = in_window_days_[comp][d];
      // If we just entered a window:
      if (in && !was_in)
      {
        hc_.Simulate(insiders_[comp]);
      }
      was_in = in; 
      for(int node = 0; node < out_trading_prob_.size(); ++node)
      {
        if ( (in && std::isnan(in_trading_prob_[node][comp])) || (!in && std::isnan(out_trading_prob_[node][comp])) )
        {
          continue;
        } 
        if (in && rnd_.get() > in_trading_prob_[node][comp])
        {
          Transaction tr;
          tr.comp_ = comp;
          tr.date_ = d;
          tr.node_ = node;
          tr.price_ = prices_[comp][d];
          // It is a singleton, so any reasonable threshold will do: 
          if (hc_.simulated_activations_[node] > 0.1)
          {
            tr.vol_ = IsBear(comp,d)?10:-10; 
          }
          else
          {
            tr.vol_ = (rnd_.get() > 0.5)?10:-10; 
          }
          tr.inside_ = true;
          result.push_back(tr);
        }
        else if (!in && rnd_.get() > out_trading_prob_[node][comp])
        {
          Transaction tr;
          tr.comp_ = comp;
          tr.date_ = d;
          tr.node_ = node;
          tr.price_ = prices_[comp][d];
          tr.vol_ = (rnd_.get() > 0.5)?10:-10; 
          tr.inside_ = false;
          result.push_back(tr);
        }
        
      }
    }
  }
  std::cerr << "Generated " << result.size() << " transactions \n";
  return result;
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
  days_begin_ = (bgtime - ref_time).days();
  days_end_ = days_begin_ + days_; 
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
        if (days_begin_ < an && an <= days_end_)
        {
          ++num_announcements_[cnum];
          announcement_days_[cnum].push_back(an-days_begin_);
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
  days_begin_ = (bgtime - ref_time).days();
  days_end_ = days_begin_ + days_; 
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
        if (date <= days_begin_ || date > days_end_)
        {
          continue;
        }
        int n_date = date-days_begin_;
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
  std::cerr << "Stored " << transacts_.size() << " transactions \n";
}

void
IndustryCascade::InitializeWindows()
{
  in_window_days_.resize(n_comp_, std::vector<bool>(days_+1,false));
  n_in_days_.resize(n_comp_,0);
  n_out_days_.resize(n_comp_,0);
  for (int comp = 0; comp < n_comp_; ++comp)
  {
    auto an_it = announcement_days_[comp].cbegin(); 
    if (an_it == announcement_days_[comp].cend())
    {
      continue;
    }
    for (int day = 1; day < days_; ++day)
    {
      if (std::isnan(prices_[comp][day]))
      {
        continue;
      }
      while (an_it != announcement_days_[comp].end() && *an_it <= day)
      {
        ++an_it;
      }
      if (an_it  == announcement_days_[comp].end() || *an_it > day+window_size_)
      {
        in_window_days_[comp][day] = false;
        ++n_out_days_[comp];
      }
      else
      {
        in_window_days_[comp][day] = true;
        ++n_in_days_[comp];
      }
    }
  }
}

void 
IndustryCascade::InitializeProbabilities()
{
  out_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));
  in_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));
  auto tr_it = transacts_.begin();
  for (int comp = 0; comp < n_comp_; ++comp)
  {
    while (tr_it != transacts_.end() && tr_it->comp_ < comp)
    {
      ++tr_it;
    }
    if (tr_it == transacts_.end() || tr_it->comp_ > comp)
    {
      continue;
    }
    for(int day = 1; day < days_; ++day)
    {
      while (tr_it != transacts_.end() && tr_it->date_ < day)
      {
        ++tr_it;
      }
      if (tr_it->date_ == day)
      {
        if (tr_it->inside_ && n_in_days_[comp] > 0.00001)
        {
          in_trading_prob_[tr_it->node_][comp] += 1.0/((double) n_in_days_[comp]);
        }
        else if (n_out_days_[comp] > 0.00001)
        {
          out_trading_prob_[tr_it->node_][comp] += 1.0/((double) n_out_days_[comp]);
        }
      }
    }
  }
}

void
IndustryCascade::InitializePrices()
{
  // Read the price index file
  foo_.ReadPriceTable();
  int n_price = 0;
  int skipped = 0;
  prices_.resize(n_comp_, std::vector<double>(days_ + window_size_, std::nan("missing")));
  for(auto pt_pair: foo_.pr_table_.pt_)
  { 
    const std::string& isin = pt_pair.first;
    if (foo_.isin_company_.find(isin) == foo_.isin_company_.end())
    {
      continue;
    }
    const auto& pt = pt_pair.second; 
    const std::string& cname = foo_.isin_company_[isin];
    const int& comp = foo_.cnames_[cname]; 
    auto pt_it = pt_pair.second.cbegin();
    for(int day = 1; day < days_ + window_size_; ++day)
    {
      while(pt_it != pt.cend() && pt_it->first - days_begin_ < day)
      {
        ++pt_it;
      }
      if (pt_it == pt.cend())
      {
        break;
      }
      if(pt_it->first > day)
      {
        continue;
      }
      prices_[comp][day] = pt_it->second;
      ++n_price;
    }
  }
  std::cerr << "price set for " << n_price << " company/day pairs, skipped " <<  skipped << " companies\n"; 
}

bool
IndustryCascade::IsBear(int comp, int day)
{
  BOOST_ASSERT(0 <= comp && comp < prices_.size());
  BOOST_ASSERT(0 <= day && day + window_size_ < prices_[comp].size());
  int w = window_size_ + day - 1;
  const double& ref_pr = prices_[comp][day];
  while (w < prices_[comp].size() && std::isnan(prices_[comp][w]))
  {
    ++w;
  }
  if (w >= prices_[comp].size())
  {
    return false; 
  }
  BOOST_ASSERT(w < prices_[comp].size());
  return (ref_pr < prices_[comp][w]);
}
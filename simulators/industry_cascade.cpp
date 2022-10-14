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



namespace simulator
{
  using namespace util;
  using namespace graphmodel;
  IndustryCascade::IndustryCascade(HiddenCascade h) : simulation_coefficient_(10), foo_(false), hc_(h)
  {
    //void
  }

  IndustryCascade::IndustryCascade(int tdate, int days): simulation_coefficient_(200), window_size_(7), date_(tdate), days_(days), n_node_(-1), n_comp_(-1), 
    def_p_(0.1), foo_(tdate), hc_(HiddenCascade(foo_.GetGraph(tdate),def_p_,def_p_,def_p_) )
  {
    using namespace boost::gregorian; 
    date ref_time(from_simple_string(REFDAY));
    int yy = VREFDATE / 10000;
    int mm = (VREFDATE - yy*10000) / 100;
    int dd = VREFDATE%100; 
    date bgtime(yy,mm,dd);
    days_begin_ = (bgtime - ref_time).days();
    days_end_ = days_begin_ + days_; 
    InitializeFoo();
    InitializeAnnouncements();
    InitializePrices();
    InitializeWindows();
    InitializeTransactions();
    InitializeProbabilities();
    // std::cerr << "Randomizer gives: " << rnd_.get() << "," << rnd_.get() << "," << rnd_.get() << "\n"; 
  }

  // Copy constructor
  IndustryCascade::IndustryCascade(const IndustryCascade& rhs):
  simulation_coefficient_(rhs.simulation_coefficient_),
  window_size_(rhs.window_size_), date_(rhs.date_), 
  days_(rhs.days_), days_begin_(rhs.days_begin_),
  days_end_(rhs.days_end_), n_node_(rhs.n_node_),
  n_comp_(rhs.n_comp_),  def_p_(rhs.def_p_), foo_(rhs.foo_), 
  hc_(rhs.hc_), rnd_(rhs.rnd_), insiders_(rhs.insiders_),
  num_announcements_(rhs.num_announcements_), announcement_days_(rhs.announcement_days_), 
  in_window_days_(rhs.in_window_days_), trading_days_(rhs.trading_days_), 
  prices_(rhs.prices_),  transacts_(rhs.transacts_),
  out_trading_prob_(rhs.out_trading_prob_), in_trading_prob_(rhs.in_trading_prob_),
  n_in_days_(rhs.n_in_days_), n_out_days_(rhs.n_out_days_)
  {
    // void
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
    std::cerr << "Windows sanity check: \n";

    for(int comp = 0; comp < n_comp_; ++comp)
    {
      int in_d = 0;
      int out_d = 0;
      for(int day = 0; day < days_; ++day)
      {
        if (in_window_days_[comp][day])
        {
          in_d++;
        }
        else
        {
          out_d++;
        }
      }
      if (n_in_days_[comp] != in_d)
      {
        std::cerr << "Company " << comp << " had " << in_d << " vs " << n_in_days_[comp] << "in days \n";
      }
      if (n_out_days_[comp] != out_d)
      {
        std::cerr << "Company " << comp << " had " << out_d << " vs " << n_out_days_[comp] << "out days \n";
      }
    }

    std::cerr << "passed. \n";
    std::cerr << "**************\n";
    std::cerr << "Probabilities sanity check: \n";
    int zero_prob_n = 0;
    for (int node = 0; node < n_node_; ++node)
    {
      bool all_zeros = true;
      for(int comp = 0; comp < n_comp_; ++comp)
      {
        //std::cerr << "node: " << node << "comp: " << comp << ", i-prob: "<< in_trading_prob_[node][comp] << "\n";
        if (in_trading_prob_[node][comp] > 1 || in_trading_prob_[node][comp] < 0)
        {
          std::cerr << "bad in probability at " << node << "," << comp << ": " << in_trading_prob_[node][comp] << "\n";
          std::cerr << "In trade days " << n_in_days_[comp] << "\n";
          throw std::logic_error("bad probability"); 
        }   
        //std::cerr << "node: " << node << "comp: " << comp << ", o-prob: "<< out_trading_prob_[node][comp] << "\n";
        if (out_trading_prob_[node][comp] > 1 || in_trading_prob_[node][comp] < 0)
        {
          std::cerr << "bad out probability at " << node << "," << comp << ": " << out_trading_prob_[node][comp] << "\n";
          std::cerr << "Out trade days " << n_out_days_[comp] << "\n";
          throw std::logic_error("bad probability"); 
        }   
        if (in_trading_prob_[node][comp] > 0.01 || out_trading_prob_[node][comp] > 0.01)
        { 
          all_zeros = false;
        }
      }
      if (all_zeros)
      {
        ++zero_prob_n;
      }
    }
    std::cerr << "Number of zero prob nodes:" << zero_prob_n << " out of " << n_node_ << "\n";
    std::cerr << "passed \n"; 
    std::cerr << "****************\n";
    

  }

  std::vector<double> 
  IndustryCascade::GetSimulationWeights()
  {
    std::vector<double> togo;
    double div = 0;
    for (int i = 0; i < static_cast<int>(insiders_.size()); ++i)
    {
      togo.push_back((double) num_announcements_[i]);
    }
    return togo; 
  }

  std::pair<double,double> 
  IndustryCascade::RunTotal()
  {
    std::vector<double> weights(n_comp_,1);
    return RunTotal(weights); 
  }

  std::pair<double,double>
  IndustryCascade::RunTotal(const std::vector<double>& weights)
  {
    BOOST_ASSERT(weights.size() == n_comp_);
    if (insiders_.empty())
    {
      throw std::logic_error("Insiders not set.");
    }    
    double totals = 0;
    double div = 0.0;
    double var = 0.0; 
    double normalizer = std::accumulate(weights.cbegin(), weights.cend(),0.0);
    // std::cerr << "Starting insider cycle\n";
    for (int i = 0; i < static_cast<int>(insiders_.size()); ++i)
    {
      // Skip the insiders that empty
      if (insiders_[i].empty())
      {
        continue;
      }
      int n = weights[i]*simulation_coefficient_;
      if (n == 0)
      {
        continue;
      }
      hc_.SetSimulationN(n);
      div += n;
      totals+= n*(hc_.Simulate(insiders_.at(i)));
      double a = weights[i]/normalizer;
      var += a*a*hc_.LastVar(); 
      BOOST_ASSERT(!std::isnan(totals));
    }
    return std::make_pair(totals/div,sqrt(var));
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
    for(int comp = 0; comp < (int) insiders_.size(); ++comp)
    {
      for(int d = 0; d < days_; ++d)
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
        for(int node = 0; node < (int) out_trading_prob_.size(); ++node)
        {
          if ( (in && std::isnan(in_trading_prob_[node][comp])) || (!in && std::isnan(out_trading_prob_[node][comp])) )
          {
            continue;
          } 
          if (in && rnd_.get() < in_trading_prob_[node][comp])
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
              tr.vol_ = (rnd_.get() < 0.5)?10:-10; 
            }
            tr.inside_ = true;
            result.push_back(tr);
          }
          else if (!in && rnd_.get() < out_trading_prob_[node][comp])
          {
            Transaction tr;
            tr.comp_ = comp;
            tr.date_ = d;
            tr.node_ = node;
            tr.price_ = prices_[comp][d];
            tr.vol_ = (rnd_.get() < 0.5)?10:-10; 
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
    inside_of_.resize(n_node_);
    const MonoGraph* graph = foo_.GetGraph(date_);
    //std::cerr << "grap ok " << graph << "\n";
    for (auto is_name: foo_.cnames_)
    {
      //DEBUG:
      //std::cerr << is_name.first << "," << is_name.second << "\n";
      int cnum = is_name.second;
      const auto& inside = graph->GetInsider(cnum);
      insiders_.at(cnum).insert(insiders_.at(cnum).end(), inside.cbegin(), inside.cend());
      for(const int insider: inside)
      {
        inside_of_.at(insider).push_back(cnum);
      }
    }
  }

  void
  IndustryCascade::InitializeAnnouncements()
  {
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
    std::cerr << "initialize transactions\n";
    std::unordered_map<int,bool> trades; 
    int unique = 0; 
    for (auto tran: foo_.tr_table_)
    {
      int node = tran.first;
      n_node_ = std::max(n_node_, node+1);
      for(auto trans_v: tran.second.pt_)
      {
        const std::string& isin = trans_v.first;
        const std::string& cname = foo_.isin_company_[isin];
        const int cnum = foo_.cnames_[cname];
        // Go through the transactions. 
        int prev_day = -1;
        for (auto trans: trans_v.second)
        {
          int date = std::get<0>(trans);
          // Skip those that are not in the data.
          if (date <= days_begin_ || date > days_end_)
          {
            continue;
          }
          int n_date = date-days_begin_;
          if (n_date != prev_day)
          {
            unique++;
          }
          prev_day = n_date; 
          Transaction tr;
          tr.node_ = node;
          tr.comp_ = cnum;
          tr.date_ = n_date;
          tr.price_ = std::get<1>(trans); 
          tr.vol_ = std::get<2>(trans);
          tr.inside_ = in_window_days_[cnum][n_date]; 
          transacts_.push_back(tr);
          trades[node] = true; 
        }
      }
    }
    std::sort(transacts_.begin(), transacts_.end());
    int num_trades = 0;
    for (auto tra: trades)
    {
      if (tra.second) {++num_trades;}
    }
    std::cerr << "Stored " << transacts_.size() << " transactions on " << unique << " nodedays, " << num_trades << "nodes active \n";
  }

  void
  IndustryCascade::InitializeWindows()
  {
    in_window_days_.resize(n_comp_, std::vector<bool>(days_+1,false));
    n_in_days_.resize(n_comp_,0);
    n_out_days_.resize(n_comp_,0);
    int n_zero_in = 0;
    int n_zero_out = 0;
    for (int comp = 0; comp < n_comp_; ++comp)
    {
      auto an_it = announcement_days_[comp].cbegin(); 
      /*
      if (an_it == announcement_days_[comp].cend())
      {
        continue;
      }
      */
      bool zero_in = true;
      bool zero_out = true;
      for (int day = 0; day < days_; ++day)
      {
        /*
        if (std::isnan(prices_[comp][day]))
        {
          continue;
        }
        */
        while (an_it != announcement_days_[comp].end() && *an_it <= day)
        {
          ++an_it;
        }
        if (an_it  == announcement_days_[comp].end() || *an_it > day+window_size_)
        {
          n_out_days_[comp]++;
          in_window_days_[comp][day] = false;
          zero_out = false;
        }
        else
        {
          n_in_days_[comp]++;
          in_window_days_[comp][day] = true;
          zero_in = false; 
        }
      }
      if (zero_in)
      {
        n_zero_in++;
      }
      if (zero_out)
      {
        n_zero_out++; 
      }
    }
    std::cerr << "Initialized windows, " << n_zero_in << " zero in-window companies, " << n_zero_out << " zero out-windows \n";
  }


  void 
  IndustryCascade::InitializeProbabilities()
  {
    std::cerr << "initializing probabilities\n";
    out_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));
    in_trading_prob_.resize(n_node_, std::vector<double>(n_comp_,0.0));

    int prev_day = -1;
    int prev_node = -1;
    int prev_comp = -1; 
    for (auto tr: transacts_)
    {
      if (tr.date_ == prev_day && tr.node_ == prev_node && tr.comp_ == prev_comp)
      {
        continue; 
      }
      if (tr.inside_)
      {
        in_trading_prob_[tr.node_][tr.comp_] += 1.0; 
      }
      else
      {
        out_trading_prob_[tr.node_][tr.comp_] += 1.0;
      }
      prev_day = tr.date_;
      prev_comp = tr.comp_;
      prev_node = tr.node_;
    }
    for (int n = 0; n < n_node_; ++n)
    {
      for (int c = 0; c < n_comp_; ++c)
      {
        if (n_out_days_[c] > 0.1)
        {
          out_trading_prob_[n][c] /= n_out_days_[c];
        }
        else
        {
          out_trading_prob_[n][c] = 0;
        }
        if (n_in_days_[c] > 0.1)
        {
          in_trading_prob_[n][c] /= n_in_days_[c];
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
        ++skipped;
        continue;
      }
      const auto& pt = pt_pair.second; 
      const std::string& cname = foo_.isin_company_[isin];
      const int& comp = foo_.cnames_[cname]; 
      auto pt_it = pt_pair.second.cbegin();
      bool found_price = false;
      for(int day = 0; day < days_ + window_size_; ++day)
      {
        while(pt_it != pt.cend() && pt_it->first - days_begin_ < day)
        {
          ++pt_it;
        }
        if (pt_it == pt.cend())
        {
          break;
        }
        if(pt_it->first - days_begin_ > day)
        {
          continue;
        }
        BOOST_ASSERT(pt_it->first - days_begin_ == day);
        found_price = true;
        prices_[comp][day] = pt_it->second;
        ++n_price;
      }
      if (!found_price)
      {
        ++skipped; 
      }
    }
    std::cerr << "price set for " << n_price << " company/day pairs, skipped " <<  skipped << " companies\n"; 
    n_price = 0;
    for (int day = 0; day < days_ + window_size_; ++day)
    {
      bool trade_day = false; 
      for (int comp = 0; comp < n_comp_; ++comp)
      {
        if (! std::isnan(prices_[comp][day]) )
        {
          trade_day = true;
          break;
        }
      }
      if (! trade_day)
      {
        continue;
      }
      for (int comp = 0; comp < n_comp_; ++comp)
      {
        if (std::isnan(prices_[comp][day]))
        {
          prices_[comp][day] = InterPolate(prices_[comp], day);
          ++n_price;
        }
      }
    }
    std::cerr << "price extrapolated for " << n_price << " company/day pairs\n"; 
  }

  bool
  IndustryCascade::IsBear(int comp, int day)
  {
    BOOST_ASSERT(0 <= comp && comp < prices_.size());
    BOOST_ASSERT(0 <= day && day + window_size_ < prices_[comp].size());
    int w = window_size_ + day - 1;
    const double& ref_pr = prices_[comp][day];
    while (w < (int) prices_[comp].size() && std::isnan(prices_[comp][w]))
    {
      ++w;
    }
    if (w >= (int) prices_[comp].size())
    {
      return false; 
    }
    BOOST_ASSERT(w < prices_[comp].size());
    return (ref_pr < prices_[comp][w]);
  }

  double InterPolate(const std::vector<double>& vec, int i)
  {
    // Precondition
    BOOST_ASSERT(vec.size() > 1);

    // Case one: It is outside and we extrapolate:
    if (i >= (int) vec.size())
    {
      int j = vec.size() - 1;
      while (j > 0 && std::isnan(vec[j]))
      {
        --j;
      }
      if (j == 0)
      {
        return vec[0];
      } 
      else
      {
        int k = j-1;
        while (k > 0 && std::isnan(vec[k]))
        {
          --k;
        }
        double unit_inc = (vec[j] - vec[k])/(j-k);
        return vec[j] + unit_inc*(i-j);
      }
    }
    int j = i-1;
    int k = i+1;
    while (j > 0 && std::isnan(vec[j]))
    {
      --j;
    } 
    while (k < (int) vec.size() && std::isnan(vec[k]))
    {
      ++k;
    }
    double unit_inc = (vec[k] - vec[j])/(k-j);
    return vec[j] + unit_inc*(i-j);
  }
  

  void 
  IndustryCascade::DeactivateFromInside(int node, int comp)
  {
    // Preconditions: Legit node, legit comp, node is insider of comp
    BOOST_ASSERT(0 <= node && node < n_node_);
    BOOST_ASSERT(0 <= comp && comp < n_comp_);
    auto pos = std::find(insiders_.at(comp).begin(), insiders_.at(comp).end(), node);
    BOOST_ASSERT(pos != insiders_.at(comp).end() );

    DisableCacheEntry to_disable(node, comp);
    disable_cache_.push_back(to_disable);
    for (auto f_it = insiders_[comp].begin(); f_it != insiders_[comp].end(); ++f_it)
    {
      if (f_it == pos)
      {
        continue;
      }
      hc_.DeactivateConnection(node,*f_it);
    }
    insiders_[comp].erase(pos);
  }

  void 
  IndustryCascade::ReactivateInside(int node, int comp)
  {
    // Precondition:
    BOOST_ASSERT(0 <= node && node < n_node_);
    BOOST_ASSERT(0 <= comp && comp < n_comp_);
    BOOST_ASSERT(std::find(insiders_.at(comp).begin(), insiders_.at(comp).end(), node) == insiders_.at(comp).end());
    auto dis_it = std::find(disable_cache_.begin(),disable_cache_.end(), DisableCacheEntry(node,comp));
    while (dis_it != disable_cache_.end())
    {
      disable_cache_.erase(dis_it);
      dis_it = std::find(disable_cache_.begin(),disable_cache_.end(), DisableCacheEntry(node,comp));
    }
    // Activate all connections:
    for (int i: insiders_[comp])
    {
      hc_.ActivateConnection(node,i);
    }
    // Re-insert into the company insiders. 
    insiders_[comp].push_back(node);
  
    // Postcondition: The entry is not in disablecache. 
    BOOST_ASSERT(dis_it == disable_cache_.end()  );
  }

  void 
  IndustryCascade::ReactivateInside(int node)
  {
    // Precondition:
    BOOST_ASSERT(0 <= node && node < n_node_);
    auto cache_it = disable_cache_.begin();  
    while (cache_it != disable_cache_.end())
    {
      if (cache_it->node_ == node)
      {
        [[maybe_unused]] int comp = cache_it->comp_;
        auto temp = cache_it;
        ++cache_it;
        disable_cache_.erase(temp);
      } 
      else
      {
        ++cache_it; 
      }
    }
  }
  void 
  IndustryCascade::ReactivateInside()
  {
    while (!disable_cache_.empty())
    {
      DisableCacheEntry entry = disable_cache_.back();
      disable_cache_.pop_back();
      ReactivateInside(entry.node_, entry.comp_);
    }
    
  }

  void 
  IndustryCascade::SetConstantProb(double p)
  {
    hc_.SetConstantProb(p);
  }

const std::vector<double>& 
IndustryCascade::GetCompCentrality()
{
  if (comp_adj_.size() != n_comp_)
  {
    InitializeCompAdj();
  }
  // Precondition
  if (comp_centrality_.size() != n_comp_)
  {
    comp_centrality_.resize(n_comp_, 0);

    std::vector<std::list<int>> P(n_comp_);
    for (int s = 0; s < n_comp_; ++s)
    {
      std::list<int> S;
      std::vector<std::list<int>> P(n_comp_);
      std::list<int> Q; 
      std::vector<double> sigma(n_comp_,0);
      std::vector<int> d(n_comp_, -1);
      d[s] = 0;
      sigma[s] = 1.0; 
      Q.push_back(s);
      while (!Q.empty())
      {
        int u = Q.front();
        Q.pop_front();
        S.push_back(u);
        for (int v: comp_adj_[u])
        {
          if (d[v] < 0)
          {
            d[v] = d[u] + 1;
            Q.push_back(v);
          }
          if (d[v] == d[u] + 1)
          { 
            sigma[v] = sigma[v] + sigma[u];
            P[v].push_back(u);
          } 
        }
      }
      std::vector<double> delta(n_comp_,0.0);
      while(!S.empty())
      {
        // How many from s to w 
        int w = S.back();
        S.pop_back();
        for (int v: P[w])
        {
          delta[v] = delta[v] + (sigma[v]/sigma[w])*(1 + delta[w]);
        }
        if(w != s)
        {
          comp_centrality_[w] = comp_centrality_[w] + delta[w];
        }
      }
    }
  }
  return comp_centrality_; 
}

void
IndustryCascade::InitializeCompAdj()
{
  BOOST_ASSERT(true);
  // Precondition.
  comp_adj_.resize(n_comp_);
  for (int i = 0; i < n_comp_-1; ++i)
  {
    std::set<int> adjs(insiders_[i].begin(), insiders_[i].end());
    for (int j = i+1; j < n_comp_; ++j)
    {
      for (int v: insiders_[j])
      {
        if (adjs.find(v) != adjs.end())
        {
          comp_adj_[i].push_back(j);
          comp_adj_[j].push_back(i);
          break;
        }
      }    
    }
  }
}





void 
IndustryCascade::GraphDiagnostics(std::ostream& out)
{
  const auto& comp_c = GetCompCentrality();
  const auto& node_c = hc_.GetNodeCentrality();
  out << "* diagnostic start *\n";
  out.flush();

  for (int i = 0; i < n_comp_; ++i)
  {
    out << "comp: " << i << " centrality: " << comp_c.at(i) << "\n";
    out.flush(); 
  } 
  out << "*** \n"; 
  for (int i = 0; i < n_node_; ++i)
  {
    out << "Node: " << i << " centrality: " << node_c.at(i) << "\n";
    out.flush(); 
  } 

  
  return;
}

const std::unordered_map<size_t,double>& 
IndustryCascade::GetNeighbourOutActivation()
{
  // Precondition
  BOOST_ASSERT(hc_.gather_statistic_ && !hc_.edge_activations_.empty());
  if (!hc_.gather_statistic_ || hc_.edge_activations_.empty())
  {
    throw std::logic_error("Hidden cascade not calculated properly, not possible to calculate Neighbour activation");
  }
  for (int node = 0; node < n_node_; ++node)
  {
    for (int comp = 0; comp < n_comp_; ++comp)
    {
      if (std::find(insiders_.at(comp).begin(),insiders_.at(comp).end(), node) == insiders_.at(comp).end())
      {
        continue;
      }
      size_t key = hc_.Key(node,comp);
      neighbour_out_activation_[key] = 0;
      for (const auto& other_node: insiders_.at(comp))
      {
        if(other_node == node)
        {
          continue;
        }
        size_t key_2 = hc_.Key(node,other_node);
        neighbour_out_activation_[key] += hc_.edge_activations_.at(key_2);
      }
    }
  }
  return neighbour_out_activation_;
}

const std::unordered_map<size_t,double>& 
IndustryCascade::GetNeighbourInActivation()
{
  // Precondition
  BOOST_ASSERT(hc_.gather_statistic_ && !hc_.edge_activations_.empty());
  if (!hc_.gather_statistic_ || hc_.edge_activations_.empty())
  {
    throw std::logic_error("Hidden cascade not calculated properly, not possible to calculate Neighbour activation");
  }
  for (int node = 0; node < n_node_; ++node)
  {
    for (int comp = 0; comp < n_comp_; ++comp)
    {
      if (std::find(insiders_.at(comp).begin(),insiders_.at(comp).end(), node) == insiders_.at(comp).end())
      {
        continue;
      }
      size_t key = hc_.Key(node,comp);
      neighbour_in_activation_[key] = 0;
      for (const auto& other_node: hc_.adj_.at(node))
      {
        if(other_node == node ||  std::find(insiders_.at(comp).begin(),insiders_.at(comp).end(), other_node) != insiders_.at(comp).end())
        {
          continue;
        }
        size_t key_2 = hc_.Key(node,other_node);
        neighbour_in_activation_[key] += hc_.edge_activations_.at(key_2);
      }
    }
  }
  return neighbour_in_activation_;
}


const std::vector<double>& 
IndustryCascade::GetNodeCentrality()
{
  return hc_.GetNodeCentrality();
}

}
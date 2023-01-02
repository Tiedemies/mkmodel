#include"hidden_cascade.hpp"

#include"../graphmodel/graphmodel.hpp"
#include"../utils/h_random.hpp"

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

namespace simulator
{
  using namespace graphmodel;
  using namespace util;
// The constructor uses uniform p, false positive and true positive probabilities
HiddenCascade::HiddenCascade(const MonoGraph* mg, const double& p, const double& fp, const double& tp)
{
    gather_statistic_ = false; 
    // We copy the structure of the monograph and ascribe constant probabilities to all links
    int i = 0;
    uni_p_ = p; 
    for (auto apair: mg->adj_)
    {
        int src = apair.first;
        for (int tgt: apair.second)
          {
              p_map_[Key(src,tgt)] = p;
              i = std::max(i, tgt);
          }
          i = std::max(i, src);
    }
    adj_.resize(i+1);
    for (auto apair: mg->adj_)
    {
        int src = apair.first;
        for (int tgt: apair.second)
        {
            adj_[src].push_back(tgt); 
        }         
    }
    false_positive_prob_.resize(i+1,fp);
    true_positive_prob_.resize(i+1,tp);
    is_disabled_.resize(i+1,false);
    simulated_activations_.resize(i+1,0.0);
    sim_n_ = 5*i;
}

HiddenCascade::~HiddenCascade()
{
    p_map_.clear();
    simulated_activations_.clear();
    adj_.clear();

    // void;
}

HiddenCascade::HiddenCascade(const HiddenCascade& rhs):
  sim_n_(rhs.sim_n_),  false_positive_prob_(rhs.false_positive_prob_),
  true_positive_prob_(rhs.true_positive_prob_), p_map_(rhs.p_map_),
  adj_(rhs.adj_),uni_p_(rhs.uni_p_),is_disabled_(rhs.is_disabled_), disabled_(rhs.disabled_) 
{
  // pass
}


size_t
HiddenCascade::Key(const int& u, const int& v) const
{
  return (size_t) (( (size_t) u << 32) | (unsigned int) v);
}


// Run a simulation. 
// Store the information to simulated activations.
double
HiddenCascade::Simulate(const std::vector<int>& inside) 
{
  BOOST_ASSERT(sim_n_ > 0);
  size_t nodes = simulated_activations_.size();
  //max_activated_ = 0;
  //min_activated_ = simulated_activations_.size();
  // std::fill(simulated_activations_.begin(),simulated_activations_.end(), 0.0);
  infected_by_sim_ = std::vector<int>(sim_n_, 0);
  double num_active = 0;
  if (gather_statistic_)
  {
    for (auto x: p_map_)
    {
      edge_activations_[x.first] = 0; 
    }
  }
  #pragma omp parallel for shared(num_active, edge_activations_)
  for (size_t i = 0; i < sim_n_; i=i+1)
  {
    std::unordered_map<size_t, double> activation_coefficient_;
    if (gather_statistic_)
    {
      for (auto x: p_map_)
      {
        activation_coefficient_[x.first] = 0; 
      }
    }
    Random rnd;
    // std::queue<double> random_numbers;
    int num_infected = 0;
    //for (int pass = 0; pass < 2; ++pass)
    //{
      std::vector<bool> is_infected(nodes+1,false);
      std::stack<int> infected;
    
      // infected.reserve(nodes/2);
      for (const int j: inside)
      {
        //deriv[j] = 0;
        if(!is_disabled_.at(j))
        {
          infected.push(j);
        } 
        is_infected.at(j) = true;
      }
      for (const int j: disabled_)
      {
        is_infected.at(j) = true; 
      }
      while(!infected.empty())
      {
        const int j = infected.top(); 
        infected.pop(); 
        for (auto l = adj_.at(j).cbegin(); l != adj_.at(j).cend(); ++l)
        {
          const int& k = *l; 
          if (is_infected.at(k))
          {
            continue;
          }
          double r_num =  rnd.get(); // (pass==0 || random_numbers.empty())?rnd.get():1-random_numbers.back();
          /*
          if (pass && !random_numbers.empty())
          {
            random_numbers.pop();
          }
          if (!pass)
          {
            random_numbers.push(r_num);
          }
          */
          if (IsSuccess(j,k,r_num))
          {
            is_infected.at(k) = true;
            infected.push(k);
            ++num_infected;
            if (gather_statistic_)
            {
              activation_coefficient_[Key(j,k)] += 1;
            }
          }
        } 
      //}
    }
    
    #pragma omp atomic update
    num_active += num_infected; 
    
    #pragma omp atomic update
    infected_by_sim_[i] += num_infected;

    if (gather_statistic_)
    {
      for (auto x: p_map_)
      {
        #pragma omp atomic update
        edge_activations_[x.first] += activation_coefficient_[x.first]/sim_n_;
      }
    }
    
    //std::cerr << "Simulated: " << num_infected << "\n";
    //min_activated_ = std::min(min_activated_, static_cast<int>(informed.size()));
    //max_activated_ = std::max(min_activated_, static_cast<int>(informed.size()));
  }
  // std::cerr << "complete\n";
  // Normalize
  /*
  #pragma omp parallel for
  for(size_t i = 0; i < simulated_activations_.size(); ++i)
  {
    simulated_activations_[i] /= sim_n_;
  }
  num_active = std::accumulate(simulated_activations_.begin(), simulated_activations_.end(), 0.0);
  BOOST_ASSERT(!std::isnan(num_active));
  */
  return num_active / sim_n_;
}

int
HiddenCascade::SimulateConst(const std::vector<int>& inside, const std::unordered_map<size_t,double>& dice, bool invert) const
{
  int num_infected = 0;
  std::vector<bool> is_infected(GetN(),false);
  std::queue<int> infected;
  for (const int j: inside)
  {
    if(disabled_.size() == 0 || !is_disabled_.at(j))
    {
      infected.push(j);
    } 
    is_infected.at(j) = true;
  }
  for (const int j: disabled_)
  {
    is_infected.at(j) = true; 
  }
  while(!infected.empty())
  {
    const int j = infected.front(); 
    infected.pop(); 
    for (auto l = adj_.at(j).cbegin(); l != adj_.at(j).cend(); ++l)
    {
      const int& k = *l; 
      if (is_infected.at(k))
      {
        continue;
      }
      double r_num = invert?1-dice.at(Key(j,k)):dice.at(Key(j,k));
      if (IsSuccess(j,k,r_num))
      {
        is_infected.at(k) = true;
        infected.push(k);
        ++num_infected;
      }
    } 
  }     
  return num_infected;
}



// Calculte a single success over a connection. 
bool
HiddenCascade::IsSuccess(const int& u, const int& v, Random& rnd) const
{
  return (p_map_.at(Key(u,v)) - rnd.get() > 0);
}

bool
HiddenCascade::IsSuccess(const int& u, const int& v, const double& r_num) const
{
  return (p_map_.at(Key(u,v)) - r_num > 0);
}


// Change the simulation's number of simulations
void 
HiddenCascade::SetSimulationN(size_t n)
{
  sim_n_ = n;
}

// Generate a trading pattern for a day, depending on whether it is inside the time window
std::vector<double> 
HiddenCascade::Generate(const std::vector<int>& inside, bool window_day)
{
  if (window_day)
  {
    SetSimulationN(1);
    Simulate(inside);
  }
}

int
HiddenCascade::GetMaxActivated() const
{
  return max_activated_;
}

int
HiddenCascade::GetMinActivated() const
{
  return min_activated_;
}

void 
HiddenCascade::DeactivateConnection(int u, int v)
{
  // Remove (u,v) and (v,u) from the adjacency lists.
  auto &z = adj_.at(u);
  auto &w = adj_.at(v);
  z.erase(std::remove(z.begin(), z.end(), v), z.end());
  w.erase(std::remove(w.begin(), w.end(), u), w.end());
}

void 
HiddenCascade::ActivateConnection(int u, int v, double p)
{
  auto p_it = p_map_.find(Key(u,v));
  // If p  is not given it is negative, and then it is assumed to be the uniform probability
  if (p < 0) 
  {
    if (p_it == p_map_.end())
    {
      p = uni_p_;
    }
  }
  // Add the element to the adjacency lists
  auto &z = adj_.at(u);
  auto &w = adj_.at(v);
  if (std::find(z.begin(), z.end(), v) == z.end())
  {
    z.push_back(v);
  }
  if (std::find(w.begin(), w.end(), u) == w.end())
  {
    w.push_back(u);
  }
  // Give them the probability unless it exists in which case you reactivate:
  if (p >= 0 && p <= 1)
  {
    p_map_[Key(u,v)] = p;
    p_map_[Key(v,u)] = p;
  }
}

void
HiddenCascade::RandomizeWeights(double mu, double sigma)
{
  // Precondition
  BOOST_ASSERT(mu > 0.0);
  BOOST_ASSERT(mu < 1.0);

  double a = 0.0;
  double b = 0.0;
  if (sigma < 0)
  {
    b = 50; 
    a = b / (1-mu);
  }
  else 
  {
    BOOST_ASSERT(sigma > 0);
    BOOST_ASSERT(sigma < mu*(1-mu));
    a = ( (1-mu) / (sigma*sigma) - 1/mu)*mu*mu;
    b = a*(1/mu - 1);
  }
  BetaRandom gen(a,b);

  /* We iterate over the adjacency list and generate probabilities from beta-distribution */
  for(int u = 0; u < adj_.size(); ++u)
  {
    for(int v: adj_.at(u))
    {
      p_map_[Key(u,v)] = gen.get();
    }
  }
}

int
HiddenCascade::GetN() const
{
  return adj_.size();
}

void 
HiddenCascade::DeactivateNode(int n)
{
  disabled_.insert(n);
  is_disabled_.at(n) = true;
}

void 
HiddenCascade::ReactivateNode(int n)
{
  disabled_.erase(n);
  is_disabled_.at(n) = false;
}

void 
HiddenCascade::SetConstantProb(double p)
{
  BOOST_ASSERT(p >= 0.0);
  BOOST_ASSERT(p <= 1.0);
  uni_p_ = p;
  for (auto& entry: p_map_)
  {
    p_map_[entry.first] = p;
  }
}

double
HiddenCascade::LastVar() const
{
  //return 0;
  
  double avg = std::accumulate(infected_by_sim_.cbegin(), infected_by_sim_.cend(), 0.0) / sim_n_; 
  double err = 0;
  for(size_t i = 0; i < sim_n_; ++i)
  {
      err += (avg - infected_by_sim_[i])*(avg - infected_by_sim_[i]) / (sim_n_- 1);
  }
  // std::cerr << "Variance :" << err << " with ste " << sqrt(err) << " and avg " << avg << "\n";
  return err;
  
}

const std::vector<double>& 
HiddenCascade::GetNodeCentrality()
{
  // Precondition
  if (node_centrality_.size() != adj_.size())
  {
    int n = (int) adj_.size();
    node_centrality_.resize(n, 0);

    std::vector<std::list<int>> P(n);
    for (int s = 0; s < n; ++s)
    {
      std::list<int> S;
      std::vector<std::list<int>> P(n);
      std::list<int> Q; 
      std::vector<double> sigma(n,0);
      std::vector<int> d(n, -1);
      d[s] = 0;
      sigma[s] = 1.0; 
      Q.push_back(s);
      while (!Q.empty())
      {
        int u = Q.front();
        Q.pop_front();
        S.push_back(u);
        const auto adj_u = adj_[u];
        if (adj_u.empty())
        {
          // no neighbours
          continue;
        }
        for (int v: adj_u)
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
      std::vector<double> delta(n,0.0);
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
          node_centrality_[w] = node_centrality_[w] + delta[w];
        }
      }
    }
  }
  return node_centrality_; 
}

}
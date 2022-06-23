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

// The constructor uses uniform p, false positive and true positive probabilities
HiddenCascade::HiddenCascade(const MonoGraph* mg, const double& p, const double& fp, const double& tp)
{
    // We copy the structure of the monograph and ascribe constant probabilities to all links
    int i = 0;
    for (auto apair: mg->adj_)
    {
        int src = apair.first;
        auto alist = apair.second; 
        adj_[src] = alist; 
        for (int tgt: alist)
        {
            p_map_[Key(src,tgt)] = p;
            i = std::max(i, tgt);
        }
        i = std::max(i, src);
    }
    false_positive_prob_.resize(i+1,fp);
    true_positive_prob_.resize(i+1,tp);
    is_disabled_.resize(i+1,false);
    simulated_activations_.resize(i+1,0.0);
    simulated_profits_.resize(i+1,0.0);
    sim_n_ = 5*i;
}

HiddenCascade::~HiddenCascade()
{
    p_map_.clear();
    simulated_activations_.clear();
    simulated_profits_.clear();
    adj_.clear();

    // void;
}
inline 
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
  size_t nodes = simulated_profits_.size();
  std::fill(simulated_activations_.begin(),simulated_activations_.end(), 0.0);
  double num_active = 0;
  #pragma omp parallel for 
  for (size_t i = 0; i < sim_n_; ++i)
  {
    std::vector<bool> is_infected(nodes+1,false);
    std::stack<int> infected;
    std::stack<int> informed;
    // infected.reserve(nodes/2);
    for (int j: inside)
    {
      //deriv[j] = 0;
      if(!is_disabled_.at(j))
      {
        infected.push(j);
      } 
      is_infected.at(j) = true;
    }
    for (int j: disabled_)
    {
      is_infected.at(j) = true; 
    }
    while(!infected.empty())
    {
      int j = infected.top(); 
      infected.pop(); 
      for (int k: adj_.at(j))
      {
        if (is_infected.at(k))
        {
          continue;
        }
        if (IsSuccess(j,k))
        {
          is_infected.at(k) = true;
          #pragma omp critical(A)
          {
            simulated_activations_.at(k) += 1.0;
            num_active += 1;
          }
          infected.push(k);
          informed.push(k);
        }
      } 
    }
  }
  // Normalize
  #pragma omp parallel for
  for(size_t i = 0; i < simulated_profits_.size(); ++i)
  {
    simulated_activations_[i] /= sim_n_;
  }
  return num_active/sim_n_;
}

// Calculte a single success over a connection. 
bool
HiddenCascade::IsSuccess(const int& u, const int& v) const
{
    return (p_map_.at(Key(u,v)) - rnd_.get() > 0);
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
 
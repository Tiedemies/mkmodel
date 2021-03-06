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

// The constructor uses uniform p, false positive and true positive probabilities
HiddenCascade::HiddenCascade(const MonoGraph* mg, const double& p, const double& fp, const double& tp)
{
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
  max_activated_ = 0;
  min_activated_ = simulated_activations_.size();
  std::fill(simulated_activations_.begin(),simulated_activations_.end(), 0.0);
  double num_active = 0;
  #pragma omp parallel for shared(simulated_activations_)
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
      informed.push(j);
    }
    for (int j: disabled_)
    {
      is_infected.at(j) = true; 
    }
    while(!infected.empty())
    {
      const int j = infected.top(); 
      infected.pop(); 
      for (auto l = adj_.at(j).begin(); l != adj_.at(j).end(); ++l)
      {
        const int& k = *l; 
        if (is_infected.at(k))
        {
          continue;
        }
        if (IsSuccess(j,k))
        {
          is_infected.at(k) = true;
          #pragma omp atomic   
          ++simulated_activations_.at(k);
         
          infected.push(k);
          informed.push(k);
        }
      } 
    }
    min_activated_ = std::min(min_activated_, static_cast<int>(informed.size()));
    max_activated_ = std::max(min_activated_, static_cast<int>(informed.size()));
  }
  // Normalize
  #pragma omp parallel for reduction(+:num_active)
  for(size_t i = 0; i < simulated_profits_.size(); ++i)
  {
    simulated_activations_[i] /= sim_n_;
    num_active += simulated_activations_[i];
  }
  return num_active;
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
  // If p  is not given it is negative, and then it is assumed to be the uniform probability
  if (p < 0) 
  {
    p = uni_p_;
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
  // Give them the probability
  p_map_[Key(u,v)] = p;
  p_map_[Key(v,u)] = p;
}


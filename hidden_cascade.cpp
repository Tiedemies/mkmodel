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


void
HiddenCascade::Simulate(const std::vector<int>& inside,bool simulate_profits) 
{
  size_t nodes = simulated_profits_.size();
  double prof_fil = simulate_profits?0.0:std::nan("n");
  std::fill(simulated_profits_.begin(),simulated_profits_.end(), prof_fil);
  std::fill(simulated_activations_.begin(),simulated_activations_.end(), 0.0);

  #pragma omp parallel for 
  for (size_t i = 0; i < sim_n_; ++i)
  {
    std::vector<bool> is_infected(nodes+1,false);
    std::queue<int> infected;
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
      int j = infected.back(); 
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
          }
          infected.push(k);
        }
      } 
    }
    if (!simulate_profits)
    {
      continue;
    }
    for(size_t i = 0; i < simulated_profits_.size(); ++i)
    {
      // Infected gets true positive result
      if (is_infected.at(i))
      {
        #pragma omp critical(B)
        { simulated_profits_[i] += (true_positive_prob_[i] > rnd_.get()?1.0:0.0); }
      }
      // Others get random false positive result 
      else if (false_positive_prob_[i] > rnd_.get())
      {
        #pragma omp critical(B)
        { simulated_profits_[i] += rnd_.get() < 0.5 ? 1.0 : -1.0; } 
      } 
    }

  }
  // Normalize
  #pragma omp parallel for
  for(size_t i = 0; i < simulated_profits_.size(); ++i)
  {
    simulated_activations_[i] /= sim_n_;
    if (simulate_profits)
    {
      simulated_profits_[i] /= sim_n_;
    }
  }
}

bool
HiddenCascade::IsSuccess(const int& u, const int& v) const
{
    return (p_map_.at(Key(u,v)) - rnd_.get() > 0);
}

void 
HiddenCascade::SetSimulationN(size_t n)
{
  sim_n_ = n;
}

/*
    public:
        HiddenCascde(const MonoGraph& mg, const double& p);
        ~HiddenCascde();
        std::vector<double> Simulate
        void Disable(int node);
    private:
        std::vector<double> q_vec_;
        std::unordered_map<std::pair<int,int>, double> p_map_;
*/

 
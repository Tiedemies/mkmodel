#include"announcement_cascade.hpp"
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

  AnnouncementCascade::AnnouncementCascade(int date, int days):ic_(IndustryCascade(date,days))
  {
    //void
  }

  AnnouncementCascade::AnnouncementCascade(const IndustryCascade& ind):ic_(ind)
  {
    //void
  }

  AnnouncementCascade::~AnnouncementCascade()
  {
    // Void
  }

  // Execute one cascade of announcements over all of the industry.
  double 
  AnnouncementCascade::RunSingleCascade(bool anti) const
  {
    double c_h = 0;
    double div = 0;
    //std::cerr << "ic insiders count: " << ic_.insiders_.size() << ", n comp: " << ic_.n_comp_ 
    // << " num announcment size: " << ic_.num_announcements_.size() <<  "\n"; 
    #pragma omp parallel for
    for(int i =0; i < ic_.n_comp_;++i)
    {
      if (i >= ic_.insiders_.size() || i >= ic_.num_announcements_.size())
      {
        continue;
      }
       
      const auto& i_vec = ic_.insiders_.at(i);
      const int n = ic_.num_announcements_.at(i);
      const double prop_div = ic_.n_node_ - i_vec.size();
      if (i_vec.empty() || n < 1)
      {
        continue;
      }
      double xbar_i = 0;
      div += anti?((double)2*n):(double) n;
      for (int l = 0; l < n; ++l)
      {
        const auto dice = GetDies();
        const int a = ic_.hc_.SimulateConst(i_vec,dice,false);
        const int a2 = anti?ic_.hc_.SimulateConst(i_vec,dice,true):0;
        xbar_i += (double)(a + a2)/prop_div; 
      }
      #pragma omp critical
      {
        c_h += xbar_i;
      }
    }
    return c_h/div; 
  }


  const std::unordered_map<size_t,double> 
  AnnouncementCascade::GetDies() const
  {
    Random die;
    std::unordered_map<size_t,double> dies; 
    for (int u = 0; u < ic_.n_node_; ++u)
    {
      const auto& aa = ic_.hc_.adj_.at(u);
      if (aa.empty())
      {
        continue;
      }
      for (int v: aa)
      {
        dies[ic_.hc_.Key(u,v)] = die.get();
      }
    }
    return dies;
  }


}
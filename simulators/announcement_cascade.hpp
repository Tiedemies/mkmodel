#ifndef A_CASCADE
#define A_CASCADE
#include "hidden_cascade.hpp"
#include "industry_cascade.hpp"
#include "../utils/defs.hpp"
#include "../utils/h_random.hpp"
#include "../utils/ioroutines.hpp"

#include "../graphmodel/graphmodel.hpp"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>
#include<list>

namespace algorithm
{
  class InfluenceMinimizer;
}

namespace simulator
{
  using namespace util;
  using namespace graphmodel;
  /* 
  * The AnnouncementCascade class is a wrapper for the Industrial cascade that runs the simulations in parallel.  
  * about announcement numbers and insider sets, and we can run the whole system through them.
  */
  class AnnouncementCascade
  {
    public:
      AnnouncementCascade() = delete; 
      AnnouncementCascade(int date, int days = 365);
      AnnouncementCascade(const IndustryCascade& ind);
      ~AnnouncementCascade();
  
      double RunSingleCascade(bool anti = false) const;
      const std::unordered_map<size_t,double> GetDies() const;

      friend class algorithm::InfluenceMinimizer; 

    private:
      IndustryCascade ic_;
  };
}

#endif
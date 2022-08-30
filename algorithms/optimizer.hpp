#ifndef OPTIM
#define OPTIM
#include "../graphmodel/graphmodel.hpp"
#include "../mkmodel/mkmodel.hpp"
#include "../utils/ioroutines.hpp"
#include "../utils/defs.hpp"
#include "../utils/h_random.hpp"
#include "../simulators/mcsim.hpp"


#include <chrono>
#include <set>
#include <vector>
#include <map>
#include <random>
namespace algorithm
{
using namespace util;
using namespace simulator;
using namespace markov;
using namespace graphmodel;
class Optimizer
{ 
public: 
  
  Optimizer(Simulator* sim, IoR* zed);
  ~Optimizer();
  std::pair<double,double> Optimize(double p = 0.5, double q = 0.5);
  // Optimize per company, find p and q
  std::pair<double,double> CompanyOptimize(int k, double p = 0.5, double q = 0.5);

  // Optimize per company, find p and q
  std::pair<double,double> CompanyOptimizeQ(int k, double p = 0.5, double q = 0.5);

  // Calculate the error for the constant model
  std::tuple<double,double,double> CalculateError(const DatePmap& dm, double p, double q, const DatePmap& deriv);

  // Calculate the error for the whole model
  std::tuple<double,double,double> CalculateError(const DatePmap& dm, double p, double q, int k, const DatePmap& deriv);

  // Calculate the error per company
  std::tuple<double,double,double> CalculateCompanyError(const DatePmap& dm, double p, double q, int k, const DatePmap& deriv);

  // Calculate the Q given a dm and p and company k and q vector
  std::tuple<double,double,double> CalculateCompanyErrorQ(const DatePmap& dm, double p, std::vector<int>& qv, double q, int k, const DatePmap& deriv);
  
  void FisherTest();
    
  
private:
  double norm(const std::vector<double>& x);
  double avg(const std::vector<double>& x);
  double dot(const std::vector<double>& x, const std::vector<double>& y);
  double sum(const std::vector<double>& x, const std::vector<double>& y);
  double qexp(const std::vector<double>& x, const std::vector<double>& y);
  double EvaluateDeviation(const int ns);
  void InitializeDM(DatePmap& dm);
  void InitializeUp();
  void CalculateInvestorR(); 

  void CalculateDM(DatePmap& dm, double p, int n, DatePmap& deriv); 
  void CalculateDM(DatePmap& dm, double p, int n, DatePmap& deriv, int k); 
  void InitializeCY();

  simulator::Simulator* sim_;
  IoR* zed_;
  AnnouncementDates andat_;
  AnnouncementDict andi_;
  NumTimeDict insiders_;
  OwnerFractionDict frac_;
  int n_companies_;
  int n_nodes_;
  std::set<int> years_; 
  std::unordered_map<int,std::vector<double>> up_;
  std::vector<double> investorR_; 
  std::unordered_map<int,std::set<int>> c_years_; 
  std::unordered_map<int,std::unordered_map<int,int>> c_year_first_;
  // This will contain the profitability p-value inside, according to fisher test
  std::unordered_map<int,std::vector<double>> fisher_p_;
  // This will contain the p-value of trading rate being uniform. We assume it is poisson. 
  std::unordered_map<int,std::vector<double>> rate_p_; 
 
};

}
#endif

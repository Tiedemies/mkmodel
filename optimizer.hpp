#ifndef OPTIM
#define OPTIM
#include "graphmodel.hpp"
#include "mkmodel.hpp"
#include "ioroutines.hpp"
#include "mcsim.hpp"
#include <map>
#include "defs.hpp"
#include <random>
#include "h_random.hpp"
#include "mcsim.hpp"
#include <chrono>
#include <set>
#include <vector>

class Optimizer
{
public: 
  Optimizer(Simulator* sim, IoR* zed);
  ~Optimizer();
  std::pair<double,double> Optimize(double p = 0.5, double q = 0.5);
  // Optimize per company, find p and q
  std::pair<double,double> CompanyOptimize(int k, double p = 0.5, double q = 0.5);

  // Calculate the error for the constant model
  std::tuple<double,double,double> CalculateError(const DatePmap& dm, double p, double q, const DatePmap& deriv);

  // Calculate the error for the whole model
  std::tuple<double,double,double> CalculateError(const DatePmap& dm, double p, double q, int k, const DatePmap& deriv);

  // Calculate the error per company
  std::tuple<double,double,double> CalculateCompanyError(const DatePmap& dm, double p, double q, int k, const DatePmap& deriv);
  
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

  Simulator* sim_;
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
};


#endif

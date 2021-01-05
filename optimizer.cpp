#include "graphmodel.hpp"
#include "mkmodel.hpp"
#include "ioroutines.hpp"
#include "optimizer.hpp"
#include "mcsim.hpp"
#include <map>
#include "defs.hpp"
#include <random>
#include "h_random.hpp"
#include "mcsim.hpp"
#include <chrono>
#include <math.h>
#define FIXED_FLOAT(x) std::fixed <<std::setprecision(2)<<(x)


Optimizer::Optimizer(Simulator* sim, IoR* zed)
{
  zed_ = zed; 
  zed->ReadCompanyDictionary();
  zed->ReadNodeDictionary();
  frac_ = zed->ReadFractions(); 
  andi_ = zed->ReadAnnouncements();
  andat_ = zed->GetDates(); 
  sim_ = sim;
  sim_->Simplify();
  n_companies_ = zed->max_com_;
  n_nodes_ = zed->max_tra_;
  years_ = zed->years_; 
}

Optimizer::~Optimizer()
{
  //void
}


// Sämplää 30 pistettä yrityksistä. TODO: muuta niin, että ne on eri yrityksistä satunnaisia pisteitä. 
double Optimizer::EvaluateDeviation(const int ns)
{
  int k = 0;
  int date = c_year_first_.at(k).at(2006);
  const auto& inside = sim_->GetInside(date,k);
  // Calculate vectors:
  std::vector<std::vector<double>> ps;
  for(int i = 0; i < 30; ++i)
  {
    ps.push_back(sim_->Simulate(inside, date, ns).first);
  } 
  // Calculare averages
  std::vector<double> average;
  for(int i = 0; i < ps[0].size(); ++i)
  {
    average.push_back(0.0);
    for(int j = 0; j < 30; ++j)
    {
      average.at(i) += ps[j][i] / 30.0;
    }
  }
  // Calculate variance
  std::vector<double> var;
  for(int i = 0; i < ps[0].size(); ++i)
  {
    var.push_back(0.0);
    for(int j = 0; j < 30; ++j)
    {
      var.at(i) += pow(average[i] - ps[j][i], 2) / 29.0;
    }
  }
  return avg(var); 
}

// Calculate the error over the whole dm
std::tuple<double,double,double> 
Optimizer::CalculateError(const DatePmap& dm, double p, double q, const DatePmap& deriv)
{
  double e = 0;
  double dq = 0;
  double dp = 0;
  std::cerr << "calculating error from dm with " << dm.size() << " entries\n";
  for (int k = 0; k<= n_companies_;++k)
  {
    if (dm.find(k) == dm.end() || up_.find(k) == up_.end())
    {
      continue;
    }
    auto D = deriv.at(k);
    for (auto& P:  dm.at(k)) 
    {
      for (int i = 0; i < P.second.size(); ++i)
      {
        const double& pp = P.second.at(i); 
        const double& dip = D.at(P.first).at(i); 
        e = e + pow((up_.at(k)).at(i) - (q*pp + 0.5*(1-pp)), 2);
        dp += 2*(up_.at(k).at(i) - (q*pp + 0.5*(1-pp)))*dip;
        dq += 2*(up_.at(k).at(i) - (q*pp + 0.5*(1-pp)))*pp;
      } 
    }
  }
  return std::make_tuple(e,dp,dq);
}

// Calculate the error for a given company. 
std::tuple<double,double,double> 
Optimizer::CalculateCompanyError(const DatePmap& dm, double p, double q, int k, 
                          const DatePmap& deriv)
{
  double e = 0;
  double dq = 0;
  double dp = 0;
  if (dm.find(k) == dm.end() || up_.find(k) == up_.end())
  {
    throw std::domain_error("company has no data.\n");
  }
  auto D = deriv.at(k);
  for (auto& P:  dm.at(k)) 
  {
    for (int i = 0; i < n_nodes_; ++i)
    {
      // Default rate is the rate over *all other* companies
      double def_rate = 0.5;
      try
      {
        def_rate = (investorR_.at(i)*n_companies_ - up_.at(k).at(i))/(n_companies_-1);
      }
      catch(const std::exception& e)
      {
        std::cerr << e.what() << " at " << i << '\n';
      }      
      const double& pp = P.second.at(i); 
      const double& dip = D.at(P.first).at(i); 
      e = e + pow((up_.at(k)).at(i) - (q*pp + def_rate*(1-pp)), 2);
      dp += 2*(up_.at(k).at(i) - (q*pp + def_rate*(1-pp)))*dip;
      dq += 2*(up_.at(k).at(i) - (q*pp + def_rate*(1-pp)))*pp;
    } 
  }
  return std::make_tuple(e,dp,dq);
}

double Optimizer::norm(const std::vector<double>& x) 
{
  double res = 0;
  for (double a: x)
  {  
    res += a*a;
  }
  return res; 
}

double Optimizer::avg(const std::vector<double>& x)
{
  double res = 0.0;
  const double n = (double) x.size();
  for(int i = 0; i < n; ++i)
  {
    res += x[i];
  }
  return res/n;
}

double Optimizer::qexp(const std::vector<double>& u, const std::vector<double>& p)
{
  const int n = u.size();
  double up = 0;
  double down = 0;
  for(int i = 0 ; i<n ; ++i)
  {
    up += (u.at(i) - 0.5*(1 - p.at(i)))*p.at(i);
    down += p.at(i)*p.at(i);
  }
  return up/down; 
}

void Optimizer::InitializeCY()
{
  for (int comp = 0; comp <= n_companies_; ++comp)
    {
      std::unordered_map<int,int> first;
      if (andi_.find(comp) == andi_.end())
      {
	      continue;
      }
      bool c_d = false;
      for (int date: andi_.at(comp))
      {
	      int year = date/10000;
	      c_years_[comp].insert(year);
	      if (first.find(year) == first.end())
	      {
	        first[year] = date;
	      }
      }
      c_year_first_[comp] = first;
    }
}


void Optimizer::InitializeDM(DatePmap& dm)
{
  for (int y: years_)
  {
    CompanyPmap cmap;
    for( int k = 0; k <= n_companies_; ++k )
    {
      Pmap foo(n_nodes_+1,0.0);
      cmap[k] = foo;
    }
    dm[y] = cmap;
  }
}

void Optimizer::InitializeUp()
{
  investorR_.clear();
  double avg = 0.0; 
  double diver = 0.0; 
  for (int i = 0; i <= n_nodes_; ++i)
  {
    auto X = frac_.find(i);
    double sum = 0.0;
    double div = 0.0;
    if (X == frac_.end())
    {
      investorR_.push_back(0.0);
      continue;
    }
    auto C = frac_.at(i); 
    for (int k = 0; k <= n_companies_; ++k)
    {
      auto Y = C.find(k);
      if (Y == C.end())
      {
	      continue;
      }	
      if (up_.find(k) == up_.end())
      {
	      std::vector<double> Kv(n_nodes_+1, 0.0);
	      up_[k] = Kv;
      }
      if(std::get<0>(C.at(k)) != 0)
      {
	      (up_.at(k)).at(i) = ( (double) std::get<1>(C.at(k))) / (double) std::get<2>(C.at(k));
        sum += up_.at(k).at(i);
        div += 1.0;
      }
      else
      {
        div += 1.0;
      }
      
    }
    investorR_.push_back(sum/div);
    avg += (sum/div);
    diver += 1.0; 
  }
  std::cerr << "avg profitability: " << avg/diver << "\n";
}

void Optimizer::CalculateDM(DatePmap& dm, double p, int n, DatePmap& deriv)
{
  sim_->SetTransProp(p); 
  #pragma omp parallel for
  for (int k = 0;k <= n_companies_; ++k)
  {
    if (zed_->cids_.find(k) == zed_->cids_.end())
    {
      continue;
    } 
    //std::cerr << "simulating company " << k << ": " << zed_->cids_.at(k) << "\n";
    CompanyPmap bar1;
    CompanyPmap bar2;
    auto c_iter = c_year_first_.find(k);
    if (c_iter == c_year_first_.end())
    {
      //std::cerr << "no announcements, skipping\n";
      continue;
    }
    auto& first = c_year_first_.at(k); 
    for (int y: c_years_.at(k))
    {
      int date = first.at(y);
      //std::cerr << "date: " << date << "\n";
      const auto& inside = sim_->GetInside(date,k);
      //std::cerr << "simulating..\n";
      std::pair<Pmap,Pmap> foo = sim_->Simulate(inside, date, n);  
      bar1[y] = foo.first;
      bar2[y] = foo.second; 
      //std::cerr << "simulated.\n";
    }
    dm[k] = bar1;
    deriv[k] = bar2; 
  }
  //std::cerr << "done.\n";
}

/*
* Calculate the datepmap with the given propagation probability, given company k.
*/
void 
Optimizer::CalculateDM(DatePmap& dm, double p, int n, DatePmap& deriv, int k)
{
  sim_->SetTransProp(p); 
  CompanyPmap bar1;
  CompanyPmap bar2;
  auto c_iter = c_year_first_.find(k);
  if (c_iter == c_year_first_.end())
  {
    std::cerr << "no announcements, skipping\n";
    return;
  }
  auto& first = c_year_first_.at(k);
  for (auto mm = first.begin(); mm != first.end(); ++mm)
  {
    int y = mm->first; 
    int date = mm->second;
    const auto& inside = sim_->GetInside(date,k);
    std::pair<Pmap,Pmap> foo = sim_->Simulate(inside, date, n);  
    bar1[y] = foo.first;
    bar2[y] = foo.second; 
  }
  dm[k] = bar1;
  deriv[k] = bar2; 
}

// The bulk optimization function. 
std::pair<double,double> 
Optimizer::Optimize(double p, double q)
{
  int n = 500; // Number of iterations per year
  // Initialize up_ by iterating over nodes (i) and company (k)
  std::cerr << "Initializing UP...\n";
  InitializeUp();
  std::cerr << "Initializing C_years.. \n";
  InitializeCY();
  // Initialize dm. This is the P's that we will have for different dates. 
  DatePmap dm;
  DatePmap deriv;
  int l = 0;
  double st = EvaluateDeviation(n);
  std::cerr << "average deviation " << st << "\n";
  while (l < 10)
  {
    std::cerr << "Calculating DM for point [" << p << "," << q << "]\n";
    CalculateDM(dm,p,n,deriv);
    auto r = CalculateError(dm, p,q,deriv);
    //std::cerr << "Error: " << std::get<0>(r) << "\n";
    p = p-0.01;
    q = q-0.01;
    ++l;
  }
  return std::make_pair(p,q); 
}

// Optimizer function for singe company
std::pair<double,double> 
Optimizer::CompanyOptimize(int k, double p, double q)
{
  int n = 500; // Number of iterations per year
  int max_iter = 100; // Maximum number of iterations. 
  // Initialize up_ by iterating over nodes (i) and company (k)
  std::cerr << "Initializing UP...\n";
  InitializeUp();
  std::cerr << "Initializing C_years.. \n";
  InitializeCY();
  // Initialize dm. This is the P's that we will have for different dates. 
  DatePmap dm;
  DatePmap deriv;

  // l is our iteration counter
  int l = 0;

  std::cerr << "Optimizing for company " << zed_->cids_.at(k) << "\n";

  // these keep track of the values
  double best_p = p;
  double best_q = q;
  double best_error = -1;
  double current_p = p;
  double current_q = q;
  double current_error = -1;
  double T = 1000;
  double scale = 0.1;
 
  double top_p = 1.0;
  double bottom_p = 0.0;
   double top_q = 1.0;
  double bottom_q = 0.0;  
  bool flip = false;

  // We initialize a new random variable for all runs, so that the seed is initialized fresh
  Random hrnd; 

  while (l < 3000)
  {
    ++l;
    CalculateDM(dm,p,n,deriv,k);
    auto r = CalculateCompanyError(dm, p, q, k, deriv);
    double err = std::get<0>(r);
    double dp = std::get<1>(r);
    double dq = std::get<2>(r);
    if (best_error < 0 || best_error > err)
    {
      best_error = err; 
      current_error = err;
      best_p = p;
      current_p = p;
      best_q = q;
      current_q = q;
      // std::cerr << "found better point at " << p << "," << q << "\n";
    } 
    p += (0.1-0.2*hrnd.get());
    q += (0.1-0.2*hrnd.get());
    if (p > 1.0)
    {
      p = 1.0;
    }
    if (q > 1.0)
    {
      q = 1.0;
    }
    if (p < 0.0)
    {
      p = 0.0;
    }
    if (q < 0.0)
    {
      q = 0.0;
    }
  }
  std::cerr << "best point found was: \n[" << best_p << " , " << best_q << "]\n";
  std::cerr << "error at the point was " << best_error << "\n";
}

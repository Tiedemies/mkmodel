#ifndef I_CASCADE
#define I_CASCADE
#include "hidden_cascade.hpp"
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

/* 
 * The IndustryCascade class is a wrapper for a hidden cascade that also includes the information 
 * about announcement numbers and insider sets, and we can run the whole system through them.
 */
class IndustryCascade
{
  public:
    IndustryCascade(HiddenCascade h);
    IndustryCascade(int date, int days = 365);
    ~IndustryCascade();
    std::vector<double> RunTotal();
    std::vector<double> RunSingle(int i);
    void RandomizeCoefs(double mu);
    void SetSimulationCoef(int n);
    NodeTransactionTable RunGeneration();

    struct Transaction
    {
      int node_;
      int comp_;
      int date_;
      double price_;
      double ref_price_; 
      double vol_;
      bool inside_;
      bool operator<(const Transaction& rhs) const;
    };
    

  private:
    void InitializeFoo();
    void InitializeAnnouncements();
    void InitializeWindows();
    void InitializeTransactions();
    void InitializeProbabilities(); 
    int simulation_coefficient_; 
    int date_; 
    int days_;
    int n_node_;
    int n_comp_; 
    double def_p_;
    IoR foo_; 
    HiddenCascade hc_;
    std::vector<std::vector<int>> insiders_;
    std::vector<int> num_announcements_;
    std::vector<std::vector<int>> announcement_days_; 
    std::vector<std::vector<bool>> in_window_days_;
    std::vector<std::vector<bool>> trading_days_; 
    std::vector<Transaction> transacts_;
    std::vector<std::vector<double>> out_trading_prob_;
    std::vector<std::vector<double>> in_trading_prob_;
};



#endif
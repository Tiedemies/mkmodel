#ifndef STAT_TEST
#define STAT_TEST
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
#include <boost/numeric/ublas/matrix.hpp>


class StatTester
{
    public:
    StatTester();
    ~StatTester();
    // set window size:
    void SetWindowSize(int);
    // Calculate profits inside and outside windows
    void CreateProfitWindows(); 
    // Calculate the numbers of inside and outside days 
    void CreateInsideDayWindows(); 
    // Do the hypergeometric tests of inside vs outside. 
    void TestHyperG(); 

    void GenerateDataMatrix(); 

    private:
    // Days in the window
    int window_size_; 

    // Io-class handle
    IoR ior_;
    // How much profit is inside
    NodeProfitMap profit_inside_;
    // How many days are inside for each ISIN
    IsinCountMap days_inside_; 
    // How much profit is outside
    NodeProfitMap profit_outside_;
    // How many days are outside for each ISIN  
    IsinCountMap days_outside_;

    //Nodewise count of how many transactions done out/in:
    NodeIsinCountMap n_profit_inside_;
    NodeIsinCountMap n_profit_outside_;

    // P-values for the node-company test
    NodeCompanyPValueMap hg_pvalues_; 
    int num_hg_tests_; 

    //Numbers:
    int num_traders_;
    int num_companies_;
    int num_transactions_; 

    // The data matrix for regression 
    boost::numeric::ublas::matrix<double> X_; 
    boost::numeric::ublas::matrix<double> y_;
};
#endif
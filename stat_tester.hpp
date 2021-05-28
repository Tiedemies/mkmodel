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


class StatTester
{
    public:
    StatTester();
    ~StatTester();
    void SetWindowSize(int);
    void CreateProfitWindows(); 
    
    private:

    int window_size_; 
    IoR ior_;
    NodeTransactionTable delayed_profits_; 

};
#endif
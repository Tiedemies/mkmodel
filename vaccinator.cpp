#include "defs.hpp"
#include "ioroutines.hpp"
#include "graphmodel.hpp"
#include "hidden_cascade.hpp"
// #include "stat_routines.hpp"
#include <chrono>
#include <iostream>
#include <cmath>

#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices_nan.txt"
#define TRANSACTFILE "table_transacitions.txt"
#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"

#define REFDATE 20101231

double st_error(std::vector<double> input)
{
    double n = static_cast<double>(input.size());
    double avg = std::accumulate(input.begin(), input.end(), 0.0) / n;
    double err = 0.0;
    for(auto k: input)
    {
        err += (avg - k)*(avg-k);
    }
    err /= (n-1);
    return err; 
}

int main()
{
    // Start timing.
    auto start = std::chrono::high_resolution_clock::now(); 
    // Initialize a metagraph.      
    /*
    IoR zed;
    zed.SetAnnouncementDirectory(ADIR);
    zed.SetCompanyDictionaryFile(CDFILE);
    zed.SetNodeDictionaryFile(NDFILE);
    zed.ReadCompanyDictionary();
    zed.ReadNodeDictionary();
    zed.ReadAnnouncements();
    */
    std::cout << "Intialize IoR and Metagraph for " << REFDATE << "\n";
    IoR foo(REFDATE);
    std::cout << "Initializing monograph for " << REFDATE << "\n";
    MonoGraph* bar = foo.GetGraph(REFDATE);
    if (bar == nullptr)
    {
        std::cerr << "Error: Null metagraph \n";
        return -1;
    }
    double p = 0.08;
    double tp = 0.1;
    double fp = 0.08;
    std::cout << "Initializing hidden cascade model \n"; 
    HiddenCascade cas(bar,p,fp,tp);
    auto stop = std::chrono::high_resolution_clock::now();   
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << "Initialization took " << duration.count()/1000 << "ms total\n";

    auto board = bar->GetInsider(6);
    std::vector<int> insiders(board.begin(), board.end());
    size_t n = 5000;
    int k = 30;
    cas.SetSimulationN(n);
    std::cout << "starting simulations. \n";
    std::vector<double> num_active;
    double sum =0.0;
    num_active.resize(k,0.0);
    for (int i = 0; i < k; ++i)
    {
        num_active[i] = cas.Simulate(insiders);
        sum += num_active[i]/k; 
    }
    auto stop2 = std::chrono::high_resolution_clock::now();   
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-stop);
    

    std::cout << "On average " << sum << " nodes were activated \n";
    std::cout << "Sigma: " << std::sqrt(st_error(num_active)) << "\n";
    std::cout <<  n*k << " simulations took " << duration2.count()/1000 << "ms total\n";
    std::cout << duration2.count()/(n*k) << " microseconds per simulation \n";
    std::cout << static_cast<double>(n*k)/(duration2.count()/1000000) << " simulations per second \n";
    return 0;
}
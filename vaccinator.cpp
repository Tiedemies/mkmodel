#include "defs.hpp"
#include "ioroutines.hpp"
#include "graphmodel.hpp"
#include "hidden_cascade.hpp"
#include <chrono>
#include <iostream>

#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices_nan.txt"
#define TRANSACTFILE "table_transacitions.txt"
#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"

#define REFDATE 20101231

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
    std::cout << "Reading Metagraph for " << REFDATE << "\n";
    MetaGraph foo(NWDIR, REFDATE);
    std::cout << "Initializing monograph for " << REFDATE << "\n";
    MonoGraph* bar = foo.GetGraph(REFDATE);
    if (bar == nullptr)
    {
        std::cerr << "Error: Null metagraph \n";
        return -1;
    }
    double p = 0.2;
    double tp = 0.1;
    double fp = 0.08;
    std::cout << "Initializing hidden cascade model \n"; 
    HiddenCascade cas(bar,p,fp,tp);
    auto stop = std::chrono::high_resolution_clock::now();   
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << "Initialization took " << duration.count()/1000 << "ms total\n";

    auto board = bar->GetInsider(5);
    std::vector<int> insiders(board.begin(), board.end());
    size_t n = 1500;
    int k = 100;
    cas.SetSimulationN(n);
    std::cout << "starting simulations. \n";
    for (int i = 0; i < k; ++i)
    {
        cas.Simulate(insiders,true);
    }
    auto stop2 = std::chrono::high_resolution_clock::now();   
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-stop);
    std::cout <<  n*k << " simulations took " << duration2.count()/1000 << "ms total\n";
    std::cout << duration2.count()/(n*k) << " microseconds per simulation \n";
    std::cout << static_cast<double>(n*k)/(duration2.count()/1000000) << "simulations per second \n";
    return 0;
}
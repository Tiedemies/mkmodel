#include "defs.hpp"
#include "ioroutines.hpp"
#include "graphmodel.hpp"
#include "hidden_cascade.hpp"
// #include "stat_routines.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>



#define REFDATE 20101231
#define PLOTFILE "plotfile.txt"

double st_error(std::vector<double> input)
{
    const int& n = static_cast<int>(input.size());
    const double& avg = std::accumulate(input.begin(), input.end(), 0.0) / n;
    double err = 0.0;
    #pragma omp parallel for reduction(+:err)
    for(int i = 0; i <= n-1;++i)
    {
        err += (avg - input[i])*(avg - input[i]);
    }
    err /= (n-1);
    return std::sqrt(err); 
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
    std::cout << "Initializing hidden cascade model \n"; 
    
    auto stop = std::chrono::high_resolution_clock::now();   
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << "Initialization took " << duration.count()/1000 << "ms total\n";

    auto board = bar->GetInsider(6);
    std::vector<int> insiders(board.begin(), board.end());
    size_t n = 500;
    int k = 30;
    std::cout << "starting simulations. \n";
    std::vector<double> num_active;
    double sum =0.0;
    num_active.resize(k,0.0);
    
    std::ofstream out;
    out.open(PLOTFILE);
    for (int delta = 0; delta < 100; ++delta)
    {
        double p = static_cast<double>(delta)/100.0 + 0.001;
        double tp = 0.3;
        double fp = 0.08;
        HiddenCascade cas(bar,p,fp,tp);
        int mina = std::numeric_limits<int>::max();
        int maxa = std::numeric_limits<int>::min();
        sum = 0.0;
        for (int i = 0; i < k; ++i)
        {
            cas.SetSimulationN(n);
            num_active[i] = cas.Simulate(insiders);
            sum += num_active[i]/k; 
            mina = std::min(mina, cas.GetMinActivated());
            maxa = std::max(maxa, cas.GetMaxActivated());
        }
        out <<  p << " " << mina << " " << maxa << " " << sum << " " << st_error(num_active) << "\n";
    }
    out.close();
    auto stop2 = std::chrono::high_resolution_clock::now();   
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-stop);
    
    return 0;
    /*
    std::cout << "On average " << sum << " nodes were activated \n";
    std::cout << "Min: " << mina << ", Max: " << maxa << "\n";
    std::cout << "Sigma: " << st_error(num_active) << "\n";
    std::cout <<  n*k << " simulations took " << duration2.count()/1000 << "ms total\n";
    std::cout << duration2.count()/(n*k) << " microseconds per simulation \n";
    std::cout << static_cast<double>(n*k)/(duration2.count()/1000000) << " simulations per second \n";
    return 0;
    */
}
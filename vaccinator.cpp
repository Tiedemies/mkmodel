#include "utils/defs.hpp"
#include "utils/ioroutines.hpp"
#include "graphmodel/graphmodel.hpp"
#include "simulators/hidden_cascade.hpp"
#include "simulators/industry_cascade.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>



#define REFDATE 20081231
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
    std::cout << "Intialize IndustryCascade for " << REFDATE << "\n";
    IndustryCascade foo(REFDATE);

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
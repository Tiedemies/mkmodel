#include "utils/defs.hpp"
#include "utils/ioroutines.hpp"
#include "graphmodel/graphmodel.hpp"
#include "simulators/hidden_cascade.hpp"
#include "simulators/industry_cascade.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>

double sum(const std::vector<double>& input)
{
    return std::accumulate(input.cbegin(), input.cend(), 0.0);
}

double avg(const std::vector<double>& input)
{
    return sum(input) / input.size();
}

double st_error(std::vector<double> input)
{
    const int& n = static_cast<int>(input.size());
    const double& avgs = avg(input);
    double err = 0.0;
    #pragma omp parallel for reduction(+:err)
    for(int i = 0; i <= n-1;++i)
    {
        err += (avgs - input[i])*(avgs - input[i]);
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
    std::cout << "Intialize IndustryCascade for " << VREFDATE << "\n";
    IndustryCascade foo(VREFDATE);
    std::cerr << "all created.\n";

    auto tvec = foo.RunTotal();
    std::cerr << "Simulations run. Total number of activations: " << sum(tvec) << "\n";
    std::cerr << "Avg: " << avg(tvec) << ", stde:" << st_error(tvec) << "\n";
    return 0;
}
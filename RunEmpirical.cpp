#include "graphmodel.hpp"
#include "mkmodel.hpp"
#include "ioroutines.hpp"
#include "mcsim.hpp"
#include "optimizer.hpp"
#include <map>
#include "defs.hpp"
#include <random>
#include "h_random.hpp"
#include <iostream>
#include <fstream>
#include <new>
#include "mcsim.hpp"
#include <chrono>



#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"

#define FFILE "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/fraction of trades/fraction of trades before announcement period/mean/fraction_of_trades_before_and_outside_announcement_period_5d.txt"



int main()
{

  //[[maybe_unused]] double p = 0.7;
  //[[maybe_unused]] double w = 1;
  //[[maybe_unused]] double wt = 1;
  //[[maybe_unused]] double q = 0.4;
  //[[maybe_unused]] double noise = 0.1;
  auto start = std::chrono::high_resolution_clock::now(); 

  // auto stop1 = std::chrono::high_resolution_clock::now();
  //[[maybe_unused]] auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop1-start); 
 
 
  IoR zed;
  zed.SetAnnouncementDirectory(ADIR);
  zed.SetCompanyDictionaryFile(CDFILE);
  zed.SetNodeDictionaryFile(NDFILE);
  zed.SetFractionsFile(FFILE);
  zed.ReadCompanyDictionary();
  zed.ReadNodeDictionary();
  zed.ReadAnnouncements();
 

  MetaGraph foo(NWDIR);
  Random bar;
  MkModel mk(bar);
  Simulator sim(&foo,mk);
  Optimizer opt(&sim, &zed);
  std::cout << "Optimizer initialized.\n"; 
  auto stop2 = std::chrono::high_resolution_clock::now();   
  auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-start);
  std::cout << "I/O took " << duration2.count()/1000 << "ms total\n";
  /*
  for (auto entry: zed.cids_)
  {
    stop2 = std::chrono::high_resolution_clock::now();
    try
    {
      opt.CompanyOptimize(entry.first, 0.5,0.5);
    }
    catch(std::exception& e)
    {
      std::cerr << e.what();
    }
    auto stop3 = std::chrono::high_resolution_clock::now(); 
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3-stop2);
    std::cout << "Optimize took " << duration3.count()/1000 << "ms for company " << entry.second << "\n";
  }
  */
  opt.CompanyOptimize(zed.cnames_.at("Nokia"), 0.5,0.5);
  auto stop3 = std::chrono::high_resolution_clock::now();     
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3-stop2);
  std::cout << "Optimize took " << duration3.count()/1000 << "ms for company \n";
  opt.CompanyOptimizeQ(zed.cnames_.at("Nokia"), 0.5,0.5);
  auto stop4 = std::chrono::high_resolution_clock::now();     
  auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4-stop3);
  std::cout << "Optimize Q took " << duration4.count()/1000 << "ms for company \n";

  return 0;

}

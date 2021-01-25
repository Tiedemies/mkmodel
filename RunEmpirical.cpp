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


// Here define the network directory, announcement directory, company dictionary file. 
#define NWDIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/output/"
#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"
#define CDFILE "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/dicts/comp_dict_unified_INTERNAL.txt"
#define NDFILE "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/dicts/nodes_dict_unified_INTERNAL.txt"

#define FFILE "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/fraction of trades/fraction of trades before announcement period/mean/fraction_of_trades_before_and_outside_announcement_period_5d.txt"



int main()
{

  double p = 0.7;
  double w = 1;
  double wt = 1;
  double q = 0.4;
  double noise = 0.1;
  auto start = std::chrono::high_resolution_clock::now(); 

  auto stop1 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop1-start); 
 
  IoR zed;
  zed.SetAnnouncementDirectory(ADIR);
  zed.SetCompanyDictionaryFile(CDFILE);
  zed.SetNodeDictionaryFile(NDFILE);
  zed.SetFractionsFile(FFILE);
  zed.ReadCompanyDictionary();
  zed.ReadNodeDictionary();
  zed.ReadAnnouncements();
 

  MetaGraph foo(NWDIR,zed.GetDates());
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
  opt.CompanyOptimize(zed.cnames_.at("Kone"), 0.5,0.5);
  auto stop3 = std::chrono::high_resolution_clock::now();     
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3-stop2);
  std::cout << "Optimize took " << duration3.count()/1000 << "ms for company \n";
  opt.CompanyOptimizeQ(zed.cnames_.at("Kone"), 0.5,0.5);
  auto stop4 = std::chrono::high_resolution_clock::now();     
  auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4-stop3);
  std::cout << "Optimize Q took " << duration4.count()/1000 << "ms for company \n";

  return 0;

}

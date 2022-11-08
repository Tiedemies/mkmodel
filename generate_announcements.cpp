#include "utils/defs.hpp"
#include "utils/ioroutines.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/date_time/gregorian/gregorian.hpp>


#define CSVANNOUNCE "announcements.csv"

int main()
{
  util::IoR iors; 
  iors.ReadCompanyDictionary();
  iors.ReadAnnouncements();
  iors.ReadAnnounceTable();
  std::ofstream out;
  out.open(CSVANNOUNCE);
  using namespace boost::gregorian; 
  date ref_time(from_simple_string(REFDAY));

  for (auto an_pair: iors.an_table_)
  {
    std::string isin = an_pair.first;
    //std::cerr << "isin: " << isin << "\n";
    if(iors.isin_company_.find(isin) == iors.isin_company_.end())
    {
      continue;
    }
    std::string cname = iors.isin_company_[isin];
    std::cerr << "company: " << cname << "\n";
    if (iors.cnames_.find(cname) == iors.cnames_.end())
    {
      continue;
    }
    int cnum = iors.cnames_.at(cname);
    std::cerr << "Number: " << cnum << ";" << an_pair.second.size() << " announcements \n";
    for (auto an: an_pair.second)
    {
      date cur_date = ref_time+days(an);
      out << cname << ";" << cnum <<  ";" << "0;" << cur_date << "\n";
    }
    out.flush();
  }
  for (auto an_pair: iors.an_table_sc_)
  {
    std::string isin = an_pair.first;
    //std::cerr << "isin: " << isin << "\n";
    if(iors.isin_company_.find(isin) == iors.isin_company_.end())
    {
      continue;
    }
    std::string cname = iors.isin_company_[isin];
    //std::cerr << "company: " << cname << "\n";
    if (iors.cnames_.find(cname) == iors.cnames_.end())
    {
      continue;
    }
    int cnum = iors.cnames_.at(cname);
    //std::cerr << "Number: " << cnum << "\n";
    for (auto an: an_pair.second)
    {
      date cur_date = ref_time+days(an);
      out << cname << ";" << cnum <<  ";" << "1;" << cur_date << "\n";
    }
    out.flush();
  }
  out.close(); 
  return 0;
}
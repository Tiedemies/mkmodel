#include "ioroutines.hpp"
#include "defs.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include <set>

#define ANFILE "/announcement_periods_with_insiders.txt"

IoR::IoR()
{
  //void, default
}

IoR::~IoR()
{
  //void default
}

void IoR::SetAnnouncementDirectory(std::string andir)
{
  andir_ = andir;
}
void IoR::SetCompanyDictionaryFile(std::string cdictfile)
{
  cdictfile_ = cdictfile;
}
void IoR::SetNodeDictionaryFile(std::string ndictfile)
{
  ndictfile_ = ndictfile;
}
void IoR::SetFractionsFile(std::string fractionsfile)
{
  fractionsfile_ = fractionsfile;
}

void IoR::ReadCompanyDictionary()
{
  if (cdictfile_.empty())
  {
    std::cerr << "Warning, company dictionary file not defined. \n";
    return;
  }
  max_com_ = 0;
  std::ifstream in;
  in.open(cdictfile_);
  int c;
  if (in.is_open())
  {
    in >> c;
    while(!in.eof() && in.good())
    {
      std::string namestr; 
      in.get();
      char x = in.get();
      while(x != '\n')
      {
	if (!iscntrl(x)) 
	  {
	    namestr.push_back(x);
	  }
	x = in.get(); 
      }
      if (c > max_com_)
      {
	max_com_ = c; 
      }
      cids_[c] = namestr;
      cnames_[namestr] = c;
      //std::cerr << "Gave company: \n" << namestr << "\n a number.\n";
      //std::cerr << namestr.size() << " characters. \n"; 
      in >> c;
    }
  }
  else
  {
    std::cerr << "Failed to open dictionary file.\n";
  }
  
}

void IoR::ReadNodeDictionary()
{
  int node;
  int owner_id;
  int n = 0;
  if (ndictfile_.empty())
  {
    std::cerr << "Warning, node dictionary file not defined. \n";
    return;
  }
  max_tra_ = 0;
  std::ifstream in;
  in.open(ndictfile_);
  while(!in.eof() && in.good())
  {  
    in >> node;
    in.get();
    in >> owner_id;
    if (!in.good())
    {
      continue;
    }
    nodedict_[node] = owner_id;
    nodeinvdict_[owner_id] = node;
    ++n;
    if (node > max_tra_)
    {
      max_tra_ = node;
    }
  }
  std::cerr << "Read " << n << " nodes\n";
  in.close();
}

std::vector<int> IoR::ReadTimes(std::istream& in)
{
  std::vector<int> togo;
  years_.clear();
  char next = in.get();
  next = in.get();
  int yy;
  int mm; 
  int dd;
  int date;
  while (next != ']' && in.good())
  {
    in >> yy;
    next = in.get();    
    //std::cerr << "read " << yy << " next character: " << next << "\n";
    in  >> mm;
    next = in.get(); 
    //std::cerr << "read " << mm << " next character: " << next << "\n";
    in >> dd;
    //std::cerr << "read " << dd << "\n";
    date = 10000*yy + 100*mm + dd;
    //std::cerr << date << " extracted \n";
    togo.push_back(date); 
    years_.insert(yy);
    while (next != ']' && next != '(')
    {
      next = in.get(); 
    }
  }
  while (next != '\n' && in.good())
    {
      next = in.get(); 
    }
  return togo; 
}
AnnouncementDict IoR::ReadAnnouncements()
{
  AnnouncementDict togo;
  if (andir_.empty())
  {
    std::cerr << "Warning, no announcement directory, null announcements\n";
    return togo;
  }
  if (cids_.empty())
  {
    std::cerr << "Warning, no company dictionary, null announcements\n";
    return togo;
  }
  std::ifstream in;
  in.open(andir_ + ANFILE);
  if (!in.is_open())
  {
    std::cerr << "Warning, announcement file not found, null announcements\n";
    return togo;
  }
  while (!in.eof() && in.good())
  {
    std::string namestr;
    char next = in.get();
    while(next != ';' && in.good())
    {
      namestr.push_back(next);
      next = in.get();
    }
    if (cnames_.find(namestr) == cnames_.end() && in.good())
    {
      std::cerr << "Warning, company " << namestr <<  " not found.\n";
      ReadTimes(in);
    }
    else if (in.good())
      {
	int i = cnames_.at(namestr);
	togo[i] = ReadTimes(in);
	std::copy(togo[i].begin(), togo[i].end(), std::inserter(an_dates_, an_dates_.end()));
      }
  }
  std::cerr << "announcements done\n";
  return togo; 

} 

AnnouncementDates 
IoR::GetDates()
{
  return an_dates_; 
}


OwnerFractionDict 
IoR::ReadFractions()
{
  OwnerFractionDict togo;
  std::ifstream in;
  in.open(fractionsfile_);
  if (!in.is_open())
  {
    std::cerr << "Warning, fractions file not found, null fractions\n";
    return togo;
  }
  std::string line; // assume there is a header line
  std::getline(in, line);
  std::set<int> nfind;
  std::set<int> find ;
  while (!in.eof() && in.good())
  {
    int owner_id;
    in >> owner_id;
    int node_id;
    if (nodeinvdict_.find(owner_id) == nodeinvdict_.end())
    {
	    node_id = owner_id;
	    nfind.insert(owner_id);
    }
    else
    {
	    node_id = nodeinvdict_.at(owner_id); 
	    find.insert(owner_id);
    }
    // std::cerr << "read owner nr: " << owner_id << " ... ";
    std::string cname; 
    in.get();
    char next = in.get();
    while (next != ';' && !in.eof() && in.good())
    {
      cname.push_back(next);
      next = in.get(); 
    } 
    //std::cerr << " got company " << cname << " ... ";
    int n_trade;
    in >> n_trade;
    in.get();
    double n_profit;
    in >> n_profit;
    //std::cerr << " got " << n_profit << " inside profit days ... ";
    in.get();
    int n_an;
    in >> n_an;
    int count_del = 0;
    while (count_del < 6 && !in.eof() && in.good())
    {
      char x = in.get();
      while(x != ';' &&  !in.eof() && in.good())
	    {
	      x = in.get(); 
	    }
      ++count_del;
    }
    int n_market_days;
    in >> n_market_days;
    //std::cerr << " got " << n_market_days << " market days ... ";
    char x = in.get();
    while(x!=';' && !in.eof() && in.good())
    {
      x = in.get();
    }
    double n_trade_outside;
    in >> n_trade_outside;
    //std::cerr << " got " << n_trade_outside << " outside days ... ";
    x = in.get();
    while(x!=';' && !in.eof() && in.good())
    {
      x = in.get();
    }
    double n_profit_outside;
    in >> n_profit_outside;
    //std::cerr << " of which " << n_trade_outside << " days are profitable. \n ";
    if (cnames_.find(cname) == cnames_.end())
    {
	    std::cerr << "Warning, company " << cname << " is missing.\n"; 
	    continue;
    }
    int cnum = cnames_.at(cname);
    if (togo.find(node_id) == togo.end())
    {
      CompanyFractions newfrac;
      togo[node_id] = newfrac;
    }
    if (n_trade_outside < n_profit_outside)
    {
      std::cerr << "Warning: Company " << cname <<  " has trader " << owner_id << " do " << n_trade_outside << " trades and " << n_profit_outside << " profit ??\n";
      throw std::logic_error("foobar");
    }
    FractionEntry foo = std::make_tuple(n_trade,n_profit,n_an,n_market_days,n_trade_outside, n_profit_outside);
    
    (togo[node_id])[cnum] = foo;
    //throw std::logic_error("test");
  }
  //std::cerr << "fractions done. " << nfind.size() << " unidentified nodes " << find.size() << " identified.\n"; 
  return togo; 
}


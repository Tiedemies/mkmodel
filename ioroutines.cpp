#include "ioroutines.hpp"
#include "defs.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include <set>
#include <new>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#define TRANS_ERRFILEHANDLE "bad_transactions.txt"


IoR::IoR()
{
  tablesdir_ = TABLEDIR; 
  atfile_ = ANFILE;
  ptfile_ = PRICEFILE;
  ttfile_ = TRANSACTFILE;
  graphdir_ = NWDIR;
  metag_ = new MetaGraph(graphdir_);
}

IoR::~IoR()
{
  //void default
}

void IoR::SetAnnouncementDirectory(const std::string& andir)
{
  andir_ = andir;
}
void IoR::SetCompanyDictionaryFile(const std::string& cdictfile)
{
  cdictfile_ = cdictfile;
}
void IoR::SetNodeDictionaryFile(const std::string& ndictfile)
{
  ndictfile_ = ndictfile;
}
void IoR::SetFractionsFile(const std::string& fractionsfile)
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

AnnouncementDict 
IoR::ReadAnnouncements()
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

void 
IoR::ReadDates()
{
  trade_days_.clear(); 
  std::cerr << "day read initializing \n";
  std::ifstream in;
  in.open(std::string(TABLEDIR) + std::string(DATEFILE));
  if (!in.is_open())
  {
    std::cerr << "Warning, date table file not found, null table set\n";
  }
  using namespace boost::gregorian; 
  // Reference date:
  date ref_time(from_simple_string(REFDAY));
  // Loop reads the line. 
  int count = 0;
  while (!in.eof() && in.good())
  {
    // First read the date. 
    std::string datestr = ReadNext(in);   
    date tt;
    try
    {
      tt = from_simple_string(datestr);
      int days = (tt-ref_time).days();
      if (days >= 0)
      {
        trade_days_.insert(days);
      }
    }
    catch(const std::exception& e)
    {
      continue;
    }
  }
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
// Tables:
void 
IoR::SetTablesDirectory(const std::string& tdir)
{
  tablesdir_ = tdir;
}

void 
IoR::SetAnnounceTableFile(const std::string& atfile)
{
  atfile_ = atfile;
}

void 
IoR::SetInsiderTableFile(const std::string& itfile)
{
  itfile_ = itfile;
}

void 
IoR::SetPriceTableFile(const std::string& ptfile)
{
  ptfile_ = ptfile;
}

void 
IoR::SetTransactionTableFile(const std::string& ttfile)
{
  ttfile_ = ttfile; 
}

// Now read everything: 
void 
IoR::ReadTables()
{
  ReadAnnounceTable();
  // ReadInsiderTable();
  ReadPriceTable();
  ReadTransactionTable();
}

// Read table from tablesfile. 
void
IoR::ReadAnnounceTable()
{
  // This is for the announcment dates.   

  // Company Dictionary needs to be read. 
  if (cids_.empty())
  {
    //std::cerr << "Warning, no company dictionary, reading now";
    //void ReadCompanyDictionary();
  }

  /// Read the tables now 
  using namespace boost::gregorian; 
  std::ifstream in;
  in.open(tablesdir_ +  atfile_);
  if (!in.is_open())
  {
    std::cerr << "Warning, announcement table file not found, null table set\n";
    return; 
  }
  // Read one line  which is the header and throw away. 
  std::string line;
  std::getline(in, line);
  date ref_time(from_simple_string(REFDAY));
  // Loop reads the line. 
  int count = 0;
  while (!in.eof() && in.good())
  {
    // First read the date. 
    std::string datestr = ReadNext(in);   
    date tt;
    try
    {
       tt = from_simple_string(datestr);
    }
    catch(const std::exception& e)
    {
      // std::cerr << e.what() << '\n';
      break;
    }
    int days = (tt-ref_time).days();
   
    // Skip NetDays;
    ReadNext(in);
    // Skip nextDay
    ReadNext(in);
    // Read Company Name
    std::string cname = ReadNext(in);
    // Read ISIN
    std::string isin = ReadNext(in);

    auto isin_it = isin_company_.find(isin);
    if (isin_it != isin_company_.end())
    {
        if (isin_it->second != cname)
        {
          std::cerr << "WARNING: duplicate name " << cname << " for isin " << isin << ", old name " << isin_it->second << "\n";
        }
    }
    else
    {
      isin_company_[isin] = cname; 
    }
    
    // Skip id; relatedto:
    ReadNext(in);
    ReadNext(in);
    std::string sched = ReadNext(in);
    if (sched == "Non-scheduled")
    {
      ++count;
      an_table_[isin].push_back(days);
    }
    SkipLine(in);
  }
  //Debug output:
  std::cerr << "Read " << count <<  " announcements for " << an_table_.size() << " isin codes \n";
  for (auto XY: an_table_)
  {
    std::sort(XY.second.begin(), XY.second.end());
  } 
}

void
IoR::ReadPriceTable()
{
 
  /// Read the tables now 
  using namespace boost::gregorian; 
  std::ifstream in;
  // trade_days_.clear(); 
  in.open(tablesdir_ +  ptfile_);
  if (!in.is_open())
  {
    std::cerr << "Warning, price table file not found, null table set\n";
  }
  // Read one line  which is the header and throw away. 
  std::string line;
  std::getline(in, line);
  // Reference date:
  date ref_time(from_simple_string(REFDAY));
  // Loop reads the line. 
  int count = 0;
  while (!in.eof() && in.good())
  {
    // First read the date. 
    std::string datestr = ReadNext(in);   
    date tt;
    try
    {
       tt = from_simple_string(datestr);
    }
    catch(const std::exception& e)
    {
      // std::cerr << e.what() << '\n';
      break;
    }
    int days = (tt-ref_time).days();
    if (days >= 0)
    {
      // trade_days_.insert(days);
    }
    // ReadIsin;
    std::string isin = ReadNext(in);
    isin_set_.insert(isin);
    std::string pricestr = ReadNext(in);
    double price = std::nan("missing");
    if (days >= 0)
    {
      try
      {
        price = std::stod(pricestr);/* code */
      }
      catch(const std::exception& e)
      {
        price = std::nan("missing");
      }
      pr_table_.AddPCompanyDayPrice(isin, days, price);
      ++count; 
    }
    SkipLine(in);
  }
  //Debug output:
  std::cerr << "Read " << count <<  " prices for " << pr_table_.size() << " isin codes.\n";
  // pr_table_.Sort();  
}

void
IoR::ReadTransactionTable()
{
  // This is for the announcment dates.   

  /// Read the tables now 
  using namespace boost::gregorian; 
  std::ifstream in;
  std::ofstream out;
  in.open(tablesdir_ +  ttfile_);
  out.open(TRANS_ERRFILEHANDLE);
  if (!in.is_open())
  {
    std::cerr << "Warning, transaction table file not found, null table set\n";
    return;
  }
  // Read one line  which is the header and throw away. 
  std::string line;
  std::getline(in, line);
  // Reference date:
  date ref_time(from_simple_string(REFDAY));
  // Loop reads the line. 
  int count = 0;
  int foo = 0;
  int non_market = 0;
  while (!in.eof() && in.good())
  {
    // First read the date. 
    std::string datestr = ReadNext(in);   
    date tt;
    try
    {
       tt = from_simple_string(datestr);
    }
    catch(const std::exception& e)
    {
      // std::cerr << e.what() << '\n';
      break;
    }
    int days = (tt-ref_time).days();
   
    // ReadIsin;
    // skip owner_id; we use node id;
    std::string owner = ReadNext(in); 
    std::string isin = ReadNext(in);
    std::string price_str = ReadNext(in);
    double price = std::stod(price_str);
    std::string volume_str = ReadNext(in);
    double volume = std::stod(volume_str);
    // ReadNext(in);
    //skip vol_price
    std::string vol_price = ReadNext(in);
    std::string temp = ReadNext(in);
    int nodeid = -1; 
    try
    {
      nodeid = std::stoi(temp);
    }
    catch(const std::exception& e)
    {
      ++foo; 
      // std::cerr << e.what() << " from " << temp << '\n';
      out << datestr << ";" << owner << ";" << isin << ";" << price_str << ";" << volume_str 
          << ";" << vol_price << ";" << temp << ";";
      SkipLine(in,out);
      continue; 
    }
    int x = std::stoi(ReadNext(in));
    // Skip non-market trades
    if (x != 1)
    {
      continue;
      ++non_market;
    }
    
    if (tr_table_.find(nodeid) == tr_table_.end())
    {
      TransactionTable novel;
      tr_table_[nodeid] = novel;
    }
    if (days >= 0)
    {
      tr_table_[nodeid].AddPCompanyDayPriceTransaction(isin, days, price, volume);
      ++count;
    } 
  }
  //Debug output:
  std::cerr << "Read " << count <<  " transactions for " << tr_table_.size() << " traders. " << foo << " bad.\n";
  std::cerr << non_market << " total non-market transactions\n"; 
}

std::string
IoR::ReadNext(std::istream& in)
{
  std::string togo;
  char next = in.get();
  while(next != ';' && next != '\n' && in.good())
  {
    togo.push_back(next);
    next = in.get();
  }
  return togo;
}

void 
IoR::SkipLine(std::istream& in)
{
  char x = in.get();;
  while ((x !='\n') && in.good())
  {
    x = in.get();
  }
}

void 
IoR::SkipLine(std::istream& in, std::ostream& out)
{
  char x = in.get();
  out.put(x);
  while ((x !='\n') && in.good())
  {
    x = in.get();
    out.put(x);
  }
}


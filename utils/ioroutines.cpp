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

namespace util
{
  using namespace graphmodel;
double 
PilkkuPisteeksi(const std::string&  luku)
{
  double cumulative = 0;
  double aftercom = 0.1;
  bool comma = false;
  for (auto sit = luku.begin(); sit != luku.end(); ++sit)
  {
    if (!comma && std::isdigit(*sit))
    {
      int i = (*sit - '0');
      cumulative = 10*cumulative + i;
    }
    else if (comma && std::isdigit(*sit))
    {
      int i = (*sit - '0');
      cumulative += aftercom*i;
      aftercom = aftercom * 0.1; 
    }
    else
    {
      comma = true;
    }
  }
  return cumulative;
}

std::set<int> 
BracketArrayToSet(const std::string& array)
{
  std::string nums = array.substr(1,array.size()-2); 
  std::stringstream input(nums);
  std::set<int> togo;
  std::string item; 
  while (std::getline(input, item, ',')) {
      togo.insert(std::stoi(item));
  }
  return togo;
}

IoR::IoR()
{
  tablesdir_ = TABLEDIR; 
  atfile_ = ANFILE;
  ptfile_ = PRICEFILE;
  ttfile_ = TRANSACTFILE;
  hhfile_ = HHTFILE;
  graphdir_ = NWDIR;
  indexfile_ = INDEXFILE;
  reasonfile_ = REASONFILE;
  cdictfile_ = CDFILE;
  std::cerr << "Creating metagraph\n";
  metag_ = new MetaGraph(graphdir_, -1);
  std::cerr << "metagraph created\n";
}

IoR::IoR(const IoR& rhs):
max_com_(rhs.max_com_), 
max_tra_(rhs.max_tra_),
years_(rhs.years_),
trade_days_(rhs.trade_days_),
isin_set_(rhs.isin_set_),
cnames_(rhs.cnames_),
cids_(rhs.cids_),
isin_company_(rhs.isin_company_),
an_table_(rhs.an_table_),
relativemap_(rhs.relativemap_),
metag_(rhs.metag_), 
nodedict_(rhs.nodedict_),
nodeinvdict_(rhs.nodeinvdict_),
dateindex_(rhs.dateindex_),
an_dates_(rhs.an_dates_),
announcements_(rhs.announcements_)
{
  // void
}

IoR::IoR(bool none)
{
  // void
}

IoR::IoR(int year)
{
  tablesdir_ = TABLEDIR;
  andir_ = ADIR; 
  atfile_ = ANFILE;
  ptfile_ = PRICEFILE;
  ttfile_ = TRANSACTFILE;
  hhfile_ = HHTFILE;
  graphdir_ = NWDIR;
  indexfile_ = INDEXFILE;
  reasonfile_ = REASONFILE;
  cdictfile_ = CDFILE;
  //std::cerr << "Creating metagraph\n";
  metag_ = new MetaGraph(graphdir_, year);
  //std::cerr << "metagraph created\n";
}

IoR::~IoR()
{
  //void default
}

void IoR::SetAnnouncementDirectory(const std::string& andir)
{
  andir_ = andir;
}

void IoR::SetIndexFile(const std::string& indexfile)
{
  indexfile_ = indexfile;
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
  //[[maybe_unused]] int count = 0;
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
  std::cerr << "days read\n";
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
  //ReadInsiderTable();
  ReadPriceTable();
  ReadTransactionTable();
  ReadHouseHoldTransactionTable();
  ReadIndex();
}

void
IoR::ReadIndex()
{
  /// Read the tables now 
  using namespace boost::gregorian; 
  std::ifstream in;
  // trade_days_.clear(); 
  std::cerr << "Reading index values\n";
  in.open(indexfile_);
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
      std::cerr << e.what() << '\n';
      continue;
    }
    int days = (tt-ref_time).days();
    if (days >= 0)
    {
      // trade_days_.insert(days);
    }
    // skip high
    ReadNext(in);
    //skip low
    ReadNext(in);
    // Median
    std::string pricestr = ReadNext(in);
    double price = std::nan("missing");
    if (days >= 0)
    {
      price = PilkkuPisteeksi(pricestr);
      pr_table_.AddPCompanyDayPrice("index", days, price);
      ++count; 
    }
    SkipLine(in);
  }
  //Debug output:
  std::cerr << "Read " << count <<  "index values \n"; 
  /*
  for (auto xx: dateindex_)
  {
    std::cerr << "day: " << xx.first << ", index value: " << xx.second << "\n";
  }
  */
}
  


// Read table from tablesfile. 
void
IoR::ReadAnnounceTable()
{
  // This is for the announcment dates.   
  std::cerr << "Reading announcements\n";
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
  int sccount = 0;
  while (!in.eof() && in.good())
  {
    // First read the date. 
    std::string datestr = ReadNext(in); 
    //std::cerr << "At date: " << datestr << "\n";  
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
    //std::cerr << "Got nam: " << cname << "\n";
    // Read ISIN
    std::string isin = ReadNext(in);
    //std::cerr << "Got isin: " << isin << "\n";


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
    std::string temp = ReadNext(in);
    //std::cerr << "skipped: " << temp;
    temp = ReadNext(in);
    //std::cerr << "skipped: " << temp;
    std::string sched = ReadNext(in);
    //std::cerr << "Got sch: " << sched << "\n";

    if (sched == "Non-scheduled")
    {
      ++count;
      an_table_[isin].push_back(days);
    }
    else if (sched == "Scheduled by date")
    {
      ++sccount;
      an_table_sc_[isin].push_back(days);
    }
    SkipLine(in);
  }
  //Debug output:
  std::cerr << "Read " << count <<  " unscheduled and " << sccount << " scheduled announcements for " << an_table_.size() << " isin codes \n";
  for (auto XY: an_table_)
  {
    std::sort(XY.second.begin(), XY.second.end());
  } 
}

void
IoR::ReadPriceTable()
{
 
  std::cerr << "Reading price table\n";
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

void
IoR::ReadHouseHoldTransactionTable()
{
  // This is for the announcment dates.   
  std::cerr << "Reading household transactions\n";
  /// Read the tables now 
  using namespace boost::gregorian; 
  std::ifstream in;
  std::ofstream out;
  in.open(tablesdir_ +  hhfile_);
  out.open(TRANS_ERRFILEHANDLE);
  if (!in.is_open())
  {
    std::cerr << "Warning, household table file not found, null table set\n";
    return;
  }
  // Read one line  which is the header and throw away. 
  std::string line;
  std::getline(in, line);
  // Reference date:
  date ref_time(from_simple_string(REFDAY));
  // Loop reads the line. 
  int count = 0;
  int nodes = 0;
  // //[[maybe_unused]] int foo = 0;
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
    int ownerid = -1;
    try
    {
      ownerid = -1*std::stoi(owner);
    }
    catch(const std::exception& e)
    {
      std::cerr << e.what() << '\n';
      continue;
    }
     
    std::string isin = ReadNext(in);
    std::string price_str = ReadNext(in);
    double price = -1;
    try
    {
      price = std::stod(price_str);
    }
    catch(const std::exception& e)
    {
      std::cerr << e.what() << '\n';
      SkipLine(in,out);
    }
    std::string volume_str = ReadNext(in);
    double volume = std::stod(volume_str);
    // ReadNext(in);
    //skip vol_price
    std::string vol_price = ReadNext(in);
    // Skip node id. 
    std::string temp = ReadNext(in);
    // Does it have a node id? 
    if (! temp.empty() ) 
    {
      ++nodes;
      SkipLine(in);
      continue;
    }

   // Market thing?
    int x = std::stoi(ReadNext(in));
    // Skip non-market trades
    if (x != 1)
    {
      ++non_market;
      continue;
    }
    
    if (tr_table_.find(ownerid) == hh_table_.end())
    {
      TransactionTable novel;
      tr_table_[ownerid] = novel;
    }
    if (days >= 0)
    {
      tr_table_[ownerid].AddPCompanyDayPriceTransaction(isin, days, price, volume);
      ++count;
    } 
  }
  //Debug output:
  std::cerr << "Read " << count <<  " household transactions for " << hh_table_.size() << " housholds, " << nodes << " skipped insiders \n";
  std::cerr << non_market << " total non-market transactions\n"; 
}

AnnouncementDates
IoR::GetDates()
{
  return an_dates_;
}

void 
IoR::ReadReasonsTable()
{
   // This is for the announcment dates.   
  std::cerr << "Reading insider labels file";
  /// Read the tables now 
  std::ifstream in;
  in.open(tablesdir_ +  reasonfile_);
  // throw away the first line; 
  std::string line;
  std::getline(in, line);
  while (!in.eof() && in.good())
  {
    // First read actual id 
    std::string id_str = ReadNext(in);
    if (in.eof())
    {
      break; 
    }
    int id = std::stoi(id_str);
    // Second read to company id
    std::string cname = ReadNext(in);
    std::string basis = ReadNext(in);
    std::string related_id = ReadNext(in);
    if (related_id.empty())
    {
      SkipLine(in);
      continue;
    }
    int node_id = std::stoi(ReadNext(in));
    std::string related_node = ReadNext(in);
    std::set<int> bas_set = BracketArrayToSet(basis);
    relativemap_[node_id] = std::stoi(related_node);
  }
}

MonoGraph* 
IoR::GetGraph(int date)
{
  return metag_->GetGraph(date);
} 

void
IoR::PrintISINTable(std::ofstream& out)
{
  for (auto isin_t: isin_company_)
  {
    std::string isin = isin_t.first;
    std::string cname = isin_t.second;
    int cnum = cnames_[cname];
    out << isin << ";" << cname << ";" << cnum << "\n";
  }
}
}
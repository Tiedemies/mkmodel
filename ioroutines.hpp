#ifndef IO_R
#define IO_R
#include "defs.hpp"
#include "graphmodel.hpp"
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<set>

class IoR
{
public: 
  IoR();
  ~IoR();
  void SetAnnouncementDirectory(const std::string& andir);
  // Tables
  void SetTablesDirectory(const std::string& tdir);
  void SetAnnounceTableFile(const std::string& atfile);
  void SetInsiderTableFile(const std::string& itfile);
  void SetPriceTableFile(const std::string& ptfile);
  void SetTransactionTableFile(const std::string& ttfile);
  void ReadTables();
  void ReadAnnounceTable();
  void ReadInsiderTable();
  void ReadPriceTable();
  void ReadTransactionTable();

  // Dictionaries
  void SetCompanyDictionaryFile(const std::string& cdictfile);
  void SetNodeDictionaryFile(const std::string& ndictfile);
  void SetFractionsFile(const std::string& fractionsfile);
  void ReadCompanyDictionary(); 
  void ReadNodeDictionary();
  OwnerFractionDict ReadFractions();
  NumTimeDict ReadInsiders();
  AnnouncementDict ReadAnnouncements();
  void ReadDates();
  int max_com_; 
  int max_tra_;
  std::set<int> years_;
  std::set<int> trade_days_; 
  std::set<std::string> isin_set_; 
  CompanyNameDict cnames_;
  CompanyNameInvDict cids_;
private:
  std::string ReadNext(std::istream& in);
  void SkipLine(std::istream& in);
  void SkipLine(std::istream& in, std::ostream& out);
  std::vector<int> ReadTimes(std::istream& in);
  std::string graphdir_; 
  std::string ndictfile_;
  std::string andir_;
  std::string tablesdir_;
  std::string itfile_;
  std::string atfile_;
  std::string ttfile_;
  std::string ptfile_; 

  MetaGraph* metag_; 

  std::string cdictfile_;
  std::string fractionsfile_;
  NodeDict nodedict_;
  NodeDict nodeinvdict_;
  AnnouncementDates an_dates_; 
  AnnouncementDict announcements_;

  // Tables:
  AnnouncementTable an_table_;
  PriceTable pr_table_; 
  NodeTransactionTable tr_table_; 

  // Transformations:
  IsinCompanyMap isin_company_; 

  friend class StatTester; 
};

#endif

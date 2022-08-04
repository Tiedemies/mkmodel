#ifndef IO_R
#define IO_R
#include "defs.hpp"
#include "../graphmodel/graphmodel.hpp"
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<set>

class IoR
{
public: 
  IoR();
  IoR(int year);
  IoR(bool none);
  ~IoR();
  void SetAnnouncementDirectory(const std::string& andir);
  // Tables
  void SetTablesDirectory(const std::string& tdir);
  void SetAnnounceTableFile(const std::string& atfile);
  void SetInsiderTableFile(const std::string& itfile);
  void SetReasonsTableFile(const std::string& rcfile);
  void SetPriceTableFile(const std::string& ptfile);
  void SetTransactionTableFile(const std::string& ttfile);
  void ReadTables();
  void ReadAnnounceTable();
  void ReadInsiderTable();
  void ReadPriceTable();
  void ReadTransactionTable();
  void ReadHouseHoldTransactionTable();
  void ReadReasonsTable();
  void ReadIndex();
  MonoGraph* GetGraph(int year);

  // Dictionaries
  void SetCompanyDictionaryFile(const std::string& cdictfile);
  void SetNodeDictionaryFile(const std::string& ndictfile);
  void SetFractionsFile(const std::string& fractionsfile);
  void SetIndexFile(const std::string& indexfile);
  void ReadCompanyDictionary(); 
  void ReadNodeDictionary();
  OwnerFractionDict ReadFractions();
  NumTimeDict ReadInsiders();
  AnnouncementDict ReadAnnouncements();
  void ReadDates();
  AnnouncementDates GetDates();
  int max_com_; 
  int max_tra_;
  std::set<int> years_;
  std::set<int> trade_days_; 
  std::set<std::string> isin_set_; 
  CompanyNameDict cnames_;
  CompanyNameInvDict cids_;
  IsinCompanyMap isin_company_; 
  AnnouncementTable an_table_;
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
  std::string hhfile_;
  std::string indexfile_;
  std::string reasonfile_; 
  MetaGraph* metag_; 

  std::string cdictfile_;
  std::string fractionsfile_;
  NodeDict nodedict_;
  NodeDict nodeinvdict_;
  IndexMap dateindex_; 
  AnnouncementDates an_dates_; 
  AnnouncementDict announcements_;

  // Tables:
  
  AnnouncementTable an_table_sc_;
  PriceTable pr_table_; 
  NodeTransactionTable tr_table_; 
  NodeTransactionTable hh_table_;

  // Transformations:
  

  friend class StatTester; 
};

#endif

#ifndef IO_R
#define IO_R
#include "defs.hpp"
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
  void SetAnnouncementDirectory(std::string andir);
  void SetCompanyDictionaryFile(std::string cdictfile);
  void SetNodeDictionaryFile(std::string ndictfile);
  void SetFractionsFile(std::string fractionsfile);
  void ReadCompanyDictionary(); 
  void ReadNodeDictionary();
  OwnerFractionDict ReadFractions();
  NumTimeDict ReadInsiders();
  AnnouncementDict ReadAnnouncements();
  AnnouncementDates GetDates();
  int max_com_; 
  int max_tra_;
  std::set<int> years_;
  CompanyNameDict cnames_;
  CompanyNameInvDict cids_;
private:
  std::vector<int> ReadTimes(std::istream& in);
  std::string ndictfile_;
  std::string andir_;
  std::string cdictfile_;
  std::string fractionsfile_;
  NodeDict nodedict_;
  NodeDict nodeinvdict_;
  AnnouncementDates an_dates_; 
  AnnouncementDict announcements_;




};

#endif

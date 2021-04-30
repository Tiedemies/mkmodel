// Type definitions for the optimizer

#ifndef DEFS
#define DEFS 

#include<unordered_map>
#include<vector>
#include<string>
#include<set>

// Number of days
#define NDAYS 5


// Dictionaries 
// Fraction entry: <#possible inside, #profit pi, #announcements, #market days, #outside trades, #outside profit trades>
typedef std::tuple<int,double,int,int,double, double> FractionEntry;
typedef std::unordered_map<int,int> NodeDict;
typedef std::unordered_map<int,std::vector<int>> NumTimeDict;
typedef std::unordered_map<int,FractionEntry> CompanyFractions;
typedef std::unordered_map<int,CompanyFractions> OwnerFractionDict;
typedef std::unordered_map<std::string,int> CompanyNameDict;
typedef std::unordered_map<int, std::string> CompanyNameInvDict; 
typedef std::unordered_map<int, std::vector<int>> AnnouncementDict; 
typedef std::set<int> AnnouncementDates; 
typedef std::vector<double> Pmap;
typedef std::unordered_map<int,Pmap> CompanyPmap; 
typedef std::unordered_map<int,CompanyPmap> DatePmap; 

// Graph
typedef std::vector<int> Alist; 
typedef std::unordered_map<int,Alist> AdjLists;

//Mkmod el 
typedef std::pair<int,int> MkTrans;
typedef std::vector<MkTrans> MkAlist;
typedef std::vector<MkAlist> MkLists;

#endif

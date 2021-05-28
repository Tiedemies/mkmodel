// Type definitions for the optimizer

#ifndef DEFS
#define DEFS 

#include<unordered_map>
#include<vector>
#include<string>
#include<set>

// Number of days
#define NDAYS 5
#define REFDAY "2001-01-01 12:00:00"

// The default directories and files. 
#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices.txt"
#define TRANSACTFILE "table_transacitions.txt"


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
typedef std::vector<int> AnnouncementVector; 
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

//tables 
typedef std::unordered_map<std::string, std::vector<int>> AnnouncementTable;
typedef std::vector<std::pair<int,double>> DatePriceVector; 
typedef std::unordered_map<std::string,DatePriceVector> InternalPriceTable; 
// Pricetable and functionality
class PriceTable
{
    private:
        InternalPriceTable pt_;
        bool sorted_;
    public:
        PriceTable();
        ~PriceTable();
        void AddPCompanyDayPrice(const std::string& cname, int day, double price);
        double GetCompanyDayPrice(const std::string& cname, int day, int offset) const;
        void Sort();
        int size();
        friend class StatTester; 
};

typedef std::unordered_map<int,PriceTable*> NodeTransactionTable;  


#endif

// Type definitions for the optimizer

#ifndef DEFS
#define DEFS 

/*** Definitions for directories and files 
 * EDIT THIS SECTION IF YOUR FILES ARE ELSEWHERE ***/ 

// The default directories and files. 
#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define DATEFILE "tradingDates.csv"
//#define PRICEFILE "table_prices.txt"
#define PRICEFILE "table_prices_nan.txt"
#define TRANSACTFILE "table_transacitions.txt"
#define HHTFILE "table_households.txt"
#define REASONFILE "table_insiderships.txt"


// Here define the network directory
#define NWDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/output/"

// Define the announcement directory
#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"

// Dictionaries 
#define CDFILE "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/dicts/comp_dict_unified_INTERNAL.txt"
#define NDFILE "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/unified_networks/dicts/nodes_dict_unified_INTERNAL.txt"

//Outputfiles
#define XINCSV "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/X_in.csv"
#define YINCSV "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/Y_in.csv"

#define XHHCSV "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/X_hh.csv"
#define YHHCSV "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/Y_hh.csv"

//index file
#define INDEXFILE "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/indexdata.csv"

/*** Basic definition of constants ***/

// Number of days
#define NDAYS 5
// Reference day is 5 days before the first announcement
#define REFDAY "2005-05-21 12:00:00"
 

// Define significance levels:
#define SIGNIFICANT 0.05
#define VERY_SIGNIFICANT 0.01
#define EXTRA_SIGNIFICANT 0.001

// What return is considered profit:
#define PROFIT_THRESHOLD 0.0

//Vaccinator begin date:
#define VREFDATE 20071231

#define GREFDATE 20051231

#include<map>
#include<unordered_map>
#include<vector>
#include<string>
#include<set>
#include<boost/assert.hpp>
#include<algorithm>
#include<math.h>
#include<numeric>


// forward declaration

namespace data
{
  class StatTester; 
}

namespace simulator
{
  class IndustryCascade; 
}

namespace util
{
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
typedef std::map<int,double> IndexMap;

typedef std::unordered_map<std::string, std::string> IsinCompanyMap;

// Graph
typedef std::vector<int> Alist; 
typedef std::unordered_map<int,Alist> AdjLists;
typedef std::set<int> Alist2;

//Mkmod el 
typedef std::pair<int,int> MkTrans;
typedef std::vector<MkTrans> MkAlist;
typedef std::vector<MkAlist> MkLists;

//tables 
typedef std::unordered_map<std::string, std::vector<int>> AnnouncementTable;
typedef std::map<int,double> DatePriceMap;
typedef std::vector<std::tuple<int,double,double>> DatePriceVolumeVector; 
typedef std::unordered_map<std::string,DatePriceMap> InternalPriceTable;
typedef std::unordered_map<std::string,DatePriceVolumeVector> InternalTransactionTable; 

// Pricetable and functionality
class PriceTable
{
  private:
    InternalPriceTable pt_;
    // bool sorted_;
  public:
    PriceTable();
    ~PriceTable();
    void AddPCompanyDayPrice(const std::string& cname, int day, double price);
    double GetCompanyDayPrice(const std::string& cname, int day, int offset) const;
    std::pair<int,double> GetFirstChangePrice(const std::string& cname, int day, int offset, const std::set<int>& dates) const;
    // void Sort();
    int size();
    friend class data::StatTester;
    friend class simulator::IndustryCascade; 
};

class TransactionTable
{
  private:
    InternalTransactionTable pt_;
    bool sorted_;
  public:
    TransactionTable();
    ~TransactionTable();
    void AddPCompanyDayPriceTransaction(const std::string& cname, int day, double price, double volume);
    void Sort();
    int size();
    friend class data::StatTester; 
    friend class simulator::IndustryCascade;
};



typedef std::unordered_map<int,TransactionTable> NodeTransactionTable;  

typedef std::unordered_map<std::string, std::vector<double>> CompanyProfitTable; 
typedef std::unordered_map<int, CompanyProfitTable> NodeProfitMap;
typedef std::unordered_map<std::string, int> IsinCountMap; 
typedef std::unordered_map<int, IsinCountMap> NodeIsinCountMap; 
typedef std::unordered_map<int, double> CompanyPValueMap;
typedef std::unordered_map<int, CompanyPValueMap> NodeCompanyPValueMap;


double sum(const std::vector<double>& input);


double avg(const std::vector<double>& input);

double w_avg(const std::vector<double>& input, const std::vector<double>& weights);

double st_error(std::vector<double> input);


class NLargest
{
  public:
    NLargest() = delete;
    NLargest(int n);
    ~NLargest();
    void Insert(int i, double v);
    const std::vector<int>& GetIndeces() const;
    const std::vector<double>& GetValues() const; 
    std::pair<int,double> PopSmallest(); 
  private:
    void BuildHeap();
    void Heapify(int cur = 0); 
    const int n_; 
    std::vector<double> values_;
    std::vector<int> indices_;
};

}

#endif

#include "stat_tester.hpp"
#include "defs.hpp"
#include<string>

// 
StatTester::StatTester()
{
    ior_.ReadTables();
    ior_.pr_table_.Sort(); 
    window_size_ = 5; 
}

// Default destructor. 
StatTester::~StatTester()
{

}

void
StatTester::SetWindowSize(int size)
{
    if (size < 1)
    {
        throw std::runtime_error("illegal window size");
    }
    window_size_ = size;
}

void 
StatTester::CreateProfitWindows()
{
    // These are the node transaction and price tables that we use to create the  delayed version. 
    const NodeTransactionTable& transacts = ior_.tr_table_;
    const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;

    for (auto nodepair: transacts)
    {
        const int& node = nodepair.first;
        const PriceTable& c_trans = *(nodepair.second);
        for (auto isin_vector_pair: c_trans.pt_)
        {
            const std::string& isin = isin_vector_pair.first;
            const DatePriceVector& trans = isin_vector_pair.second;
            auto ans_it = ans.find(isin);
            {
                if(ans_it == ans.end())
                {
                    continue;
                }
            }    
            const auto& ansvector = ans_it->second; 
            #pragma omp parallel for
            for (unsigned int i = 0; i < trans.size(); ++i)
            {
                const int& date = trans[i].first;
                const double& price = trans[i].second; 
                const double& refprice = pricetable.GetCompanyDayPrice(isin,date,window_size_);
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                
            }
        }
    }
}
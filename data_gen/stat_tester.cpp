#include "stat_tester.hpp"
#include "../utils/defs.hpp"
#include<string>
#include<set>
#include<iostream>
#include <boost/math/distributions/hypergeometric.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"



// 
StatTester::StatTester()
{
    ior_.ReadDates();
    ior_.ReadTables();
    // ior_.pr_table_.Sort(); 
    ior_.SetCompanyDictionaryFile(CDFILE);
    ior_.SetNodeDictionaryFile(NDFILE);
    ior_.ReadCompanyDictionary();
    ior_.ReadNodeDictionary();
    std::cerr << ior_.trade_days_.size(); 
    profit_window_size_ = 5;
    profit_window_size_2_ = 21;
    inside_window_size_ = 5;
    inside_window_size_2_ = 10;

}

// Default destructor. 
StatTester::~StatTester()
{
    //void
}

void
StatTester::SetWindowSize(int size)
{
    if (size < 1)
    {
        throw std::runtime_error("illegal window size");
    }
    profit_window_size_ = size;
}
void 
StatTester::CreateInsideDayWindows()
{
    // we do this for every isin
    for(auto isin: ior_.isin_set_)
    {
        std::string comp = ior_.isin_company_[isin];
        if (days_inside_.find(comp) == days_inside_.end()) 
        {
            days_inside_[comp] = 0;
            days_outside_[comp] = 0;
        }
        auto ans_it = ior_.an_table_.find(isin);
        if (ans_it == ior_.an_table_.end())
        {
            continue;
        }
        const auto& ans_vector = ans_it->second;
        // iterate over days. 
        for(auto day: ior_.trade_days_)
        {
            auto ans_v_it = std::lower_bound(ans_vector.begin(), ans_vector.end(), day);
            // inside:
            if (*ans_v_it > day && *ans_v_it < day + inside_window_size_)
            {
                ++days_inside_[comp];
            }
            else
            {
                ++days_outside_[comp];
            }
            
        }
    }

}

void 
StatTester::CreateProfitWindows()
{
    // These are the node transaction and price tables that we use to create the  delayed version. 
    // Also: Calculate the numbers and sets for the major cluster. 
    const NodeTransactionTable& transacts = ior_.tr_table_;
    const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;
    num_traders_ = ior_.nodedict_.size();
    num_companies_ = ior_.cnames_.size();
    num_transactions_ = 0;
    int n_inside_p = 0;
    int n_outside_p = 0; 
    int n_reversed = 0;
    // Have to do this otherwise
    for (auto nodepair: transacts)
    {
        const int node = nodepair.first;
        const TransactionTable& c_trans = nodepair.second;
        if (profit_inside_.find(node) == profit_inside_.end())
        {
            CompanyProfitTable novel1;
            CompanyProfitTable novel2;
            profit_inside_[node] = novel1;
            profit_outside_[node] = novel2;
            CompanyProfitTable novelA;
            CompanyProfitTable novelB;
            profit_inside_2_[node] = novelA;
            profit_outside_2_[node] = novelB;
            IsinCountMap novel3;
            IsinCountMap novel4;
            IsinCountMap novel5;
            IsinCountMap novel6;
            n_profit_inside_[node] = novel3;
            n_profit_outside_[node] = novel4;
            n_profit_inside_2_[node] = novel5;
            n_profit_outside_2_[node] = novel6;
        }
        for (auto isin_vector_pair: c_trans.pt_)
        {
            const std::string& isin = isin_vector_pair.first;
            const std::string& comp = ior_.isin_company_[isin];
            const DatePriceVolumeVector& trans = isin_vector_pair.second;
            auto ans_it = ans.find(isin);
            {
                if(ans_it == ans.end())
                {
                    continue;
                }
            } 
            //std::cerr << comp << "\n";
            // ansvector is the announcement dates. 
            if (profit_outside_[node].find(comp) == profit_outside_[node].end())
            {
                std::vector<double> novel1;
                std::vector<double> novel2;
                std::vector<double> novel3;
                std::vector<double> novel4;
                profit_outside_[node][comp] = novel1;
                profit_inside_[node][comp] = novel2;
                profit_outside_2_[node][comp] = novel3;
                profit_inside_2_[node][comp] = novel4;
                n_profit_inside_[node][comp] = 0;
                n_profit_outside_[node][comp] = 0;
            }    
            const auto& ansvector = ans_it->second; 
            for (unsigned int i = 0; i < trans.size(); ++i)
            {
                const int date = std::get<0>(trans[i]);
                const double price = std::get<1>(trans[i]);
                const double volume = std::get<2>(trans[i]);
                std::pair<int, double> refpair;
                std::pair<int, double> refpair2;
                double refprice;
                double refprice2; 
                try
                {
                    refpair = pricetable.GetFirstChangePrice(isin,date,profit_window_size_, ior_.trade_days_); 
                    refprice = refpair.second; 
                    refpair2 = pricetable.GetFirstChangePrice(isin,date,profit_window_size_2_, ior_.trade_days_);
                    refprice2 = refpair2.second;
                }
                catch(const std::exception& e)
                {
                    boost::gregorian::date realday = boost::gregorian::from_simple_string(REFDAY) + boost::gregorian::days(date);
                    std::cerr << "isin:" << isin << ", date:" << realday << ", price:" << price << "\n";
                    std::cerr << e.what() << '\n';
                    continue;
                }
                
                // The next announcement. 

                double ret = refprice/price;
                double ret2 = refprice2/price;
                int sgn = volume > 0.0?1:-1;
                double lret = sgn*log(ret);
                double lret2 = sgn*log(ret2);
            
                //debug
                if (refpair.first > refpair2.first && !std::isnan(refprice2))
                {
                    
                    boost::gregorian::date realday = boost::gregorian::from_simple_string(REFDAY) + boost::gregorian::days(date);
                    boost::gregorian::date realday1 = realday + boost::gregorian::days(refpair.first);
                    boost::gregorian::date realday2 = realday + boost::gregorian::days(refpair2.first); 
                    
                    std::cerr << "isin:" << isin << ", date:" << realday << ", price:" << price << ", refday:" << realday1 << ", refprice:" << 
                        refprice << ", refday 2:" << realday2 << ", refprice2:" << refprice2 << "\n"; 
                    std::cerr << refpair.first << "," << refpair2.first << "\n";
                    
                    // throw std::logic_error("dates reveresed");
                    ++n_reversed;
                }

                ++num_transactions_;
                // std::cerr << "Return: " << ret;
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                // if(std::isnan(lret) || std::isnan(lret2))
                // Is is inside?
                if (ans_v_it != ansvector.end() && *ans_v_it > date && *ans_v_it <= date + inside_window_size_)
                {
                    //std::cerr << "inside;";
                    // Its inside
                    if(!std::isnan(lret))
                    {
                        profit_inside_[node][comp].push_back(lret);
                    }
                    if(!std::isnan(lret2))
                    {
                        profit_inside_2_[node][comp].push_back(lret2);
                    }
                    if (lret > PROFIT_THRESHOLD && !std::isnan(lret))
                    {
                        //std::cerr << "profit;";
                        ++n_profit_inside_[node][comp];
                    }
                    if (lret2 > PROFIT_THRESHOLD && !std::isnan(lret2))
                    {
                        //std::cerr << "profit;";
                        ++n_profit_inside_2_[node][comp];
                    }
                    ++n_inside_p;
                }
                else
                {
                    //std::cerr << "outside;";
                    if (lret > PROFIT_THRESHOLD)
                    {
                        //std::cerr << "profit;";
                        ++n_profit_outside_[node][comp];
                    }
                    if (lret2 > PROFIT_THRESHOLD)
                    {
                        //std::cerr << "profit;";
                        ++n_profit_outside_2_[node][comp];
                    }
                    if(!std::isnan(lret))
                    {
                        profit_outside_[node][comp].push_back(lret);
                    }
                    if(!std::isnan(lret2))
                    {
                        profit_outside_2_[node][comp].push_back(lret2);
                    }
                    ++n_outside_p;
                }
                //std::cerr << "\n";
            }
        }
    }
    //std::cerr << "Inside transactions: " << n_inside_p << ", Outside transactions: " << n_outside_p << "\n";
    //std::cerr << "reversed days: " << n_reversed << "\n"; 
    for (auto x: ior_.hh_table_)
    {
        num_transactions_+=x.second.size();
    }
}

void
StatTester::TestHyperG()
{ 
    num_hg_tests_ = 0; 
    for (auto i_p: profit_inside_)
    {
        // Investor i
        int i = i_p.first;
        CompanyPValueMap imap;
        hg_pvalues_[i] = imap;  
        for (auto k_p: i_p.second)
        {
            std::string cname = k_p.first;
            //std::cerr << "cname: "  << cname; 
            int k = ior_.cnames_[cname];
            //std::cerr << "  k: " << k << "\n";
            int n_in = k_p.second.size();
            int n_out = profit_outside_[i][cname].size();
            // If there are no inside or outside transactions, then nothing should be done. 
            if (n_in == 0 || n_out == 0)
            {
                continue; 
            }
            int n_pin = n_profit_inside_[i][cname];
            int n_pout = n_profit_outside_[i][cname];
            double p = 1.0;
            // std::cerr << n_in << "," << n_pin << " --- " << n_out + n_in << "," << n_pout + n_pin << "\n";
            try
            {
                boost::math::hypergeometric_distribution<double> hg_dist(n_pin + n_pout, n_in, n_in+n_out);
                p =  (1.0 - boost::math::cdf<double>(hg_dist, n_pin)) + boost::math::pdf<double>(hg_dist, n_pin); 
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                std::cerr << n_in << "," << n_pin << " --- " << n_out + n_in << "," << n_pout + n_pin << "\n";
                throw e;
            }
            hg_pvalues_[i][k] = p;
            if (std::isnan(p))
            {
                std::cerr << "Company:" << ior_.cnames_[cname] << ", Investor:" << i << ", inside:" << n_in << ", outside: " << n_out
                << ", profit in:" << n_pin << ", profit out:" << n_pout << ", p-value:" << p << "\n";
            }
            ++num_hg_tests_;
        }
    }

    for (auto i_p: profit_inside_2_)
    {
        // Investor i
        int i = i_p.first;
        CompanyPValueMap imap;
        hg_pvalues_2_[i] = imap;  
        for (auto k_p: i_p.second)
        {
            std::string cname = k_p.first;
            //std::cerr << "cname: "  << cname; 
            int k = ior_.cnames_[cname];
            //std::cerr << "  k: " << k << "\n";
            int n_in = k_p.second.size();
            int n_out = profit_outside_2_[i][cname].size();
            // If there are no inside or outside transactions, then nothing should be done. 
            if (n_in == 0 || n_out == 0)
            {
                continue; 
            }
            int n_pin = n_profit_inside_2_[i][cname];
            int n_pout = n_profit_outside_2_[i][cname];
            double p = 1.0;
            // std::cerr << n_in << "," << n_pin << " --- " << n_out + n_in << "," << n_pout + n_pin << "\n";
            try
            {
                boost::math::hypergeometric_distribution<double> hg_dist(n_pin + n_pout, n_in, n_in+n_out);
                p =  (1.0 - boost::math::cdf<double>(hg_dist, n_pin)) + boost::math::pdf<double>(hg_dist, n_pin); 
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                std::cerr << n_in << "," << n_pin << " --- " << n_out + n_in << "," << n_pout + n_pin << "\n";
                throw e;
            }
            hg_pvalues_2_[i][k] = p;
            ++num_hg_tests_;
        }
    }
}


void
StatTester::GenerateCSV()
{
    std::ofstream outI;
    outI.open(XINCSV);
    std::ofstream outH;
    outH.open(XHHCSV);
    outI << std::setprecision(std::numeric_limits<double>::digits10+2);
    outH << std::setprecision(std::numeric_limits<double>::digits10+2);
    for (int i = 0; i < static_cast<int>(X_.size1()); ++i)
    {
        bool is_hh = (X_(i,0) < 0);
        std::ofstream& out1 = is_hh ? outH : outI;
        for (int j = 0; j < static_cast<int>(X_.size2()); ++j)
        {
            if (is_hh && 2 <= j && j <= 5)
            {
                continue;
            }
            if (j != 10)
            {
                if (std::isnan(X_(i,j)))
                {
                    out1 << "NaN" << ",";
                }
                else
                {
                    out1 << X_(i,j) << ",";
                }
            }
            else
            {
                out1 << (int) std::round(X_(i,j)) << ",";
            }
            
        }
        out1 << "\n";
    }
    outH.close();
    outI.close();

    std::ofstream outYIN;
    outYIN.open(YINCSV);
    std::ofstream outYHH;
    outYHH.open(YHHCSV);
    for (int i = 0; i < static_cast<int>(y_.size1()); ++i)
    {
        bool is_hh = (X_(i,0) < 0);
        std::ofstream & out2 = is_hh?outYHH:outYIN;
        for (int j = 0; j < static_cast<int>(y_.size2()); ++j)
        {
            out2 << y_(i,j) << ",";
        }
        out2 << "\n";
    }
    outYIN.close();
    outYHH.close();
}

void StatTester::PrintHGTest()
{
    // Bonferoni corrected values:
    double bsig = SIGNIFICANT;
    double bvsig = VERY_SIGNIFICANT;
    double besig = EXTRA_SIGNIFICANT;
    int nsig = 0;
    int nvsig = 0;
    int nesig = 0;
    int nn = 0;
    for (auto x: hg_pvalues_)
    {
        // int node = x.first;
        for (auto y: x.second)
        {
            double p = y.second;
            if (p < bsig)
            {
                ++nsig;
            }
            if (p < bvsig)
            {
                ++nvsig;
            }
            if (p < besig)
            {
                ++nesig;
            }
        }
        ++nn;
    }
    std::cerr << "Total number of HG tests " << nn << "\n";
    std::cerr << "Hypergeometric test results, cumulative, profit window size " << profit_window_size_ <<  " days \n";
    std::cerr << "p < " << SIGNIFICANT << ": " << nsig << "\n";
    std::cerr << "p < " << VERY_SIGNIFICANT << ": " << nvsig << "\n";
    std::cerr << "p < " << EXTRA_SIGNIFICANT << ": " << nesig << "\n";
    nsig = 0;
    nvsig = 0;
    nesig = 0;
    for (auto x: hg_pvalues_2_)
    {
        // int node = x.first;
        for (auto y: x.second)
        {
            double p = y.second;
            if (p < bsig)
            {
                ++nsig;
            }
            if (p < bvsig)
            {
                ++nvsig;
            }
            if (p < besig)
            {
                ++nesig;
            }
        }
    }
    std::cerr << "Hypergeometric test results, cumulative, profit window size " << profit_window_size_2_ << " days \n";
    std::cerr << "p < " << SIGNIFICANT << ": " << nsig << "\n";
    std::cerr << "p < " << VERY_SIGNIFICANT << ": " << nvsig << "\n";
    std::cerr << "p < " << EXTRA_SIGNIFICANT << ": " << nesig << "\n";
}
 
void
StatTester::GenerateSmallDataMatrix()
{
       // These are the node transaction and price tables that we use to create the  delayed version. 
    const NodeTransactionTable& transacts = ior_.tr_table_;  
    const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;
    const AnnouncementTable& ans_sc_ = ior_.an_table_sc_;


    // LEGEND X: #0:id #1:company #2:distance #3:centrality #4:normalized_degree #5:degree 
    //           #6:in_window (non-scheduled) #7:in_window_2 (non-scheduled) #8:in_window (scheduled) #9: in_window_2 (scheduled)

    //         #10: date #11:volume*price  #12:business_day #13:past_inside #14:future_inside #15:How many boards #16: how many inside
    //         # 17: Insiders actually exist
    int p_length = 18;  
    X_.resize(num_transactions_,p_length, 0.0);

    // LEGEND Y: #0: Return in window 1 #1: actual trading days for window 1, #2 Return in window 2, #3 actual trading days in window 2
    //  #4: Market return in window 1 #5: actual trading days for window 1, #6 Market Return in window 2, #7 actual trading days in window 2

    y_.resize(num_transactions_,8,0.0);
    std::cerr << "generating " << num_transactions_ << " x " << p_length << " matrix \n";   
    int row = 0;

    bool done = false;
    int extras = 0;
    for (auto nodepair: transacts)
    {    
        const int node = nodepair.first;
        int i = node;
        const TransactionTable& c_trans = nodepair.second;
        for (auto isin_vector_pair: c_trans.pt_)
        {
            const std::string& isin = isin_vector_pair.first;
            const std::string& comp = ior_.isin_company_[isin];
            // std::cerr << "Company: " << comp << "\n";
            // K is for company. 
            auto comp_it = ior_.cnames_.find(comp);
            if (comp_it == ior_.cnames_.end())
            {
                std::cerr << "Warning: company " << comp << " not in dictionary\n"; 
                //throw std::runtime_error("Illegal company");
            }
            int k = ior_.cnames_[comp];
            const DatePriceVolumeVector& trans = isin_vector_pair.second;
            auto ans_it = ans.find(isin);
            auto ans_sc_it = ans_sc_.find(isin);
            bool nonscheduled_exist = (ans_it != ans.end());
            bool scheduled_exist = (ans_sc_it != ans_sc_.end());
        
            if(!(nonscheduled_exist || scheduled_exist))
            {
                std::cerr << "Warning, announcements for " << isin << " not found\n";
                continue;
            }
            std::vector<int> dummy;
             
            const auto& ansvector = nonscheduled_exist?ans_it->second:dummy;
            const auto& ansvector_sc = scheduled_exist?ans_sc_it->second:dummy;
            // std::cerr << "starting transactions:\n";
            for (unsigned int j = 0; j < trans.size(); ++j)
            {
                if (done)
                {
                    ++extras;
                    continue;
                }
                const int date = std::get<0>(trans[j]);
                int business_day = 1;
                if (ior_.trade_days_.find(date) == ior_.trade_days_.end())
                {
                    business_day = 0;
                }
                // std::cerr << date << "\n";
                const double price = std::get<1>(trans[j]);
                const double volume = std::get<2>(trans[j]); 
                // double refprice = pricetable.GetFirstChangePrice(isin,date,profit_window_size_);
                auto refpair = pricetable.GetFirstChangePrice(isin,date,profit_window_size_, ior_.trade_days_);
                double refprice = refpair.second;
                auto refpair2 = pricetable.GetFirstChangePrice(isin,date,profit_window_size_2_, ior_.trade_days_);
                double refprice2 = refpair2.second;

                
                
                
               
                double ret = refprice/price;
                double ret2 = refprice2/price;

                const double index_price = pricetable.GetCompanyDayPrice("index", date, 0);
                double m_refprice = pricetable.GetCompanyDayPrice("index", date+refpair.first,0);
                double m_refprice2 = pricetable.GetCompanyDayPrice("index", date+refpair2.first,0); 
                // std::cerr << index_price << " index vs " << m_refprice << " reference \n";
                double market_ret = m_refprice/index_price;
                double market_ret2 = m_refprice2/index_price;
        
               
                int sgn = volume > 0.0?1:-1;
                
                y_(row,0) = log(ret)*sgn;
                y_(row,1) = refpair.first;
                y_(row,2) = log(ret2)*sgn;
                y_(row,3) = refpair2.first;
                y_(row,4) = log(market_ret)*sgn;
                y_(row,5) = refpair.first;
                y_(row,6) = log(market_ret2)*sgn;
                y_(row,7) = refpair2.first;

                double dist = std::nan("household");
                double centr = std::nan("household");
                double n_deg = std::nan("household");
                double deg =  std::nan("household");
                boost::gregorian::date ref_time(boost::gregorian::from_simple_string(REFDAY));
                boost::gregorian::date truetime = ref_time + boost::gregorian::days(date);
                // Only take the beginning of the month 
                int fixed_time = 10000*truetime.year() + 100*truetime.month(); // + truetime.day();
                // std::cerr << truetime << "\n";
                auto gg = ior_.metag_->GetGraph(fixed_time);

                int insiders_exist = 1;

                // The following code is only done for insiders.    
                if (i >= 0)
                { 
                    //std::cerr << "got graph. \n";
                    try
                    {
                        dist = (double) gg->GetDistance(k,i);
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << "distance" << '\n';
                        throw e;
                    }                
                    if (dist < 0)
                    {
                        if (dist < -1.9)
                        {
                            insiders_exist = 0;
                        }
                        dist = std::nan("missing");
                    }
                    try
                    {
                        centr = gg->GetCentrality(i);
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << "Centrality" << '\n';
                        throw e;
                    }
                    try
                    {
                        n_deg = gg->NormalDegree(i);
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << "normal degree" << '\n';
                        throw e;
                    }
                    try
                    {
                        deg = gg->Degree(i);
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << "raw degree" << '\n';
                        throw e;
                    }

                    if (!std::isnan(dist) && deg == 0)
                    {
                        throw std::logic_error("distance exists but zero degree");
                    }
                }
                
                X_(row, 0) = i;
                X_(row, 1) = k;
                X_(row, 2) = dist;
                X_(row, 3) = centr;
                X_(row, 4) = n_deg;
                X_(row, 5) = deg;

                X_(row, 6) = 0;
                X_(row, 7) = 0;

                X_(row, 8) = 0;
                X_(row, 9) = 0;

                // The next announcements after date 

                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                auto ans_sc_v_it = ansvector_sc.begin();
                try
                {
                    /* code */                 
                    ans_sc_v_it = std::lower_bound(ansvector_sc.begin(), ansvector_sc.end(), date);
                }                
                catch(const std::exception& e)
                {
                    std::cerr << "It was the lowerbound.\n";
                    std::cerr << e.what() << '\n';
                }


                // Count the number of non-scheduled announcements in the two time windows
                while (ans_v_it != ansvector.end() && *ans_v_it >= date)
                {
                    if(*ans_v_it <= date + inside_window_size_)
                    {
                        ++X_(row, 6);
                    }
                    if(*ans_v_it <= date + inside_window_size_2_)
                    {
                        ++X_(row, 7);
                    }
                    else
                    {
                        break;
                    } 
                    ++ans_v_it; 
                }
                try
                {
                
                // count the number of scheduled announcements
                while (ans_sc_v_it != ansvector_sc.end() && *ans_sc_v_it >= date)
                {
                    if(*ans_sc_v_it <= date + inside_window_size_)
                    {
                        ++X_(row, 8);
                    }
                    if(*ans_sc_v_it <= date + inside_window_size_2_)
                    {
                        ++X_(row, 9);
                    }
                    else
                    {
                        break;
                    } 
                    ++ans_sc_v_it; 
                }
                     /* code */
                }
                catch(const std::exception& e)
                {
                    std::cerr << "It was the iteration";
                    std::cerr << e.what() << '\n';
                }
               
                //std::cerr << "announcements ok\n";

                X_(row, 10) = 10000*truetime.year() + 100*truetime.month() + truetime.day();
                X_(row, 11) = volume*price;
                // on business day
                X_(row, 12) = business_day;
                //std::cerr << "bdays done.\n";
                auto g_past = ior_.metag_->GetGraph(fixed_time - 10000);
                auto g_fut = ior_.metag_->GetGraph(fixed_time + 10000);
                //std::cerr << "got metagraphs.\n";
                if (!g_past)
                {
                    X_(row, 13) = std::nan("past missing");
                }
                else
                {
                    auto inset = g_past->GetInsider(k);
                    if (inset.find(i) != inset.end())
                    {
                        X_(row,13) = 1;
                    }
                }
                //std::cerr << "gpast ok.\n";

                if (!g_fut)
                {
                    X_(row, 14) = std::nan("future missing");
                }
                else
                {
                    auto inset = g_fut->GetInsider(k);
                    if (inset.find(i) != inset.end())
                    {
                        X_(row,14) = 1;
                    } 
                }
                //std::cerr << "gfut ok.\n";
                X_(row,15) = gg->GetBoardOf(i).size();
                X_(row,16) = gg->GetInsiderOf(i).size();
                if (X_(row,15) > X_(row,16))
                {
                    std::cerr << "something bad\n";
                }
                X_(row,17) = insiders_exist;
                ++row;
                if (row >= num_transactions_)
                {
                    std::cerr << "Warning, too many transactions!\n";
                    done = true;
                }
                //std::cerr << "this transaction done\n";
            }
        }
    }
    std::cerr << "Generated. Skipped " << extras << " transactions. \n";
}

void 
StatTester::DoGraphTests()
{
       // These are the node transaction and price tables that we use to create the  delayed version. 
    // const NodeTransactionTable& transacts = ior_.tr_table_;
    // const PriceTable& pricetable = ior_.pr_table_;
    // const AnnouncementTable& ans = ior_.an_table_;
            
    auto gg = ior_.metag_->GetGraph(0);
    std::cerr << "graph read.\n";
    if (!gg)
    {
        throw std::runtime_error("null graph");
    }
    for (auto cpair: ior_.cids_)
    {
        int c = cpair.first;
        std::vector<int> corrupted;
        for (auto foo: hg_pvalues_)
        {
            double hp = nan("none");
            int i = foo.first;    
            auto aux_it = foo.second.find(c);
            if (aux_it != foo.second.end())
            {
                hp = aux_it->second;
            }
            if (!std::isnan(hp) && hp < VERY_SIGNIFICANT)
            {
                corrupted.push_back(i);
            }
        }
        // std::cerr << "corrupted vector created, size " << corrupted.size() << "\n";
        if (corrupted.empty())
        {
            continue;
        }
        auto max_neigbours = gg->GetMaxComp(corrupted);
        std::cerr << "Company " << cpair.second << " corrupted " << corrupted.size() << " max component " << max_neigbours.size() << "\n";
    }
    
}

void StatTester::PrintIsinMap()
{
    for (auto isin_it: ior_.isin_company_)
    {
        auto isin = isin_it.first;
        auto cname = isin_it.second;
        auto id = ior_.cnames_[cname];
        std::cout << isin << ", " << cname << ", " << id << std::endl;
    }
}


void 
StatTester::TestGraphIntegrity()
{
    const NodeTransactionTable& transacts = ior_.tr_table_;
    // const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;
    for (auto nodepair: transacts)
    {
        const int node = nodepair.first;
        int i = node;
        const TransactionTable& c_trans = nodepair.second;
        for (auto isin_vector_pair: c_trans.pt_)
        {
            const std::string& isin = isin_vector_pair.first;
            const std::string& comp = ior_.isin_company_[isin];
            auto comp_it = ior_.cnames_.find(comp);
            if (comp_it == ior_.cnames_.end())
            {
                std::cerr << "company " << node << " not in dictionary\n"; 
                continue;
                // throw std::runtime_error("Illegal company");
            }
            int k = ior_.cnames_[comp];
            const DatePriceVolumeVector& trans = isin_vector_pair.second;
            auto ans_it = ans.find(isin);
        
            if(ans_it == ans.end())
            {
                continue;
            }
             
            // const auto& ansvector = ans_it->second; 
            for (unsigned int j = 0; j < trans.size(); ++j)
            {
                const int date = std::get<0>(trans[j]);
                /* int business_day = 1;
                if (ior_.trade_days_.find(date) == ior_.trade_days_.end())
                {
                    business_day = 0;
                }*/
                boost::gregorian::date ref_time(boost::gregorian::from_simple_string(REFDAY));
                boost::gregorian::date truetime = ref_time + boost::gregorian::days(date);
                // Only take the beginning of the year network. 
                int fixed_time = 10000*truetime.year() + 100*truetime.month();// + truetime.day();
                // std::cerr << truetime << "\n";
                auto gg = ior_.metag_->GetGraph(fixed_time);
                //std::cerr << "got graph. \n";
                double d = (double) gg->GetDistance(k,i);  
                double c = gg->GetCentrality(i);
                double p = gg->NormalDegree(i);
                double dd = gg->Degree(i);
                if (p == 0 && d > 0)
                {
                    std::cerr << "Distance :" << d << " but adjacent: " << p << "\n";
                    std::cerr << "Time stamp for graph: " << fixed_time; 
                    throw std::logic_error("Graph mismatch");
                }
                std::cerr << "Distance: " << d << " adjacent: " << dd << " Centrality: " << c << "\n";
            }
        }
    }
}


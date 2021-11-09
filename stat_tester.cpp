#include "stat_tester.hpp"
#include "defs.hpp"
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
    std::cerr << "Inside transactions: " << n_inside_p << ", Outside transactions: " << n_outside_p << "\n";
    std::cerr << "reversed days: " << n_reversed << "\n"; 
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
    std::ofstream out1;
    out1.open(XCSV);
    out1 << std::setprecision(std::numeric_limits<double>::digits10+2);
    for (int i = 0; i < X_.size1(); ++i)
    {
        for (int j = 0; j < X_.size2(); ++j)
        {
            if (j != 7)
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
    out1.close();

    std::ofstream out2;
    out2.open(YCSV);
    for (int i = 0; i < y_.size1(); ++i)
    {
        for (int j = 0; j < y_.size2(); ++j)
        {
            out2 << y_(i,j) << ",";
        }
        out2 << "\n";
    }
    out2.close();
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
        int node = x.first;
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
        int node = x.first;
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
   
    // Investor ID, Company ID, Distance, Centrality, Normalized Degree, IsInside, hg p-value, date
    int p_length = 10;  
    X_.resize(num_transactions_,p_length, 0.0);
    y_.resize(num_transactions_,4,0.0);
    std::cerr << "generating " << num_transactions_ << " x " << p_length << " matrix \n";   
    int row = 0;
    for (auto nodepair: transacts)
    {
        const int node = nodepair.first;
        /*
        auto node_it = ior_.nodeinvdict_.find(node); 
        // Find the index of the node:
        if (node_it == ior_.nodeinvdict_.end())
        {
            std::cerr << "node " << node << " not in dictionary\n"; 
            throw std::runtime_error("Illegal node");
        }
        */
        int i = node;
        const TransactionTable& c_trans = nodepair.second;
        for (auto isin_vector_pair: c_trans.pt_)
        {
            const std::string& isin = isin_vector_pair.first;
            const std::string& comp = ior_.isin_company_[isin];
            // K is for company. 
            auto comp_it = ior_.cnames_.find(comp);
            if (comp_it == ior_.cnames_.end())
            {
                std::cerr << "company " << node << " not in dictionary\n"; 
                throw std::runtime_error("Illegal company");
            }
            int k = ior_.cnames_[comp];
            const DatePriceVolumeVector& trans = isin_vector_pair.second;
            auto ans_it = ans.find(isin);
        
            if(ans_it == ans.end())
            {
                continue;
            }
             
            const auto& ansvector = ans_it->second; 
            for (unsigned int j = 0; j < trans.size(); ++j)
            {
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
                /*
                if (refprice < 0.000001)
                {
                    refprice = price;
                }
                if (refprice2 < 0.000001)
                {
                    refprice2 = price;
                }
                */
                // The next announcement. 
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                double ret = refprice/price;
                double ret2 = refprice2/price;
               
                int sgn = volume > 0.0?1:-1;
                //std::cerr << "adding: " << ret << " to row " << row << "\n";
                // We need a fixed stat. 
            
                
                boost::gregorian::date ref_time(boost::gregorian::from_simple_string(REFDAY));
                boost::gregorian::date truetime = ref_time + boost::gregorian::days(date);
                // Only take the beginning of the year network. 
                int fixed_time = 10000*truetime.year(); // + 100*truetime.month() + truetime.day();
                // std::cerr << truetime << "\n";
                auto gg = ior_.metag_->GetGraph(fixed_time);
                //std::cerr << "got graph. \n";
                double d = (double) gg->GetDistance(k,i);  
                if (d < 0)
                {
                    d = std::nan("missing");
                }
                double c = gg->GetCentrality(i);
                double p = gg->PageRank(i);
                // std::cerr << "Centrality " << c << "\n";
                y_(row,0) = log(ret)*sgn;
                y_(row,1) = refpair.first;
                y_(row,2) = log(ret2)*sgn;
                y_(row,3) = refpair2.first;
                
                //std::cerr << "row " << row << " node " << i << " company number " << k << " offset " << offset << " distance " << d << "\n";
                X_(row, 0) = i;
                X_(row, 1) = k;
                X_(row, 2) = d;
                X_(row, 3) = c;
                X_(row, 4) = p;

                if (ans_v_it != ansvector.end() && *ans_v_it >= date && *ans_v_it <= date + inside_window_size_)
                {
                    X_(row, 5) = 1;
                }
                double hp = std::nan("miss");
                auto hg_it = hg_pvalues_.find(i);
                if (hg_it != hg_pvalues_.end())
                {
                    auto aux_it = hg_it->second.find(k);
                    if (aux_it != hg_it->second.end())
                    {
                        hp = aux_it->second;
                    }
                }
                X_(row, 6) = hp;
                X_(row, 7) = 10000*truetime.year() + 100*truetime.month() + truetime.day();
                // std::cerr << 10000*truetime.year() + 100*truetime.month() + truetime.day() << "\n";
                // buy/sell
                X_(row, 8) = volume*price;
                // on business day
                X_(row, 9) = business_day;

                //else
                //{
                //}
                //std::cerr << "\n";
                ++row;
                if (i == 11 && k == 50)
                {
                   std::cerr << "row " << row << " node " << i << " company number " << k << " normdegree " << p << " distance " << d << "\n";
                }
            }
        }
    }
}

void 
StatTester::DoGraphTests()
{
       // These are the node transaction and price tables that we use to create the  delayed version. 
    const NodeTransactionTable& transacts = ior_.tr_table_;
    const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;
            
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
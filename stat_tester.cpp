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
    ior_.ReadTables();
    ior_.pr_table_.Sort(); 
    ior_.SetCompanyDictionaryFile(CDFILE);
    ior_.SetNodeDictionaryFile(NDFILE);
    ior_.ReadCompanyDictionary();
    ior_.ReadNodeDictionary();
    window_size_ = 5; 
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
    window_size_ = size;
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
            if (*ans_v_it > day && *ans_v_it < day + window_size_)
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
            IsinCountMap novel3;
            IsinCountMap novel4;
            n_profit_inside_[node] = novel3;
            n_profit_outside_[node] = novel4;
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
                (profit_outside_[node])[comp] = novel1;
                (profit_inside_[node])[comp] = novel2;
                n_profit_inside_[node][comp] = 0;
                n_profit_outside_[node][comp] = 0;
            }    
            const auto& ansvector = ans_it->second; 
            for (unsigned int i = 0; i < trans.size(); ++i)
            {
                const int date = std::get<0>(trans[i]);
                const double price = std::get<1>(trans[i]);
                const double volume = std::get<2>(trans[i]); 
                double refprice = pricetable.GetCompanyDayPrice(isin,date,window_size_);
                double refprice2 = pricetable.GetFirstChangePrice(isin,date,window_size_);
                // The next announcement. 

                double profit = volume*(refprice-price);
                
                if (refprice < 0.000001)
                {
                    refprice = price;
                }
                if (refprice2 < 0.000001)
                {
                    refprice2 = price;
                }
                // The next announcement. 
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                double ret = refprice/price;
                if( fabs(refprice2/price - 1) > fabs(ret -1))
                {
                    ret = refprice2/price;
                }

                int sgn = volume > 0.0?1:-1;
                double lret = sgn*log(ret);

                ++num_transactions_;
                // std::cerr << "Return: " << ret;
                // Is is inside?
                if (ans_v_it != ansvector.end() && *ans_v_it > date && *ans_v_it <= date + window_size_)
                {
                    //std::cerr << "inside;";
                    // Its inside
                    profit_inside_[node][comp].push_back(ret);
                    if (lret > PROFIT_THRESHOLD)
                    {
                        //std::cerr << "profit;";
                        ++n_profit_inside_[node][comp];
                    }
                }
                else
                {
                    //std::cerr << "outside;";
                    if (lret > PROFIT_THRESHOLD)
                    {
                        //std::cerr << "profit;";
                        ++n_profit_outside_[node][comp];
                    }
                    profit_outside_[node][comp].push_back(ret);
                }
                //std::cerr << "\n";
            }
        }
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
            ++num_hg_tests_;
        }
    }
}

void 
StatTester::GenerateDataMatrix()
{
       // These are the node transaction and price tables that we use to create the  delayed version. 
    const NodeTransactionTable& transacts = ior_.tr_table_;
    const PriceTable& pricetable = ior_.pr_table_;
    const AnnouncementTable& ans = ior_.an_table_;
   
    int offset = num_companies_ + num_traders_ + 2;
    int p_length = 2*offset;  
    X_.resize(num_transactions_,p_length, 0);
    y_.resize(num_transactions_,1);
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
            {
                if(ans_it == ans.end())
                {
                    continue;
                }
            } 
            const auto& ansvector = ans_it->second; 
            for (unsigned int j = 0; j < trans.size(); ++j)
            {
                const int date = std::get<0>(trans[j]);
                const double price = std::get<1>(trans[j]);
                const double volume = std::get<2>(trans[j]); 
                double refprice = pricetable.GetCompanyDayPrice(isin,date,window_size_);
                double refprice2 = pricetable.GetFirstChangePrice(isin,date,window_size_);
                if (refprice < 0.000001)
                {
                    refprice = price;
                }
                if (refprice2 < 0.000001)
                {
                    refprice2 = price;
                }
                // The next announcement. 
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                double ret = refprice/price;
                if( fabs(refprice2/price - 1) > fabs(ret -1))
                {
                    ret = refprice2/price;
                }
                // The next announcement. 
                int sgn = volume > 0.0?1:-1;
                //std::cerr << "adding: " << ret << " to row " << row << "\n";
                y_(row,0) = log(ret)*sgn;

                // We need a fixed stat.
                
                boost::gregorian::date ref_time(boost::gregorian::from_simple_string(REFDAY));
                boost::gregorian::date truetime = ref_time + boost::gregorian::days(date);
                // Only take the beginning of the year network. 
                int fixed_time = 10000*truetime.year(); // + 100*truetime.month() + truetime.day();
                
                auto gg = ior_.metag_->GetGraph(fixed_time);
                //std::cerr << "got graph. \n";
                int d = gg->GetDistance(k,i);   
                double c = gg->GetCentrality(i);
                // std::cerr << "Centrality " << c << "\n";

                //std::cerr << "row " << row << " node " << i << " company number " << k << " offset " << offset << " distance " << d << "\n";
                if (ans_v_it != ansvector.end() && *ans_v_it >= date && *ans_v_it <= date + window_size_)
                {
                    // Its inside
                    //std::cerr << "inside.\n";
                    X_(row, offset + i) = 1;
                    //std::cerr << "first ok.";
                    //std::cerr << offset + num_traders_ + k << " vs " << p_length << "\n"; 
                    //std::cerr << k << " vs " << num_companies_ << " and "  << i << " vs " << num_traders_ << "\n";
                    X_(row, offset + num_traders_ + k) = 1;
                    //std::cerr << "seconf ok.";
                    X_(row, offset + num_traders_ + num_companies_) = d;
                    //std::cerr << "third ok\n";
                    X_(row, offset + num_traders_ + num_companies_ + 1) = c;
                }
                else
                {
                    //std::cerr << "outside;\n";
                    X_(row, i) = 1;
                    //std::cerr << "first ok.";
                    X_(row, num_traders_ + k) = 1;
                    //std::cerr << "seconf ok.";
                    X_(row, num_traders_ + num_companies_) = d;
                    //std::cerr << "third ok\n";
                    X_(row, num_traders_ + num_companies_+ 1) = c;
                }
                //std::cerr << "\n";
                ++row;
            }
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
                out1 << X_(i,j) << ",";
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
    double bsig = SIGNIFICANT/num_hg_tests_;
    double bvsig = VERY_SIGNIFICANT/num_hg_tests_;
    double besig = EXTRA_SIGNIFICANT/num_hg_tests_;
    int nsig = 0;
    int nvsig = 0;
    int nesig = 0;
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
    }
    std::cerr << "Hypergeometric test results, cumulative, Bonferoni corrected values: \n";
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
    int p_length = 9;  
    X_.resize(num_transactions_,p_length, 0.0);
    y_.resize(num_transactions_,1);
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
                const double price = std::get<1>(trans[j]);
                const double volume = std::get<2>(trans[j]); 
                double refprice = pricetable.GetFirstChangePrice(isin,date,window_size_);
                double refprice2 = pricetable.GetCompanyDayPrice(isin,date,window_size_);
                if (refprice < 0.000001)
                {
                    refprice = price;
                }
                if (refprice2 < 0.000001)
                {
                    refprice2 = price;
                }
                // The next announcement. 
                auto ans_v_it = std::lower_bound(ansvector.begin(), ansvector.end(), date);
                double ret = refprice/price;
                if( fabs(refprice2/price - 1) > fabs(ret -1))
                {
                    ret = refprice2/price;
                }
                int sgn = volume > 0.0?1:-1;
                //std::cerr << "adding: " << ret << " to row " << row << "\n";
                // We need a fixed stat.
                
                boost::gregorian::date ref_time(boost::gregorian::from_simple_string(REFDAY));
                boost::gregorian::date truetime = ref_time + boost::gregorian::days(date);
                // Only take the beginning of the year network. 
                int fixed_time = 10000*truetime.year(); // + 100*truetime.month() + truetime.day();
                
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
                //std::cerr << "row " << row << " node " << i << " company number " << k << " offset " << offset << " distance " << d << "\n";
                X_(row, 0) = i;
                X_(row, 1) = k;
                X_(row, 2) = d;
                X_(row, 3) = c;
                X_(row, 4) = p;

                if (ans_v_it != ansvector.end() && *ans_v_it >= date && *ans_v_it <= date + window_size_)
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
                X_(row, 8) = sgn;

                //else
                //{
                //}
                //std::cerr << "\n";
                ++row;
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
            
    auto gg = ior_.metag_->GetGraph(20050101);
    for (auto cpair: ior_.cids_)
    {
        int c = cpair.first;
        std::vector<int> corrupted;
        for (auto ncpair: hg_pvalues_)
        {
            auto p_it = ncpair.second.find(c);
            if(p_it == ncpair.second.end())
            {
                continue;
            }
            int in = ncpair.first;
            if (hg_pvalues_[in][c] < VERY_SIGNIFICANT)
            {
                corrupted.push_back(in);
            }
        }
    }
    
}
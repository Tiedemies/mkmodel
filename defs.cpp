#include "defs.hpp"
#include <algorithm>
#include <string>

PriceTable::PriceTable()
{
    // Default
    sorted_ = false;
    return;
}
PriceTable::~PriceTable()
{
    // void
}

// Add a company day price
void 
PriceTable::AddPCompanyDayPrice(const std::string& cname, int day, double price)
{
    pt_[cname].push_back(std::make_pair(day,price));
}

double 
PriceTable::GetCompanyDayPrice(const std::string& cname, int day, int offset) const
{
    if (!sorted_)
    {
        throw std::runtime_error("unsorted pricetables"); 
    }
    auto pt_it = pt_.find(cname);
    if(pt_it == pt_.end())
    {
        return nan("notfound");
    }
    auto lower = std::lower_bound(pt_it->second.begin(), pt_it->second.end(), std::make_pair(day, 0.0),std::less_equal<std::pair<int,double>>());
    while (lower != pt_it->second.end() && (*lower).first < day+offset)
    {
        ++lower;
    }
    return lower->second; 
}

void 
PriceTable::Sort()
{
    if(sorted_)
    {
        return;
    }
    for (auto Z: pt_)
    {
        std::sort(Z.second.begin(),Z.second.end(),std::less<std::pair<int,double>>()); 
    }
    sorted_ = true; 
}

int 
PriceTable::size()
{
    return (int) pt_.size(); 
}


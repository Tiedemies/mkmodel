#include "defs.hpp"
#include <algorithm>
#include <string>
#include <iostream>

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
    auto lower = std::lower_bound(pt_it->second.begin(), pt_it->second.end(), std::make_pair(day + offset, 0.0),std::less_equal<std::pair<int,double>>());
    while (lower != pt_it->second.end() && std::isnan(lower->second))
    {
        ++lower;
    }
    return lower->second; 
}

double 
PriceTable::GetFirstChangePrice(const std::string& cname, int day, int offset) const
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
    auto lower = std::upper_bound(pt_it->second.begin(), pt_it->second.end(), std::make_pair(day, 0.0),std::less_equal<std::pair<int,double>>());
    double ref = lower->second; 
    bool changed = false; 
    int num = lower->first;
    while ((!changed || num-day < offset))
    {
        ++lower;
        num = lower->first; 
        if (lower != pt_it->second.end())
        {
            changed = (fabs(lower->second - ref) > 0.00001) && !std::isnan(lower->second);
        }
        else return(ref);
    }
    // std::cerr << "Date: " << day << " refday " << lower->first << " price " << lower->second << " original " << ref << "\n";
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


// Transaction table

TransactionTable::TransactionTable()
{
    // Default
    sorted_ = false;
    return;
}
TransactionTable::~TransactionTable()
{
    // void
}

// Add a company day price
void 
TransactionTable::AddPCompanyDayPriceTransaction(const std::string& cname, int day, double price, double volume)
{
    pt_[cname].push_back(std::make_tuple(day,price, volume));
}

void 
TransactionTable::Sort()
{
    if(sorted_)
    {
        return;
    }
    for (auto Z: pt_)
    {
        std::sort(Z.second.begin(),Z.second.end(),std::less<std::tuple<int,double,double>>()); 
    }
    sorted_ = true; 
}

int 
TransactionTable::size()
{
    return (int) pt_.size(); 
}


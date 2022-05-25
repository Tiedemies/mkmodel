#include "defs.hpp"
#include <algorithm>
#include <string>
#include <iostream>

PriceTable::PriceTable()
{
    // Default
}
PriceTable::~PriceTable()
{
    // void
}

// Add a company day price
void 
PriceTable::AddPCompanyDayPrice(const std::string& cname, int day, double price)
{
    pt_[cname][day] = price;
}

double 
PriceTable::GetCompanyDayPrice(const std::string& cname, int day, int offset) const
{
    auto pt_it = pt_.find(cname);
    if(pt_it == pt_.end())
    {
        return std::nan("notfound");
    }
    auto d_pair = pt_it->second.find(day);
    if (d_pair == pt_it->second.end())
    {
        return std::nan("notfound");
    }
    else return d_pair->second; 
}


std::pair<int,double> 
PriceTable::GetFirstChangePrice(const std::string& cname, int day, int offset, const std::set<int>& dates) const
{
    auto pt_it = pt_.find(cname);
    if(pt_it == pt_.end())
    {
        return std::make_pair(-1, std::nan("notfound"));
    }
    auto lower = pt_it->second.lower_bound(day);
    if (lower == pt_it->second.end())
    {
        std::make_pair(-1, std::nan("notfound"));
    }
    auto tradingday = dates.lower_bound(day);
    for (int i=0; i < offset; ++i)
    {   
        if (tradingday != dates.end())
        {
            ++tradingday;
        }
    }
    int offset2 = offset; 
    if (tradingday != dates.end())
    {
        offset2 = (*tradingday - day); 
    }
    else
    {
        return std::make_pair(-1, std::nan("notfound"));
    }
    
    // int num = lower->first;
    while ((lower->first-day < offset2) || std::isnan(lower->second))
    {
        ++lower;
        if (lower == pt_it->second.end())
        {
            return std::make_pair(-1, std::nan("missing"));
        }
    }
    // std::cerr << "Date: " << day << " refday " << lower->first << " price " << lower->second << " original " << ref << "\n";
    return std::make_pair(lower->first-day ,lower->second); 
}


/*
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
*/

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


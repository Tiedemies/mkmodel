#include "defs.hpp"
#include <algorithm>
#include <string>
#include <iostream>
#include <cmath>

namespace util
{
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


NLargest::~NLargest()
{
  // void
}

NLargest::NLargest(int n): n_(n)
{
  // void
}

void 
NLargest::Insert(int i, double v)
{
  if (values_.size() < n_)
  {
    values_.push_back(v);
    indices_.push_back(i);
    if (values_.size() == n_)
    {
      BuildHeap();
    }
    return;
  }
  // IsHeap
  if (values_[0] > v)
  {
    return;
  }
  values_[0] = v;
  indices_[0] = i;
  Heapify();
}

const std::vector<int>& 
NLargest::GetIndeces() const
{
  return indices_;
}

const std::vector<double>& 
NLargest::GetValues() const
{
  return values_;
} 

void
NLargest::Heapify(int cur)
{
  int n = values_.size();
  while (true)
  {
    int smallest = cur; 
    int left = 2*cur + 1;
    int right = 2*cur + 2;  
    if (left >= n) 
    {
      return;
    }
    if (values_[left] < values_[cur])
    {
      smallest = left;
    }
    if (right < n && values_[right] < values_[smallest])
    {
      smallest = right; 
    }    
    if (smallest == cur)
    {
      return;
    }
    int temp = indices_[cur];
    indices_[cur] = indices_[smallest];
    indices_[smallest] = temp;
    double temp2 = values_[cur];
    values_[cur] = values_[smallest];
    values_[smallest] = temp2;
    cur = smallest; 
  }
}

void
NLargest::BuildHeap()
{
  for (int i = (int) values_.size()/2; i >= 0; --i)
  {
    Heapify(i);
  }
}

std::pair<int,double>
NLargest::PopSmallest()
{
  double temp_value = values_.back();
  int temp_index = indices_.back();
  double ret_value = values_[0];
  int ret_index = indices_[0];
  values_[0] = temp_value;
  indices_[0] = temp_index;
  values_.pop_back();
  indices_.pop_back();
  Heapify();
}

}
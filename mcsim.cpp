#include"mcsim.hpp"
#include"mkmodel.hpp"
#include"defs.hpp"
#include<vector>
#include<list>
#include<iostream>

Simulator::Simulator(MetaGraph* mg, MkModel mk)
{
  mgraph_ = mg;
  mkproto_ = mk;
  n_markovs_ = 1; 
  int n = mg->GetNumber();
  markovs_.resize(n+1,&mkproto_);
  // std::cerr << "Debug: Instantiated simulator with " << n << " markov copies\n";
}
  
Simulator::Simulator(MetaGraph* mg, std::vector<MkModel*> mk)
{
  mgraph_ = mg;
  n_markovs_ = mk.size();
  markovs_ = mk;
  if (n_markovs_ != mg->GetNumber())
  {
    throw std::logic_error("different number of nodes and markovs"); 
  }
} 
 
Simulator::~Simulator() 
{
  //void
}
void Simulator::SetTransProps(std::vector<double> props)
{
  if (props.size() != markovs_.size())
  {
    throw std::logic_error("different number of trans props and markovs");
  }
  int i = 0;
  for (MkModel* mk: markovs_)
  {
    mk->setTransmitProb(props[i]);
    ++i;
  }

}

int Simulator::YearFirst(int year)
{
  return yfirstmap_[year]; 
}

void Simulator::SetTransProp(double prop)
{
  if (n_markovs_ != 1)
  {
    throw std::logic_error("not singleton model for setTrans");
  }
  mkproto_.setTransmitProb(prop);
}

void Simulator::SetTradeProps(std::vector<double> props)
{
  if (props.size() != markovs_.size())
  {
    throw std::logic_error("different number of trade props and markovs");
  }
  int i = 0;
  for (MkModel* mk: markovs_)
  {
    mk->setTradeProb(props[i]);
    ++i;
  }

}

std::vector<int> 
Simulator::GetInside(int date, int k)
{
  const MonoGraph* M = mgraph_->GetGraph(date);
  //std::cerr << "got graph." << M << "\n"; 
  return M->GetInsider(k);  
}

void Simulator::Simplify()
{
  auto dates = mgraph_->GetDates();
  // std::cerr << "simplified " << dates.size() << " dates \n"; 
  int year = 0;
  int refdate = 0;
  for (int date: dates)
  {
    if (year != date/10000)
    {
      year = date/10000;
      refdate = date;
      yfirstmap_[year] = date; 
      //std::cerr << "date : " << date << " year: " << year << "\n"; 
    }
    yearmap_[date] = refdate;
  }
  // std::cerr << "years: " << yfirstmap_.size() << ", dates: " << yearmap_.size() ; 
}

void Simulator::SetTradeProp(double prop)
{
  if (n_markovs_ != 1)
  {
    throw std::logic_error("not singleton model for trade");
  }
  mkproto_.setTradeProb(prop);
}
  
std::pair<std::vector<double>, std::vector<double>> 
Simulator::Simulate(const std::vector<int>& inside, int date, int n) const
{
  const MonoGraph* graph = mgraph_->GetGraph(date);
  int nodes = graph->GetNumber();
  std::vector<double> result(nodes+1,0); 
  std::vector<double> margin(nodes+1,0); 
  #pragma omp parallel for 
  for (int i = 0; i < n; ++i)
  {
    std::vector<bool> is_infected(nodes+1,false);
    std::list<int> infected;
    for (auto j: inside)
    {
      infected.push_back(j); 
      is_infected.at(j) = true;
    }
    while(!infected.empty())
    {
      int j = infected.back(); 
      infected.pop_front(); 
      for (int k: graph->GetNeighbours(j))
      {
        if (!is_infected.at(k) && markovs_.at(j)->isTransmitting(0))
        {
          is_infected.at(k) = true;
          result.at(k) += 1.0/n;
          infected.push_back(k);
        }
      } 
    }
  }
  return std::make_pair(result,margin); 
}

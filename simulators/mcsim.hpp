#ifndef MC_SIM
#define MC_SIM

#include "../graphmodel/graphmodel.hpp"
#include "../mkmodel/mkmodel.hpp"
#include "../utils/defs.hpp"


#include<vector> 

class Simulator
{
public:
  Simulator(MetaGraph* mg, MkModel mk);
  Simulator(MetaGraph* mg, std::vector<MkModel*> mk);
  ~Simulator();
  void SetTransProps(std::vector<double> props);
  void SetTransProp(double prop);
  void SetTradeProps(std::vector<double> props);
  void SetTradeProp(double props);
  void Simplify();
  int YearFirst(int year);
  std::vector<int> GetInside(int date, int k);
  std::pair<std::vector<double>, std::vector<double>> Simulate(const std::vector<int>& inside, int date, int n) const; 
private:
 
  MetaGraph* mgraph_;
  std::vector<MkModel*> markovs_;
  MkModel mkproto_;
  int n_markovs_; 
  std::unordered_map<int,int> yearmap_;
  std::unordered_map<int,int> yfirstmap_;
};


#endif

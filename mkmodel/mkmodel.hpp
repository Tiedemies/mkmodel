#ifndef MK_MODEL
#define MK_MODEL

#include "../utils/defs.hpp"
#include "../utils/h_random.hpp"

#include<vector>

class MkModel
{
public:
  MkModel();
  MkModel(Random& rnd);
  MkModel(const MkModel &rhs);
  ~MkModel();
  bool isTrading(int i);
  bool isTransmitting(int i);
  void setTradeProb(double p);
  void setTransmitProb(double p); 
  double getTradeProp();
  double getTransmitProp();
  int Move(int i); 
private:
  MkLists alist_;
  std::vector<bool> tradestates_;
  std::vector<bool> transmitstates_;
  double transmitprop_;
  double tradeprop_; 
  Random rnd_;

};

#endif

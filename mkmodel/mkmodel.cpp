#include "mkmodel.hpp"

#include "../utils/h_random.hpp"
#include "../utils/defs.hpp"
#include <chrono>
#include <vector>

MkModel::MkModel()
{
  //void
}

MkModel::MkModel(Random& rnd)
{
  tradestates_ = {true};
  transmitstates_ = {true};
  transmitprop_ = 0.5;
  tradeprop_ = 0.5;
  rnd_ = rnd; 
}

// Copy constructor. 
MkModel::MkModel(const MkModel& rhs)
{
  tradestates_ = rhs.tradestates_;
  transmitstates_ = rhs.transmitstates_;
  transmitprop_ = rhs.transmitprop_;
  tradeprop_ = rhs.tradeprop_;
  rnd_ = rhs.rnd_; 
}

MkModel::~MkModel()
{
  //void
}


bool MkModel::isTrading(int i)
{
  if (tradestates_[i%tradestates_.size()])
  {
    return (rnd_.get() - tradeprop_ > 0);
  }
  return false;
}


bool MkModel::isTransmitting(int i)
{
  if (transmitstates_[i%tradestates_.size()])
  {
    return (rnd_.get() - transmitprop_ > 0);
  }
  return false;
}

double MkModel::getTransmitProp()
{
  return transmitprop_; 
}

void MkModel::setTradeProb(double p)
{
  tradeprop_ = p;
}

void MkModel::setTransmitProb(double p)
{
  transmitprop_ = p;
}

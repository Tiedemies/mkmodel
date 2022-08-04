#include<random>
#include "h_random.hpp"

Random::Random()
{
  dis_ = std::uniform_real_distribution<double>(0,1);
}

Random::~Random()
{
  //void
}

double Random::get()
{
  return dis_(gen_);
}

NormalRandom::NormalRandom(double mu, double sigma)
{
  dis_ = std::normal_distribution<double>(mu,sigma);
}

NormalRandom::~NormalRandom()
{
  //void
}

double NormalRandom::get()
{
  return dis_(gen_);
}

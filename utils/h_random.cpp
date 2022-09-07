#include<random>
#include "h_random.hpp"
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/beta.hpp>

namespace util
{
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

BetaRandom::BetaRandom(double a, double b)
{
  BOOST_ASSERT(a > 0);
  BOOST_ASSERT(b > 0);
  dist_ = BDist(a,b);
  rand_gen_ = RNDGenerator(2018);
}

BetaRandom::~BetaRandom()
{
  //void
}

double BetaRandom::get()
{
    return dist_(rand_gen_);
}
}
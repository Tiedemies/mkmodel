#include<random>
#include "h_random.hpp"
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/beta.hpp>
#include <chrono>

namespace util
{
Random::Random()
{
  dis_ = std::uniform_real_distribution<double>(0,1);
  gen_.seed(std::chrono::system_clock::now().time_since_epoch().count());
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
  dis1_ = std::uniform_real_distribution<double>(0,1);
}

BetaRandom::~BetaRandom()
{
  //void
}

double BetaRandom::get()
{
    return boost::math::quantile(dist_, dis1_(rand_gen_));
}
}
#ifndef H_RANDOM
#define H_RANDOM

#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/beta.hpp>

namespace util{
// Global random generator
class Random
{
public:
  Random();
  ~Random();
  double get();
private:
  std::default_random_engine gen_;
  std::uniform_real_distribution<double> dis_; 
};

class NormalRandom
{
  public:
    NormalRandom(double mu=0, double sigma=1);
    ~NormalRandom();
    double get();
  private:
    std::default_random_engine gen_;
    std::normal_distribution<double> dis_;
};

class BetaRandom
{
  typedef boost::random::mt19937 RNDGenerator;
  typedef boost::math::beta_distribution<> BDist;
  public:
    BetaRandom() = delete;
    BetaRandom(double a, double b);
    ~BetaRandom();
    double get();
  private:
    BDist dist_;
    RNDGenerator rand_gen_;
    std::uniform_real_distribution<double> dis1_;
};
}
#endif


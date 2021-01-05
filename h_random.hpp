#ifndef H_RANDOM
#define H_RANDOM

#include<random>

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

#endif


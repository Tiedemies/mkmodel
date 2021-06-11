
#ifndef GRAPH
#define GRAPH
#include "defs.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>

class MonoGraph
{
  typedef std::unordered_map<int,int> DistMap;
  typedef std::unordered_map<int,DistMap> CompanyDistanceMap;
public:
  MonoGraph(std::ifstream& in);
  ~MonoGraph();
  int GetNumber() const;
  // Get neighbours
  const Alist& GetNeighbours(int i) const;
  void ReadInsiders(std::ifstream& in);
  std::vector<int> GetInsider(int k) const; 
  // Calculate the distance of insiders to company. 
  const std::unordered_map<int,int>& GetDistances(int comp);
  int GetDistance(int comp, int node); 
private:
  // Adjacency list
  AdjLists adj_;
  int number_;
  std::unordered_map<int,std::vector<int>> insiderdict_;
  CompanyDistanceMap dists_;  
};

class MetaGraph
{

public:
  MetaGraph(std::string dir);
  MetaGraph(std::string dir, const AnnouncementDates& an);
  ~MetaGraph();
  MonoGraph* GetGraph(int date) const;
  int GetNumber() const;
  std::vector<int> GetDates() const;
private:
  std::map<int,MonoGraph*> graphs_;
  int number_;
};

#endif

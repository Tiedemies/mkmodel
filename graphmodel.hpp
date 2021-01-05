
#ifndef GRAPH
#define GRAPH
#include "defs.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<vector>

class MonoGraph
{
public:
  MonoGraph(std::ifstream& in);
  ~MonoGraph();
  int GetNumber() const;
  // Get neighbours
  const Alist& GetNeighbours(int i) const;
  void ReadInsiders(std::ifstream& in);
  std::vector<int> GetInsider(int k) const; 
private:
  // Adjacency list
  AdjLists adj_;
  int number_;
  std::unordered_map<int,std::vector<int>> insiderdict_;
 
};

class MetaGraph
{
public:
  MetaGraph(std::string dir);
  MetaGraph(std::string dir, const AnnouncementDates& an);
  ~MetaGraph();
  const MonoGraph* GetGraph(int date) const;
  int GetNumber() const;
  std::vector<int> GetDates() const;
private:
  std::unordered_map<int,MonoGraph*> graphs_;
  int number_; 
};

#endif

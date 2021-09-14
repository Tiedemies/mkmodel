
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
  
  /* 
   * Distance between node and the inner circle of companyt comp. 
   * Calculated with breadt-first search. 
   */ 
  int GetDistance(int comp, int node); 

  /*
   * Betweenness centrality, calculated with the Brandes- algorithm (2001)
   * Then it is normalized with in- and out degree. 
   */
  double GetCentrality(int node);

  /*
   * PageRank- centrality
   */
  double PageRank(int node, int comp=0);

private:
  // Adjacency list
  AdjLists adj_;

  /*normalize centrality*/
  void NormalizeCentrality(); 

  int number_;
  int num_edges_; 
  std::unordered_map<int,std::vector<int>> insiderdict_;
  CompanyDistanceMap dists_;  
  std::vector<double> centrality_;
  std::vector<double> P_;
  int max_distance_;  
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
 
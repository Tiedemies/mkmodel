
#ifndef GRAPH
#define GRAPH
#include "../utils/defs.hpp"
#include<iostream>
#include<fstream>
#include<unordered_map>
#include<map>
#include<vector>
#include<set>

namespace graphmodel
{
  using namespace util;
class MonoGraph
{
  typedef std::unordered_map<int,int> DistMap;
  typedef std::unordered_map<int,DistMap> CompanyDistanceMap;
public:
  MonoGraph(std::ifstream& in);
  ~MonoGraph();
  int GetNumber() const;
  // Get neighbours
  const Alist2& GetNeighbours(int i) const;
  void ReadInsiders(std::ifstream& in);
  void ReadBoardMembers(std::ifstream& in);
  std::set<int> GetInsider(int k) const; 
  std::set<int> GetInsiderOf(int i);
  std::set<int> GetBoard(int k); 
  std::set<int> GetBoardOf(int i); 
  std::set<int> GetMaxComp(std::vector<int> input) const;
  // Calculate the distance of insiders to company. 
  const std::unordered_map<int,int>& GetDistances(int comp);
  const Alist2 def_; 
  
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
   * Normalized Degree
   */
  double NormalDegree(int node);

    /*
   * Normalized Degree
   */
  double Degree(int node);

  std::unordered_map<int,Alist2> adj_;

private:
  // Adjacency list

  /*normalize centrality*/
  void NormalizeCentrality(); 

  int number_;
  int num_edges_; 
  // Company --> Insiders
  std::unordered_map<int,std::set<int>> insiderdict_;
  // Insider --> Companies
  std::unordered_map<int,std::set<int>> insider_of_; 
  std::unordered_map<int,std::set<int>> boarddict_;
  std::unordered_map<int,std::set<int>> board_of_;
  CompanyDistanceMap dists_;  
  std::vector<double> centrality_;
  std::vector<double> P_;
  bool p_calculated_; 
  int max_distance_;  
};

class MetaGraph
{

public:
  MetaGraph(std::string dir, int date = -1);
  MetaGraph(std::string dir, const AnnouncementDates& an);
  ~MetaGraph();
  MonoGraph* GetGraph(int date) const;
  int GetNumber() const;
  std::vector<int> GetDates() const;
  bool CheckValidity(int node); 
private:
  std::map<int,MonoGraph*> graphs_;
  int number_;
  std::set<int> validated_nodes_; 
};
}
#endif
 
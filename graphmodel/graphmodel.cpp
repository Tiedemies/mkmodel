/*
* Graphmodel binary files
*/

#include "graphmodel.hpp"
#include "../utils/defs.hpp"
#include<dirent.h>
#include<new>
#include<algorithm>
#include<list>
#include<iostream>
#include<sstream>

#define FHANDLE "/egde_list.txt"
#define CHANDLE "/company_dict_insiders.txt"
#define BHANDLE "/company_dict_board_members.txt"

namespace graphmodel
{
  using namespace util;
MonoGraph::MonoGraph(std::ifstream& in):def_({}),p_calculated_(false)
{
  // Read the graph in
  std::string c;
  int i = 0;
  number_ = 0;
  num_edges_ = 0;
  if (in.is_open())
  {
    in >> c;
    while(!in.eof())
    {
      ++i;
      std::stringstream foo(c);
      int s_node;
      char sep;
      int t_node;
      foo >> s_node >> sep >> t_node;
      auto a_it = adj_.find(s_node);
      if (a_it == adj_.end())
      {
	      adj_[s_node] = Alist2({t_node});
        ++num_edges_;
      }
      else
      {
	      adj_[s_node].insert(t_node);
        ++num_edges_; 
      }
      if (adj_.find(t_node) == adj_.end())
      {
        adj_[t_node] = Alist2({s_node});
      }
      else
      {
        adj_[t_node].insert(s_node);
      }
      
      if (s_node > number_) 
      {
        number_ = s_node;
      }
      if (t_node > number_) 
      {
        number_ = t_node; 
      }
      in >> c;
    }
    //std::cerr << adj_.size() << " ";
  }
  else
  {
    std::cerr << "Failed to open file.\n";
  }
  max_distance_ = 0;
}
MonoGraph::~MonoGraph()
{
  //void
}

void MonoGraph::ReadInsiders(std::ifstream& in)
{
  int company;
  int insider;
  char delimit1;
  char delimit2;
  if (!in.is_open())
  {
    std::cerr << "file not open!! \n";
  }
  in >> company >> delimit1 >> delimit2;
  while (!in.eof() && in.is_open())
  {
    while (in >> insider >> delimit2)
    {
      insiderdict_[company].insert(insider);
      insider_of_[insider].insert(company);
      if (delimit2 == ']')
	    {
	      break;
	    }
    }
    in >> company >> delimit1 >> delimit2;
  }
  /*
  for (auto insider: insider_of_)
  {
    std::cerr << "insider " << insider.first << " had " << insider.second.size() << " companies \n";
  } 
  */
}

void MonoGraph::ReadBoardMembers(std::ifstream& in)
{
  int company;
  int insider;
  char delimit1;
  char delimit2;
  if (!in.is_open())
  {
    std::cerr << "Board file not open!! \n";
  }
  in >> company >> delimit1 >> delimit2;
  while (!in.eof() && in.is_open())
  {
    while (in >> insider >> delimit2)
    {
      boarddict_[company].insert(insider);
      board_of_[insider].insert(company);
      if (delimit2 == ']')
	    {
	      break;
	    }
    }
    in >> company >> delimit1 >> delimit2;
  }
  /*
  for (auto insider: board_of_)
  {
    std::cerr << "board member " << insider.first << " had " << insider.second.size() << " companies \n";
  } 
  */
}


const std::set<int> MonoGraph::GetInsider(int k) const
{
  auto foo = insiderdict_.find(k);
  if (foo == insiderdict_.end())
  {
    std::set<int> bar;
    return bar;
  }
  return insiderdict_.at(k);
}

const std::set<int> MonoGraph::GetInsiderOf(int i) const
{
  // std::cerr << "insider dictionary size: " << insiderdict_.size() << "\n";
  return insider_of_.at(i); 
}

const std::set<int> MonoGraph::GetBoard(int k)
{
  return boarddict_[k];
}

const std::set<int> MonoGraph::GetBoardOf(int i) 
{
  // std::cerr << "insider dictionary size: " << insiderdict_.size() << "\n";
  return board_of_[i]; 
}


const Alist2& MonoGraph::GetNeighbours(int i) const
{
  //std::cerr << "n\n";
  auto togo = adj_.find(i);
  if (togo == adj_.end())
  {
    return def_;
  }
  else
  {
    return togo->second;
  }
}

int MonoGraph::GetNumber() const
{
  return number_;
}

int 
MonoGraph::GetDistance(int comp, int node)
{
  NodeDict d;
  try
  {
    d = GetDistances(comp);
  }
  catch(const std::exception& e)
  {
    return -2;
  }
  

  if (d.find(node) == d.end())
  {
    //std::cerr << "distance for " << node << " in " << comp << "-1\n";
    return -1; //2*max_distance_;
  }
  //std::cerr << "distance for " << node << " in " << comp << ": " << d.at(node) << "\n";
  return d.at(node); 
}

/* Betweenness cntrality*/
double
MonoGraph::GetCentrality(int node)
{
  // This gives the betweenness-centrality
  //std::cerr << "cen. called\n";
  if (node < 0)
  {
    throw std::runtime_error("Negative node index");
  }
  if (centrality_.size() > 0 && static_cast<int>(centrality_.size()) > node)
  {
    //std::cerr << centrality_.size() << " and " << node << "\n"; 
    return centrality_.at(node);
  }
  //std::cerr << "cen. executing\n";
  if (centrality_.size() > 0)
  {
    throw std::runtime_error("node value too large");
  }
  // We calculate betweenness-centrality using Brandes- Algorithm
  // Journal of Math. Sociology 25(2001), 2, pp 163-177
  const int n = number_ + 1;
  //std::cerr << "n " << n << "\n";
  centrality_.resize(n,0.0);
  std::vector<std::list<int>> P(n);
  for (int s = 0; s < n; ++s)
  {
    std::list<int> S;
    std::vector<std::list<int>> P(n);
    std::list<int> Q; 
    std::vector<double> sigma(n,0);
    std::vector<int> d(n, -1);
    d[s] = 0;
    sigma[s] = 1.0; 
    Q.push_back(s);
    while (!Q.empty())
    {
      int u = Q.front();
      Q.pop_front();
      S.push_back(u);
      auto adj_it = adj_.find(u);
      if (adj_it == adj_.end())
      {
        // no neighbours
        continue;
      }
      for (int v: adj_it->second)
      {
        if (d[v] < 0)
        {
          d[v] = d[u] + 1;
          Q.push_back(v);
        }
        if (d[v] == d[u] + 1)
        { 
          sigma[v] = sigma[v] + sigma[u];
          P[v].push_back(u);
        } 
      }
    }
    std::vector<double> delta(n,0.0);
    while(!S.empty())
    {
      // How many from s to w 
      int w = S.back();
      S.pop_back();
      for (int v: P[w])
      {
        delta[v] = delta[v] + (sigma[v]/sigma[w])*(1 + delta[w]);
      }
      if(w != s)
      {
        centrality_[w] = centrality_[w] + delta[w];
      }
    }
  }  
  NormalizeCentrality();
  // std::cerr << " centrality calculated at least once\n";
  return centrality_[node];
}

void 
MonoGraph::NormalizeCentrality()
{
  double min = centrality_[0];
  double max = -1;
  for (int i = 0; i < number_ + 1; ++i)
  {
    if (centrality_[i] > max)
    {
      max = centrality_[i];
    }
    if (min > centrality_[i])
    {
      min = centrality_[i];
    }
  }
  for (int i = 0; i < number_ +1; ++i)
  {
    double c = centrality_[i];
    c = (c - min)/(max - min);
    centrality_[i] = c; 
  }
}


// Breadth first distances from company c
const std::unordered_map<int,int>& MonoGraph::GetDistances(int c)
{
  if (dists_.find(c) != dists_.end())
  {
    return dists_[c];
  }
  // std::cerr << "calculating distances for " << c << "\n";
  // Make new. 
  std::unordered_map<int,int>& d = dists_[c];
  // No insiders, return the empty. 
  if (insiderdict_.find(c) == insiderdict_.end())
  {
    std::cerr << "no insiders for " << c << "\n";
    throw std::runtime_error("no insiders");
  }
  std::list<int> q;
  // first initialize:
  for (int j: insiderdict_.at(c))
  {
    d[j] = 0;
    q.push_back(j);
  }
  while (!q.empty())
  {
    int node = q.front();
    int dist = d[node];
    q.pop_front();
    auto adj_list_it = adj_.find(node);
    for (int next: adj_list_it->second)
    {
      auto dist_it = d.find(next);
      if(dist_it != d.end())
      {
        continue;
      }
      d[next] = dist + 1;
      q.push_back(next);
    }
    if (dist >= max_distance_)
    {
      max_distance_ = dist + 1;
    }
  }
  return d;
}
/*
 * Degree of node
 */
double 
MonoGraph::Degree(int node)
{
  return (double) adj_[node].size();
}


/*
 * Normalized Degree of node:
 */
double 
MonoGraph::NormalDegree(int node)
{

  if (p_calculated_)
  {
    double kk = P_.at(node);
    return kk;
  } 
  P_.resize(number_+1, 0.0);
  double max = 0.0;
  for(int i = 0; i < number_+1; ++i)
  {    
    double c = (double) GetNeighbours(i).size();
    P_[i] = c;
    if (c > max)
    {
      max = c;
    }
  }
  for (int i = 0; i < number_+1; ++i)
  {
    P_[i] /= max;
  }
  p_calculated_ = true;
  return P_[node];
}


MetaGraph::MetaGraph(std::string dir, int rdate)
{
  DIR* dirp = opendir(dir.c_str());
  std::string fname = FHANDLE; 
  struct dirent * dp;
  number_ = 0;
  while((dp = readdir(dirp)) != NULL)
  {
    std::string filename = dir+dp->d_name+fname;
    std::stringstream s(dp->d_name);
    std::string insiders = dir+dp->d_name+CHANDLE; 
    std::string boards = dir+dp->d_name+BHANDLE; 
    int date;
    s >> date;
    if (rdate > 0 && date != rdate)
    {
      continue; 
    }
    std::ifstream in;
    in.open(filename);
    if (!in.is_open())
    {
      continue;
    }
    auto ss = new MonoGraph(in);
    in.close();
    graphs_[date] = ss;
    if (ss->GetNumber() > number_)
    {
      number_ = ss->GetNumber(); 
    }
    std::ifstream in2; 
    // std::cerr << "reading insiderfile: " << insiders << "\n"; 
    in2.open(insiders);
    // std::cerr << "trying to read file: " << insiders << "\n";
    ss->ReadInsiders(in2);
    in2.close();
    std::ifstream in3; 
    // std::cerr << "reading insiderfile: " << insiders << "\n"; 
    in3.open(boards);
    // std::cerr << "trying to read file: " << insiders << "\n";
    ss->ReadBoardMembers(in3);
    in3.close();
  }  

}
// Read only specific dates
MetaGraph::MetaGraph(std::string dir, const AnnouncementDates& an)
{
  DIR* dirp = opendir(dir.c_str());
  std::string fname = FHANDLE; 
  struct dirent * dp;
  number_ = 0;
  while((dp = readdir(dirp)) != NULL)
  {
    std::string filename = dir+dp->d_name+fname;
    std::stringstream s(dp->d_name);
    std::string insiders = dir+dp->d_name+CHANDLE; 
    int date;
    s >> date;
    std::ifstream in;
    if (an.find(date) == an.end())
    {
      continue;
    }
    in.open(filename);
    if (!in.is_open())
    {
      continue;
    }
    auto ss = new MonoGraph(in);
    in.close();
    graphs_[date] = ss;
    if (ss->GetNumber() > number_)
    {
      number_ = ss->GetNumber(); 
    }
    std::ifstream in2;
    in2.open(insiders);
    // std::cerr << "reading file " << insiders << "\n"; 
    ss->ReadInsiders(in2);
    in2.close();
  }  

}

MetaGraph::~MetaGraph()
{
  // void
}

MonoGraph* MetaGraph::GetGraph(int date) const
{
  if (graphs_.size() == 1)
  {
    return graphs_.begin()->second;
  }
  auto it = graphs_.lower_bound(date);
  if (it == graphs_.end())
  {
    return 0;
  }
  auto it2 = graphs_.upper_bound(date);
  if (it2 == graphs_.end())
  {
    return 0;
  }
  return it->second; 
}

int MetaGraph::GetNumber() const
{
  return number_; 
}

std::vector<int>
MetaGraph::GetDates() const
{
  std::vector<int> togo;
  for (auto x: graphs_)
  {
    togo.push_back(x.first);
  }
  std::sort(togo.begin(),togo.end());
  return togo; 
}

std::set<int> 
MonoGraph::GetMaxComp(std::vector<int> set) const
{
  //std::cerr << "foo\n";
  std::vector<std::set<int>> comps;
  std::set<int> foo(set.begin(),set.end());
  for (int i: set)
  {
    //std::cerr << "i:" << i << "\n";
    bool ff = false;
    for (std::set<int> found: comps)
    {
      if (found.find(i) != found.end())
      {
        ff = true;
        break;
      }
    }
    if (ff)
    {
      continue; 
    }
    std::list<int> S;
    std::set<int> F;
    S.push_back(i);
    F.insert(i);
    //std::cerr << "starting...";
    while (!S.empty())
    {
      int u = S.back();
      S.pop_back();
      for (int v: GetNeighbours(u))
      {
        if (foo.find(v) == foo.end() || F.find(v) != F.end())
        {
          continue;
        }
        S.push_back(v);
        F.insert(v);
      }
    }
    comps.push_back(F);
  }
  int i = 0;
  int j = 0;
  int k = 0;
  for (auto ss: comps)
  {
    if (static_cast<int>(ss.size()) > i)
    {
      j = k;
      i = ss.size();
    }
    ++k;
  }
  return comps[j];
}

bool 
MetaGraph::CheckValidity(int node)
{
  if (validated_nodes_.empty())
  {
    for (auto gg_pair: graphs_)
    {
      const MonoGraph* gg = gg_pair.second;
      for (auto a_pair: gg->adj_)
      {
        if(!a_pair.second.empty())
        {
          validated_nodes_.insert(a_pair.first);
        }
      }
    }
  }
  return (validated_nodes_.find(node) != validated_nodes_.end());
}
}
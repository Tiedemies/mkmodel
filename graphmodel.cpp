#include<iostream>
#include<sstream>
#include "defs.hpp"
#include "graphmodel.hpp"
#include<dirent.h>
#include<new>
#include<algorithm>
#include<list>

#define FHANDLE "/egde_list.txt"
#define CHANDLE "/company_dict_insiders.txt"


MonoGraph::MonoGraph(std::ifstream& in)
{
  // Read the graph in
  std::string c;
  int i = 0;
  adj_[0] = {};
  number_ = 0;
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
	      adj_[s_node] = Alist{t_node};
      }
      else
      {
	      adj_[s_node].push_back(t_node); 
      }
      if (s_node > number_) number_ = s_node;
      if (t_node > number_) number_ = t_node; 

      in >> c;
    }
    //std::cerr << adj_.size() << " ";
  }
  else
  {
    std::cerr << "Failed to open file.\n";
  }
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
      insiderdict_[company].push_back(insider);
      if (delimit2 == ']')
	    {
	      break;
	    }
    }
    in >> company >> delimit1 >> delimit2;
  } 
}

std::vector<int> MonoGraph::GetInsider(int k) const
{
  // std::cerr << "insider dictionary size: " << insiderdict_.size() << "\n";
  return insiderdict_.at(k); 
}


const Alist& MonoGraph::GetNeighbours(int i) const
{
  auto togo = adj_.find(i);
  if (togo == adj_.end())
  {
    return adj_.at(0);
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
  const auto d = GetDistances(comp);

  if (d.find(node) == d.end())
  {
    return 0;
  }
  return d.at(node); 
}

// Breadth first:
const std::unordered_map<int,int>& MonoGraph::GetDistances(int c)
{
  if (dists_.find(c) != dists_.end())
  {
    return dists_[c];
  }
  std::unordered_map<int,int>& d = dists_[c];
  if (insiderdict_.find(c) == insiderdict_.end())
  {
    return d; 
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
    if (adj_list_it == adj_.end())
    {
      continue;
    }
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
  }
  return d;
}


MetaGraph::MetaGraph(std::string dir)
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
    ss->ReadInsiders(in2);
    in2.close();
  }  

}

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
  auto it = graphs_.lower_bound(date);
  if (it == graphs_.end())
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

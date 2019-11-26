#pragma once

namespace Grid {
namespace QCD {



double linear_interpolation(double x, double x_lower, double y_lower, double x_upper, double y_upper) {
	return ( (x - x_lower) * y_upper +  (x_upper - x) * y_lower ) / (x_upper - x_lower);
}

int my_smod(int x, int L) {
  if(x<=L/2) return x;
  else return x - L;
}

std::vector<int> my_smod(const std::vector<int> &x, const std::vector<int> &L) {
  std::vector<int> ret(x.size());
  for(int i=0; i<x.size(); ++i) ret[i] = my_smod(x[i], L[i]);
  return ret;
}


Coordinate my_smod(const Coordinate &x, const Coordinate &L) {
  Coordinate ret(x.size());
  for(int i=0; i<x.size(); ++i) ret[i] = my_smod(x[i], L[i]);
  return ret;
}

template<class T>
double len(const std::vector<T> &vec){
  double ret = 0.;
  for(auto x: vec) ret += x * x;
  return std::sqrt(ret);
}

// void localIndexToLocalGlobalCoor(GridBase *grid, int ss, std::vector<int> &lcoor, std::vector<int> &gcoor) {
void localIndexToLocalGlobalCoor(GridBase *grid, int ss, Coordinate &lcoor, Coordinate &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret._grid->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  Coordinate processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}


// distance between two points with periodic boundary condition
int distance(int t1, int t2, int T) {
  int tmp = std::abs(t1 - t2);
  if(tmp <= T/2) return tmp;
  else return T-tmp;
}

int rightPoint(int t_base, int t, int T) { // after shift t_base to 0, determine wheter t is on the right or t_base is on the right
  int tmp = t - t_base; // shift t_base to 0
  if(tmp > 0 && tmp <=T/2) return t;
  else return t_base; // tmp < 0 || tmp > T/2
}

int leftPoint(int t_base, int t, int T) { // after shift t_base to 0, determine wheter t is on the left or t_base is on the left
  int tmp = t - t_base; // shift t_base to 0
  if(tmp > 0 && tmp <=T/2) return t_base;
  else return t; // tmp < 0 || tmp > T/2
}

}}

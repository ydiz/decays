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
  // ss is local index; parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++)
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

inline Coordinate operator-(const Coordinate &v1, const Coordinate &v2)
{
  const int N = v1.size();
  Coordinate ret(N);
  for(size_t i=0; i<N; ++i) ret[i] = v1[i] - v2[i];
  return ret;
}

template<class T>
Lattice<T> get_reflection(const Lattice<T> &lat) {
  Lattice<T> new_lat(lat.Grid());

  // pull corresponding piece from another node.
  // std::vector<int> pcoor;
  Coordinate pcoor;
  lat.Grid()->ProcessorCoorFromRank(lat.Grid()->ThisRank(), pcoor);
  // std::vector<int> new_pcoor = lat.Grid()->_processors - std::vector<int>{1,1,1,1} - pcoor;
  Coordinate new_pcoor = lat.Grid()->_processors - Coordinate(std::vector<int>{1,1,1,1}) - pcoor;
  int partner = lat.Grid()->RankFromProcessorCoor(new_pcoor);

  // int bytes = sizeof(T) * lat._odata.size();
  long int bytes = sizeof(T) * lat.View().size(); // if use int, there will be integer overflow for 48I ensemble
  MPI_Request recv_request, send_request, requests_array[2];

  long int cur_idx = 0;
  while(bytes>0) {            // the "count" parameter in MPI_Irecv and MPI_Isend are "int", so cannot be larger than INT_MAX
    int batch_size = std::min(bytes, (long int)INT_MAX);
    std::cout << "batch_size: " << batch_size << std::endl;
    assert(batch_size <= INT_MAX);
    MPI_Irecv((char *)&new_lat.View()[0] + cur_idx, batch_size, MPI_BYTE, partner, 0, MPI_COMM_WORLD, &recv_request); // non-bloking communication
    MPI_Isend((char *)&lat.View()[0] + cur_idx, batch_size, MPI_BYTE, partner, 0, MPI_COMM_WORLD, &send_request);

    bytes -= batch_size;
    cur_idx += batch_size;
    requests_array[0] = recv_request;  requests_array[1] = send_request;
    MPI_Waitall(2, requests_array, MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // reflection inside node
  Lattice<T> ret(lat.Grid());
  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++) {
    typename T::scalar_object m;
    // std::vector<int> lcoor(4);
    Coordinate lcoor(4);
    lat.Grid()->LocalIndexToLocalCoor(ss, lcoor);

    // std::vector<int> new_lcoor(4);
    Coordinate new_lcoor(4);
    new_lcoor = lat.Grid()->_ldimensions - Coordinate(std::vector<int>{1,1,1,1}) - lcoor;
    peekLocalSite(m, new_lat, new_lcoor);
    pokeLocalSite(m, ret, lcoor);
  }

  // after cshift, we would get the reflection we want. Note that under reflection, 0 -> 0, i -> L - i when i!=0
  for(int mu=0; mu<4; ++mu) ret = Cshift(ret, mu, -1);

  return ret;
}

// // test reflection
// Lattice<iScalar<iScalar<iVector<vComplex, 4>>>> coor(grid);
// coor=zero;
// for(int d=0;d<4;d++){
//   LatticeComplex t(grid);
//   LatticeCoordinate(t,d);
//   pokeColour(coor, t, d);
// }
// cout << get_reflection(coor) << endl;






}}

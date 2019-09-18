#pragma once

#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

inline std::vector<int> operator-(const std::vector<int> &v1, const std::vector<int> &v2)
{
  const int N = v1.size();
  std::vector<int> ret(N);
  for(size_t i=0; i<N; ++i) ret[i] = v1[i] - v2[i];
  return ret;
}

template<class T>
Lattice<T> get_reflection(const Lattice<T> &lat) {
  Lattice<T> new_lat(lat._grid);

  // pull corresponding piece from another node.
  std::vector<int> pcoor;
  lat._grid->ProcessorCoorFromRank(lat._grid->ThisRank(), pcoor);
  std::vector<int> new_pcoor = lat._grid->_processors - std::vector<int>{1,1,1,1} - pcoor;
  int partner = lat._grid->RankFromProcessorCoor(new_pcoor);

  int bytes = sizeof(T) * lat._odata.size();
  MPI_Request recv_request, send_request, requests_array[2];
  MPI_Irecv(new_lat._odata.data(), bytes, MPI_BYTE, partner, 0, MPI_COMM_WORLD, &recv_request); // non-bloking communication
  MPI_Isend(lat._odata.data(), bytes, MPI_BYTE, partner, 0, MPI_COMM_WORLD, &send_request);
  requests_array[0] = recv_request;  requests_array[1] = send_request;
  MPI_Waitall(2, requests_array, MPI_STATUSES_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);

  // reflection inside node
  Lattice<T> ret(lat._grid);
  parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    typename T::scalar_object m;
    std::vector<int> lcoor(4);
    lat._grid->LocalIndexToLocalCoor(ss, lcoor);

    std::vector<int> new_lcoor(4);
    new_lcoor = lat._grid->_ldimensions - std::vector<int>{1,1,1,1} - lcoor;
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


// The J -- pion loop in disconnected diagram
LatticeLoop loops_contraction(const std::vector<LatticePropagator> &wall_props, int t_min) {
  LatticeLoop ret(wall_props[0]._grid);

  Gamma gamma5(Gamma::Algebra::Gamma5);
  // Gamma::gmu[0], Gamma::gmu[1], Gamma::gmu[2], Gamma::gmu[3];// GammaX, GammaY, GammaZ, GammaT

  parallel_for(int ss=0; ss<ret._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(ret._grid, ss, lcoor, gcoor);

    // cheng's tsep
    int Tsize = ret._grid->_fdimensions[Tdir];
    int t_wall = gcoor[3] - t_min; 
    if(t_wall < 0) t_wall += Tsize; // if wall is on the left side of the current

    typename LatticePropagator::vector_object::scalar_object wall_to_v;//, v_to_wall;
    typename LatticeLoop::vector_object::scalar_object ret_site;

    peekLocalSite(wall_to_v, wall_props[t_wall], lcoor);

    ret_site = 0.;
    for(int nu=0; nu<4; ++nu) ret_site()()(nu) = trace(gamma5 * Gamma::gmu[nu] * wall_to_v * adj(wall_to_v));

    pokeLocalSite(ret_site, ret, lcoor);

  }
  return ret;
}


void get_loop_exp_minus(LatticeComplex &lat, int t_min, double Mpi) {
  int T = lat._grid->_fdimensions[Tdir];

  std::vector<double> exps(T);
  for(int t=0; t<T; ++t) {
    int xt = (t <= T/2) ? t : t - T;
    exps[t] = std::exp(Mpi * (t_min - 0.5 * xt));
  }

  parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    typename LatticeComplex::vector_object::scalar_object m;
    m()()() = Complex(exps[gcoor[3]], 0.);
    pokeLocalSite(m, lat, lcoor);
  }
}



void get_loop_exp_plus(LatticeComplex &lat, int t_min, double Mpi) {
  int T = lat._grid->_fdimensions[Tdir];

  std::vector<double> exps(T);
  for(int t=0; t<T; ++t) {
    int xt = (t <= T/2) ? t : t - T;
    exps[t] = std::exp(Mpi * (t_min + 0.5 * xt));
  }

  parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    typename LatticeComplex::vector_object::scalar_object m;
    m()()() = Complex(exps[gcoor[3]], 0.);
    pokeLocalSite(m, lat, lcoor);
  }
}


// BOth loop 1 and loop 2 should be purely imaginary
LatticePGG three_point_loop(const LatticeLoop &loop1, const LatticeLoop &loop2, int t_min, double Mpi) {

  FFT theFFT((GridCartesian *)loop1._grid);
  LatticeLoop loop1_fft(loop1._grid), loop2_fft(loop1._grid);
  LatticePGG ret(loop1._grid), ret_u(loop1._grid), ret_v(loop1._grid), ret_u_fft(loop1._grid),ret_v_fft(loop1._grid);
  LatticeComplex t(loop1._grid);

  LatticeLoop tmp = imag(loop1);
  theFFT.FFT_all_dim(loop1_fft, tmp, FFT::forward);
  tmp = imag(loop2);
  theFFT.FFT_all_dim(loop2_fft, tmp, FFT::forward);

  parallel_for(int ss=0; ss<ret._grid->oSites(); ss++) {
    for(int mu=0; mu<4; ++mu)
      for(int nu=0; nu<4; ++nu) {
        ret_u_fft[ss]()()(mu, nu) = adj(loop1_fft[ss]()()(mu)) * loop2_fft[ss]()()(nu);
        ret_v_fft[ss]()()(mu, nu) = adj(loop2_fft[ss]()()(mu)) * loop1_fft[ss]()()(nu);
      }
  }
  theFFT.FFT_all_dim(ret_u, ret_u_fft, FFT::backward);
  theFFT.FFT_all_dim(ret_v, ret_v_fft, FFT::backward);

  get_loop_exp_minus(t, t_min, Mpi);
  ret_u = ret_u * t;

  get_loop_exp_plus(t, t_min, Mpi);
  ret_v = ret_v * t;

  ret = ret_u + ret_v;
  const std::vector<int> &gc = ret._grid->_fdimensions;
  double V = gc[0] * gc[1] * gc[2] * gc[3];
  ret = ret * (- 1. / V);  // Note the -1 here. at the begining, we took imag(loop1), i.e. multiplied it by -i
  return ret;

}

void calculate_disc(LatticePGG &disc, const std::string &ensemble, int traj, int t_min, double Mpi) {

  LatticeLoop loop1(disc._grid);
  read_loop(loop1, loop_path_24ID(traj));

  int T = disc._grid->_fdimensions[Tdir];
  std::vector<LatticePropagator> wall_props(T, disc._grid);
  read_wall_src_props(ensemble, traj, wall_props);
  LatticeLoop loop2 = loops_contraction(wall_props, t_min);
  // // writeScidac(loop2, "loop2_1370.lat");

  // LatticeLoop loop2(grid);
  // readScidac(loop2, "loop2_1370.lat");

  disc = three_point_loop(loop1, loop2, t_min, Mpi);
}

}}

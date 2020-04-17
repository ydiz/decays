#pragma once

#include <qlat/qlat.h>
#include <Grid/Grid.h>
#include <headers/headers.h>

namespace Grid {
namespace QCD {


void read_wall_src_props(const std::string &ensemble, int traj, std::vector<LatticePropagator> &wall_props) {
  using namespace qlat;
  assert(wall_props.size() != 0);
  std::string gauge_transform_path;
  std::string (*wall_path)(int, int);

  if(ensemble == "24ID") {
    gauge_transform_path = gauge_transform_path_24D(traj);
    wall_path = wall_path_l_24ID;
  }
  else if(ensemble == "32ID") {
    gauge_transform_path = gauge_transform_path_32ID(traj);
    wall_path = wall_path_ud_32ID;
  }
  else assert(0);
  assert(dirExists(wall_path(traj, 0)));
  assert(dirExists(gauge_transform_path));

  // read gauge transformation
  std::cout << "Load Gauge Transform And Get Inv: " <<  gauge_transform_path << std::endl;
  GaugeTransform qlat_gtinv;
  {
    GaugeTransform qlat_gt;
    dist_read_field(qlat_gt, gauge_transform_path);
    to_from_big_endian_64(get_data(qlat_gt)); 
    gt_inverse(qlat_gtinv, qlat_gt);
  }
  LatticeColourMatrix gt(wall_props[0].Grid());
  grid_convert(gt, qlat_gtinv);

  std::cout << "reading wall source propagators and applying gauge transformations" << std::endl;
  for(int t=0; t<wall_props[0].Grid()->_fdimensions[Tdir]; ++t) {
    // read_qlat_propagator(wall_props[t], wall_subdirs[t]);
    read_qlat_propagator(wall_props[t], wall_path(traj, t));
    wall_props[t] = gt * wall_props[t];
  }
}






inline std::vector<int> operator-(const std::vector<int> &v1, const std::vector<int> &v2)
{
  const int N = v1.size();
  std::vector<int> ret(N);
  for(size_t i=0; i<N; ++i) ret[i] = v1[i] - v2[i];
  return ret;
}



// The J -- pion loop in disconnected diagram
LatticeLoop loops_contraction(const std::vector<LatticePropagator> &wall_props, int t_min) {
  LatticeLoop ret(wall_props[0].Grid());

  Gamma gamma5(Gamma::Algebra::Gamma5);
  // Gamma::gmu[0], Gamma::gmu[1], Gamma::gmu[2], Gamma::gmu[3];// GammaX, GammaY, GammaZ, GammaT

  parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++){

    // std::vector<int> lcoor, gcoor;
    // localIndexToLocalGlobalCoor(ret._grid, ss, lcoor, gcoor);

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(ret.Grid(), ss, lcoor, gcoor);
    // cheng's tsep
    int Tsize = ret.Grid()->_fdimensions[Tdir];
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
  int T = lat.Grid()->_fdimensions[Tdir];

  std::vector<double> exps(T);
  for(int t=0; t<T; ++t) {
    int xt = (t <= T/2) ? t : t - T;
    exps[t] = std::exp(Mpi * (t_min - 0.5 * xt));
  }

  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    // std::vector<int> lcoor, gcoor;
    // localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
    typename LatticeComplex::vector_object::scalar_object m;
    m()()() = Complex(exps[gcoor[3]], 0.);
    pokeLocalSite(m, lat, lcoor);
  }
}



void get_loop_exp_plus(LatticeComplex &lat, int t_min, double Mpi) {
  int T = lat.Grid()->_fdimensions[Tdir];

  std::vector<double> exps(T);
  for(int t=0; t<T; ++t) {
    int xt = (t <= T/2) ? t : t - T;
    exps[t] = std::exp(Mpi * (t_min + 0.5 * xt));
  }

  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    // std::vector<int> lcoor, gcoor;
    // localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
    typename LatticeComplex::vector_object::scalar_object m;
    m()()() = Complex(exps[gcoor[3]], 0.);
    pokeLocalSite(m, lat, lcoor);
  }
}


// BOth loop 1 and loop 2 should be purely imaginary
LatticePGG three_point_loop(const LatticeLoop &loop1, const LatticeLoop &loop2, int t_min, double Mpi) {

  FFT theFFT((GridCartesian *)loop1.Grid());
  LatticeLoop loop1_fft(loop1.Grid()), loop2_fft(loop1.Grid());
  LatticePGG ret(loop1.Grid()), ret_u(loop1.Grid()), ret_v(loop1.Grid()), ret_u_fft(loop1.Grid()),ret_v_fft(loop1.Grid());
  LatticeComplex t(loop1.Grid());

  LatticeLoop tmp = imag(loop1);   // I am taking two iamg(.) here, and multiplying -1 at the end of this function
  theFFT.FFT_all_dim(loop1_fft, tmp, FFT::forward);
  tmp = imag(loop2);
  theFFT.FFT_all_dim(loop2_fft, tmp, FFT::forward);

  parallel_for(int ss=0; ss<ret.Grid()->oSites(); ss++) {
    for(int mu=0; mu<4; ++mu)
      for(int nu=0; nu<4; ++nu) {
        auto loop1_v = loop1_fft.View();
        auto loop2_v = loop2_fft.View();
        auto ret_u_view = ret_u_fft.View();
        auto ret_v_view = ret_v_fft.View();

        // ret_u_fft[ss]()()(mu, nu) = adj(loop1_fft[ss]()()(mu)) * loop2_fft[ss]()()(nu);
        // ret_v_fft[ss]()()(mu, nu) = adj(loop2_fft[ss]()()(mu)) * loop1_fft[ss]()()(nu);
        ret_u_view[ss]()()(mu, nu) = adj(loop1_v[ss]()()(mu)) * loop2_v[ss]()()(nu);
        ret_v_view[ss]()()(mu, nu) = adj(loop2_v[ss]()()(mu)) * loop1_v[ss]()()(nu);
      }
  }
  theFFT.FFT_all_dim(ret_u, ret_u_fft, FFT::backward);
  theFFT.FFT_all_dim(ret_v, ret_v_fft, FFT::backward);

  get_loop_exp_minus(t, t_min, Mpi);
  ret_u = ret_u * t;

  get_loop_exp_plus(t, t_min, Mpi);
  ret_v = ret_v * t;

  ret = ret_u + ret_v;
  const Coordinate &gc = ret.Grid()->_fdimensions;
  double V = gc[0] * gc[1] * gc[2] * gc[3];
  ret = ret * (- 1. / V);  // Note the -1 here. at the begining, we took imag(loop1), i.e. multiplied it by -i
  return ret;

}

void calculate_disc(LatticePGG &disc, const std::string &ensemble, int traj, int t_min, double Mpi) {

  LatticeLoop loop1(disc.Grid());
  read_loop(loop1, loop_path_24ID(traj));

  int T = disc.Grid()->_fdimensions[Tdir];
  std::vector<LatticePropagator> wall_props(T, disc.Grid());
  read_wall_src_props(ensemble, traj, wall_props);
  LatticeLoop loop2 = loops_contraction(wall_props, t_min);
  // // writeScidac(loop2, "loop2_1370.lat");

  // LatticeLoop loop2(grid);
  // readScidac(loop2, "loop2_1370.lat");

  disc = three_point_loop(loop1, loop2, t_min, Mpi);
}

}}

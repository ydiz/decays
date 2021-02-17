#pragma once


// #include "/sdcc/u/ydzhao/A2AGrid/io.h"

namespace Grid {


class Env {
public:


  GridCartesian *grid;
  std::string ensemble;
  std::string out_prefix;
  std::vector<int> lat_size;


  LatticeColourMatrix *g_test;

  int traj;

  double M_K; // mass of kaon in lattice unit
  double M_pi; // mass of pion in lattice unit
  double N_K;  
  double N_pi;
  double Z_V;

  int N_pt_src; // number of point sources to average over; -1 means use all 

  std::vector<std::vector<int>> xgs_l;
  std::vector<std::vector<int>> xgs_s;

  void setup_traj(int _traj);
  LatticePropagator get_point(const std::vector<int> &src, char quark) const;
  std::vector<LatticePropagator> get_wall(char quark, bool useCoulombSink = false) const;
  // LatticePropagator get_wall(int t, char quark) const;
  LatticePropagator get_Lxx() const;
  LatticePropagator get_sequential(const std::vector<int> &src) const;
  LatticeColourMatrix get_gaugeTransform() const;
  // std::vector<LatticeFermionD> get_a2a(char vw) const;

  // void toCoulombVW(const LatticeColourMatrix &gt, std::vector<LatticeFermionD> &in) const;
  LatticePropagator toCoulombSink(const LatticeColourMatrix &gt, const LatticePropagator &in) const;

  Env(const std::string &_ensemble);

  std::string point_path(char quark) const;
  std::string point_path(const std::vector<int> &src, char quark) const;
  std::string sequential_path(const std::vector<int> &src) const;
  std::string wall_path(int t, char quark) const;
  std::string gauge_transform_path() const;
  std::string Lxx_path() const;

  std::vector<std::vector<int>> get_xgs(char quark);

};






Env::Env(const std::string &_ensemble) // cannot initialize grid in initializer list; it depends on lat_size
{

  ensemble = _ensemble;

  if(ensemble=="FreeField_8nt8") {
    lat_size = {8,8,8,8};
    M_pi = 0.14;
    N_pi = 50;

    M_K = 0.5;
    N_K = 50;

    Z_V = 0.7;

    out_prefix = "./xxxxxxxxx";
  }
  else assert(0);

  grid = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi()); 
  grid->show_decomposition();

  // assert(dirExists(out_prefix));


}



// Free field test // The source can only be at (0,0,0,0)
std::vector<std::vector<int>> Env::get_xgs(char quark) {
  // return {{0,0,0,0}, {1,2,3,4}, {6,4,7,2}};
  // return {{0,0,0,0}};
  return {{0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,0,3}, {0,0,0,4}, {0,0,0,5}, {0,0,0,6}, {0,0,0,7}};
}

void Env::setup_traj(int _traj) {
  using namespace std;
  traj = _traj;
  cout << "[traj: " << to_string(traj) << "]" << endl;
  // get point source locations

  xgs_l = get_xgs('l');
  xgs_s = get_xgs('s');
}


LatticePropagator Env::get_Lxx() const {
  LatticePropagator point_prop = get_point({0,0,0,0}, 'l');
  LatticePropagatorSite site;
  peekSite(site, point_prop, Coordinate({0,0,0,0}));

  LatticePropagator Lxx(grid);
  Lxx = site;

  return Lxx;
}


LatticeColourMatrix Env::get_gaugeTransform() const {
  LatticeColourMatrix gt(grid);
  gt = 1.0;
  return gt;
}



// std::vector<LatticeFermionD> Env::get_a2a(char vw) const {
//   int nl = 2000, nh = 768;
//   // int nl = 1, nh = 1;
//   std::string prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/timeDilutedA2AVectors";
//
//   std::vector<LatticeFermionD> low(nl, grid), high(nh, grid);
//   print_memory();
//
//   if(vw == 'v') {
//     A2AVectorsIo::read(low, prefix + "/vl", false, traj);
//     A2AVectorsIo::read(high, prefix + "/vh", false, traj);
//   }
//   else if(vw == 'w') {
//     A2AVectorsIo::read(low, prefix + "/wl", false, traj);
//     A2AVectorsIo::read(high, prefix + "/wh", false, traj);
//   }
//   else assert(0);
//
//   low.insert(low.end(), high.begin(), high.end());          high.clear();
//   print_memory();
//   return low;
// }



///////////////////////////////
// paths
//////////////////////////////


std::vector<LatticePropagator> Env::get_wall(char quark, bool useCoulombSink/* = false */) const {
  int T = grid->_fdimensions[Tdir];
  std::vector<LatticePropagator> wall_props(T, grid);  // sometimes this fails, do not know why.

  std::string prefix;
  if(quark == 'l')  prefix = "/home/ydzhao/cuth/decays/propagators/wall_s/free_field_test/free_field_props/wall_l/";   
  else if(quark == 's') prefix = "/home/ydzhao/cuth/decays/propagators/wall_s/free_field_test/free_field_props/wall_s/";

  for(int t=0; t<T; ++t) {
    std::string fname = prefix + "/" + std::to_string(t);
    readScidac_prop_f2d(wall_props[t], fname);
  }

  return wall_props;
}

LatticePropagator Env::get_point(const std::vector<int> &src, char quark) const {


  LatticePropagator point_prop(grid);

  std::string fname;
  if(quark == 'l')  fname="/home/ydzhao/cuth/decays/propagators/point_s/free_field_test/free_field_props/point_l/" + coor2str(src);   
  else if(quark == 's') fname="/home/ydzhao/cuth/decays/propagators/point_s/free_field_test/free_field_props/point_s/" + coor2str(src);
  readScidac_prop_f2d(point_prop, fname);

  return point_prop;
}


LatticePropagator Env::get_sequential(const std::vector<int> &src) const {

  std::string fname = "/home/ydzhao/cuth/decays/propagators/sequential/free_field_test/free_field_props/" + coor2str(src);
  LatticePropagator seq_prop(grid);
  readScidac_prop_f2d(seq_prop, fname);

  // for(int mu=0; mu<4; ++mu) {
  //   if(src[mu]!=0) seq_prop = Cshift(seq_prop, mu, -src[mu]);
  // }
  return seq_prop;
}



// void Env::toCoulombVW(const LatticeColourMatrix &gt, std::vector<LatticeFermionD> &in) const {
//   for(LatticeFermionD &x: in) x = gt * x;
// }

LatticePropagator Env::toCoulombSink(const LatticeColourMatrix &gt, const LatticePropagator &in) const {
  LatticePropagator out(in.Grid());
  out = gt * in;

  return out;
}


std::string Env::point_path(char quark) const {
  return "xxxxxxxxxxxxx";
}


std::string Env::point_path(const std::vector<int> &src, char quark) const {
  return "xxxxxxxxxxxxx";
}

std::string Env::sequential_path(const std::vector<int> &src) const {
  return "xxxxxxxxxxxxx";
}


std::string Env::wall_path(int t, char quark) const {
  return "xxxxxxxxxxxxx";
}

std::string Env::gauge_transform_path() const {
  return "xxxxxxxxxxxxx";
}

std::string Env::Lxx_path() const {
  return "xxxxxxxxxxxxx";
}










}






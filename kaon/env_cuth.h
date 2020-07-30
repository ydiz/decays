#pragma once


// #include "/sdcc/u/ydzhao/A2AGrid/io.h"

namespace Grid {
namespace QCD {


class Env {
public:
  // std::vector<int> lat_size;
  GridCartesian *grid;
  std::string ensemble;
  std::string out_prefix;
  std::vector<int> lat_size;

  int traj;

  double M_K; // mass of kaon in lattice unit
  double M_pi; // mass of pion in lattice unit
  double N_K;  
  double N_pi;
  double Z_V;

  // int T_wall;
  // int T_u;
  int N_pt_src; // number of point sources to average over; -1 means use all 

  // // type II
  // int T_wall_typeII;
  // int num_pairs;

  std::vector<std::vector<int>> xgs_l;
  std::vector<std::vector<int>> xgs_s;

#ifdef USE_SPARSE
  qlat::FieldSelection fsel;
#endif

  void setup_traj(int _traj);
  LatticePropagator get_point(const std::vector<int> &src, char quark) const;
  std::vector<LatticePropagator> get_wall(char quark, bool useCoulombSink = false) const;
  LatticePropagator get_wall(int t, char quark) const;
  LatticePropagator get_Lxx() const;
  LatticeColourMatrix get_gaugeTransform() const;
  std::vector<LatticeFermionD> get_a2a(char vw) const;

  void toCoulombVW(const LatticeColourMatrix &gt, std::vector<LatticeFermionD> &in) const;
  LatticePropagator toCoulombSink(const LatticeColourMatrix &gt, const LatticePropagator &in) const;

  Env(const std::string &_ensemble);

  std::string point_path(const std::vector<int> &src, char quark) const;
  std::string wall_path(int t, char quark) const;
  std::string gauge_transform_path() const;
  std::string Lxx_path() const;

  std::vector<std::vector<int>> get_xgs(char quark);

};

std::vector<std::vector<int>> Env::get_xgs(char quark) {
  std::vector<std::vector<int>> xgs;
  std::string path = point_path({}, quark);
  // std::cout << path << std::endl;

  DIR *dir;
  dir = opendir(path.c_str());
  assert(dir!=NULL); // make sure directory exists
  struct dirent *entry;
  while ((entry = readdir(dir)) != NULL) {
    // printf ("%s\n", entry->d_name);
    std::string fname = entry->d_name;
    if(isdigit(fname[0])) xgs.push_back(CSL2coor(fname)); // This is also dir "." and ".."
  }
  closedir(dir);
  return xgs;
}

Env::Env(const std::string &_ensemble) {

  ensemble = _ensemble;

  if(ensemble=="24ID") {
    lat_size = {24, 24, 24, 64};
    M_pi = 0.13975;
    N_pi = 51.561594;

    M_K = 0.504154;
    N_K = 53.3291; // FIXME: update this  // Using 1 traj, tsep=6, I got 53.531806
    std::cout << "FIXME: Using the N_h calculated from only one trajtory " << std::endl;

    Z_V = 0.72672;

    out_prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/results/";
  }
  else assert(0);

  grid = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  grid->show_decomposition();

  // assert(dirExists(out_prefix));
}

void Env::setup_traj(int _traj) {
  using namespace std;
  traj = _traj;
  cout << "[traj: " << to_string(traj) << "]" << endl;
  // get point source locations

#ifdef USE_MY_PROPAGATOR
  xgs_l = get_xgs('l');
  xgs_s = get_xgs('s');
#endif


#ifdef USE_SPARSE
  string fname_l = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/zyd_point_sources_l.txt";
  string fname_s = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/zyd_point_sources_s.txt";
  xgs_l = coor_from_file(fname_l);   
  xgs_s = coor_from_file(fname_s);

  const long spatial_vol = lat_size[0] * lat_size[1] * lat_size[2];
  const long n_per_tslice = spatial_vol / 16;
  string f1 = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/f-rank.field";

  fsel.init();
  qlat::read_field_selection(fsel, f1, n_per_tslice);
#endif
}


#ifdef USE_SPARSE
std::vector<LatticePropagator> Env::get_wall(char quark, bool useCoulombSink/* = false */) const {
  using namespace qlat;
  using namespace std;

  int T = grid->_fdimensions[Tdir];
  vector<LatticePropagator> wall_props(T, grid);  // sometimes this fails, do not know why.

  string type;
  if(quark == 'l') type = "0";
  else if (quark == 's') type = "1";

  for(int t=0; t<T; ++t) {
    // readScidac_prop_f2d(wall_props[t], wall_path(t, quark));
    string fname = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/prop-wall-src/tslice=" + to_string(t) + " ; type=" + type + " ; accuracy=1.sfield";
    std::cout << GridLogMessage << fname << std::endl;
    assert(dirExists(fname));

    Propagator4d qlat_prop;
    read_selected_field_double_from_float(qlat_prop, fname, fsel);

    grid_convert(wall_props[t], qlat_prop);
  }

  if(useCoulombSink) {}   // For the stored Coulomb propagators, sink is in Coulomb gauge
  else {
    qlat::GaugeTransform qlat_gt;
    string gt_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/gauge-transform.field";
    std::cout << GridLogMessage << gt_path << std::endl;
    assert(dirExists(gt_path));
    read_field_double(qlat_gt, gt_path);

    qlat::GaugeTransform qlat_gtinv;
    qlat::gt_inverse(qlat_gtinv, qlat_gt); // FIXME: make sure that I should use qlat_gtinv instead of qlat_gt

    LatticeColourMatrix gt(grid);
    grid_convert(gt, qlat_gtinv);

    for(int t=0; t<T; ++t) wall_props[t] = gt * wall_props[t];
  }
  return wall_props;
}

LatticePropagator Env::get_point(const std::vector<int> &src, char quark) const {
  using namespace qlat;
  using namespace std;

  LatticePropagator point_prop(grid);

  string fname = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/propagators/results=" + to_string(traj) + "/prop-point-src/xg=(" + coor2CSL(src)  +") ; type=0 ; accuracy=0.sfield";
  std::cout << GridLogMessage << fname << std::endl;
  assert(dirExists(fname));
  Propagator4d qlat_prop;
  read_selected_field_double_from_float(qlat_prop, fname, fsel);

  qlat::print_qlat_field(qlat_prop);
  assert(0);
  grid_convert(point_prop, qlat_prop);

  print_grid_field_site(point_prop, {10,0,0,0});
  return point_prop;
}
#endif





LatticePropagator Env::get_Lxx() const {
  LatticePropagator Lxx(grid);
  readScidac(Lxx, Lxx_path());
  return Lxx;
}


LatticeColourMatrix Env::get_gaugeTransform() const {
  LatticeColourMatrix gt(grid);
  readScidac(gt, gauge_transform_path());
  return gt;
}



std::vector<LatticeFermionD> Env::get_a2a(char vw) const {
  // int nl = 2000, nh = 768;
  // // int nl = 1, nh = 1;
  // std::string prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID";
  //
  // std::vector<LatticeFermionD> low(nl, grid), high(nh, grid);
  // print_memory();
  //
  // if(vw == 'v') {
  //   A2AVectorsIo::read(low, prefix + "/MADWF_A2AVector_vl", false, traj);
  //   A2AVectorsIo::read(high, prefix + "/MADWF_A2AVector_vh", false, traj);
  // }
  // else if(vw == 'w') {
  //   A2AVectorsIo::read(low, prefix + "/MADWF_A2AVector_wl", false, traj);
  //   A2AVectorsIo::read(high, prefix + "/MADWF_A2AVector_wh", false, traj);
  // }
  // else assert(0);
  //
  // low.insert(low.end(), high.begin(), high.end());          high.clear();
  // print_memory();
  // return low;
}

///////////////////////////////
// paths
//////////////////////////////


#ifdef USE_MY_PROPAGATOR
std::vector<LatticePropagator> Env::get_wall(char quark, bool useCoulombSink/* = false */) const {
  int T = grid->_fdimensions[Tdir];
  std::cout << "before allocating vector of wall source propagators" << std::endl;
  print_memory();
  std::vector<LatticePropagator> wall_props(T, grid);  // sometimes this fails, do not know why.
  std::cout << "after allocating vector of wall source propagators" << std::endl;
  print_memory();
  for(int t=0; t<T; ++t) {
    readScidac_prop_f2d(wall_props[t], wall_path(t, quark));
  }


  if(useCoulombSink) {
    LatticeColourMatrix gt(grid);
    readScidac(gt, gauge_transform_path());
    for(int t=0; t<T; ++t) wall_props[t] = gt * wall_props[t];
  }

  return wall_props;
}

LatticePropagator Env::get_point(const std::vector<int> &src, char quark) const {
  LatticePropagator point_prop(grid);
  readScidac_prop_f2d(point_prop, point_path(src, quark));
  return point_prop;
}
#endif   // end of #ifdef USE_MY_PROPAGATOR

LatticePropagator Env::get_wall(int t, char quark) const {
  LatticePropagator wall_prop(grid);
  readScidac_prop_f2d(wall_prop, wall_path(t, quark));
  return wall_prop;
}


void Env::toCoulombVW(const LatticeColourMatrix &gt, std::vector<LatticeFermionD> &in) const {
  for(LatticeFermionD &x: in) x = gt * x;
}

LatticePropagator Env::toCoulombSink(const LatticeColourMatrix &gt, const LatticePropagator &in) const {
  LatticePropagator out(in.Grid());
  out = gt * in;

  return out;
}



std::string Env::point_path(const std::vector<int> &src, char quark) const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/point_" + std::string(1, quark) + "/" + std::to_string(traj) + "/" + coor2CSL(src);
  else assert(0);
  std::cout << "reading from " << path << std::endl;
  assert(dirExists(path));
  return path;
}

std::string Env::wall_path(int t, char quark) const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/wall_" + std::string(1, quark) + "/"  + std::to_string(traj) + "/" + std::to_string(t);
  else assert(0);
  std::cout << "reading from " << path << std::endl;
  assert(dirExists(path));
  return path;
}

std::string Env::gauge_transform_path() const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/gauge_transform/"  + std::to_string(traj) ;
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::Lxx_path() const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/Lxx/Lxx." + std::to_string(traj); // calculated with A2A propagator without time dilution
  else assert(0);
  assert(dirExists(path));
  return path;
}


}}

#pragma once

namespace Grid {
namespace QCD {

// g5, gmu, gmu5, gL
/////cannot use Gamma::gmu to initialize global variables. Gamma::gmu is uninitialzed at this point. Do not know why.

Gamma g5(Gamma::Algebra::Gamma5);
std::array<const Gamma, 4> gmu = {Gamma(Gamma::Algebra::GammaX), Gamma(Gamma::Algebra::GammaY), Gamma(Gamma::Algebra::GammaZ), Gamma(Gamma::Algebra::GammaT)};
std::array<const Gamma, 4> gmu5 = { // gamma_mu * gamma5
                                  Gamma(Gamma::Algebra::GammaXGamma5),
                                  Gamma(Gamma::Algebra::GammaYGamma5),
                                  Gamma(Gamma::Algebra::GammaZGamma5),
                                  Gamma(Gamma::Algebra::GammaTGamma5)};

// Gamma-left matrices gL_mu = g_mu*(1 - g5)
std::array<const GammaL, 4> gL {GammaL(gmu[0]), GammaL(gmu[1]), GammaL(gmu[2]), GammaL(gmu[3])};


class Env {
public:
  // std::vector<int> lat_size;
  GridCartesian *grid;
  std::string ensemble;
  int traj;
  double M_K; // mass of kaon in lattice unit
  double wilson_c1; // wilson coefficient for q1
  double wilson_c2; // wilson coefficient for q2

  std::vector<std::vector<int>> xgs_l;
  std::vector<std::vector<int>> xgs_s;

  void setup_traj(int _traj);
  LatticePropagator get_point_l(const std::vector<int> &src) const;
  LatticePropagator get_point_s(const std::vector<int> &src) const;
  std::vector<LatticePropagator> get_wall_l() const;
  std::vector<LatticePropagator> get_wall_s() const;

  Env(const std::vector<int> &_lat, const std::string &_ensemble);

private:
  std::string point_path_l() const;
  std::string point_path_s() const;
  std::string wall_path_l(int t) const;
  std::string wall_path_s(int t) const;
  std::string gauge_transform_path() const;

  void read_wall_src_props(std::vector<LatticePropagator> &wall_props, char quark) const;

  std::map<std::vector<int>, std::string> point_subdirs_l;
  std::map<std::vector<int>, std::string> point_subdirs_s;
};

Env::Env(const std::vector<int> &_lat, const std::string &_ensemble) {
  grid = SpaceTimeGrid::makeFourDimGrid(_lat, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  ensemble = _ensemble;

  if(ensemble=="24ID") {
    M_K = 0.504154;
    wilson_c1 = -0.3735505346; // FIXME: need to update wilson coefficient // this number is Daiqian thesis for mu=1.53GeV
    wilson_c2 = 1.189590707;
  }
  else assert(0);
}

void Env::setup_traj(int _traj) {
  traj = _traj;
  std::cout << "[traj: " << std::to_string(traj) << "]" << std::endl;
  // get point source locations
  get_xgs(point_path_l(), xgs_l, point_subdirs_l, 'l');
  get_xgs(point_path_s(), xgs_s, point_subdirs_s, 's');
}

LatticePropagator Env::get_point_l(const std::vector<int> &src) const {
  LatticePropagator point_prop(grid);
  read_qlat_propagator_no_dist(point_prop, point_subdirs_l.at(src));
  return point_prop;
}

LatticePropagator Env::get_point_s(const std::vector<int> &src) const {
  LatticePropagator point_prop(grid);
  read_qlat_propagator_no_dist(point_prop, point_subdirs_s.at(src));
  return point_prop;
}

std::vector<LatticePropagator> Env::get_wall_l() const {
  std::vector<LatticePropagator> wall_props(grid->_fdimensions[Tdir], grid);
  read_wall_src_props(wall_props, 'l');
  return wall_props;
}


std::vector<LatticePropagator> Env::get_wall_s() const {
  std::vector<LatticePropagator> wall_props(grid->_fdimensions[Tdir], grid);
  read_wall_src_props(wall_props, 's');
  return wall_props;
}


////////////////////////////
//
//


void Env::read_wall_src_props(std::vector<LatticePropagator> &wall_props, char quark) const {
  using namespace qlat;
  assert(wall_props.size() != 0);

  // read gauge transformation
  GaugeTransform qlat_gtinv;
  {
    GaugeTransform qlat_gt;
    dist_read_field(qlat_gt, gauge_transform_path());
    to_from_big_endian_64(get_data(qlat_gt)); 
    gt_inverse(qlat_gtinv, qlat_gt);
  }
  LatticeColourMatrix gt(wall_props[0]._grid);
  grid_convert(gt, qlat_gtinv);

  std::cout << "reading wall source propagators and applying gauge transformations" << std::endl;
  for(int t=0; t<wall_props[0]._grid->_fdimensions[Tdir]; ++t) {
    if(quark=='l') read_qlat_propagator(wall_props[t], wall_path_l(t));
    else if(quark=='s') read_qlat_propagator(wall_props[t], wall_path_s(t));
    else assert(0);

    wall_props[t] = gt * wall_props[t];
  }
}


///////////////////////////////
// paths
//////////////////////////////

std::string Env::point_path_l() const {
  std::string path;
  if(ensemble=="24ID") path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::point_path_s() const {
  return point_path_l(); // the same dir path
}

std::string Env::wall_path_l(int t) const {
  std::string path;
  if(ensemble=="24ID") path ="/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::wall_path_s(int t) const {
  std::string path;
  if(ensemble=="24ID") path ="/home/ljin/application/Public/Qlat-CPS-cc/jobs/wall-src-strange/results/24D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator/strange ; t=" + std::to_string(t);
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::gauge_transform_path() const {
  std::string path;
  if(ensemble=="24ID") path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
  else assert(0);
  assert(dirExists(path));
  return path;
}




}}

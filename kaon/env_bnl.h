#pragma once

namespace Grid {
namespace QCD {


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
  LatticePropagator get_point(const std::vector<int> &src, char quark) const;
  std::vector<LatticePropagator> get_wall(char quark) const;
  LatticePropagator get_wall(int t, char quark) const;

  Env(const std::vector<int> &_lat, const std::string &_ensemble);

private:
  std::string point_path(const std::vector<int> &src, char quark) const;
  std::string wall_path(int t, char quark) const;
  std::string gauge_transform_path() const;

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
  xgs_l = get_xgs('l');
  xgs_s = get_xgs('s');
}

std::vector<LatticePropagator> Env::get_wall(char quark) const {
  int T = grid->_fdimensions[Tdir];
  std::vector<LatticePropagator> wall_props(T, grid);
  for(int t=0; t<T; ++t) {
    readScidac_prop_f2d(wall_props[t], wall_path(t, quark));
  }
  return wall_props;
}

LatticePropagator Env::get_wall(int t, char quark) const {
  LatticePropagator wall_prop(grid);
  readScidac_prop_f2d(wall_prop, wall_path(t, quark));
  return wall_prop;
}


LatticePropagator Env::get_point(const std::vector<int> &src, char quark) const {
  LatticePropagator point_prop(grid);
  readScidac_prop_f2d(point_prop, point_path(src, quark));
  return point_prop;
}



///////////////////////////////
// paths
//////////////////////////////

std::string Env::point_path(const std::vector<int> &src, char quark) const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/point_" + std::string(1, quark) + "/" + std::to_string(traj) + "/" + coor2CSL(src);
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::wall_path(int t, char quark) const {
  std::string path;
  if(ensemble=="24ID") path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/wall_" + std::string(1, quark) + "/"  + std::to_string(traj) + "/" + std::to_string(t);
  else assert(0);
  assert(dirExists(path));
  return path;
}

std::string Env::gauge_transform_path() const {
  std::string path;
  if(ensemble=="24ID") path = "TBD";
  else assert(0);
  assert(dirExists(path));
  return path;
}


}}

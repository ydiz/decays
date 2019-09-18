#pragma once

#include <qlat/qlat.h>
#include <qlat/grid.h>
#include <dirent.h>
#include <sys/stat.h>
#include <headers/headers.h>

#include "connected.h"

#define M_PION 0.139474

namespace Grid {
namespace QCD{




void writeSites(const std::vector<std::vector<int>> &xgs, const std::vector<typename LatticePropagator::vector_object::scalar_object> &sites, const std::string &filename) {
  using namespace std;
  assert(xgs.size()==sites.size());
  std::string header = "HEADER_START\nNUM="+ std::to_string(xgs.size()) + "\nHEADER_END\n";
  std::ofstream fout(filename, std::ofstream::binary);
  fout << header;

  int xg_size = 4 * sizeof(xgs[0][0]);
  for(const auto &xg: xgs) fout.write((char *)xg.data(), xg_size);

  fout.write((char *)sites.data(), sites.size() * sizeof(sites[0]));
  fout.close();

}



void readSites(std::vector<std::vector<int>> &xgs, std::vector<typename LatticePropagator::vector_object::scalar_object> &sites, const std::string &filename) {
  using namespace std;
  std::ifstream fin(filename, std::ifstream::binary);
  int num;
  std::string tmp;
  std::getline(fin, tmp);
  assert(tmp == "HEADER_START");
  while(tmp != "HEADER_END") {
    if(tmp.substr(0,3)=="NUM") num = std::stoi(tmp.substr(4));
    std::getline(fin, tmp);
  }
  std::cout << "NUM: " << num << std::endl;

  xgs.resize(num); sites.resize(num);
  int xg_size = 4 * sizeof(xgs[0][0]);
  for(auto &xg: xgs) {
    xg.resize(4);
    fin.read((char *)xg.data(), xg_size);
  }

  fin.read((char *)sites.data(), sites.size() * sizeof(sites[0]));
  fin.close();
}


iSpinColourMatrix<Complex> getPointPropSite(const std::string &path, const std::vector<int>& xg) {
  qlat::Propagator4d qlat_prop;
  dist_read_field_double_from_float(qlat_prop, path);

  const qlat::Geometry& geo = qlat_prop.geo;
  const qlat::Coordinate lx = geo.coordinate_l_from_g(qlat::Coordinate(xg[0], xg[1], xg[2], xg[3]));
  qlat::WilsonMatrix qlat_site;
  if (geo.is_local(lx)) {
    qlat_site = qlat_prop.get_elem(lx);
  } else {
    set_zero(qlat_site);
  }
  glb_sum_double(qlat_site);

  iSpinColourMatrix<Complex> grid_site;
  assert(sizeof(qlat_site) == sizeof(grid_site));

  Complex *p_qlat = (Complex *)&qlat_site; // T is either ComplexF or ComplexD

  for(int row=0; row<12; ++row)
    for(int column=0; column<12; ++column) {
      int grid_spin_row = row/3;
      int grid_color_row = row%3;
      int grid_spin_column = column/3;
      int grid_color_column = column%3;
      grid_site()(grid_spin_row, grid_spin_column)(grid_color_row, grid_color_column) = *(p_qlat + 12*row + column);
    }

  return grid_site;
}


void generatePointLoop(const std::string &ensemble, int traj) {
  std::vector<int> gcoor;
  std::string point_src_path;
  if(ensemble == "32ID") {
    point_src_path = point_path_32ID(traj);
    gcoor = {32,32,32,64};
  }
  else assert(0);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  std::vector<std::vector<int>> xgs;
  std::map<std::vector<int>, std::string> point_subdirs;
  get_xgs(point_src_path, xgs, point_subdirs);
  std::cout << "NUM: " << xgs.size() << std::endl;

  // xgs.erase(xgs.begin()+2, xgs.end());//for test
  // int i=0;
  std::vector<typename LatticePropagator::vector_object::scalar_object> sites;
  for(const auto &xg: xgs) {

    std::cout << GridLogMessage << "xg of point src: " << xg << std::endl;

    // LatticePropagator point_prop(grid);
    // read_qlat_propagator(point_prop, point_subdirs[xg]);
    // typename LatticePropagator::vector_object::scalar_object site;
    // peekSite(site, point_prop, xg);

    typename LatticePropagator::vector_object::scalar_object site;
    site = getPointPropSite(point_subdirs[xg], xg);
    // std::cout << site << std::endl;

    sites.push_back(site);
    // ++i;
    // if(i==2) break;
  }
  
  if(grid->ThisRank()==0) writeSites(xgs, sites,  "/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/pointPropSites/pointPropSites." + std::to_string(traj));
}




void get_loop2(std::vector<LatticeLoop> &loop2s, const std::string &ensemble, int traj, double Mpi) {
  int T = loop2s[0]._grid->_fdimensions[Tdir];
  std::vector<LatticePropagator> wall_props(T, loop2s[0]._grid);
  read_wall_src_props(ensemble, traj, wall_props);

  Gamma gamma5(Gamma::Algebra::Gamma5);
  for(int xg_t=0; xg_t<T; ++xg_t) {
    parallel_for(int ss=0; ss<loop2s[0]._grid->lSites(); ss++){

      std::vector<int> lcoor, gcoor;
      localIndexToLocalGlobalCoor(loop2s[0]._grid, ss, lcoor, gcoor);

      int t_min = 10;
      int t_wall = leftPoint(xg_t, gcoor[3], T) - t_min; 
      if(t_wall < 0) t_wall += T; // if wall is on the left side of the current
      int t_sep = distance(xg_t, t_wall, T);

      typename LatticePropagator::vector_object::scalar_object wall_to_v;//, v_to_wall;
      typename LatticeLoop::vector_object::scalar_object loop2_site;

      peekLocalSite(wall_to_v, wall_props[t_wall], lcoor);

      // FIXME: can be optimized: 1. calculate wall * wall first 2. store exps
      for(int nu=0; nu<4; ++nu) loop2_site()()(nu) = trace(gamma5 * Gamma::gmu[nu] * wall_to_v * adj(wall_to_v));
      loop2_site = loop2_site * std::exp(Mpi * t_sep);
      pokeLocalSite(loop2_site, loop2s[xg_t], lcoor);
    }
  }
}



// LatticePGG three_point_disc2(const std::vector<LatticePropagator> &wall_props, const LatticePropagator &point_prop, const std::vector<int> &xg) {
void three_point_disc2(LatticePGG &ret, const std::string &ensemble, int traj, double Mpi) {

  ret = zero;

  std::vector<std::vector<int>> xgs;
  std::vector<typename LatticePropagator::vector_object::scalar_object> sites;
  readSites(xgs, sites, "/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/pointPropSites/pointPropSites." + std::to_string(traj));

  int T = ret._grid->_fdimensions[Tdir];
  std::vector<LatticeLoop> loop2s(T, ret._grid);
  get_loop2(loop2s, ensemble, traj, Mpi);

  LatticePGG pgg(ret._grid);

  for(int i=0; i<xgs.size(); ++i) {
    std::cout << GridLogMessage << "xg of point src: " << xgs[i] << std::endl;

    typename LatticeLoop::vector_object::scalar_object loop1_xg = zero;
    for(int mu=0; mu<4; ++mu) {
      loop1_xg()()(mu) = trace(Gamma::gmu[mu] * sites[i]);
    }
    // std::cout << loop_xg << std::endl; // should be purely imaginary

    const int t = xgs[i][Tdir];
    parallel_for(int ss=0; ss<ret._grid->oSites(); ss++) {
      for(int mu=0; mu<4; ++mu)
        for(int nu=0; nu<4; ++nu) {
          pgg[ss]()()(mu, nu) = loop1_xg()()(mu) * loop2s[t][ss]()()(nu);
        }
    }

    for(int mu=0; mu<4; ++mu) pgg = Cshift(pgg, mu, xgs[i][mu]);
    ret += pgg;
  }

  ret = ret * (1. / xgs.size());
}

// if(grid->ThisRank()==2) {
//   cout.clear();
// cout << xgs << endl;
// cout << sites << endl;
// std::cout << "Finished"  << std::endl;
// }





}}

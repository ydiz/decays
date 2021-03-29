#pragma once

#include <alg/a2a/a2a_fields.h>

CPS_START_NAMESPACE

void zyd_metadata(const std::string &path, int &Ls, int &N_evec) {
  using namespace std;

  Ls = -1; N_evec = -1;
  ifstream f(path + "/metadata.txt");
  string str;
  while(getline(f, str)) {
    if(str.substr(0,4)=="s[4]") Ls = stoi(str.substr( str.find('=') + 1 ));
    if(str.substr(0,4)=="neig") N_evec = stoi(str.substr( str.find('=') + 1 ));
  }
  f.close();
  std::cout << "ZMobius Ls: " << Ls << std::endl;
  std::cout << "N_evec: " << N_evec << std::endl;

  if(Ls == -1 || N_evec == -1) {
    std::cout << "Not found s[4] or neig" << std::endl;
    assert(0);
  }
}

// void zyd_read_compressed(const std::string &evec_dir, GnoneFgridZmobius &lat, std::vector<double> &evals, std::vector<Grid::LatticeFermionF> &evecs) {
void zyd_read_compressed(const std::string &evec_dir, Grid::GridCartesian *FGridF_s, Grid::GridRedBlackCartesian *FrbGridF_s, typename A2ApoliciesDoubleAutoAlloc::FgridGFclass &lat, std::vector<double> &evals, std::vector<Grid::LatticeFermionF> &evecs) {

  int ZMobius_Ls, N_evec;
  zyd_metadata(evec_dir, ZMobius_Ls, N_evec); // read N_evec and ZMobius_Ls from metadata.txt

  int original_SnodeSites = GJP.SnodeSites();
  GJP.SnodeSites(ZMobius_Ls);   // Must set Ls to ZMobius Ls, will set it back after reading eigenvectors.
  assert(GJP.Snodes()==1);

  const int FsiteSize = 2 * 3 * 4 * ZMobius_Ls; // lat.FsiteSize() // FsiteSize = Re/Im * color * spin * ZMobius Ls
  std::cout << __func__ << " FsiteSize(): "<< FsiteSize << std::endl;
  size_t evec_size = (size_t) (GJP.VolNodeSites () / 2) * FsiteSize;

  const char *evec_name = "light_evec";
  EigenCacheGrid<Grid::ZMobiusFermionF::FermionField> *ecache = new EigenCacheGrid<Grid::ZMobiusFermionF::FermionField>(evec_name);
  ecache->alloc(N_evec, evec_size, sizeof(float)); // compressed eigenvectors must be saved in float
  std::cout << __func__ << " ecache->alloc() finished " << std::endl;
  ecache->read_compressed(evec_dir.c_str());
  std::cout << __func__ << " ecache->read_compressed() finished " << std::endl;

  // copy evals
  evals = ecache->evals;

  // // Convert evecs to Grid LatticeFermionF
  // Grid::GridCartesian *FGridF_s = Grid::SpaceTimeGrid::makeFiveDimGrid(ZMobius_Ls, lat.getUGridF());
  // Grid::GridRedBlackCartesian *FrbGridF_s = Grid::SpaceTimeGrid::makeFiveDimRedBlackGrid(ZMobius_Ls, lat.getUGridF());

  evecs.assign(N_evec, FrbGridF_s); // LatticeFermionF == Grid::QCD::ZMobiusFermionF::FermionField
  Grid::LatticeFermionF grid_f(FGridF_s);  // Note: must be "F", float
  for(int i = 0; i < N_evec; i++){
    cps::Vector *cps_vec = (cps::Vector*) ecache->evecs[i];
    lat.ImpexFermion<float, Grid::LatticeFermionF, Grid::SpinColourVectorF>(cps_vec, grid_f, 1, cps::Odd, NULL);
    Grid::pickCheckerboard(Grid::Odd, evecs[i], grid_f); // full sites -> half sites
    // // TBD: free ecache->evecs[i]
    // p.s. cannot use lat.ImportFermion
  }
  std::cout << __func__ << " Finished converting to Grid LatticeFermion" << std::endl;

  // Finished
  GJP.SnodeSites(original_SnodeSites);   
  ecache->dealloc();
  std::cout << __func__ << " Finished " << std::endl;
}


CPS_END_NAMESPACE




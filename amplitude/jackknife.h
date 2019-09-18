#pragma once

#include <Grid/Grid.h>
// #include "lep_para.h"
#include "lep.h"
#include "lep_CUBA3d.h"
#include "imaginary_part.h"
#include "form_factor.h"

namespace Grid{
namespace QCD{

RealD mean(const std::vector<RealD>& data)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/RealD(N);
}

RealD jack_mean(const std::vector<RealD>& data, int sample)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/RealD(N-1);
}

RealD jack_std(const std::vector<RealD>& jacks, RealD mean)
{
    int N = jacks.size();
    RealD std(0.0);
    for(int i=0; i<N; ++i){ std += std::pow(jacks[i]-mean, 2.0); }
    return std::sqrt(RealD(N-1)/RealD(N)*std);
}

std::vector<RealD> jack_stats(const std::vector<RealD>& data)
{
  int N = data.size();
  std::vector<RealD> jack_samples(N);
  std::vector<RealD> jack_stats(2);

  jack_stats[0] = mean(data);
  for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
  jack_stats[1] = jack_std(jack_samples, jack_stats[0]);

  std::cout << "jackknife average: " << jack_stats[0] << std::endl;
  std::cout << "jackknife error: " << jack_stats[1] << std::endl;

  return jack_stats;
}


struct Jack_para {

  // hadronic part
  std::string ensemble;
  double M_h; // mass of pion/kaon in lattice unit
  double N_h; // normalization of wall src operator. N_h = <0 | pi(0) | pi>
  double Z_V;
  double hadron_coeff;
  double BR_coeff; // decay rate = BR_coeff * |M|^2; BR does not contain lattice quantities. (use physical electron and pion mass)

  // leptonic part
  std::string target;
  std::string file_p3;
  std::string file_p1;
  double lep_coeff;
  double leptonic_space_limit;
  double leptonic_time_limit;

  std::vector<int> lat_size;
  std::vector<int> traj_skip; // skip these trajectories
  int traj_start, traj_end, traj_sep, traj_num;
  int time_cutoff_start, time_cutoff_end, time_cutoff_num;

  void get_leptonic(LatticePGG &lat);
  void get_three_point(LatticePGG &three_point, int traj);
  std::vector<double> get_result_with_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic);
  
};

void Jack_para::get_leptonic(LatticePGG &lat) {
  if(target == "real_CUBA3d" || target == "imag_CUBA3d") get_leptonic_CUBA3d(file_p1, file_p3, lat, leptonic_space_limit, leptonic_time_limit);
  else if(target == "real") Grid::QCD::get_leptonic(file_p3, lat, leptonic_space_limit, leptonic_time_limit);
  else if(target == "imag_analytic") imag_part(lat, M_h);
  else if(target == "form_factor") form_factor_integrand(lat, M_h); // technically this is not leptonic part; but for convience I put it in get_leptonic function
  else assert(0);
}

// return <Jmu(-w/2) Jnu(w/2) |pi> / hadron_coeff
void Jack_para::get_three_point(LatticePGG &three_point, int traj) {

	static LatticeComplex pp(three_point._grid); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp, M_h); // translational factor 
    pp_initialzed = true;
  }

  if(ensemble == "Pion_32ID") {
    std::string file = three_point_32ID(traj); 
    read_cheng_PGG(three_point, file);
    three_point = imag(three_point) * pp;
  }
  else if(ensemble == "Pion_32IDF") {
    std::string file = three_point_path_32IDF(traj);
    read_luchang_PGG(three_point, file); // FIXME: change this after cheng generated his three point functions
    three_point = imag(three_point) * pp;
  }
  else if(ensemble == "Pion_48I") {
    std::string file = three_point_48I(traj);
    read_luchang_PGG(three_point, file); // FIXME: change this after cheng generated his three point functions
    three_point = imag(three_point) * pp;
  }
  else if(ensemble == "Pion_24ID") {
    std::string file = three_point_24ID(traj);
    read_cheng_PGG(three_point, file);
    three_point = imag(three_point) * pp;
  }
  else if(ensemble == "Pion_24ID_disc") {
    std::string file = three_point_disc_24ID(traj);
    readScidac(three_point, file);
    // three_point = imag(three_point) * pp;
  }
  else if(ensemble == "Pion_32ID_disc2") {
    std::string file = three_point_disc2_32ID(traj);
    readScidac(three_point, file);
    three_point = real(three_point) * pp;
  }
  else assert(0);
  if(target=="form_factor") three_point = three_point / pp;
}

std::vector<double> Jack_para::get_result_with_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic) {
  if(target=="form_factor") return form_factor(three_point, leptonic, hadron_coeff, M_h);
  else if(target == "real" || target == "real_CUBA3d" || target=="imag_analytic" || target == "imag_CUBA3d") {
    return calculate_decay_amplitude_cutoff(three_point, leptonic, lep_coeff, hadron_coeff);
  }
  else assert(0);
}








}}

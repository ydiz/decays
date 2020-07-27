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
  std::vector<std::string> targets;
  std::string file_p3;
  std::string file_p1;
  double leptonic_space_limit;
  double leptonic_time_limit;

  std::string cutoff_type;

  std::string output_prefix;

  std::vector<int> lat_size;
  std::vector<int> traj_skip; // skip these trajectories
  int traj_start, traj_end, traj_sep, traj_num;

  LatticePGG get_leptonic(GridCartesian *grid, const std::string &target);
  double get_lep_coef(const std::string &target);
  void get_three_point(LatticePGG &three_point, int traj);
  std::vector<double> get_result_with_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic, const std::string &target);
  
};

LatticePGG Jack_para::get_leptonic(GridCartesian * grid, const std::string &target) {
    LatticePGG lat(grid);
  if(target == "real_CUBA3d" || target == "imag_CUBA3d") get_leptonic_CUBA3d(file_p1, file_p3, lat, leptonic_space_limit, leptonic_time_limit);
  else if(target == "real") Grid::QCD::get_leptonic(file_p3, lat, leptonic_space_limit, leptonic_time_limit);
  else if(target == "imag_analytic") imag_part(lat, M_h);
  else if(target == "form_factor") form_factor_integrand(lat, M_h); // technically this is not leptonic part; but for convience I put it in get_leptonic function
  else assert(0);

  // For form_factor, the three point function should be <J(0) J(x)|pi >, while the get_three_point function returns <Jmu(-w/2) Jnu(w/2) |pi>; so I add an exponential factor to the leptonic part
  if(target=="form_factor") {
	  LatticeComplex pp(grid); 
    get_translational_factor(pp, M_h); // translational factor 
    lat = lat / pp; 
  }

  return lat;
}

double Jack_para::get_lep_coef(const std::string &target) {
  // me is in eV; all other masses are in lattice unit
  double lep_coef;
  double me = 511000; // Unit is eV
  if(target=="real") lep_coef = 2. / (M_PI) / 137. / 137. * me;
  else if(target=="real_CUBA3d") lep_coef = 1. / (2 * M_PI) / 137. / 137. * me;
  else if(target=="imag_analytic") {
    double Mpi = 135000000;  
    double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi)); //FIXME: change this for kaon -> mu+ mu-
    lep_coef = me / M_h * M_PI / 137. / 137. * (1. / beta * std::log((1 + beta) / (1 - beta)));
  }
  else if(target=="imag_CUBA3d") lep_coef = 1. / 2. / 137. / 137. * me;
  else if(target=="form_factor") {
    lep_coef = 2. / std::pow(M_h, 4);

    if(ensemble.substr(0,4)=="Pion") { // For pion, calculate F / F_ABJ instead of F itself
      double Fpi = 93. * (M_h / 135.); // Fpi = 93MeV. Convert Fpi to lattice unit;  // F_\pi = 0.092424
      double F_ABJ = 1. / (4. * M_PI * M_PI) / Fpi;
      lep_coef /= F_ABJ;
    }
  }
  else assert(0);

  std::cout << "Leptonic Coeff: " << lep_coef << std::endl;
  return lep_coef;
}



// always return <Jmu(-w/2) Jnu(w/2) |pi> / hadron_coeff
void Jack_para::get_three_point(LatticePGG &three_point, int traj) {

	static LatticeComplex pp(three_point.Grid()); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp, M_h); // translational factor 
    pp_initialzed = true;
  }


  if(ensemble == "Pion_24ID" || ensemble == "Pion_32ID" || ensemble == "Pion_32IDF" || ensemble == "Pion_48I" || ensemble == "Pion_48I_pqpm" || ensemble=="Pion_64I") {
    bool useCheng = false, useLuchang = true;  // use luchang's propagator // the same as (cheng's + fission) * 0.5
    // bool useCheng = true, useLuchang = false; // use cheng's propagator
    assert((useCheng && useLuchang) == false);

    if(useCheng) {
      std::string file = three_point_path(traj, ensemble, "decay_cheng");
      read_cheng_PGG(three_point, file);
      three_point = imag(three_point) * pp;
    }

    if(useLuchang) {

      std::string file_decay = three_point_path(traj, ensemble, "decay");
      if(ensemble=="Pion_64I" || ensemble=="Pion_48I_pqpm") { // for 64I and new 48I luchange did not add type1 and type2 to final "decay"

        read_luchang_PGG(three_point, file_decay+"_type_1");  // Luchang's type II contraction is the complex conjugate of type I
        three_point = 2. * real(three_point);
        if(ensemble=="Pion_64I") {
          three_point = three_point * 16.;
        std::cout << "Luchang missed a factor of 16, so I am multiplying the three point function by 16 here. Remove this once Luchange fixed this bug." << std::endl;
        }
      }
      else read_luchang_PGG(three_point, file_decay);

      LatticePGG three_point_fission(three_point.Grid());
      std::string file_fission = three_point_path(traj, ensemble, "fission");
      if(ensemble=="Pion_64I"|| ensemble=="Pion_48I_pqpm") {

        read_luchang_PGG(three_point_fission, file_fission + "_type_1");
        three_point_fission = 2. * real(three_point_fission);

        if(ensemble=="Pion_64I") {
          three_point_fission = three_point_fission * 16.;
          std::cout << "Luchang missed a factor of 16, so I am multiplying the three point function by 16 here. Remove this once Luchange fixed this bug." << std::endl;
        }
      }
      else read_luchang_PGG(three_point_fission, file_fission);

      three_point = 0.5 * (three_point + get_reflection(three_point_fission)); // average over "decay" and "fission"

      static LatticeComplex luchang_exp(three_point.Grid());
      static bool luchange_exp_initialized = false;
      if(!luchange_exp_initialized) {
        std::map<std::string, double> tmins {{"Pion_24ID", 10.}, {"Pion_32ID", 10.}, {"Pion_32IDF", 14.}, 
                                                {"Pion_48I", 16.}, {"Pion_48I_pqpm", 16.}, {"Pion_64I", 22.}};
        get_luchang_exp_factor(luchang_exp, M_h, tmins.at(ensemble)); // can be optimized; do not calculate every time
        luchange_exp_initialized = true;
      }

      three_point = three_point * luchang_exp; // multiply it by exp(Mpi * t_pi)
      three_point = real(three_point) * pp;
      // print_grid_field_site(three_point, {1,1,1,1});
      // print_grid_field_site(three_point, {20,20,20,5});
    }
  }
  else if(ensemble == "Pion_24ID_disc") {
    std::string file = three_point_disc_24ID(traj);
    readScidac(three_point, file);
  }
  else if(ensemble == "Pion_32ID_disc2") {
    std::string file = three_point_disc2_32ID(traj);
    readScidac(three_point, file);
    three_point = real(three_point) * pp;
  }  
  // ////////////////// Kaon //////////////////
  // else if(ensemble == "Kaon_24ID") { // kaon four point function should be (1) a real function, (2) H(0, x); not shift to H(-w/2, w/2)
  //   std::string file = Kaon_four_point_24ID(traj);
  //   readScidac(three_point, file);
  //   three_point = real(three_point) * pp;
  // }
  else assert(0);
}

std::vector<double> Jack_para::get_result_with_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic, const std::string &target) {
  double lep_coeff = get_lep_coef(target);
  if(target=="form_factor") return form_factor(three_point, leptonic, hadron_coeff, lep_coeff, cutoff_type);
  else if(target == "real" || target == "real_CUBA3d" || target=="imag_analytic" || target == "imag_CUBA3d") {
    return calculate_decay_amplitude_cutoff(three_point, leptonic, lep_coeff, hadron_coeff, cutoff_type);
  }
  else assert(0);
}








}}

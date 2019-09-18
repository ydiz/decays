// g++ -I/home/ydzhao/cuth/install/boost/include -L/home/ydzhao/cuth/install/boost/lib GF_para.cc -lboost_program_options
#pragma once

#include <stdlib.h>
#include <boost/program_options.hpp>
#include "jackknife.h"

namespace po = boost::program_options;

namespace Grid {
namespace QCD {


void cmdOptionIntVector(const std::string &str,std::vector<int> &vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  int i;
  while (ss >> i){
    vec.push_back(i);
    if(std::ispunct(ss.peek()))
      ss.ignore();
  }
  return;
}

void init_para(int argc, char **argv, Jack_para &para)
{
  po::options_description desc("jackknife options");
  desc.add_options()("help", "help message")
                    ("ensemble", po::value<std::string>(&para.ensemble))
                    // ("lat_size", po::value<std::string>())
                    // ("traj_start", po::value<int>(&para.traj_start)->default_value(1000))
                    // ("traj_end", po::value<int>(&para.traj_end)->default_value(1000))
                    // ("traj_sep", po::value<int>(&para.traj_sep)->default_value(10))
                    // ("traj_skip", po::value<std::string>()->default_value(""))
                    ("time_cutoff_start", po::value<int>(&para.time_cutoff_start)->default_value(1))
                    ("time_cutoff_end", po::value<int>(&para.time_cutoff_end)->default_value(16))
                    ("target", po::value<std::string>(&para.target)->default_value(""))
                    ("file_p3", po::value<std::string>(&para.file_p3)->default_value(""))
                    ("file_p1", po::value<std::string>(&para.file_p1)->default_value(""))
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options // command line options have higher priority
  po::store(po::parse_config_file<char>("jack_init.ini", desc), vm);
  po::notify(vm);

  /////////////////////////////////////////////

  // cmdOptionIntVector(vm["lat_size"].as<std::string>(), para.lat_size);
  // cmdOptionIntVector(vm["traj_skip"].as<std::string>(), para.traj_skip);

  if(para.ensemble.substr(0,4)=="Pion") {
      double me = 511000;
      double Mpi = 135000000;
      double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
      double Gamma_coeff = 2.0 * beta / (16 * M_PI * Mpi); // the first factor 2.0 comes from adding two possible polarizations
      double Gamma_photons = 7.82;
      para.BR_coeff = Gamma_coeff / Gamma_photons;
  }
  else if(para.ensemble.substr(0,4)=="Kaon") {
    //TODO
    para.BR_coeff = 0.;
  }
  else assert(0);
  assert(para.BR_coeff != 0.);


  // hadronic part
  if(para.ensemble == "Pion_32ID") {
    para.M_h = 0.139474;
    para.N_h = 52.089753;
    para.Z_V = 0.7260;

    para.lat_size = {32, 32, 32, 64};
    para.traj_start = 980;
    para.traj_end = 1370;
    para.traj_sep = 10;
    para.traj_skip = {};
  }
  else if(para.ensemble == "Pion_32ID_disc2") {
    para.M_h = 0.139474;
    para.N_h = 52.089753;
    para.Z_V = 0.7260;

    para.lat_size = {32, 32, 32, 64};
    para.traj_start = 1080;
    para.traj_end = 1370;
    para.traj_sep = 10;
    para.traj_skip = {1200, 1210};
  }
  else if(para.ensemble == "Pion_24ID") {
    para.M_h = 0.13975;
    para.N_h = 51.561594;
    para.Z_V = 0.7267;

    para.lat_size = {24, 24, 24, 64};
    para.traj_start = 2260;
    para.traj_end = 2640;
    para.traj_sep = 10;
    para.traj_skip = {2360, 2520, 2540, 2580};
  }
  else if(para.ensemble == "Pion_24ID_disc") {
    para.M_h = 0.13975;
    para.N_h = 51.561594;
    para.Z_V = 0.7267;

    para.lat_size = {24, 24, 24, 64};
    para.traj_start = 1000;
    para.traj_end = 2290;
    para.traj_sep = 10;
    para.traj_skip = {1020, 1060, 1100, 1340};
  }
  else if(para.ensemble == "Pion_32IDF") {
    para.M_h = 0.10468;
    para.N_h = 69.268015;
    para.Z_V = 0.68339;

    para.lat_size = {32, 32, 32, 64};
    para.traj_start = 270;
    para.traj_end = 430;
    para.traj_sep = 10;
    para.traj_skip = {};
  }
  else if(para.ensemble == "Pion_48I") {
    para.M_h = 0.08049;
    para.N_h = 85.866659;
    para.Z_V = 0.71076;

    para.lat_size = {48, 48, 48, 96};
    para.traj_start = 1290;
    para.traj_end = 1730;
    para.traj_sep = 20;
    para.traj_skip = {1450, 1470, 1490};
  }

  else assert(0);

  para.hadron_coeff = 1./ (3 * std::sqrt(2)) * para.Z_V * para.Z_V * 2. * para.M_h / para.N_h;

  // trajectoies
  for(int t: para.traj_skip) assert(t > para.traj_start && t < para.traj_end);
  para.traj_num = (para.traj_end - para.traj_start) / para.traj_sep + 1 - para.traj_skip.size();

  // reading leptonic part
  para.leptonic_space_limit = para.lat_size[0] / 2; // for reading leptonic part // this should match the CUBA calculation of leptonic part
  // para.leptonic_time_limit = 16; // for reading leptonic part // this should match the CUBA calculation of leptonic part
  para.leptonic_time_limit = para.lat_size[3] / 4; // for reading leptonic part // this should match the CUBA calculation of leptonic part
  para.time_cutoff_num = para.time_cutoff_end - para.time_cutoff_start + 1;
  assert(para.time_cutoff_end <= para.leptonic_time_limit); // If target requires reading integrals, larger time cutoff does not make any sense since leptonic part is 0.

  // leptonic part
  // me is in eV; all other masses are in lattice unit
  double me = 511000; // Unit is eV
  if(para.target=="real") para.lep_coeff = 2. / (M_PI) / 137. / 137. * me;
  else if(para.target=="real_CUBA3d") para.lep_coeff = 1. / (2 * M_PI) / 137. / 137. * me;
  else if(para.target=="imag_analytic") {
    double Mpi = 135000000;  
    double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi)); //FIXME: change this for kaon -> mu+ mu-
    para.lep_coeff = me / para.M_h * M_PI / 137. / 137. * (1. / beta * std::log((1 + beta) / (1 - beta)));
  }
  else if(para.target=="imag_CUBA3d") para.lep_coeff = 1. / 2. / 137. / 137. * me;
  else if(para.target=="form_factor") ;
  else assert(0);


  ///////////////////////////////////////

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

  assert(para.lat_size.size()==4);

  /////////////////////////////////////
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "ensemble: " << para.ensemble << std::endl;
  std::cout << "M_h: " << para.M_h << std::endl;
  std::cout << "N_h: " << para.N_h << std::endl;
  std::cout << "Z_V: " << para.Z_V << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "target: " << para.target << std::endl;
  std::cout << "file_p3: " << para.file_p3 << std::endl;
  std::cout << "file_p1: " << para.file_p1 << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "lat_size: " << para.lat_size << std::endl;
  std::cout << "traj_start: " << para.traj_start << std::endl;
  std::cout << "traj_end: " << para.traj_end << std::endl;
  std::cout << "traj_sep: " << para.traj_sep << std::endl;
  std::cout << "traj_skip: " << para.traj_skip << std::endl;
  std::cout << "traj_num: " << para.traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "time_cutoff_start: " << para.time_cutoff_start << std::endl;
  std::cout << "time_cutoff_end: " << para.time_cutoff_end << std::endl;
  std::cout << "time_cutoff_num: " << para.time_cutoff_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

}



}}

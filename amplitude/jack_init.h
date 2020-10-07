// g++ -I/home/ydzhao/cuth/install/boost/include -L/home/ydzhao/cuth/install/boost/lib GF_para.cc -lboost_program_options
#pragma once

#include <stdlib.h>
#include <boost/program_options.hpp>
#include "jackknife.h"
#include <map>
#include <sstream>

namespace po = boost::program_options;

namespace Grid {
namespace QCD {

template<typename T>
void cmdOptionToVector(const std::string &str,std::vector<T> &vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  T i;
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
                    ("targets", po::value<std::string>()->default_value(""))
                    ("file_p3", po::value<std::string>(&para.file_p3)->default_value(""))
                    ("file_p1", po::value<std::string>(&para.file_p1)->default_value(""))
                    ("cutoff_type", po::value<std::string>(&para.cutoff_type)->default_value("time or 4D"))
                    ("QuickTest", po::value<bool>()->default_value(false))
                    ("output_prefix", po::value<std::string>(&para.output_prefix)->default_value("."))
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options // command line options have higher priority
  po::store(po::parse_config_file<char>("jack_init.ini", desc), vm);
  po::notify(vm);

  /////////////////////////////////////////////
  cmdOptionToVector(vm["targets"].as<std::string>(), para.targets);

  // cmdOptionIntVector(vm["lat_size"].as<std::string>(), para.lat_size);
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
    double m_mu = 105658000;
    double M_K = 497611000;
    double beta = std::sqrt(1 - 4*m_mu*m_mu / (M_K*M_K));
    double Gamma_coeff = 2.0 * beta / (16 * M_PI * M_K); // the first factor 2.0 comes from adding two possible polarizations
    double Gamma_photons = 7.037567e-12;
    para.BR_coeff = Gamma_coeff / Gamma_photons;
  }
  else assert(0);
  assert(para.BR_coeff != 0.);


  // hadronic part
  if(para.ensemble == "Pion_24ID") {
    para.M_h = 0.13975;
    para.N_h = 51.561594;
    para.Z_V = 0.7267;

    para.lat_size = {24, 24, 24, 64};
    para.traj_start = 1010;
    para.traj_end = 2640;
    para.traj_sep = 10;
    para.traj_skip = {1020,1060,1100,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500,1510,1520,1530,1540,1550,1560,1570,1580,1590,1600,1610,1620,1630,1640,1650,1660,1670,1680,1690,1700,1710,1720,1730,1740,1750,1760,1770,1780,1790,1800,1810,1820,1830,1840,1850,1860,1870,1880,1890,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100,2110,2120,2130,2140,2150,2160,2170,2180,2190,2200,2210,2220,2230,2240,2250,2360,2520,2540,2580};
  }
  else if(para.ensemble == "Pion_32ID") {
    para.M_h = 0.139474;
    para.N_h = 52.089753;
    para.Z_V = 0.7260;

    para.lat_size = {32, 32, 32, 64};
    para.traj_start = 690;
    para.traj_end = 1370;
    para.traj_sep = 10;
    para.traj_skip = {770,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990};
  }
  else if(para.ensemble == "Pion_32IDF") {
    para.M_h = 0.10468;
    para.N_h = 69.268015;
    para.Z_V = 0.68339;

    para.lat_size = {32, 32, 32, 64};
    para.traj_start = 500;
    para.traj_end = 1270;
    para.traj_sep = 10;
    para.traj_skip = {650,670,690,700,750,790,800,810,1100,1110,1120,1130,1140,1220,1230,1240,1250};
  }
  else if(para.ensemble == "Pion_48I") {
    para.M_h = 0.08049;
    para.N_h = 85.866659;
    para.Z_V = 0.71076;

    para.traj_start = 990; // the old three point functions that Luchang accidentally deleted
    para.traj_end = 1850;
    para.traj_sep = 20;
    para.traj_skip = {1050,1070,1150,1170,1190,1230,1250,1270,1450,1470,1490,1810,1830};
  }
  else if(para.ensemble == "Pion_48I_pqpm") { // partially_quenched_pion_mass
    para.M_h = 0.0783812;
    para.N_h = 87.1353496;
    para.Z_V = 0.71076;

    para.lat_size = {48, 48, 48, 96};
    para.traj_start = 630;
    para.traj_end = 1250;
    para.traj_sep = 10;
    para.traj_skip = {720, 760, 800, 840, 860, 880, 900, 920, 940, 960, 980, 990, 1000, 1010, 1020, 1030, 1040, 1060, 1080,
                      1090, 1100, 1110, 1120, 1130, 1140, 1160, 1180, 1200, 1210, 1220, 1240};
  }
  else if(para.ensemble == "Pion_64I") {
    para.M_h = 0.057328; // this is different from the 2014 paper. This is the pion mass used to calculate propagator, which is different from the light quark mass in generating the ensemble which is in the 2014 paper.  
    para.N_h = 107.01; 
    para.Z_V = 0.74293;

    para.lat_size = {64, 64, 64, 128};
    para.traj_start = 1290;
    para.traj_end = 2490;
    para.traj_sep = 20;
    para.traj_skip = {1350,1550,1590,1630,1770,1890,1930,1970,2010,2050,2090,2130};
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
  ///////////////////////////// Kaon /////////////////////
  else if(para.ensemble == "Kaon_24ID") {
    para.M_h = 0.50365;
    para.N_h = 51.561594;
    para.Z_V = 0.72672;

    para.lat_size = {24, 24, 24, 64};
    // para.traj_start = 2300;
    // para.traj_end = 2300;
    para.traj_start = 2160;
    para.traj_end = 2250;
    para.traj_sep = 10;
    para.traj_skip = {};
  }
  else assert(0);

  // hadron_coeff
  if(para.ensemble.substr(0,4)=="Pion") {
    para.hadron_coeff = 1./ (3 * std::sqrt(2)) * para.Z_V * para.Z_V * 2. * para.M_h / para.N_h;
  }
  // Note, for sBar_d diagram, G_F is also already multiplied
  else if(para.ensemble.substr(0,4)=="Kaon") {
    // double Vud = 0.97446, Vus = 0.22452;
    // double G_F = 1.1663787e-5 * std::pow(para.M_h / 0.497611, -2); // G_F in lattice unit // unit of G_F is GeV^-2. Thus,  G_{F, lat} / G_F = (M_{K,lat} / M_K)^-2
    // para.hadron_coeff = (G_F * Vud * Vus/ std::sqrt(2)) * para.Z_V * para.Z_V * 2. * para.M_h / para.N_h;
    para.hadron_coeff = para.Z_V * para.Z_V * 2. * para.M_h / para.N_h;
  }
  else assert(0);

  std::cout << "Hadronic Coeff: " << para.hadron_coeff << std::endl;


  if(vm["QuickTest"].as<bool>()) {
    para.traj_end = para.traj_start;
    para.traj_skip.clear();
  }


  // trajectoies
  for(int t: para.traj_skip) assert(t > para.traj_start && t < para.traj_end);
  para.traj_num = (para.traj_end - para.traj_start) / para.traj_sep + 1 - para.traj_skip.size();

  // reading leptonic part
  para.leptonic_space_limit = para.lat_size[0] / 2; // for reading leptonic part // this should match the CUBA calculation of leptonic part
  para.leptonic_time_limit = para.lat_size[3] / 4; // for reading leptonic part // this should match the CUBA calculation of leptonic part




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
  std::cout << "targets: " << para.targets << std::endl;
  std::cout << "file_p3: " << para.file_p3 << std::endl;
  std::cout << "file_p1: " << para.file_p1 << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "cutoff_type: " << para.cutoff_type << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "output_prefix: " << para.output_prefix << std::endl;
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "lat_size: " << para.lat_size << std::endl;
  std::cout << "traj_start: " << para.traj_start << std::endl;
  std::cout << "traj_end: " << para.traj_end << std::endl;
  std::cout << "traj_sep: " << para.traj_sep << std::endl;
  std::cout << "traj_skip: " << para.traj_skip << std::endl;
  std::cout << "traj_num: " << para.traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

}



}}

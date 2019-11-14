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

void init_para(int argc, char **argv, Env &env)
{
  po::options_description desc("kaon options");
  desc.add_options()("help", "help message")
                    ("T_wall", po::value<int>(&para.T_wall)->default_value(12))
                    ("T_u", po::value<int>(&para.time_cutoff_end)->default_value(6))
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options // command line options have higher priority
  po::store(po::parse_config_file<char>("kaon_init.ini", desc), vm);
  po::notify(vm);

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

  /////////////////////////////////////
  std::cout << std::string(20, '*') << std::endl;
  std::cout << "T_wall: " << env.T_wall << std::endl;
  std::cout << "T_u: " << env.T_u << std::endl;
  std::cout << std::string(20, '*') << std::endl;
}



}}
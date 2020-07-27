#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include <headers/headers.h>


// user declared OpenMP reduction for C++ vectors of a specific type:

#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))




namespace Grid{
namespace QCD{

std::vector<double> mult_HL_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic, const std::string &cutoff_type) {

  int T = hadronic.Grid()->_fdimensions[Tdir];

  LatticeComplex tmp(hadronic.Grid());
  tmp = 0.;

  auto tmp_v = tmp.View();
  auto h_v = hadronic.View();
  auto l_v = leptonic.View();
  parallel_for(int ss=0; ss<tmp.Grid()->oSites(); ++ss){
    tmp_v[ss]()()() = h_v[ss](0, 1)()() * l_v[ss](0, 1)()() + h_v[ss](1, 0)()() * l_v[ss](1, 0)()(); 
    tmp_v[ss]()()() += h_v[ss](0, 2)()() * l_v[ss](0, 2)()() + h_v[ss](2, 0)()() * l_v[ss](2, 0)()(); 
    tmp_v[ss]()()() += h_v[ss](2, 1)()() * l_v[ss](2, 1)()() + h_v[ss](1, 2)()() * l_v[ss](1, 2)()(); 
  }


  std::vector<double> ret_real(T / 2 + 1, 0); // [0, T/2] // Finally, in jackknife.cc, we are only using [0, T/4] 

  if(cutoff_type == "time") {
    std::vector<iSinglet<Complex>> ret(T);
    sliceSum(tmp, ret, Tdir);

    ret_real[0] = ret[0]()()().real();
    ret_real[T/2] = ret[T/2]()()().real();
    for(int i=1; i<ret_real.size()-1; ++i) ret_real[i] = ret[i]()()().real() + ret[T-i]()()().real(); // add t and -t
  }
  else if(cutoff_type == "4D") {
    // parallel_for(int ss=0; ss<hadronic.Grid()->lSites(); ss++){
    // #pragma omp parallel for reduction(vec_double_plus : res)  # must use reduction to prevent data race
    for(int ss=0; ss<hadronic.Grid()->lSites(); ss++){
      Coordinate lcoor, gcoor;
      localIndexToLocalGlobalCoor(hadronic.Grid(), ss, lcoor, gcoor);

      gcoor = my_smod(gcoor, hadronic.Grid()->_fdimensions);
      double d = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2] + gcoor[3]*gcoor[3]);

      if(int(d) <= T/2) {
        typename LatticeComplex::vector_object::scalar_object m;
        peekLocalSite(m, tmp, lcoor); // tmp = hadronic * leptonic
        ret_real[int(d)] += m()()().real();
      }
    }
    hadronic.Grid()->GlobalSumVector(ret_real.data(), ret_real.size()); // Sum over Nodes
  }
  else assert(0);


  std::vector<double> ret_cumulative(T / 2 + 1);
  ret_cumulative[0] = ret_real[0];
  for(int i=1; i<ret_real.size(); ++i) ret_cumulative[i] = ret_real[i] + ret_cumulative[i-1]; // result with cutoff
  return ret_cumulative;
}



// // L_{mu nu}(w) = lepton_coeff * leptonic
// // H_{mu nu}(w) = <0| Jmu(w/2) Jnu(-w/2) |pi> = hadron_coeff * three piont function
std::vector<double> calculate_decay_amplitude_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic, double lepton_coeff, double hadron_coeff, const std::string &cutoff_type) { 
	std::vector<double> ret = mult_HL_cutoff(three_point, leptonic, cutoff_type);

  std::vector<double> amplitude_M(ret.size());
  for(int i=0; i<ret.size(); ++i) amplitude_M[i] = hadron_coeff * lepton_coeff * ret[i];

  return amplitude_M; // return real part of branching ratio
}



}}


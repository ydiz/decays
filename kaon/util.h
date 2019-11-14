#pragma once

#include <sstream>
#include <iterator>

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


using LatticeKGG = Lattice<iMatrix<iScalar<iScalar<vComplex> >, 4>>;
using LatticeLorentzColour = Lattice<iMatrix<iScalar<iMatrix<vComplex, 3> >, 4>>;

template<typename vtype>
inline iColourMatrix<vtype> traceS(const iSpinColourMatrix<vtype> &p) {
  iColourMatrix<vtype> rst;
  rst()() = p()(0,0);
  for(int mu=1; mu<4; ++mu) rst()() += p()(mu,mu);
  return rst;
}

LatticeColourMatrix traceS(const LatticePropagator &p) {
  LatticeColourMatrix rst(p.Grid());
  auto rst_v = rst.View();
  auto p_v = p.View();
  parallel_for(int ss=0; ss<p.Grid()->oSites(); ++ss)
    rst_v[ss] = traceS(p_v[ss]);
  return rst;
}

inline iMatrix<iScalar<iScalar<vComplex>>, 4> traceC(const iMatrix<iScalar<iMatrix<vComplex, 3> >, 4> &p) {
  iMatrix<iScalar<iScalar<vComplex>>, 4> rst;
  rst = Zero();
  for(int mu=0; mu<4; ++mu)
    for(int nu=0; nu<4; ++nu)
      for(int c=0; c<3; ++c) {
        rst(mu, nu)()() += p(mu, nu)()(c, c);
      }

  return rst;
}


LatticeKGG traceC(const LatticeLorentzColour &p) {
  LatticeKGG rst(p.Grid());
  auto rst_v = rst.View();
  auto p_v = p.View();
  parallel_for(int ss=0; ss<p.Grid()->oSites(); ++ss) {
    rst_v[ss] = traceC(p_v[ss]);
  }
  return rst;
}


std::string coor2CSL(const std::vector<int> &x) {
  using namespace std;
  stringstream ss;
  copy(x.begin(), x.end()-1, ostream_iterator<int>(ss, ","));
  if(!x.empty()) ss << x.back();
  return ss.str();
}


std::vector<int> CSL2coor(const std::string& s) {
  using namespace std;
  stringstream ss(s);
  vector<int> rst;
  int i;
  while(ss>>i) {
    rst.push_back(i);
    ss.ignore();
  }
  return rst;
}


int abs_distance(int t1, int t2, int T) { // |t1 - t2|, result is betwen (-T/2, T/2]
  int rst = std::abs(t1 - t2);
  if(rst>T/2) rst = T - rst;
  return rst;
}

int left_distance(int t1, int t2, int T) { // the distance from t1 to t2 if you can only move to the left 
  int tmp = t1 - t2;
  if(tmp >= 0) return tmp;
  else return T+tmp;
}

int right_distance(int t1, int t2, int T) { // the distance from t1 to t2 if you can only move to the right
  int rst = T - left_distance(t1, t2, T);
  if(rst==T) return 0;
  else return rst;
}

void zero_mask(LatticeComplex &lat, int t_wall, int right_min = 5) {
  int T = lat.Grid()->_fdimensions[3];
  std::vector<int> vec(T);
  for(int t=0; t<T; ++t) {
    if(left_distance(t_wall, t, T) < T/2 || right_distance(t_wall, t, T) < right_min) vec[t] = 0.;
    else vec[t] = 1.;
  }
  
	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    // std::vector<int> lcoor, gcoor;
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    int t = gcoor[3];
		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(vec[t], 0.);
		pokeLocalSite(m, lat, lcoor);
  }
}



// return exp(coeff * t_x), where tx is in [-T/2, T/2]
void exp_lat(LatticeComplex &lat, double coeff) {

	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    // std::vector<int> lcoor, gcoor;
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		double val;
		int xt = qlat::smod(gcoor[Tdir], lat.Grid()->_fdimensions[Tdir]);
    val = std::exp(coeff * xt); // translation factor

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
		pokeLocalSite(m, lat, lcoor);
	}
}

// return exp(M_K * (u_0 - t_x))
void expUMinusTwall(LatticeComplex &lat, int t_wall, int t_x, int T_u, double M_K) {

  int T = lat.Grid()->_fdimensions[3];
	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		typename LatticeComplex::vector_object::scalar_object m;

		double val;
		int xt = gcoor[Tdir];
    int dist_to_Hw = abs_distance(xt, t_x, T);
    // if(left_dist>T/2) m()()() = Complex(0., 0.);
    // if(left_dist>T/2 - 8) m()()() = Complex(0., 0.); // FIXME: 
    // if(left_dist>T/2 - 16) m()()() = Complex(0., 0.); // FIXME: 
    // if(left_dist>T/2 - 20) m()()() = Complex(0., 0.); // FIXME: 
    // if(left_dist>T/2 - 24) m()()() = Complex(0., 0.); // FIXME: 
    // if(left_dist>T/2 - 18) m()()() = Complex(0., 0.); // FIXME: 
    if(dist_to_Hw > T_u) m()()() = Complex(0., 0.); // FIXME: 
    else {
      int left_dist = left_distance(xt, t_wall, T);
      val = std::exp(M_K * left_dist); // translation factor
      m()()() = Complex(val, 0.);
    }

		pokeLocalSite(m, lat, lcoor);
	}
}


template<class vobj> 
void pokeLorentz(Lattice<vobj> &lhs, const Lattice<decltype(peekIndex<LorentzIndex>(vobj(),0,0))> & rhs, int i, int j)
{
  PokeIndex<LorentzIndex>(lhs,rhs,i,j);
}


}}

// template<typename T>
// void CrossCorrelation_typeI(const Lattice<T> &lat1, const Lattice<T> &lat2) {
//   FFT theFFT((GridCartesian *)lat1._grid);
//   Lattice<T> loop1_fft(loop1._grid), loop2_fft(loop1._grid);
//   LatticePGG ret(loop1._grid), ret_u(loop1._grid), ret_v(loop1._grid), ret_u_fft(loop1._grid),ret_v_fft(loop1._grid);
//   LatticeComplex t(loop1._grid);
//
//   LatticeLoop tmp = imag(loop1);
//   theFFT.FFT_all_dim(loop1_fft, tmp, FFT::forward);
//   tmp = imag(loop2);
//   theFFT.FFT_all_dim(loop2_fft, tmp, FFT::forward);
//
//   parallel_for(int ss=0; ss<ret._grid->oSites(); ss++) {
//     for(int mu=0; mu<4; ++mu)
//       for(int nu=0; nu<4; ++nu) {
//         ret_u_fft[ss]()()(mu, nu) = adj(loop1_fft[ss]()()(mu)) * loop2_fft[ss]()()(nu);
//         ret_v_fft[ss]()()(mu, nu) = adj(loop2_fft[ss]()()(mu)) * loop1_fft[ss]()()(nu);
//       }
//   }
//   theFFT.FFT_all_dim(ret_u, ret_u_fft, FFT::backward);
//   theFFT.FFT_all_dim(ret_v, ret_v_fft, FFT::backward);
//
//   get_loop_exp_minus(t, t_min, Mpi);
//   ret_u = ret_u * t;
//
//   get_loop_exp_plus(t, t_min, Mpi);
//   ret_v = ret_v * t;
//
//   ret = ret_u + ret_v;
//   const std::vector<int> &gc = ret._grid->_fdimensions;
//   double V = gc[0] * gc[1] * gc[2] * gc[3];
//   ret = ret * (- 1. / V);  // Note the -1 here. at the begining, we took imag(loop1), i.e. multiplied it by -i
//   return ret;
// }

  // FFT theFFT((GridCartesian *)loop1._grid);
  // LatticeLoop loop1_fft(loop1._grid), loop2_fft(loop1._grid);
  // LatticePGG ret(loop1._grid), ret_u(loop1._grid), ret_v(loop1._grid), ret_u_fft(loop1._grid),ret_v_fft(loop1._grid);
  // LatticeComplex t(loop1._grid);
  //
  // LatticeLoop tmp = imag(loop1);
  // theFFT.FFT_all_dim(loop1_fft, tmp, FFT::forward);
  // tmp = imag(loop2);
  // theFFT.FFT_all_dim(loop2_fft, tmp, FFT::forward);
  //
  // parallel_for(int ss=0; ss<ret._grid->oSites(); ss++) {
  //   for(int mu=0; mu<4; ++mu)
  //     for(int nu=0; nu<4; ++nu) {
  //       ret_u_fft[ss]()()(mu, nu) = adj(loop1_fft[ss]()()(mu)) * loop2_fft[ss]()()(nu);
  //       ret_v_fft[ss]()()(mu, nu) = adj(loop2_fft[ss]()()(mu)) * loop1_fft[ss]()()(nu);
  //     }
  // }
  // theFFT.FFT_all_dim(ret_u, ret_u_fft, FFT::backward);
  // theFFT.FFT_all_dim(ret_v, ret_v_fft, FFT::backward);
  //
  // get_loop_exp_minus(t, t_min, Mpi);
  // ret_u = ret_u * t;
  //
  // get_loop_exp_plus(t, t_min, Mpi);
  // ret_v = ret_v * t;
  //
  // ret = ret_u + ret_v;
  // const std::vector<int> &gc = ret._grid->_fdimensions;
  // double V = gc[0] * gc[1] * gc[2] * gc[3];
  // ret = ret * (- 1. / V);  // Note the -1 here. at the begining, we took imag(loop1), i.e. multiplied it by -i
  // return ret;
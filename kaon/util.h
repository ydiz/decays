#pragma once

namespace Grid {
namespace QCD {

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
  int T = lat._grid->_fdimensions[3];
  std::vector<int> vec(T);
  for(int t=0; t<T; ++t) {
    if(left_distance(t_wall, t, T) < T/2 || right_distance(t_wall, t, T) < right_min) vec[t] = 0.;
    else vec[t] = 1.;
  }
  
	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    int t = gcoor[3];
		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(vec[t], 0.);
		pokeLocalSite(m, lat, lcoor);
  }
}



// return exp(coeff * t_x), where tx is in [-T/2, T/2]
void exp_lat(LatticeComplex &lat, double coeff) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

		double val;
		int xt = qlat::smod(gcoor[Tdir], lat._grid->_fdimensions[Tdir]);
    val = std::exp(coeff * xt); // translation factor

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
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

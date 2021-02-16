#pragma once

#include <sstream>
#include <iterator>

namespace Grid {

template<class vobj> auto peekLorentz(const Lattice<vobj> &rhs, int i, int j) -> decltype(PeekIndex<LorentzIndex>(rhs,0,0))
{
    return PeekIndex<LorentzIndex>(rhs, i, j);
}


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



// EM factor on half lattice; if r_t < 0 or r_t > max_uv_sep, Euv = 0; If r_t == 0, do not change Euv; If 0 < r_t <= max_uv_sep, multiply Euv by 2.
// The factor 2 / M_K^4 is included.
void EM_factor_half_lattice(LatticeKGG &lat, const std::vector<int> &v, double M_K, int max_uv_sep) { // max sep is maximum allowed |u0-v0|

  lat = Zero();

  const int T = lat.Grid()->_fdimensions[3];

  autoView(lat_v, lat, CpuWrite);
  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++) {

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat.Grid()->_fdimensions);  // smod: if x > T/2, map x to x-T

    if(gcoor[Tdir] < 0 || gcoor[Tdir] > max_uv_sep) continue; // Non-zero only for 0 <= t <= max_uv_sep
    else {  

      double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

      double val;
      if(w==0) val = 0.;
      else {
        double t = M_K * 0.5 * w;
        double sin_t, cos_t;
        sincos(t, &sin_t, &cos_t);
        val = ( - cos_t * M_K * w + 2 * sin_t) / w / w / w * std::exp(M_K * 0.5 * gcoor[Tdir]);
      }

      if(gcoor[Tdir] > 0) val *= 2.;  // if t > 0, multiply it by 2.

      typename LatticeKGG::vector_object::scalar_object m;
      m = Zero();
      m(0, 1)()() = Complex(val * gcoor[Zdir], 0); 
      m(0, 2)()() = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
      m(1, 2)()() = Complex(val * gcoor[Xdir], 0); 
      m(1, 0)()() = - m(0, 1)()();
      m(2, 0)()() = - m(0, 2)()();
      m(2, 1)()() = - m(1, 2)()();

      // pokeLocalSite(m, lat, lcoor);
      pokeLocalSite(m, lat_v, lcoor);
    }
  }

  lat *= (2. / M_K / M_K / M_K / M_K);


  for(int mu=0; mu<4; ++mu) {  // shift center to v
    if(v[mu] != 0) lat = Cshift(lat, mu, -v[mu]);
  }

  // std::cout << lat << std::endl;
  // exit(0);
}


struct Sum_Interval {
  int lower_bound, upper_bound, T;

  Sum_Interval(int _lower_bound, int _upper_bound, int _T) : lower_bound(_lower_bound), upper_bound(_upper_bound), T(_T) {
    assert(lower_bound >=0 && upper_bound>=0 && upper_bound >= lower_bound); // both lower_bound and upper_bound can be greater than T
  };

  template <class LatticeType>
  typename LatticeType::vector_object::scalar_object operator()(const LatticeType &lat) const { 
    using Site = typename LatticeType::vector_object::scalar_object;
    std::vector<Site> lat_sliceSum;
    sliceSum(lat, lat_sliceSum, Tdir);

    Site rst = Zero();
    for(int t=lower_bound; t<=upper_bound; ++t) {
      rst += lat_sliceSum[t % T];
    }
    return rst;
  }
};


struct Sum_Interval_TimeSlice {
  int lower_bound, upper_bound, T;

  Sum_Interval_TimeSlice(int _lower_bound, int _upper_bound, int _T) : lower_bound(_lower_bound), upper_bound(_upper_bound), T(_T) {
    assert(lower_bound >=0 && upper_bound>=0 && upper_bound >= lower_bound); // both lower_bound and upper_bound can be greater than T
  };

  template <class LatticeType>
  std::vector<typename LatticeType::vector_object::scalar_object> operator()(const LatticeType &lat) const { 
    using Site = typename LatticeType::vector_object::scalar_object;
    std::vector<Site> lat_sliceSum;
    sliceSum(lat, lat_sliceSum, Tdir);

    std::vector<Site> rst(upper_bound - lower_bound + 1);
    for(int t=lower_bound; t<=upper_bound; ++t) {
      // rst += lat_sliceSum[t % T];
      rst[t - lower_bound] = lat_sliceSum[t % T];   // rst[0] is the time slice sum at lower found; rst.back() is at upper bound
    }
    return rst;
  }
};





template<typename vtype>
inline iColourMatrix<vtype> traceS(const iSpinColourMatrix<vtype> &p) {
  iColourMatrix<vtype> rst;
  rst()() = p()(0,0);
  for(int mu=1; mu<4; ++mu) rst()() += p()(mu,mu);
  return rst;
}

LatticeColourMatrix traceS(const LatticePropagator &p) {
  LatticeColourMatrix rst(p.Grid());
  // auto rst_v = rst.View();
  // auto p_v = p.View();
  autoView(rst_v, rst, CpuWrite);
  autoView(p_v, p, CpuRead);
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
  // auto rst_v = rst.View();
  // auto p_v = p.View();
  autoView(rst_v, rst, CpuWrite);
  autoView(p_v, p, CpuRead);
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




std::vector<int> str2coor(const std::string& s) {
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

std::string coor2str(const std::vector<int> &x) { // [1,2,3,4] -> 1_2_3_4
  using namespace std;
  stringstream ss;
  copy(x.begin(), x.end()-1, ostream_iterator<int>(ss, "_"));
  if(!x.empty()) ss << x.back();
  return ss.str();
}


// void zero_mask(LatticeComplex &lat, int t_wall, int right_min = 5) {
//   int T = lat.Grid()->_fdimensions[3];
//   std::vector<int> vec(T);
//   for(int t=0; t<T; ++t) {
//     if(left_distance(t_wall, t, T) < T/2 || right_distance(t_wall, t, T) < right_min) vec[t] = 0.;
//     else vec[t] = 1.;
//   }
//   
// 	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
//     // std::vector<int> lcoor, gcoor;
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
//
//     int t = gcoor[3];
// 		typename LatticeComplex::vector_object::scalar_object m;
// 		m()()() = Complex(vec[t], 0.);
// 		pokeLocalSite(m, lat, lcoor);
//   }
// }
//
//
//
// // return exp(coeff * t_x), where tx is in [-T/2, T/2]
// void exp_lat(LatticeComplex &lat, double coeff) {
//
// 	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){
//     // std::vector<int> lcoor, gcoor;
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
//
// 		double val;
// 		int xt = my_smod(gcoor[Tdir], lat.Grid()->_fdimensions[Tdir]);
//     val = std::exp(coeff * xt); // translation factor
//
// 		typename LatticeComplex::vector_object::scalar_object m;
// 		m()()() = Complex(val, 0.);
// 		pokeLocalSite(m, lat, lcoor);
// 	}
// }



template<class vobj> 
void pokeLorentz(Lattice<vobj> &lhs, const Lattice<decltype(peekIndex<LorentzIndex>(vobj(),0,0))> & rhs, int i, int j)
{
  PokeIndex<LorentzIndex>(lhs,rhs,i,j);
}


void conjugateU(LatticePropagator &prop) {
  // To change the P[U] to P[U^*], we use P[U^*] = (C * gamma_5) P[U]^* (C * gamma_5)^{-1} = (C * gamma_5) P[U]^* gamma_5 C^{-1}
  prop = - ((gmu[2] * gmu[4] * g5) * conjugate(prop) * (g5 * gmu[2] * gmu[4]));
}





std::vector<std::vector<int>> my_get_xgs(const std::string &path, bool allowPathNotExists=false) {  // point source prop file name has the format 1_2_3_4
  std::vector<std::vector<int>> xgs;

  DIR *dir;
  dir = opendir(path.c_str());
  std::cout << "path: " << path << std::endl;

  if(!allowPathNotExists) assert(dir!=NULL); // make sure directory exists

  if(dir != NULL) {
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {     
      std::string subdir_name = std::string(entry->d_name);
      if(!(subdir_name[0] >= '0' && subdir_name[0] <= '9')) continue;  // there can be "." and ".."
      std::vector<int> xg = str2coor(subdir_name); 
      xgs.push_back(xg);
    }
    closedir(dir);
  }

  std::cout << "Number of point sources: " << xgs.size() << std::endl;


  // readdir does not gaurantee the order of files
  // make sure the order of points are the same on all nodes
  std::sort(xgs.begin(), xgs.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  });

  return xgs;
}












}

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

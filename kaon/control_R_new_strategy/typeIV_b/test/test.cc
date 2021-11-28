#include "kaon/kaon.h"
using namespace std;
using namespace Grid;

template <class T>
Lattice<T> my_3D_forward_FFT_3D(const Lattice<T> &lat) {

  GridBase *grid = lat.Grid();
  Lattice<T> lat_fft(grid);
  FFT theFFT((GridCartesian *)grid);
  Coordinate mask({1,1,1,0});  // 3D FFT
  theFFT.FFT_dim_mask(lat_fft, lat, mask, FFT::forward);
  return lat_fft;
}

template <class T>
Lattice<T> my_3D_backward_FFT_3D(const Lattice<T> &lat) {

  GridBase *grid = lat.Grid();
  Lattice<T> lat_fft(grid);
  FFT theFFT((GridCartesian *)grid);
  Coordinate mask({1,1,1,0});  // 3D FFT
  theFFT.FFT_dim_mask(lat_fft, lat, mask, FFT::backward);
  return lat_fft;
}

// generate table3d[xt][|vec{x}|][vt]; initialize it to 0
vector<vector<vector<Complex>>> control_R_initialize_table3d(const Coordinate &fdims) { 

    const int T = fdims[3];
    int X = fdims[0], Y = fdims[1], Z = fdims[2];
    int num_R = 1 + int(calc_3d_vec_norm(X/2, Y/2, Z/2));

    vector<vector<vector<Complex>>> table3d(T);  
    for(int xt=0; xt<T; ++xt) {
      table3d[xt].resize(num_R);
      for(int R=0; R<num_R; ++R) {
        table3d[xt][R].resize(T, 0.);
      }
    }
    return table3d;
}


// table3d[xt][R][vt] = sum_{vec{v}} trace( f(xt, vec{x}=vec{v}+R) * g(vt, vec{v}) )
vector<vector<vector<Complex>>> sum_by_xt_R_vt(const LatticePropagator &fx, const LatticePropagator &gv) {
  LatticePropagator fx_fft = my_3D_forward_FFT_3D(fx);  // Fourier[f]
  LatticePropagator gv_conj = conjugate(gv);
  LatticePropagator gv_fft = conjugate(my_3D_forward_FFT_3D(gv_conj));  // Fourier[f^*]^*

  vector<vector<vector<Complex>>> table3d = control_R_initialize_table3d(fx.Grid()->_fdimensions); // table3d[xt][|vec{x}|][vt] // resize and set initial values to 0

  int num_T = table3d.size();
  for(int delta_t=0; delta_t<table3d.size(); ++delta_t) { // delta_t :  vt = xt + delta_t
    LatticePropagator tmp = fx_fft * Cshift(gv_fft, Tdir, delta_t); // f(x) * g(v = x+delta_t)
    tmp = my_3D_backward_FFT_3D(tmp);
    LatticeComplex amplitude = trace(tmp);

    int num_R = table3d[0].size();
    vector<vector<Complex>> amplitude_by_xt_R = sum_T_R(amplitude, {0,0,0,0});
    for(int xt=0; xt<num_T; ++xt) {
      int vt = (xt + delta_t) % num_T;
      for(int R=0; R<num_R; ++R) {
        table3d[xt][R][vt] = amplitude_by_xt_R[xt][R];
      }
    }
  }

  return table3d;
}


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  const int T = grid->_fdimensions[3];

  // LatticeComplex f(grid), g(grid);
  LatticePropagator f(grid), g(grid);
  f = 1.0; g = 1.0;
  vector<vector<vector<Complex>>> table3d = sum_by_xt_R_vt(f, g);

  std::cout << table3d << std::endl;


  std::cout << GridLogMessage << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}

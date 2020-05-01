#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// |t_u - t_x| and |t_v - t_x| must be <= tsep - tsep2
void typeIII_set0(LatticePropagator &lat, int t_x, int tsep, int tsep2) {
  const int T = lat.Grid()->_fdimensions[3];
  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    if(distance(t_x, gcoor[3], T) > tsep - tsep2) { // set to zero if |t_u - t_x| or |t_v - t_x| > tsep - tsep2
      LatticePropagatorSite tmp; tmp = Zero();
      pokeLocalSite(tmp, lat, lcoor);
    }; 
  });
}

// return exp(M_K * (v_0 - tK)), for any v_0 that satisfies |v_0 - tx| <= tsep - tsep2
LatticeComplex typeIII_exp(GridCartesian *grid, int t_x, int tsep, int tsep2, double M_K) {
/*
tsep, tsep2: window size. Keep only points between [x - (tsep-tsep2), x + (tsep-tsep2)]
M_K: kaon mass on the lattice
*/
  LatticeComplex lat(grid);
  int T = lat.Grid()->_fdimensions[3];

  int tK = t_x - tsep;
  if(tK < 0) tK += T;

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		LatticeComplexSite m;

    int dist = distance(gcoor[3], t_x, T); // |tv - tx|
    if(dist > tsep - tsep2) m = Zero();
    else {
      int left_dist = left_distance(gcoor[3], tK, T); // v_0 - tK // the distance from v0 to tK if you can only move to the left
      double val = std::exp(M_K * left_dist); 
      m()()() = Complex(val, 0.);
    }

		pokeLocalSite(m, lat, lcoor);
  });

  return lat;
}






std::vector<int> gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  // FIXME: change those parameters
  int tsep = 16;
  int tsep2 = 6;

  Env env(gcoor, "24ID");
  // init_para(argc, argv, env);
  env.N_pt_src = 1;  // FIXME: keep only one point

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    if(env.N_pt_src != -1) env.xgs_s.erase(env.xgs_s.begin() + env.N_pt_src, env.xgs_s.end());

    // std::vector<LatticePropagator> wl = env.get_wall('l');
    // std::vector<LatticePropagator> ws = env.get_wall('s');

    //FIXME: For test
    std::vector<LatticePropagator> wl(T, env.grid);
    std::vector<LatticePropagator> ws(T, env.grid);

    LatticeKGG rst_D3Q1_allsrc(env.grid), rst_D3Q2_allsrc(env.grid), rst_sBar_d_D3_allsrc(env.grid); 
    vector<LatticeKGG*> rst_vec_allsrc= {&rst_D3Q1_allsrc, &rst_D3Q2_allsrc, &rst_sBar_d_D3_allsrc}; 
    for(auto rst: rst_vec_allsrc) *rst = Zero();

    for(const auto &x: env.xgs_s) {

      LatticePropagator pl = env.get_point(x, 'l'); // pl = L(x, v) or L(u, v)
      LatticePropagator ps = env.get_point(x, 's'); // ps = H(x, v) or H(u, v)
      LatticePropagatorSite Lxx; // L(x, x)
      peekSite(Lxx, pl, x);

      int tK = x[3] - tsep;
      if(tK < 0) tK += T;

      //FIXME: For test
      wl[tK] = env.get_wall(tK, 'l');
      ws[tK] = env.get_wall(tK, 's');

      vector<LatticePropagator> Fux(4, env.grid); // F_mu(u, x)
      for(int mu=0; mu<4; ++mu) Fux[mu] = adj(pl) * gmu5[mu] * wl[tK];

      vector<LatticePropagator> Gvx(4, env.grid); // G_nu(v, x)
      for(int nu=0; nu<4; ++nu) Gvx[nu] = adj(ws[tK]) * gmu5[nu] * ps; 
      
      print_grid_field_site(Fux[0], {1,1,1,1});
      print_grid_field_site(Gvx[0], {1,1,1,1});


      for(int mu=0; mu<4; ++mu) typeIII_set0(Fux[mu], x[3], tsep, tsep2); // set to zero if |t_u - t_x| > tsep - tsep2
      for(int nu=0; nu<4; ++nu) typeIII_set0(Gvx[nu], x[3], tsep, tsep2);

      LatticeComplex exp_factor = typeIII_exp(env.grid, x[3], tsep, tsep2, env.M_K);
      for(int nu=0; nu<4; ++nu) Gvx[nu] = Gvx[nu] * exp_factor;           // G_nu(v, x) *= exp(M_k * (v_0 - t_K))


      vector<vector<LatticePropagator>> FuGv(4, vector<LatticePropagator>(4, env.grid)); // \sum_v F(v+r,x) G(v, x) 
      std::cout << GridLogMessage << "before fft" << std::endl;
      FFT theFFT((GridCartesian *)env.grid);
      std::vector<LatticePropagator> Fux_fft(4, env.grid), Gvx_conj_fft(4, env.grid); 
      for(int mu=0; mu<4; ++mu) {
        theFFT.FFT_all_dim(Fux_fft[mu], Fux[mu], FFT::forward);
        LatticePropagator Gvx_conj = conjugate(Gvx[mu]);
        theFFT.FFT_all_dim(Gvx_conj_fft[mu], Gvx_conj, FFT::forward); // Fourier transform G(v,x)^*
      }
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          LatticePropagator tmp_fft = Fux_fft[mu] * conjugate(Gvx_conj_fft[nu]);
          theFFT.FFT_all_dim(FuGv[mu][nu], tmp_fft, FFT::backward);
        }
      }
      std::cout << GridLogMessage << "after fft" << std::endl;

      for(int mu=0; mu<4; ++mu) 
        for(int nu=0; nu<4; ++nu) 
          FuGv[mu][nu] = FuGv[mu][nu] - adj(FuGv[mu][nu]);  // !!! add contribution from Kbar and Qbar

      LatticeKGG rst_D3Q1(env.grid), rst_D3Q2(env.grid), rst_sBar_d_D3(env.grid); 
      vector<LatticeKGG*> rst_vec = {&rst_D3Q1, &rst_D3Q2, &rst_sBar_d_D3};
      for(auto rst: rst_vec) *rst = Zero();

      // rst_D3_Q1
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          LatticeComplex rst_mu_nu(env.grid);
          rst_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) {
            rst_mu_nu += trace(gL[rho] * Lxx) * trace(gL[rho] * FuGv[mu][nu]);
          }
          pokeLorentz(rst_D3Q1, rst_mu_nu, mu, nu);
        }
      }

      // rst_D3_Q2
      LatticePropagatorSite tmp;
      tmp = Zero();
      for(int rho=0; rho<4; ++rho) tmp += gL[rho] * Lxx * gL[rho];
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          LatticeComplex rst_mu_nu(env.grid);
          rst_mu_nu = trace(tmp * FuGv[mu][nu]);
          pokeLorentz(rst_D3Q2, rst_mu_nu, mu, nu);
        }
      }

      // rst_sBar_d_D3
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          LatticeComplex rst_mu_nu(env.grid);
          rst_mu_nu = trace(g5 * FuGv[mu][nu]);
          pokeLorentz(rst_sBar_d_D3, rst_mu_nu, mu, nu);
        }
      }

      for(int i=0; i<rst_vec.size(); ++i)  *rst_vec_allsrc[i] += *rst_vec[i];


    } // end of point source loop

    for(auto rst: rst_vec_allsrc) *rst = *rst * (1. / double(env.N_pt_src));

    writeScidac(rst_D3Q1_allsrc, env.out_prefix + "/typeIII/D3Q1." + to_string(traj));
    writeScidac(rst_D3Q2_allsrc, env.out_prefix + "/typeIII/D3Q2." + to_string(traj));
    writeScidac(rst_sBar_d_D3_allsrc, env.out_prefix + "/JJ_sBar_d_K/D3." + to_string(traj));

  } // end of traj loop


  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}

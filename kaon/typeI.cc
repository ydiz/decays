
// On 8 nodes, need ~9h for one trajectory (500 point sources)

#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// std::vector<int> gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

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

  Env env("24ID");
  // init_para(argc, argv, env);
  // env.N_pt_src = 1;  // FIXME: keep only one point
  env.N_pt_src = -1;  // Use all points

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);


    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    LatticeKGG rst_D1Q1_allsrc(env.grid), rst_D1Q2_allsrc(env.grid), rst_D2Q1_allsrc(env.grid), rst_D2Q2_allsrc(env.grid); 
    vector<LatticeKGG*> rst_vec_allsrc= {&rst_D1Q1_allsrc, &rst_D1Q2_allsrc, &rst_D2Q1_allsrc, &rst_D2Q2_allsrc}; 
    for(auto rst: rst_vec_allsrc) *rst = Zero();

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_s.resize(env.N_pt_src);
    for(const auto &x: env.xgs_s) {
      ++num_pt_src;

      LatticePropagator pl = env.get_point(x, 'l'); // pl = L(x, v) or L(u, v)
      LatticePropagator ps = env.get_point(x, 's'); // ps = H(x, v) or H(u, v)

      int tK = x[3] - tsep;
      if(tK < 0) tK += T;

      LatticePropagatorSite wl_x, ws_x; // L(x, tK), H(x, tK)
      peekSite(wl_x, wl[tK], x);
      peekSite(ws_x, ws[tK], x);

      vector<vector<LatticeColourMatrix>> Fux(4, vector<LatticeColourMatrix>(4, env.grid)); // F_{mu,rho}(u, x)
      for(int mu=0; mu<4; ++mu) {
        LatticePropagator tmp = adj(pl) * gmu5[mu] * pl;
        for(int rho=0; rho<4; ++rho) {
          Fux[mu][rho] = traceS(gL[rho] * tmp); 
        }
      }

      vector<vector<LatticeColourMatrix>> Gvx_D1(4, vector<LatticeColourMatrix>(4, env.grid)); // G_{nu,rho}(v, x)
      vector<vector<LatticeColourMatrix>> Gvx_D2(4, vector<LatticeColourMatrix>(4, env.grid)); 
      for(int nu=0; nu<4; ++nu) {
        LatticePropagator tmp_D1 = adj(pl) * gmu5[nu] * wl[tK] * adj(ws_x);
        LatticePropagator tmp_D2 = wl_x * adj(ws[tK]) * gmu5[nu] * ps;
        tmp_D1 = tmp_D1 - adj(tmp_D1);              //  !!! add contribution from Kbar and Qbar
        tmp_D2 = tmp_D2 - adj(tmp_D2);              //  !!! add contribution from Kbar and Qbar
        for(int rho=0; rho<4; ++rho) {
          Gvx_D1[nu][rho] = traceS(gL[rho] * tmp_D1); 
          Gvx_D2[nu][rho] = traceS(gL[rho] * tmp_D2); 
        }
      }

      LatticeComplex exp_factor = convolution_exp(env.grid, x[3], tsep, tsep2, env.M_K);
      for(int mu=0; mu<4; ++mu) {
        for(int rho=0; rho<4; ++rho) {
          convolution_set0(Fux[mu][rho], x[3], tsep, tsep2); // set to zero if |t_u - t_x| > tsep - tsep2
          convolution_set0(Gvx_D1[mu][rho], x[3], tsep, tsep2);
          convolution_set0(Gvx_D2[mu][rho], x[3], tsep, tsep2);

          Gvx_D1[mu][rho] = Gvx_D1[mu][rho] * exp_factor;  // G(v, x) *= exp(M_k * (v_0 - t_K))
          Gvx_D2[mu][rho] = Gvx_D2[mu][rho] * exp_factor;
        }
      }

      LatticeKGG rst_D1Q1(env.grid), rst_D1Q2(env.grid), rst_D2Q1(env.grid), rst_D2Q2(env.grid); 
      vector<LatticeKGG*> rst_vec = {&rst_D1Q1, &rst_D1Q2, &rst_D2Q1, &rst_D2Q2};
      for(auto rst: rst_vec) *rst = Zero();

      std::cout << GridLogMessage << "before fft" << std::endl;
      FFT theFFT((GridCartesian *)env.grid);
      vector<vector<LatticeColourMatrix>> Fux_fft(4, vector<LatticeColourMatrix>(4, env.grid)); 
      vector<vector<LatticeColourMatrix>> Gvx_D1_cj_fft_cj(4, vector<LatticeColourMatrix>(4, env.grid)); //cj: complex conjugate
      vector<vector<LatticeColourMatrix>> Gvx_D2_cj_fft_cj(4, vector<LatticeColourMatrix>(4, env.grid)); 

      for(int mu=0; mu<4; ++mu) {
        for(int rho=0; rho<4; ++rho) {
          theFFT.FFT_all_dim(Fux_fft[mu][rho], Fux[mu][rho], FFT::forward);

          LatticeColourMatrix tmp_D1 = conjugate(Gvx_D1[mu][rho]);
          theFFT.FFT_all_dim(tmp_D1, tmp_D1, FFT::forward); 
          Gvx_D1_cj_fft_cj[mu][rho] = conjugate(tmp_D1);  //  Fourier{G(v,x)^*}^*

          LatticeColourMatrix tmp_D2 = conjugate(Gvx_D2[mu][rho]);
          theFFT.FFT_all_dim(tmp_D2, tmp_D2, FFT::forward); 
          Gvx_D2_cj_fft_cj[mu][rho] = conjugate(tmp_D2);
        }
      }

      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          // D1 Q1
          LatticeComplex rst_D1Q1_mu_nu(env.grid); rst_D1Q1_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) rst_D1Q1_mu_nu += trace(Fux_fft[mu][rho]) * trace(Gvx_D1_cj_fft_cj[nu][rho]);
          pokeLorentz(rst_D1Q1, rst_D1Q1_mu_nu, mu, nu);

          // D1 Q2
          LatticeComplex rst_D1Q2_mu_nu(env.grid); rst_D1Q2_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) rst_D1Q2_mu_nu += trace(Fux_fft[mu][rho] * Gvx_D1_cj_fft_cj[nu][rho]);
          pokeLorentz(rst_D1Q2, rst_D1Q2_mu_nu, mu, nu);

          // D2 Q1
          LatticeComplex rst_D2Q1_mu_nu(env.grid); rst_D2Q1_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) rst_D2Q1_mu_nu += trace(Fux_fft[mu][rho]) * trace(Gvx_D2_cj_fft_cj[nu][rho]);
          pokeLorentz(rst_D2Q1, rst_D2Q1_mu_nu, mu, nu);

          // D2 Q2
          LatticeComplex rst_D2Q2_mu_nu(env.grid); rst_D2Q2_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) rst_D2Q2_mu_nu += trace(Fux_fft[mu][rho] * Gvx_D2_cj_fft_cj[nu][rho]);
          pokeLorentz(rst_D2Q2, rst_D2Q2_mu_nu, mu, nu);

        }
      }
      theFFT.FFT_all_dim(rst_D1Q1, rst_D1Q1, FFT::backward);
      theFFT.FFT_all_dim(rst_D1Q2, rst_D1Q2, FFT::backward);
      theFFT.FFT_all_dim(rst_D2Q1, rst_D2Q1, FFT::backward);
      theFFT.FFT_all_dim(rst_D2Q2, rst_D2Q2, FFT::backward);
      std::cout << GridLogMessage << "after fft" << std::endl;

      for(int i=0; i<rst_vec.size(); ++i)  *rst_vec_allsrc[i] += *rst_vec[i];

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    for(auto rst: rst_vec_allsrc) *rst = *rst * (1. / double(num_pt_src));  

    writeScidac(rst_D1Q1_allsrc, env.out_prefix + "/typeI/D1Q1." + to_string(traj));
    writeScidac(rst_D1Q2_allsrc, env.out_prefix + "/typeI/D1Q2." + to_string(traj));
    writeScidac(rst_D2Q1_allsrc, env.out_prefix + "/typeI/D2Q1." + to_string(traj));
    writeScidac(rst_D2Q2_allsrc, env.out_prefix + "/typeI/D2Q2." + to_string(traj));

  } // end of traj loop


  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}


// On 8 nodes, Needs 10h for one trajectory  (traj 2300) 

#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  std::cout << "To run get_reflection function, must use 1 process per node; do not know why" << std::endl;

  Grid_init(&argc, &argv);

  int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  // int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
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
  // env.N_pt_src = 1;  // FIXME: keep only one point
  env.N_pt_src = -1;  // Use all points

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    LatticeKGG rst_D3Q1_allsrc(env.grid), rst_D3Q2_allsrc(env.grid), rst_sBar_d_D3_allsrc(env.grid); 
    vector<LatticeKGG*> rst_vec_allsrc= {&rst_D3Q1_allsrc, &rst_D3Q2_allsrc, &rst_sBar_d_D3_allsrc}; 
    for(auto rst: rst_vec_allsrc) *rst = Zero();

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_s.resize(env.N_pt_src);
    for(const auto &x: env.xgs_s) {
      ++num_pt_src;

      LatticePropagator pl = env.get_point(x, 'l'); // pl = L(x, v) or L(u, v)
      LatticePropagator ps = env.get_point(x, 's'); // ps = H(x, v) or H(u, v)
      LatticePropagatorSite Lxx; // L(x, x)
      peekSite(Lxx, pl, x);

      int tK = x[3] - tsep;
      if(tK < 0) tK += T;

      vector<LatticePropagator> Fux(4, env.grid); // F_mu(u, x)
      for(int mu=0; mu<4; ++mu) Fux[mu] = adj(pl) * gmu5[mu] * wl[tK];

      vector<LatticePropagator> Gvx(4, env.grid); // G_nu(v, x)
      for(int nu=0; nu<4; ++nu) Gvx[nu] = adj(ws[tK]) * gmu5[nu] * ps; 

      for(int mu=0; mu<4; ++mu) convolution_set0(Fux[mu], x[3], tsep, tsep2); // set to zero if |t_u - t_x| > tsep - tsep2
      for(int nu=0; nu<4; ++nu) convolution_set0(Gvx[nu], x[3], tsep, tsep2);

      LatticeComplex exp_factor = convolution_exp(env.grid, x[3], tsep, tsep2, env.M_K);
      for(int nu=0; nu<4; ++nu) Gvx[nu] = Gvx[nu] * exp_factor;           // G_nu(v, x) *= exp(M_k * (v_0 - t_K))


      LatticeKGG rst_D3Q1(env.grid), rst_D3Q2(env.grid), rst_sBar_d_D3(env.grid); 
      vector<LatticeKGG*> rst_vec = {&rst_D3Q1, &rst_D3Q2, &rst_sBar_d_D3};
      for(auto rst: rst_vec) *rst = Zero();

      std::cout << GridLogMessage << "before fft" << std::endl;
      FFT theFFT((GridCartesian *)env.grid);
      std::vector<LatticePropagator> Fux_fft(4, env.grid), Gvx_conj_fft(4, env.grid); 

      for(int mu=0; mu<4; ++mu) {
        theFFT.FFT_all_dim(Fux_fft[mu], Fux[mu], FFT::forward);
        LatticePropagator Gvx_conj = conjugate(Gvx[mu]);
        theFFT.FFT_all_dim(Gvx_conj_fft[mu], Gvx_conj, FFT::forward); // Fourier transform G(v,x)^*
      }
      std::cout << GridLogMessage << "After forward FFT" << std::endl;

      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          std::cout << GridLogMessage << mu << " " << nu << std::endl;
          LatticePropagator FuGv = Fux_fft[mu] * conjugate(Gvx_conj_fft[nu]);

          // sBar d
          LatticeComplex rst_sBar_d_mu_nu(env.grid);
          rst_sBar_d_mu_nu = trace(g5 * FuGv);
          pokeLorentz(rst_sBar_d_D3, rst_sBar_d_mu_nu, mu, nu);

          // Must not include Kbar in the sBar_d diagram
          LatticePropagator FuGv_withKbar(env.grid);
          FuGv_withKbar = FuGv - adj(get_reflection(FuGv));  // !!! add contribution from Kbar and Qbar

          // Q1
          LatticeComplex rst_Q1_mu_nu(env.grid);
          rst_Q1_mu_nu = Zero();
          for(int rho=0; rho<4; ++rho) rst_Q1_mu_nu += trace(gL[rho] * Lxx) * trace(gL[rho] * FuGv_withKbar);
          pokeLorentz(rst_D3Q1, rst_Q1_mu_nu, mu, nu);

          // Q2
          LatticePropagatorSite gL_Lxx_gL; gL_Lxx_gL = Zero();
          for(int rho=0; rho<4; ++rho) gL_Lxx_gL += gL[rho] * Lxx * gL[rho];

          LatticeComplex rst_Q2_mu_nu(env.grid);
          rst_Q2_mu_nu = trace(gL_Lxx_gL * FuGv_withKbar);
          pokeLorentz(rst_D3Q2, rst_Q2_mu_nu, mu, nu);

        }
      }
      theFFT.FFT_all_dim(rst_sBar_d_D3, rst_sBar_d_D3, FFT::backward);
      theFFT.FFT_all_dim(rst_D3Q1, rst_D3Q1, FFT::backward);
      theFFT.FFT_all_dim(rst_D3Q2, rst_D3Q2, FFT::backward);
      std::cout << GridLogMessage << "After backward fft" << std::endl;

      for(int i=0; i<rst_vec.size(); ++i)  *rst_vec_allsrc[i] += *rst_vec[i];

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    for(auto rst: rst_vec_allsrc) *rst = *rst * (1. / double(num_pt_src));

    writeScidac(rst_D3Q1_allsrc, env.out_prefix + "/typeIII/D3Q1." + to_string(traj));
    writeScidac(rst_D3Q2_allsrc, env.out_prefix + "/typeIII/D3Q2." + to_string(traj));
    writeScidac(rst_sBar_d_D3_allsrc, env.out_prefix + "/JJ_sBar_d_K/D3." + to_string(traj));

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}

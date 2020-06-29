#pragma once


#include "kaon.h"

namespace Grid {
namespace QCD {

LatticeKGG typeI(Env &env, int t_min) {

  LatticeKGG rst_q1_avg(env.grid), rst_q2_avg(env.grid);
  rst_q1_avg = zero; rst_q2_avg = zero;

  std::vector<LatticePropagator> wl = env.get_wall_l();
  std::vector<LatticePropagator> ws = env.get_wall_s();
  // std::vector<LatticePropagator> wl(64, env.grid);
  // std::vector<LatticePropagator> ws(64, env.grid);

  LatticeComplex exp_factor(env.grid); // for calculating cross correlation
  exp_lat(exp_factor, env.M_K);

  env.xgs_s.erase(env.xgs_s.begin()+1, env.xgs_s.end()); // FIXME: keep only one point
  for(const auto &x: env.xgs_s) {
    std::cout << "point src x: " << x << std::endl;

    // the first loop J-Hw
    LatticePropagator pl = env.get_point_l(x); // pl = L(u, x)

    LatticeLorentzColour loop1(env.grid);
    for(int mu=0; mu<4; ++mu)
      for(int rho=0; rho<4; ++rho) {
        LatticeColourMatrix tmp(env.grid);
        tmp = traceS(gL[rho] * adj(pl) * gmu5[mu] * pl);
        pokeLorentz(loop1, tmp, mu, rho);
      }


    print_grid_field_site(loop1, {1,1,1,1});

    // the second loop J-Hw-K
    LatticeLorentzColour loop2_1(env.grid), loop2_2(env.grid);
    LatticePropagator ps = env.get_point_s(x);

    int T = env.grid->_fdimensions[3];
    int t_wall = x[3] - t_min;
    if(t_wall < 0) t_wall += T; // if wall is on the left side of the current

    typename LatticePropagator::vector_object::scalar_object wall_to_x_l, wall_to_x_s;
    peekSite(wall_to_x_l, wl[t_wall], x);
    peekSite(wall_to_x_s, ws[t_wall], x);

    for(int rho=0; rho<4; ++rho)
      for(int nu=0; nu<4; ++nu) {
        LatticePropagator tmp(env.grid);
        LatticeColourMatrix tmp2(env.grid);

        tmp = adj(pl) * gmu5[nu] * wl[t_wall] * adj(wall_to_x_s);  // can be optimized by changing order and storing adj
        tmp = tmp - adj(tmp); // Q1_bar

    print_grid_field_site(tmp, {1,1,1,1});

        tmp2 = traceS(gL[rho] * tmp);
        pokeLorentz(loop2_1, tmp2, rho, nu);

        tmp = wall_to_x_l * adj(ws[t_wall]) * gmu5[nu] * ps;
        tmp = tmp - adj(tmp); // Q1_bar
    print_grid_field_site(tmp, {1,1,1,1});
    // print_grid_field_site(tmp, {1,2,3,4});
        tmp2 = traceS(gL[rho] * tmp); 
        pokeLorentz(loop2_2, tmp2, rho, nu);
      }

    LatticeLorentzColour loop2(env.grid);
    loop2 = loop2_1 - loop2_2;


    // cross correlation
    loop1 = loop1 * exp_factor;

    // before doing cross correlation, set sites on the left of the wall and close to the wall to be zero
    LatticeComplex mask(env.grid);
    zero_mask(mask, t_wall);
    loop1 = loop1 * mask;
    loop2 = loop2 * mask;


    std::cout << GridLogMessage << "before fft" << std::endl;
    FFT theFFT((GridCartesian *)loop1._grid);
    LatticeLorentzColour loop1_fft(loop1._grid), loop2_fft(loop1._grid);
    theFFT.FFT_all_dim(loop1_fft, loop1, FFT::forward);
    theFFT.FFT_all_dim(loop2_fft, loop2, FFT::forward);

    LatticeKGG rst_q1_fft(env.grid), rst_q2_fft(env.grid);
    rst_q1_fft = traceC(conjugate(loop1_fft)) * traceC(loop2_fft);
    rst_q2_fft = traceC(conjugate(loop1_fft) * loop2_fft);
    //   rst_q1 = traceC(loop1) * traceC(loop2);
    //   rst_q2 = traceC(loop1 * loop2);

    LatticeKGG rst_q1(env.grid), rst_q2(env.grid);
    theFFT.FFT_all_dim(rst_q1, rst_q1_fft, FFT::backward);
    theFFT.FFT_all_dim(rst_q2, rst_q2_fft, FFT::backward);
    // do not need to divide it by V; cancel by the V factor coming from summation over x (position of Hw)

    std::cout << GridLogMessage << "after fft" << std::endl;

    print_grid_field_site(rst_q1, {1,1,1,1});
    print_grid_field_site(rst_q2, {1,2,3,4});

    rst_q1 = rst_q1 * std::exp(- env.M_K * t_wall);
    rst_q2 = rst_q2 * std::exp(- env.M_K * t_wall);
    rst_q1_avg += rst_q1;
    rst_q2_avg += rst_q2;
  }


  LatticeKGG rst(env.grid);
  rst = (env.wilson_c1 * rst_q1_avg + env.wilson_c2 * rst_q2_avg) * (1. / (double) env.xgs_s.size());
  rst = rst * (- 4. / 9.); // Do not forget the coefficient

  std::cout << GridLogMessage << std::endl;
  return rst;
}


}}

#include <math.h>
#include <gsl/gsl_integration.h>
#include <cuba.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>


#include "cuba_wrapper.h"

// w_x from -L_LIMIT to +L_LIMIT
// w_0 from -T_LIMIT to +T_LIMIT
// // For now, my L_LIMIT is always L / 2, T_LIMIT is T / 4

// ================ ensemble parameters ===================

// const std::string ensemble = "64I";
const std::string ensemble = "48I_pqdm";


const std::string target = "Pion";

// ================ integration parameters ===================
const double lower = 0.00001;

// upper should depend on T, the lattice sites in time direction
const double upper = 10;    // For 64I, Differences of using 10 and 20 are negligible, except for very large w and w0, where the difference can be ~0.5%; but the hadronic part is supposed to be small and suppress the leptonic part. so I will choose 10
// const double upper = 20; 
const double epsrel = 1e-8;  // The difference between 10e-8 and 10e-6 is small; For 64I, at site w=1 w0=0 where the difference is largest, the difference of final result is ~0.3%.
// const double epsrel = 1e-6; 

// const double upper = 50; // wrong!!! Too large; When using Cuhre, result would be wrong for large w0
// const double epsrel = 1e-4; // wrong!!! Too small!!! Especially for w0=0, where the results depends to a large extent on the cancellation of three integrals 








//================= code ===============================


class Ensemble_info {
public:
  double M_pi;
  std::vector<int> fcoor;

  Ensemble_info(double _M_pi, const std::vector<int> &_fcoor) : M_pi(_M_pi), fcoor(_fcoor) {} 
};


std::map<std::string, Ensemble_info> ensemble_info { {"24ID", Ensemble_info(0.13975, std::vector<int>{24, 24, 24, 64})},
                                                {"32ID", Ensemble_info(0.139474, {32, 32, 32, 64})},
                                                {"32IDF", Ensemble_info(0.10468, {32, 32, 32, 64})},
                                                {"48I", Ensemble_info(0.08049, {48, 48, 48, 96})},
                                                {"48I_pqdm", Ensemble_info(0.0783812, {48, 48, 48, 96})},
                                                {"64I", Ensemble_info(0.057328, {64, 64, 64, 128})}
};


double get_beta(const std::string &tar) {
  const double me = 0.511;
  const double Mpi = 135;
  const double m_mu = 105.658;
  const double M_KL = 497.611;
  if(tar=="Pion") return std::sqrt(1 - 4*me*me / (Mpi*Mpi));
  else if(tar=="Kaon") return std::sqrt(1 - 4*m_mu*m_mu / (M_KL*M_KL));
  else assert(0);
}

// p.s. integrate_qawc will add a term (p - pole) to the denominator; DO NOT write this factor in the function
double f1(double p, void *params) {
  std::vector<double> paras = *(std::vector<double> *)params;
  double w = paras[0];
  double w0 = paras[1];
  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  return - 0.5 * std::exp(-p * w0) * (cos_pw - sin_pw / (p * w)); // f should not contain the pole term if using qaws
}

double f2(double p, void *params) {
  std::vector<double> paras = *(std::vector<double> *)params;
  double w = paras[0];
  double w0 = paras[1];
  double M_pi = paras[2];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  return std::exp(-p * w0) / (M_pi + 2 * p) * (cos_pw - sin_pw / (p * w));
}

double f3(double p, void *params) {
  std::vector<double> paras = *(std::vector<double> *)params;
  double w = paras[0];
  double w0 = paras[1];
  double M_pi = paras[2];

	double E_pe = std::sqrt(p*p + (M_pi/2.) * (M_pi/2.) - 2 * p * pe * c);

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  return  std::exp(- E_pe * w0) / E_pe * (cos_pw - sin_pw / (p * w));
}


template<class T>
void integrate_qags(T f, std::vector<double> params, double lower, double upper, double epsabs, double epsrel, gsl_integration_workspace *w, double &result, double &error) {
  gsl_function F;
  F.function = f;
  F.params = &params;
  gsl_integration_qags(&F, lower, upper, epsabs, epsrel, 1000, w, &result, &error);
}

// calculate principal part integral
template<class T>
void integrate_qawc(T f, std::vector<double> params, double lower, double upper, double pole, double epsabs, double epsrel, gsl_integration_workspace *w, double &result, double &error) {
  gsl_function F;
  F.function = f;
  F.params = &params;
  gsl_integration_qawc(&F, lower, upper, pole, epsabs, epsrel, 1000, w, &result, &error);
}

struct Func {

  double M_pi; 
	double pe; // magnitude of electron's momentum. 
  double p_interval; // integration interval
  std::vector<double> paras; // [w, w0]

  Func(double beta, double M_pi);
  std::vector<double> operator()(const std::vector<double>& v) const;
};


Func::Func(double beta, double _M_pi) {
  M_pi = _M_pi;
  // pe = 0.5 * M_PION * beta;
  pe = 0.5 * _M_pi * beta;
}

std::vector<double> Func::operator()(const std::vector<double>& v) const{

  std::vector<double> ans(1);

	double p = v[0] * p_interval; // p: [0, p_interval]
	double c = 2. * v[1] - 1.; // cos: [-1, 1]
	double Epe = std::sqrt(p*p + (M_pi/2.) * (M_pi/2.) - 2 * p * pe * c);

  double w = paras[0];
  double w0 = paras[1];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  ans[0] =  std::exp(-Epe * w0) / (Epe * (- M_pi * M_pi + 4. * pe * pe * c * c)) * (cos_pw - sin_pw / (p * w));

  ans[0] *= 2. * p_interval;

  return ans;
}

// template<class T>
// double integrate_CUBA(T func, double eps_rel, double &estimated_err,int &fail)
// {
//   std::vector<double> integral, error, prob;
//   int nregions, neval;
//
//   // int flag = 3; // verbose level: 0-3;
//   int flag = 0; // verbose level: 0-3;
//   int ndim = 2; // dimension of integral, i.e. dimension of x
//   int ncomp = 1; // dimension of return value (value of the integral), i.e. dimension of y
//
//
//   // Cuhre is much faster than Divonne, and scales better w.r.t. eps_rel; but the estimated error in Cuhre may not be reliable
// 	integrateCuhre(integral, error, prob, nregions, neval, fail, ndim, ncomp, func, 0., eps_rel, flag);
// 	// integrateDivonne(integral, error, prob, nregions, neval, fail, ndim, ncomp, func, 0., eps_rel, flag); // Got segmentation fault; do not know why
//
//   estimated_err = error[0];
//   return integral[0];
// }


gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

int main (void)
{
  double M_pi = ensemble_info.at(ensemble).M_pi;
  std::vector<int> fcoor = ensemble_info.at(ensemble).fcoor;

  int w_max = int(fcoor[0] / 2 * std::sqrt(3)) + 1;
  int w0_max = fcoor[3] / 4;

  double beta = get_beta(target);
  double log_beta = std::log((1 + beta) / (1 - beta));

  std::cout << "M_H (in lattice unit): " << M_pi << std::endl;
  std::cout << "w: [" << 0 << ", " << w_max << "]" << std::endl;
  std::cout << "w0: [" << 0 << ", " << w0_max << "]" << std::endl;
  std::cout << "beta: " << beta << std::endl;
  std::cout << std::string(30, '*') << std::endl;
  std::cout << "Upper limit of integral: " << upper << std::endl;
  std::cout << "Allowed relative error: " << epsrel << std::endl;
  std::cout << std::string(30, '*') << std::endl;
  
  Func f3(beta, M_pi);
  f3.p_interval = upper;


  // /////////////////// test////////////////////
  // int w = 1, w0 = 21;
  // double ret3;
  // int fail;
  // f3.paras = {double(w), double(w0)};
  //
  // // for(double epsrel: {1e-4}) {
  // // for(double epsrel: {1e-4, 1e-5, 1e-6, 1e-7}) {
  // for(double epsrel: {1e-4, 1e-5}) {
  //   ret3 = integrate_CUBA(f3, epsrel, fail);
  //   std::cout << "epsrel: " << epsrel << " \t ret3: " << ret3 << std::endl;
  // }
  // // std::cout << "CUBA fail: " << fail << std::endl;
  // assert(fail==0);
  //
  // assert(0);
  // /////////////////// end of test////////////////////


  for(int w=0; w<=w_max; ++w) 
    for(int w0=0; w0<=w0_max; ++w0) {
      std::cout << "w=" << w << " w0=" << w0 << std::endl;

      double ret, ret1, ret2, ret3, est_err1, est_err2, est_err3;
      if(w==0) ret =0.;
      else {
        std::vector<double> params1 = {double(w), double(w0)};
        integrate_qawc(f1, params1, lower, upper, M_pi/2., 0, epsrel, workspace, ret1, est_err1);

        std::vector<double> params2 = {double(w), double(w0), M_pi};
        integrate_qags(f2, params2, lower, upper, 0, epsrel, workspace, ret2, est_err2);

        std::vector<double> params3 = {double(w), double(w0), M_pi};
        integrate_qags(f3, params3, lower, upper, 0, epsrel, workspace, ret3, est_err3);
        // int fail;
        // f3.paras = {double(w), double(w0)};
        // ret3 = integrate_CUBA(f3, epsrel, est_err3, fail);
        // assert(fail==0);

        double pe = 0.5 * _M_pi * beta; // magnitude of momentum of emitting electrons
        double common_coef =  1/ w / w * M_PI / (pe * M_pi) * log_beta;
        double coef1 = - std::exp(0.5 * M_pi * w0) * common_coef;
        double coef2 = std::exp(-0.5 * M_pi * w0) * common_coef;
        double coef3 = - common_coef;

        // std::cout << "c1:\t" << coef1 << "\t I1: " << ret1 << "\t estimated error: " << est_err1 << std::endl;
        // std::cout << "c2:\t" << coef2 << "\t I2: " << ret2 << "\t estimated error: " << est_err2 << std::endl;
        // std::cout << "c3:\t" << coef3 << "\t I3: " << ret3 << "\t estimated error: " << est_err3 << std::endl;

        ret = coef1 * ret1 + coef2 * ret2 + coef3 * ret3;
      }
      std::cout << std::setprecision(10) << "integral = " << ret  << std::endl;
  }

  gsl_integration_workspace_free(workspace);

  return 0;
}

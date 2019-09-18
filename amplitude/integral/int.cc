#include <math.h>
#include <gsl/gsl_integration.h>
#include <cuba.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "cuba_wrapper.h"

// // For different ensembles, only need to change M_PION, L_LIMIT and target
// w_x from -L_LIMIT to +L_LIMIT
// w_0 from -T_LIMIT to +T_LIMIT
// target determines the value of beta

// // For now, my L_LIMIT is always L / 2;

// ================= physical parameters ========================

// #define M_PION 0.13975 // 24ID
// #define M_PION 0.139474 // 32ID
// #define M_PION 0.10468 // 32IDF
#define M_PION 0.08049 // 48I


// #define L_LIMIT 12
// #define L_LIMIT 16  
#define L_LIMIT 24 

// #define T_LIMIT 16
#define T_LIMIT 24

const std::string target("Pion");

// ================ integration parameters ===================
const double upper = 50; // = 30;
const double epsrel = 1e-6; // = 1e-4;

//================= code ===============================

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

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  return std::exp(-p * w0) / (M_PION + 2 * p) * (cos_pw - sin_pw / (p * w));
}


template<class T>
void integrate_qags(T f, std::vector<double> params, double lower, double upper, double epsabs, double epsrel, gsl_integration_workspace *w, double &result, double &error) {
  gsl_function F;
  F.function = f;
  F.params = &params;
  gsl_integration_qags (&F, lower, upper, epsabs, epsrel, 1000, w, &result, &error);
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

  double p_interval;
	double pe; // magnitude of electron's momentum. 
  std::vector<double> paras; // [w, w0]

  Func(double beta);
  std::vector<double> operator()(const std::vector<double>& v) const;
};


Func::Func(double beta) {
  pe = 0.5 * M_PION * beta;
}

std::vector<double> Func::operator()(const std::vector<double>& v) const{

  std::vector<double> ans(1);

	double p = v[0] * p_interval; // p: [0, p_interval]
	double c = 2. * v[1] - 1.; // cos: [-1, 1]
	double Epe = std::sqrt(p*p + (M_PION/2.)*(M_PION/2.) - 2 * p * pe * c);

  double w = paras[0];
  double w0 = paras[1];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  ans[0] =  std::exp(-Epe * w0) / (Epe * (- M_PION*M_PION + 4. * pe * pe * c * c)) * (cos_pw - sin_pw / (p * w));

  ans[0] *= 2. * p_interval;

  return ans;
}

template<class T>
double integrate_CUBA(T func, double eps_rel)
{
  std::vector<double> integral, error, prob;
  int nregions, neval, fail;

	integrateCuhre(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., eps_rel);
  return integral[0];
}


gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

int main (void)
{
  int w_max = int(L_LIMIT * std::sqrt(3)) + 1;
  int w0_max = T_LIMIT;
  double beta = get_beta(target);
  std::cout << "M_H (in lattice unit): " << M_PION << std::endl;
  std::cout << "w: [" << -w_max << ", " << w_max << "]" << std::endl;
  std::cout << "w0: [" << -w0_max << ", " << w0_max << "]" << std::endl;
  std::cout << "beta: " << beta << std::endl;
  std::cout << std::string(30, '*') << std::endl;
  std::cout << "Upper limit of integral: " << upper << std::endl;
  std::cout << "Allowed relative error: " << epsrel << std::endl;
  std::cout << std::string(30, '*') << std::endl;
  
  double ret1, ret2, ret3, error;
  // double lower = 0.00001, upper = 30, epsrel = 1e-4;
  // double lower = 0.00001, upper = 30, epsrel = 1e-6;
  // double lower = 0.00001, upper = 50, epsrel = 1e-6;
  double lower = 0.00001;

  double log_beta = std::log((1 + beta) / (1 - beta));

  Func f3(beta);
  f3.p_interval = upper;

  for(int w=0; w<=w_max; ++w) 
    for(int w0=0; w0<=w0_max; ++w0) {
      double ret;
      std::vector<double> params {double(w), double(w0)};
      std::cout << "w=" << w << " w0=" << w0 << std::endl;

      f3.paras = params;

      if(w==0) ret =0.;
      else {
        integrate_qawc(f1, params, lower, upper, M_PION/2., 0, epsrel, workspace, ret1, error);
        integrate_qags(f2, params, lower, upper, 0, epsrel, workspace, ret2, error);
        ret3 = integrate_CUBA(f3, epsrel);

        ret = - std::exp(0.5 * M_PION * w0) / w / w * M_PI / (f3.pe * M_PION) * log_beta * ret1 
              + std::exp(-0.5 * M_PION * w0) / w / w * M_PI / (f3.pe * M_PION) * log_beta * ret2
              + 2.0 * M_PI / w / w * ret3;
      }
      std::cout << std::setprecision(10) << "integral = " << ret  << std::endl;
  }

  gsl_integration_workspace_free(workspace);

  return 0;
}

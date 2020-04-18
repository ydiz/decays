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

const std::string ensemble = "64I";
// const std::string ensemble = "48I_pqdm";
const std::string target = "Pion";

// ================ integration parameters ===================
const double lower = 0.00001;
const double upper = 50; // = 30;
// const double epsrel = 1e-6; // = 1e-4;
const double epsrel = 1e-4; // = 1e-4;










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

  // return std::exp(-p * w0) / (M_PION + 2 * p) * (cos_pw - sin_pw / (p * w));
  return std::exp(-p * w0) / (M_pi + 2 * p) * (cos_pw - sin_pw / (p * w));
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
	// double Epe = std::sqrt(p*p + (M_PION/2.)*(M_PION/2.) - 2 * p * pe * c);
	double Epe = std::sqrt(p*p + (M_pi/2.) * (M_pi/2.) - 2 * p * pe * c);

  double w = paras[0];
  double w0 = paras[1];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  // ans[0] =  std::exp(-Epe * w0) / (Epe * (- M_PION*M_PION + 4. * pe * pe * c * c)) * (cos_pw - sin_pw / (p * w));
  ans[0] =  std::exp(-Epe * w0) / (Epe * (- M_pi * M_pi + 4. * pe * pe * c * c)) * (cos_pw - sin_pw / (p * w));

  ans[0] *= 2. * p_interval;

  return ans;
}

template<class T>
double integrate_CUBA(T func, double eps_rel, int &fail)
{
  std::vector<double> integral, error, prob;
  int nregions, neval;

	integrateCuhre(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., eps_rel);
  return integral[0];
}


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


  /////////////////// test////////////////////
  int w = 1, w0 = 21;
  double ret3;
  int fail;
  f3.paras = {double(w), double(w0)};

  for(double epsrel: {1e-4, 1e-5, 1e-6, 1e-7}) {
  ret3 = integrate_CUBA(f3, epsrel, fail);
  std::cout << "epsrel: " << epsrel << " \t ret3: " << ret3 << std::endl;
  }
  // std::cout << "CUBA fail: " << fail << std::endl;
  assert(fail==0);

  assert(0);
  /////////////////// end of test////////////////////


  for(int w=0; w<=w_max; ++w) 
    for(int w0=0; w0<=w0_max; ++w0) {
      std::cout << "w=" << w << " w0=" << w0 << std::endl;

      double ret, ret1, ret2, ret3, error;
      if(w==0) ret =0.;
      else {
        std::vector<double> params1 {double(w), double(w0)};
        integrate_qawc(f1, params1, lower, upper, M_pi/2., 0, epsrel, workspace, ret1, error);
        // std::cout << "qawc estimated error: " << error << std::endl;

        std::vector<double> params2 {double(w), double(w0), M_pi};
        integrate_qags(f2, params2, lower, upper, 0, epsrel, workspace, ret2, error);
        // std::cout << "qsgs estimated error: " << error << std::endl;

        int fail;
        f3.paras = {double(w), double(w0)};
        ret3 = integrate_CUBA(f3, epsrel, fail);
        // std::cout << "CUBA fail: " << fail << std::endl;
        assert(fail==0);

        double coef1 = - std::exp(0.5 * M_pi * w0) / w / w * M_PI / (f3.pe * M_pi) * log_beta;
        double coef2 = std::exp(-0.5 * M_pi * w0) / w / w * M_PI / (f3.pe * M_pi) * log_beta;
        double coef3 = 2.0 * M_PI / w / w;
        std::cout << "c1:\t" << coef1 << '\t' << ret1   << std::endl;
        std::cout << "c2:\t" << coef2 << '\t' << ret2   << std::endl;
        std::cout << "c3:\t" << coef3 << '\t' << ret3   << std::endl;
        ret = coef1 * ret1 
              + coef2 * ret2
              + coef3 * ret3;
      }
      std::cout << std::setprecision(10) << "integral = " << ret  << std::endl;
  }

  gsl_integration_workspace_free(workspace);

  return 0;
}

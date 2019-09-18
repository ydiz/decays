#pragma once

#include <cuba.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <math.h>
#include <chrono>

struct Timer {
	using clock = std::chrono::system_clock;
	clock::time_point start_time;
	Timer() {start_time = clock::now();}
} cur_time;

std::ostream& operator<<(std::ostream &out, const Timer &time) {
	using namespace std::chrono;
	auto t = system_clock::now();
	out << (duration<float>{t- time.start_time}).count();
	return out;
}

template <class F>
inline int cubaFunction(const int* ndim, const double x[],
  const int* ncomp, double f[], void* userdata)
{
  std::vector<double> vx(*ndim, 0.0);
  std::memcpy(vx.data(), x, *ndim * sizeof(double));
  std::vector<double> vf((*((const F*)userdata))(vx));
  assert(vf.size() == (size_t)*ncomp);
  std::memcpy(f, vf.data(), *ncomp * sizeof(double));
  return 0;
}

// ?? epsrel is not reached but fail is 0
template <class F>
void integrateDivonne(
    std::vector<double>& integral, std::vector<double>& error, std::vector<double>& prob,
    int& nregions, int& neval, int& fail,
    const int ndim, const int ncomp, const F& f,
    const double epsabs = 0.0, const double epsrel = 1.0e-5,
    const int flags = 0, const int seed = 23,
    const int mineval = 128, const int maxeval = 16 * 1024 * 1024 * 4,
    // const int key1 = 7, const int key2 = 7, const int key3 = 1, const int maxpass = 5,
    const int key1 = 47, const int key2 = 1, const int key3 = 1, const int maxpass = 5,
    const double border = 0.0, const double maxchisq = 10.0, const double mindeviation = 0.25,
    const double ngiven = 0, int ldxgiven = 0, double xgiven[] = NULL,
    const int nextra = 0, peakfinder_t peakfinder = NULL,
    const char* statefile = NULL, void* spin = NULL)
{
  if (0 == ldxgiven) {
    ldxgiven = ndim;
  }
  integral.resize(ncomp);
  error.resize(ncomp);
  prob.resize(ncomp);

  cubacores(0, 0); // this will make program faster
  Divonne(ndim, ncomp, cubaFunction<F>, (void*)&f, 1,
      epsrel, epsabs, flags, seed, mineval, maxeval,
      key1, key2, key3, maxpass, border, maxchisq, mindeviation,
      ngiven, ldxgiven, xgiven, nextra, peakfinder, statefile, spin,
      &nregions, &neval, &fail,
      integral.data(), error.data(), prob.data());

	std::cout << cur_time  << "s	eval = " << neval << " fail = " << fail<< std::endl;
	for(int i=0; i<ncomp; ++i) {
		std::cout << cur_time << "s	ncomp " << i << ": integral = " << integral[i] 
							<< " error = " << error[i] << " prob = " << prob[i] << std::endl;
	}

	std::cout << std::string(20, '-') << std::endl;
}

template <class F>
void integrateCuhre(
    std::vector<double>& integral, std::vector<double>& error, std::vector<double>& prob,
    int& nregions, int& neval, int& fail,
    const int ndim, const int ncomp, const F& f,
    const double epsabs = 0.0, const double epsrel = 1.0e-5,
    const int flags = 0, const int mineval = 128, const int maxeval = 16 * 1024 * 1024 * 4,
    const int key = 7, const char* statefile = NULL, void* spin = NULL)
{
  integral.resize(ncomp);
  error.resize(ncomp);
  prob.resize(ncomp);

  cubacores(0, 0); // this will make program faster
  Cuhre(ndim, ncomp, cubaFunction<F>, (void*)&f, 1,
      epsrel, epsabs, flags, mineval, maxeval,
      key, statefile, spin,
      &nregions, &neval, &fail,
      integral.data(), error.data(), prob.data());

	// // std::cout << cur_time  << "s	eval = " << neval << " fail = " << fail<< std::endl;
	// for(int i=0; i<ncomp; ++i) {
	// 	std::cout << cur_time << "s	ncomp " << i << ": integral = " << integral[i] 
	// 						<< " error = " << error[i] << " prob = " << prob[i] << std::endl;
	// }
  //
	// std::cout << std::string(20, '-') << std::endl;
}


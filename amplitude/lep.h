#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include <headers/headers.h>

// #include <constants_macro.h>
// #include "io.h"
// #include <pGG.h>
// #include "utils.h"

namespace Grid {
namespace QCD {

// IMPORTANT: gcoor should have already be "smod".
double get_integral_site(const std::vector<int> &gcoor, const std::vector<std::vector<double>> &data, const std::vector<int> &gdim, int space_limit, int time_limit) {
  
  int x = gcoor[0], y = gcoor[1], z = gcoor[2], t = gcoor[3];

	if( (x <= space_limit && x >= -space_limit) && (y <= space_limit && y >= -space_limit) && (z <= space_limit && z >= -space_limit) && (t<=time_limit && t >= -time_limit)) {
		double r = std::sqrt(x*x + y*y + z*z); 	
		int floor = int(r), ceiling = int(r) + 1;
		double ret = linear_interpolation(r, floor, data[floor][std::abs(t)], ceiling, data[ceiling][std::abs(t)]);

    return ret;
	}
	else return 0.;
}


void read_integrals(const std::string &filename, int space_limit, int time_limit, std::vector<std::vector<double>> &data) 
{

	std::vector<int> shape {int(space_limit * std::sqrt(3)) + 2, time_limit + 1};
	// assert(shape[0]==29 && shape[1]==17);	
  
  data.resize(shape[0]);
	for(int i = 0; i < shape[0]; i++)
	{
			data[i].resize(shape[1]);
	}

	std::ifstream f(filename);
	
	std::cout << "before loop " << filename << std::endl;
	for(int i=0; i<shape[0]; ++i)
		for(int j=0; j<shape[1]; ++j)
				f >> data[i][j];

	assert(!f.eof()); // make sure there aren't too few numbers
	double tmp;
	f>>tmp;
	assert(f.eof()); // make sure there aren't too many numbers
	
	f.close();
	std::cout << "Done reading integral from " << filename << std::endl;
}


void get_leptonic(const std::string &filename, LatticePGG &lat, int space_limit, int time_limit) {
	
	std::vector<std::vector<double>> data;
	read_integrals(filename, space_limit, time_limit, data);

	parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++){

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat.Grid()->_fdimensions);

    std::vector<int> gcoor_v(gcoor.begin(), gcoor.end());
    std::vector<int> fdims_v(lat.Grid()->_fdimensions.begin(), lat.Grid()->_fdimensions.end());
		double val = get_integral_site(gcoor_v, data, fdims_v, space_limit, time_limit); 

		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m(0, 1)()() = Complex(val * gcoor[Zdir], 0); 
		m(0, 2)()() = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
		m(1, 2)()() = Complex(val * gcoor[Xdir], 0); 
		m(1, 0)()() = - m(0, 1)()();
		m(2, 0)()() = - m(0, 2)()();
		m(2, 1)()() = - m(1, 2)()();

		// m()()(0, 1) = Complex(val * gcoor[Zdir], 0); 
		// m()()(0, 2) = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
		// m()()(1, 2) = Complex(val * gcoor[Xdir], 0); 
		// m()()(1, 0) = - m()()(0, 1);
		// m()()(2, 0) = - m()()(0, 2);
		// m()()(2, 1) = - m()()(1, 2);
		pokeLocalSite(m, lat, lcoor);
	}
}

}}

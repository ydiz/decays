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

using Data = std::vector<std::vector<std::vector<double>>>;
using Data_p1 = std::vector<std::vector<std::vector<std::vector<double>>>>;

// double get_p3_site(const std::vector<int> &gcoor, double data[13][17][17], const std::vector<int> &gdim) {
double get_integral_site_p3(const std::vector<int> &gcoor, const Data &data, const std::vector<int> &gdim, int space_limit, int time_limit) {

	int x = qlat::smod(gcoor[0], gdim[0]);
	int y = qlat::smod(gcoor[1], gdim[1]);
	int z = qlat::smod(gcoor[2], gdim[2]);
	int t = qlat::smod(gcoor[3], gdim[3]);

	if( (x <= space_limit && x >= -space_limit) && (y <= space_limit && y >= -space_limit) && (z <= space_limit && z >= -space_limit) && (t<=time_limit && t >= -time_limit)) {
		double r = std::sqrt(x*x + y*y); 	
		int floor = int(r), ceiling = int(r) + 1;
		double ret = linear_interpolation(r, floor, data[floor][std::abs(z)][std::abs(t)], ceiling, data[ceiling][std::abs(z)][std::abs(t)]);
    if(z<0) ret = -ret; // p3 integral is odd in z, and even in x,y,t

    return ret;
	}
	else return 0.;
}

double get_integral_site_p1(const std::vector<int> &gcoor, const Data_p1 &data, const std::vector<int> &gdim, int space_limit, int time_limit) {

	int x = qlat::smod(gcoor[0], gdim[0]);
	int y = qlat::smod(gcoor[1], gdim[1]);
	int z = qlat::smod(gcoor[2], gdim[2]);
	int t = qlat::smod(gcoor[3], gdim[3]);

	if( (x <= space_limit && x >= -space_limit) && (y <= space_limit && y >= -space_limit) && (z <= space_limit && z >= -space_limit) && (t<=time_limit && t >= -time_limit)) {
    double ret = data[std::abs(x)][std::abs(y)][std::abs(z)][std::abs(t)];
    if(x<0) ret = -ret; // p1 integral is odd in x, and even in y,z,t
    return ret;
	}
	else return 0.;
}

// double data[13][17][17]; //r[0, 12], z[-8, 8], t [0, 16]
void read_integrals_p3(const std::string &filename, int space_limit, int time_limit, Data &data) 
{

	std::vector<int> shape {int(space_limit * std::sqrt(2)) + 2, space_limit + 1, time_limit + 1};
	assert(shape[0]==24 && shape[1]==17 && shape[2]==17);	
  
  data.resize(shape[0]);
	for(int i = 0; i < shape[0]; i++)
	{
			data[i].resize(shape[1]);
			for(int j = 0; j < shape[1]; j++)
					data[i][j].resize(shape[2]);
	}

	std::ifstream f(filename);
	
	std::cout << "before loop " << filename << std::endl;
	for(int i=0; i<shape[0]; ++i)
		for(int j=0; j<shape[1]; ++j)
			for(int k=0; k<shape[2]; ++k) {
				f >> data[i][j][k];
			}

	assert(!f.eof()); // make sure there aren't too few numbers
	double tmp;
	f>>tmp;
	assert(f.eof()); // make sure there aren't too many numbers
	
	f.close();
	std::cout << "Done reading integral from " << filename << std::endl;
}


void read_integrals_p1(const std::string &filename, int space_limit, int time_limit, Data_p1 &data) 
{
	std::vector<int> shape {space_limit + 1, space_limit + 1, space_limit + 1, time_limit + 1}; // [r, z, t]
	assert(shape[0]==17 && shape[1]==17 && shape[2]==17 && shape[3]==17);	
  
  data.resize(shape[0]);
	for(int i = 0; i < shape[0]; i++)
	{
			data[i].resize(shape[1]);
			for(int j = 0; j < shape[1]; j++) {
					data[i][j].resize(shape[2]);
          for(int k = 0; k < shape[2]; k++) 
            data[i][j][k].resize(shape[3]);
      }
	}

	std::ifstream f(filename);
	
	std::cout << "before loop " << filename << std::endl;
	for(int x=0; x<shape[0]; ++x)
		for(int y=0; y<shape[1]; ++y)
			for(int z=0; z<shape[2]; ++z) 
        for(int t=0; t<shape[3]; ++t) {
          f >> data[x][y][z][t];
        }

	assert(!f.eof()); // make sure there aren't too few numbers
	double tmp;
	f>>tmp;
	assert(f.eof()); // make sure there aren't too many numbers
	
	f.close();
	std::cout << "Done reading integral from " << filename << std::endl;
}


void get_leptonic_CUBA3d(const std::string &filename_p1, const std::string &filename_p3, LatticePGG &lat, int space_limit, int time_limit) {
	
	Data data_p3;
	Data_p1 data_p1;
	read_integrals_p3(filename_p3, space_limit, time_limit, data_p3);
	read_integrals_p1(filename_p1, space_limit, time_limit, data_p1);

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

		double val_p3 = get_integral_site_p3(gcoor, data_p3, lat._grid->_fdimensions, space_limit, time_limit); 
		double val_p1 = get_integral_site_p1(gcoor, data_p1, lat._grid->_fdimensions, space_limit, time_limit); 

    std::vector<int> gcoor_p2 = gcoor;
    std::swap(gcoor_p2[0], gcoor_p2[1]); // To calculate the integral with p2, just swap x and y coordinates in the integral with p1;
		double val_p2 = get_integral_site_p1(gcoor_p2, data_p1, lat._grid->_fdimensions, space_limit, time_limit); 
		
		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m()()(0, 1) = Complex(val_p3, 0); 
		m()()(0, 2) = Complex(-val_p2, 0); 
		m()()(1, 2) = Complex(val_p1, 0); 
		m()()(1, 0) = - m()()(0, 1);
		m()()(2, 0) = - m()()(0, 2);
		m()()(2, 1) = - m()()(1, 2);

		// if( (gcoor[0] > 8 &&gcoor[0] <24) || (gcoor[1] > 8 &&gcoor[1] <24) || (gcoor[2] > 8 &&gcoor[2] <24)|| (gcoor[3] > 8 &&gcoor[3] <56)) m = 0.; 
		pokeLocalSite(m, lat, lcoor);
	}
}

}}

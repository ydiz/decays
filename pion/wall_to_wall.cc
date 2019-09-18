#include <qlat/qlat.h>
#include <qlat/grid.h>

#include <headers/headers.h>
// #include "qlat_wrapper/qlat_wrapper.h"
// #include "constants_macro.h"

#include <dirent.h>
#include <sys/stat.h>

namespace Grid {
namespace QCD {

std::vector<double> kaon_corr(int traj, GridBase *grid) {
  int T = grid->_fdimensions[Tdir];
	LatticePropagator prop_d(grid);
	LatticePropagator prop_s(grid);

  std::vector<double> corr(T, 0.);

  for(int t=0; t<T; ++t) {

    std::string path_d = wall_path_ud_32D(traj, t);
    std::string path_s = wall_path_strange_32D(traj, t);
    read_qlat_propagator(prop_d, path_d);
    read_qlat_propagator(prop_s, path_s);

    std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum_d;
    sliceSum(prop_d, slice_sum_d, Tdir);
    std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum_s;
    sliceSum(prop_s, slice_sum_s, Tdir);

    for(int i=0; i<T; ++i) {
      int sep = (i < t) ? i - t + T : i - t;
      double tmp = TensorRemove(trace(slice_sum_d[i] * adj(slice_sum_s[i]))).real();
      corr[sep] += tmp;
    }
  }

  for(auto &x: corr) x /= T;
  std::cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< std::endl;
  std::cout << corr << std::endl;

  return corr;
}


std::vector<double> pion_corr(int traj, GridBase *grid) {
  int T = grid->_fdimensions[Tdir];
	LatticePropagator prop(grid);

  std::vector<double> corr(T, 0.);

  for(int t=0; t<T; ++t) {

    // std::string path = wall_path_ud_32D(traj, t);
    // std::string path = wall_path_24ID(traj, t);
    std::string path = wall_path_ud_48I(traj, t);
    read_qlat_propagator(prop, path);

    std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum;
    sliceSum(prop, slice_sum, Tdir);

    for(int i=0; i<T; ++i) { // slice_sum[i] = P(i, t)
      int sep = (i < t) ? i - t + T : i - t;
      double tmp = TensorRemove(trace(slice_sum[i] * adj(slice_sum[i]))).real();
      corr[sep] += tmp;
    }
  }

  for(auto &x: corr) x /= T;
  std::cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< std::endl;
  std::cout << corr << std::endl;

  return corr;
}


}}


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


// std::vector<int> gcoor({32, 32, 32, 64});
// std::vector<int> gcoor({24, 24, 24, 64});
std::vector<int> gcoor({48, 48, 48, 96});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor); 
	LatticePropagator prop(grid);

	// int traj_start = 200, traj_end = 430, traj_sep = 10; // for 32IDF
	// int traj_start = 1200, traj_end = 1200, traj_sep = 10; 
	// int traj_start = 2370, traj_end = 2510, traj_sep = 10; // for 24ID
	int traj_start = 1490, traj_end = 1730, traj_sep = 20; // for 48I  (We also have traj 1290-1430)
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

	std::cout << std::string(20, '*') << std::endl;
	std::cout << "traj_start: " << traj_start << std::endl;
	std::cout << "traj_end: " << traj_end << std::endl;
	std::cout << "traj_sep: " << traj_sep << std::endl;
	std::cout << "traj_num: " << traj_num << std::endl;
	std::cout << std::string(20, '*') << std::endl;

	std::vector<double> average_corr(gcoor[3], 0.);

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    std::vector<double> corr = pion_corr(traj, grid);
    // std::vector<double> corr = kaon_corr(traj, grid);
		for(int i=0; i<corr.size(); ++i) average_corr[i] += corr[i];
	}

	for(auto &x: average_corr) x /= traj_num;
	cout << std::string(30, '*') << endl;
	cout << "average wall to wall correlator over "<< traj_num << " trajectoies:" << endl;
	cout << average_corr << endl;

  end();

  return 0;
}

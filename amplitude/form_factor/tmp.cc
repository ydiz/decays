#include <Grid/Grid.h>
#include <qlat/qlat.h>
#include <dirent.h>
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"
// #include "calculate_avg.h"
#include "form_factor_scalar_H.h"


#include "jackknife.h"



using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

  int traj = 1200;

  double Z_V = 0.72672;
  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * 2 * M_PION / (9428.492 / std::sqrt(32*32*32.));


  std::vector<std::vector<double>> jackknife_results(17, std::vector<double>(28));

    int traj_start = 1110, traj_end = 1370;
    for(int traj=traj_start; traj<traj_end; traj+=10) {

      std::cout << "traj: " << traj << std::endl;

      LatticePGG three_point(grid);
      std::string file = three_point_exact_path(traj);
      read_cheng_PGG(three_point, file);

      LatticePGG hadronic(three_point._grid);
      hadronic = imag(three_point); 
      LatticeComplex H(three_point._grid);
      H = toScalar(hadronic);
      H = hadron_coeff * H;

      LatticeComplex weight(grid);
      form_factor_integrand_scalar(weight);
      LatticeComplex prod(grid);
      prod = H * weight;

      int T = hadronic._grid->_fdimensions[Tdir];
      std::vector<iSinglet<Complex>> ret(T);
      sliceSum(prod, ret, Tdir);

      std::vector<double> ret_real(T / 2 + 1);
      ret_real[0] = ret[0]()()().real();
      ret_real[T/2] = ret[T/2]()()().real();
      for(int i=1; i<ret_real.size()-1; ++i) ret_real[i] = ret[i]()()().real() + ret[T-i]()()().real(); // add t and -t

      cout << "sum of time slice + sum of t and -t: " << ret_real << endl;

      std::vector<double> ret_cumulative(T / 2 + 1);
      ret_cumulative[0] = ret_real[0];
      for(int i=1; i<ret_real.size(); ++i) ret_cumulative[i] = ret_real[i] + ret_cumulative[i-1]; // result with cutoff
      cout << "cumulative: " << ret_cumulative << endl;

      double other_coeff = 8 * M_PI * M_PI * 0.092424 / M_PION / M_PION / M_PION / M_PION;
      std::vector<double> F(ret.size());
      for(int i=0; i<ret.size(); ++i) F[i] = hadron_coeff * other_coeff * ret_cumulative[i];

      for(int time_cutoff = 0; time_cutoff <= 16; ++time_cutoff) {
        int t_idx = time_cutoff;
        jackknife_results[t_idx][(traj - 1110)/10] = F[time_cutoff]; 
      }
    }
    
    // for(int time_cutoff = para.time_cutoff_start; time_cutoff <= para.time_cutoff_end; ++time_cutoff) {
    for(int time_cutoff = 0; time_cutoff <= 16; ++time_cutoff) {

    std::cout << std::string(20, '*') << std::endl;
    cout << "time cutoff: " << time_cutoff << endl;
    cout << "jackknife samples: " << endl;
    // cout << jackknife_results[time_cutoff - para.time_cutoff_start] << endl;
    cout << jackknife_results[time_cutoff] << endl;

    // jack_stats(jackknife_results[time_cutoff - para.time_cutoff_start]); // jackknife average and error
    jack_stats(jackknife_results[time_cutoff]); // jackknife average and error
  }
 

	Grid_finalize();
  return 0;
}

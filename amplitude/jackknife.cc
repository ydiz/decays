#include <Grid/Grid.h>
#include <dirent.h>
#include "amplitude.h"
#include "jackknife.h"
#include "jack_init.h"
#include "lep.h"
#include "lep_CUBA3d.h"

#include "imaginary_part.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;
using namespace qlat;

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  Grid::Coordinate mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, qlat::Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  Jack_para para;
  init_para(argc, argv, para);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(para.lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  int T = grid->_fdimensions[Tdir];
  assert(T%4 ==0);

  LatticePGG leptonic(grid);
  para.get_leptonic(leptonic); // generate leptonic part

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  // We alwasy calculate the time cutoff range [0, T/4]. The maximum is T/4 because: 1) for real part, leptonic factor is calculated upto T/4.  2) If time cutoff is too large, EM could effectively be on the wrong side of wall.
  std::vector<std::vector<double>> jackknife_results(T/4+1, std::vector<double>(para.traj_num));

  int skipped = 0; // number of trajectories already skipped
  for(int traj = para.traj_start; traj <= para.traj_end; traj += para.traj_sep) {

    if(std::find(para.traj_skip.begin(), para.traj_skip.end(), traj) != para.traj_skip.end()) {
      ++skipped;
      continue;
    }

    LatticePGG three_point(grid);
    para.get_three_point(three_point, traj);

    std::vector<double> cutoffs = para.get_result_with_cutoff(three_point, leptonic);

    int traj_idx = (traj - para.traj_start) / para.traj_sep - skipped;
    for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {
      jackknife_results[time_cutoff][traj_idx] = cutoffs[time_cutoff]; 
    }
  }

  // ======================================================================

  // The amplitude are summed over in [0, time_cutoff]
  for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {

    std::cout << std::string(20, '*') << std::endl;
    cout << "time cutoff: " << time_cutoff << endl;
    cout << "jackknife samples: " << endl;
    cout << jackknife_results[time_cutoff] << endl;

    std::vector<RealD> jack = jack_stats(jackknife_results[time_cutoff]); // jackknife average and error
    if(para.target != "form_factor") {
      std::cout << "Relative Branching Ratio average: " << para.BR_coeff * jack[0] * jack[0] << std::endl;
      std::cout << "Relative Branching Ratio error: " << std::abs(2. * para.BR_coeff * jack[0] * jack[1]) << std::endl;
    }
    std::cout << std::string(20, '*') << std::endl;
  }

  Grid_finalize();
  return 0;
}

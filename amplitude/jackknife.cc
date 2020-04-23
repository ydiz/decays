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

  map<string, LatticePGG> leptonic; // map: target -> leptonic part
  for(const string &target: para.targets) {
    leptonic.insert(make_pair(target, para.get_leptonic(grid, target)));
  }

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  // We alwasy calculate the time cutoff range [0, T/4]. The maximum is T/4 because: 1) for real part, leptonic factor is calculated upto T/4.  2) If time cutoff is too large, EM could effectively be on the wrong side of wall.
  map<string, vector<vector<double>>> jk_rst; // map: target -> jk[t_sep][traj]
  for(const string &target: para.targets) {
    jk_rst[target].assign(T/4+1, vector<double>(para.traj_num));
  }

  int skipped = 0; // number of trajectories already skipped
  for(int traj = para.traj_start; traj <= para.traj_end; traj += para.traj_sep) {

    if(find(para.traj_skip.begin(), para.traj_skip.end(), traj) != para.traj_skip.end()) {
      ++skipped;
      continue;
    }

    LatticePGG three_point(grid);
    para.get_three_point(three_point, traj);

    for(const string &target: para.targets) {
      vector<double> cutoffs = para.get_result_with_cutoff(three_point, leptonic.at(target), target);

      int traj_idx = (traj - para.traj_start) / para.traj_sep - skipped;
      for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {
        jk_rst.at(target)[time_cutoff][traj_idx] = cutoffs[time_cutoff]; 
      }
    }
  }

  // ======================================================================

  for(const string &target: para.targets) {

    string fname = para.output_prefix + "/" + target + (para.cutoff_type=="4D" ? "_4D_cutoff" : "") + ".txt";
    ofstream f(fname);
    f << "target: " << target << endl;
    // The amplitude are summed over in [0, time_cutoff]
    for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {

      f << string(20, '*') << endl;
      f << "time cutoff: " << time_cutoff << endl;
      f << "jackknife samples: " << endl;
      f << jk_rst.at(target)[time_cutoff] << endl;

      vector<RealD> jack = jack_stats(jk_rst.at(target)[time_cutoff]); // jackknife average and error
      f << "jackknife average: " << jack[0] << std::endl;
      f << "jackknife error: " << jack[1] << std::endl;
      if(target != "form_factor") {
        f << "Relative Branching Ratio average: " << para.BR_coeff * jack[0] * jack[0] << endl;
        f << "Relative Branching Ratio error: " << abs(2. * para.BR_coeff * jack[0] * jack[1]) << endl;
      }
      f << string(20, '*') << endl;
    }
    f.close();
  }

  Grid_finalize();
  return 0;
}

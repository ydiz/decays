#include <Grid/Grid.h>
#include "../amplitude/amplitude.h"
#include "../amplitude/jackknife.h"
#include "../amplitude/jack_init.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  Grid::Coordinate mpi_coor = GridDefaultMpi();

  std::vector<std::string> diagrams = {
                                       "JJ_sBar_d_K/D1", "JJ_sBar_d_K/D2", "JJ_sBar_d_K/D3", 
                                       "typeI/D1Q1", "typeI/D1Q2", "typeI/D2Q1", "typeI/D2Q2",
                                       "typeII/Q1", "typeII/Q2",
                                       "typeIII/D1Q1K", "typeIII/D1Q1Kbar", "typeIII/D1Q2K", "typeIII/D1Q2Kbar", 
                                       "typeIII/D2Q1K", "typeIII/D2Q1Kbar", "typeIII/D2Q2K", "typeIII/D2Q2Kbar", 
                                       "typeIII/D3Q1", "typeIII/D3Q2"
                                      };

  Jack_para para;
  init_para(argc, argv, para);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(para.lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  int T = grid->_fdimensions[Tdir];

  LatticeKGG leptonic(grid);  // Do not use para.get_leptonic, which automatically adds translational exp factor
  form_factor_integrand(leptonic, para.M_h);

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  // We alwasy calculate the time cutoff range [0, T/4]. The maximum is T/4 because: 1) for real part, leptonic factor is calculated upto T/4.  2) If time cutoff is too large, EM current could effectively be on the wrong side of wall.
  map<string, vector<vector<double>>> jk_rst; // map: diagram -> jk[t_cutoff][traj]
  for(const string &diagram: diagrams) {
    jk_rst[diagram].assign(T/4+1, vector<double>(para.traj_num));
  }

  int skipped = 0; // number of trajectories already skipped
  for(int traj = para.traj_start; traj <= para.traj_end; traj += para.traj_sep) {

    if(find(para.traj_skip.begin(), para.traj_skip.end(), traj) != para.traj_skip.end()) {
      ++skipped;
      continue;
    }

    for(const string &diagram: diagrams) {

      // For form_factor, the three point function should be <J(0) J(x)|Hw|K >
      LatticeKGG hadronic(grid); // Do not use para.get_three_point, which automatically adds translational exp factor
      std::string fname = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/results/" + diagram + "." + std::to_string(traj);
      readScidac(hadronic, fname);

      vector<double> cutoffs = para.get_result_with_cutoff(hadronic, leptonic, para.targets[0]);

      int traj_idx = (traj - para.traj_start) / para.traj_sep - skipped;
      for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {
        jk_rst.at(diagram)[time_cutoff][traj_idx] = cutoffs[time_cutoff]; 
      }
    }
  }

  // ======================================================================
  for(const string &diagram: diagrams) {
    std::cout << diagram << ": ";
    std::cout << jk_rst.at(diagram) << std::endl;
  }

  // for(const string &target: para.targets) {
  //
  //   string fname = para.output_prefix + "/" + target + (para.cutoff_type=="4D" ? "_4D_cutoff" : "") + ".txt";
  //   ofstream f(fname);
  //   f << "target: " << target << endl;
  //   // The amplitude are summed over in [0, time_cutoff]
  //   for(int time_cutoff = 0; time_cutoff <= T/4; ++time_cutoff) {
  //
  //     f << string(20, '*') << endl;
  //     f << "time cutoff: " << time_cutoff << endl;
  //     f << "jackknife samples: " << endl;
  //     f << jk_rst.at(target)[time_cutoff] << endl;
  //
  //     vector<RealD> jack = jack_stats(jk_rst.at(target)[time_cutoff]); // jackknife average and error
  //     f << "jackknife average: " << jack[0] << std::endl;
  //     f << "jackknife error: " << jack[1] << std::endl;
  //     if(target != "form_factor") {
  //       f << "Relative Branching Ratio average: " << para.BR_coeff * jack[0] * jack[0] << endl;
  //       f << "Relative Branching Ratio error: " << abs(2. * para.BR_coeff * jack[0] * jack[1]) << endl;
  //     }
  //     f << string(20, '*') << endl;
  //   }
  //   f.close();
  // }

  Grid_finalize();
  return 0;
}

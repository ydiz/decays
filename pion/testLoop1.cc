#include "connected.h"
#include "disc.h"

// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});


using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;

RealD mean(const std::vector<RealD>& data)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/RealD(N);
}

RealD jack_mean(const std::vector<RealD>& data, int sample)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/RealD(N-1);
}

RealD jack_std(const std::vector<RealD>& jacks, RealD mean)
{
    int N = jacks.size();
    RealD std(0.0);
    for(int i=0; i<N; ++i){ std += std::pow(jacks[i]-mean, 2.0); }
    return std::sqrt(RealD(N-1)/RealD(N)*std);
}

std::vector<RealD> jack_stats(const std::vector<RealD>& data)
{
  int N = data.size();
  std::vector<RealD> jack_samples(N);
  std::vector<RealD> jack_stats(2);

  jack_stats[0] = mean(data);
  for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
  jack_stats[1] = jack_std(jack_samples, jack_stats[0]);

  std::cout << "jackknife average: " << jack_stats[0] << std::endl;
  std::cout << "jackknife error: " << jack_stats[1] << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  return jack_stats;
}




int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

  int traj_start = 1000, traj_end = 2290, traj_sep = 10; // need to have both loop and wall src propagators 
  std::vector<int> traj_skip {1020, 1060, 1100, 1340};
  // int traj_start = 2000, traj_end = 2000, traj_sep = 10; // need to have both loop and wall src propagators 
  // std::vector<int> traj_skip {};
  int traj_num = (traj_end - traj_start) / traj_sep + 1 - traj_skip.size();

  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_skip: " << traj_skip << std::endl;

  int V = 1;
  for(int x: gcoor) V *= x;

  vector<double> r_avgs, i_avgs;
  for(int traj = traj_start; traj<=traj_end; traj+=traj_sep) {
    if(std::find(traj_skip.begin(), traj_skip.end(), traj) != traj_skip.end()) continue;
    std::cout << "traj: " << traj << std::endl;

    LatticeLoop loop1(grid);
    read_loop(loop1, loop_path_24ID(traj));

    std::vector<LatticeComplex> nums(4, grid);
    for(int i=0; i<4; ++i) nums[i] = peekColour(loop1, i);

    double real_avg = 0., imag_avg = 0.;
    for(int i=0; i<4; ++i) {
      auto t = sum(nums[i]); 
      real_avg += t()()().real();
      imag_avg += t()()().imag();
    }
    real_avg /= 4. * V;
    imag_avg /= 4. * V;

    r_avgs.push_back(real_avg);
    i_avgs.push_back(imag_avg);
  }

  std::cout << "Real Part" << std::endl;
  jack_stats(r_avgs);
  std::cout << "Imaginary Part" << std::endl;
  jack_stats(i_avgs);

  end();

  return 0;
} 
    // auto t = sum(loop1)()();
    // Complex avg = ( t(0) + t(1) + t(2) + t(3) ) / (4. * V);
    //
    // LatticeLoop loop1_real(grid);
    // loop1_real = real(loop1);
    //
    // iScalar<iScalar<iVector<Complex, 4>>> tmp;
    // for(auto &x: tmp) x = Complex(avg.real(), 0);
    //
    // LatticeLoop tt(grid);
    // tt = loop1_real - tmp;
    // // tt = loop1_real;
    // tt = pow(loop1_real, 2);
    // auto sigma = sqrt(sum(tt) * (1. / V));

  // rst = rst * (1. / double(traj_num));
  // std::cout << "sum of rst is: " << sum(rst) << std::endl;
  // std::cout << rst  << std::endl;



    // double real_var = 0., imag_var = 0.;
    // for(int i=0; i<4; ++i) {
    //   LatticeComplex tmp_real(grid);
    //   tmp_real = real(nums[i]);
    //   tmp_real = tmp_real - real_avg;
    //   tmp_real = tmp_real * tmp_real;
    //   auto t = sum(tmp_real); 
    //   real_var += t()()().real();
    //
    //   LatticeComplex tmp_imag(grid);
    //   tmp_imag = imag(nums[i]);
    //   tmp_imag = tmp_imag - imag_avg;
    //   tmp_imag = tmp_imag * tmp_imag;
    //   t = sum(tmp_imag); 
    //   imag_var += t()()().real();
    //   // imag_avg += t()()().imag();
    // }
    // double real_std = std::sqrt(real_var / (4. *V));
    // double imag_std = std::sqrt(imag_var / (4. *V));
    //
    // std::cout << "real part: avg: " << real_avg << " std: " << real_std << std::endl;
    // std::cout << "imag part: avg: " << imag_avg << " std: " << imag_std << std::endl;

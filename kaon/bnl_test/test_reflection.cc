
// On 8 nodes, Needs 16h for one trajectory (500 piont sources) // ~110s per point source

#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  Env env(gcoor, "24ID");


  // for(int i=0; i<64; ++i) {
  //   // int rank = env.grid->ThisRank();
  //   int rank = i; 
  //   Coordinate pcoor;
  //   env.grid->ProcessorCoorFromRank(rank, pcoor);
  //
  //   Coordinate new_pcoor = env.grid->_processors - Coordinate(std::vector<int>{1,1,1,1}) - pcoor;
  //   int partner = env.grid->RankFromProcessorCoor(new_pcoor);
  //
  //   std::cout << "rank: " << rank << "  pcoor: " << pcoor << " partner rank: " << partner << " partner pcoor: " << new_pcoor << std::endl;
  //   // std::cout << "Partner pcoor" << new_pcoor << std::endl;
  //   // std::cout << "Partner rank" << partner << std::endl;
  // }


  LatticeComplex lat(env.grid);
  get_reflection(lat);

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}

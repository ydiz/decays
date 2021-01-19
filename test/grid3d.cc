// 3-dim grid is Okay to use; see  Grid/tests/core/Test_where.cc 

#include <Grid/Grid.h>


using namespace std;
using namespace Grid;


template<class vobj>
void SumOverSlice(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim, int orthog) // adapted from ExtracSlice
{
  typedef typename vobj::scalar_object sobj;

  GridBase *lg = lowDim.Grid();
  GridBase *hg = higherDim.Grid();
  int nl = lg->_ndimension;
  int nh = hg->_ndimension;

  assert(nl+1 == nh);
  assert(orthog < nh);
  assert(orthog >= 0);
  assert(hg->_processors[orthog] == 1);

  int dl; dl = 0;
  for(int d=0; d<nh; d++){
    if ( d != orthog ) {
      assert(lg->_processors[dl]  == hg->_processors[d]);
      assert(lg->_ldimensions[dl] == hg->_ldimensions[d]);
      dl++;
    }
  }

  // the above should guarantee that the operations are local
  autoView(lowDimv, lowDim, CpuWrite);
  autoView(higherDimv, higherDim, CpuRead);
  thread_for(idx, lg->lSites(), {
      Coordinate lcoor(nl);
      lg->LocalIndexToLocalCoor(idx, lcoor);

      Coordinate hcoor(nh);
      int ddl=0;
      for(int d=0; d<nh; d++){
        if ( d!=orthog ) { 
          hcoor[d] = lcoor[ddl++];
        }
      }

      sobj s, s_sum = Zero();
      for(int i=0; i<hg->_fdimensions[orthog]; ++i) {
        hcoor[orthog] = i;
        peekLocalSite(s, higherDimv, hcoor);
        s_sum += s;
      }
      pokeLocalSite(s_sum, lowDimv, lcoor);
  });


}



int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);

  // Coordinate fdims4d({8, 8, 8, 16});
  Coordinate fdims4d({8, 8, 8, 16});
  Coordinate fdims3d({8, 8, 8});

  Coordinate simd4d = GridDefaultSimd(4, vComplexD::Nsimd());
  Coordinate simd3d = GridDefaultSimd(3, vComplexD::Nsimd());

  Coordinate mpi4d = GridDefaultMpi();
  assert(mpi4d[3] == 1);
  Coordinate mpi3d = mpi4d;
  mpi3d.resize(3);

  GridCartesian *grid4d = SpaceTimeGrid::makeFourDimGrid(fdims4d, simd4d, mpi4d);
  GridCartesian *grid3d = SpaceTimeGrid::makeFourDimGrid(fdims3d, simd3d, mpi3d);

  grid4d->show_decomposition();
  grid3d->show_decomposition();

  //////////////////////////////////////////////////////////////////////////

  // Test SumOverSlice  // Sum over time direction // 4d grid -> 3d grid
  Lattice<iVector<vComplex, 4>> coor(grid4d);  // peekLocalSite and pokeLocalSite only work with vComplex; do not work with vInteger
  for(int mu=0; mu<4; mu++){
    Lattice<iScalar<vComplex>> tmp(grid4d);    LatticeCoordinate(tmp, mu);
    PokeIndex<LorentzIndex>(coor, tmp, mu);
  }
  // std::cout << coor << std::endl;

  Lattice<iVector<vComplex, 4>> rst(grid3d); 
  // ExtractSlice(rst, coor, 2, Tp);
  SumOverSlice(rst, coor, Tp);  // Sum over time direction // 4d grid -> 3d grid
  std::cout << rst << std::endl;

  //////////////////////////////////////////////////////////////////////////

  // // Test FFT on three dimensional grid
  // FFT theFFT((GridCartesian *)grid3d);
  // Lattice<iVector<vComplex, 4>> rst_fft(grid3d); 
  // theFFT.FFT_all_dim(rst_fft, rst, FFT::forward);
  // std::cout << rst_fft << std::endl;



  Grid_finalize();
}



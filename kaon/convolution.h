#pragma once

namespace Grid {
namespace QCD {

// // |t_u - t_x| and |t_v - t_x| must be <= tsep - tsep2
// void convolution_set0(LatticePropagator &lat, int t_x, int tsep, int tsep2) {
//   const int T = lat.Grid()->_fdimensions[3];
//   thread_for(ss, lat.Grid()->lSites(), {
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
//
//     if(distance(t_x, gcoor[3], T) > tsep - tsep2) { // set to zero if |t_u - t_x| or |t_v - t_x| > tsep - tsep2
//       LatticePropagatorSite tmp; tmp = Zero();
//       pokeLocalSite(tmp, lat, lcoor);
//     }; 
//   });
// }


// |t_u - t_x| and |t_v - t_x| must be <= tsep - tsep2
template <typename U>
void convolution_set0(Lattice<U> &lat, int t_x, int tsep, int tsep2) {
  const int T = lat.Grid()->_fdimensions[3];
  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    if(distance(t_x, gcoor[3], T) > tsep - tsep2) { // set to zero if |t_u - t_x| or |t_v - t_x| > tsep - tsep2
      typename U::scalar_object tmp; tmp = Zero();
      pokeLocalSite(tmp, lat, lcoor);
    }; 
  });
}




// return exp(M_K * (v_0 - tK)), for any v_0 that satisfies |v_0 - tx| <= tsep - tsep2
LatticeComplex convolution_exp(GridCartesian *grid, int t_x, int tsep, int tsep2, double M_K) {
/*
tsep, tsep2: window size. Keep only points between [x - (tsep-tsep2), x + (tsep-tsep2)]
M_K: kaon mass on the lattice
*/
  LatticeComplex lat(grid);
  int T = lat.Grid()->_fdimensions[3];

  int tK = t_x - tsep;
  if(tK < 0) tK += T;

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		LatticeComplexSite m;

    int dist = distance(gcoor[3], t_x, T); // |tv - tx|
    if(dist > tsep - tsep2) m = Zero();
    else {
      int left_dist = left_distance(gcoor[3], tK, T); // v_0 - tK // the distance from v0 to tK if you can only move to the left
      double val = std::exp(M_K * left_dist); 
      m()()() = Complex(val, 0.);
    }

		pokeLocalSite(m, lat, lcoor);
  });

  return lat;
}








}}

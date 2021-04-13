#pragma once


namespace Grid {

// <eta sBar d K> diagrams
  
LatticeComplex sBar_d_T1D1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticeComplex rst = trace(g5 * wl_tEta_tK * adj(ws_K) * g5 * wl_eta);
  rst = 2. * real(rst); // Add contribution of K0 bar
  return rst;
}

LatticeComplex sBar_d_T1D2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticeComplex rst = trace(g5 * adj(ws_eta) * g5 * wl_K * adj(ws_tEta_tK));
  rst = 2. * real(rst); // Add contribution of K0 bar
  return rst;
}

LatticeComplex sBar_d_T2D1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticeComplex tmp = trace(g5 * wl_K * adj(ws_K));
  tmp = 2. * real(tmp); // Add contribution of K0 bar
  LatticeComplex rst = trace(g5 * wl_tEta_tEta) * tmp;
  return rst;
}
    
LatticeComplex sBar_d_T2D2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticeComplex tmp = trace(g5 * wl_K * adj(ws_K));
  tmp = 2. * real(tmp); // Add contribution of K0 bar
  LatticeComplex rst = trace(g5 * ws_tEta_tEta) * tmp;
  return rst;
}
    
    
///////////////////////////////////////////////////////////
// <eta Hw K> diagrams

LatticeComplex Hw_T1D1Q1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp_eta = wl_eta * adj(wl_eta);
  LatticePropagator tmp_K = wl_K * adj(ws_K);
  tmp_K = tmp_K + adj(tmp_K); // Add contribution of K0 bar 
  
  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * tmp_eta) * trace(gL[mu] * tmp_K);
  return rst;
}

LatticeComplex Hw_T1D1Q2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp_eta = wl_eta * adj(wl_eta);
  LatticePropagator tmp_K = wl_K * adj(ws_K);
  tmp_K = tmp_K + adj(tmp_K); // Add contribution of K0 bar 

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * tmp_eta * gL[mu] * tmp_K);
  return rst;
}

LatticeComplex Hw_T2D1Q1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_eta * g5 * wl_tEta_tK * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
  return rst;
}
    
LatticeComplex Hw_T2D1Q2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_eta * g5 * wl_tEta_tK * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * Lxx * gL[mu] * tmp);
  return rst;
}
    
LatticeComplex Hw_T2D2Q1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp =  wl_K * adj(ws_tEta_tK) * g5 * adj(ws_eta);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * Lxx) * trace(gL[mu] * tmp);
  return rst;
}
    
LatticeComplex Hw_T2D2Q2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_K * adj(ws_tEta_tK) * g5 * adj(ws_eta);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(gL[mu] * Lxx * gL[mu] * tmp);
  return rst;
}
    
LatticeComplex Hw_T3D1Q1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_K * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(g5 * wl_tEta_tEta ) * trace(gL[mu] * Lxx) * trace(gL[mu] * tmp );
  return rst;
}
    
LatticeComplex Hw_T3D1Q2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_K * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(g5 * wl_tEta_tEta ) * trace(gL[mu] * Lxx * gL[mu] * tmp );
  return rst;
}
    
LatticeComplex Hw_T3D2Q1(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_K * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(g5 * ws_tEta_tEta) * trace(gL[mu] * Lxx) * trace(gL[mu] * tmp );
  return rst;
}
    
LatticeComplex Hw_T3D2Q2(const LatticePropagator &Lxx, const LatticePropagator &wl_K, const LatticePropagator &ws_K, const LatticePropagator &wl_eta, const LatticePropagator &ws_eta, const LatticePropagatorSite &wl_tEta_tK, const LatticePropagatorSite &ws_tEta_tK, const LatticePropagatorSite &wl_tEta_tEta, const LatticePropagatorSite &ws_tEta_tEta) {
  LatticePropagator tmp = wl_K * adj(ws_K);
  tmp = tmp + adj(tmp); // Add contribution of K0 bar

  LatticeComplex rst(wl_K.Grid());  rst = Zero();
  for(int mu=0; mu<4; ++mu) rst += trace(g5 * ws_tEta_tEta) * trace(gL[mu] * Lxx * gL[mu] * tmp);
  return rst;
}


}

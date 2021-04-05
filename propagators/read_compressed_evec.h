#pragma once
#include <Grid/Grid.h>

namespace Grid {


static void cev_get_read_geometry(const GridBase* _grid,const std::vector<int>& cnodes,
				  std::map<int, std::vector<int> >& slots, 
				  std::vector<int>& slot_lvol,
				  std::vector<int>& lvol,
				  int64_t& slot_lsites,int& ntotal) {
  
  int _nd = (int)cnodes.size();
  std::vector<int> nodes = cnodes;
  
  slots.clear();
  slot_lvol.clear();
  lvol.clear();
  
  int i;
  ntotal = 1;
  int64_t lsites = 1;
  slot_lsites = 1;
  for (i=0;i<_nd;i++) {
    assert(_grid->_fdimensions[i] % nodes[i] == 0);
    slot_lvol.push_back(_grid->_fdimensions[i] / nodes[i]);
    lvol.push_back(_grid->_fdimensions[i] / _grid->_processors[i]);
    lsites *= lvol.back();
    slot_lsites *= slot_lvol.back();
    ntotal *= nodes[i];
  }
  
  std::vector<int> lcoor, gcoor, scoor;
  lcoor.resize(_nd); gcoor.resize(_nd);  scoor.resize(_nd);
  
  // create mapping of indices to slots
  for (int lidx = 0; lidx < lsites; lidx++) {
    Lexicographic::CoorFromIndex(lcoor,lidx,lvol);
    for (int i=0;i<_nd;i++) {
      gcoor[i] = lcoor[i] + _grid->_processor_coor[i]*lvol[i];
      scoor[i] = gcoor[i] / slot_lvol[i];
    }
    int slot;
    Lexicographic::IndexFromCoor(scoor,slot,nodes);
    auto sl = slots.find(slot);
    if (sl == slots.end())
      slots[slot] = std::vector<int>();
    slots[slot].push_back(lidx);
  }
}

static uint32_t cev_crc32(unsigned char* data, int64_t len, uint32_t previousCrc32 = 0) {
  
  // crc32 of zlib was incorrect for very large sizes, so do it block-wise
  uint32_t crc = previousCrc32;
  off_t blk = 0;
  off_t step = 1024*1024*1024;
  while (len > step) {
    crc = crc32(crc,&data[blk],step);
    blk += step;
    len -= step;
  }
  
  crc = crc32(crc,&data[blk],len);
  return crc;
  
}

template<typename Field>
class cev_BlockedGrid {
public:
  GridBase* _grid;
  typedef typename Field::scalar_type Coeff_t;
  typedef typename Field::vector_type vCoeff_t;
  
  Coordinate _bs; // block size
  Coordinate _nb; // number of blocks
  Coordinate _l;  // local dimensions irrespective of cb
  Coordinate _l_cb;  // local dimensions of checkerboarded vector
  Coordinate _l_cb_o;  // local dimensions of inner checkerboarded vector
  Coordinate _bs_cb; // block size in checkerboarded vector
  Coordinate _nb_o; // number of blocks of simd o-sites
  int _nd, _blocks, _cf_size, _cf_block_size, _cf_o_block_size, _o_blocks, _block_sites;
  
  cev_BlockedGrid(GridBase* grid, const std::vector<int>& block_size) :
    _grid(grid), _bs(block_size), _nd((int)_bs.size()), 
    _nb(block_size), _l(block_size), _l_cb(block_size), _nb_o(block_size),
    _l_cb_o(block_size), _bs_cb(block_size) {

    _blocks = 1;
    _o_blocks = 1;
    _l = grid->FullDimensions();
    _l_cb = grid->LocalDimensions();
    _l_cb_o = grid->_rdimensions;
    
    _cf_size = 1;
    _block_sites = 1;
    for (int i=0;i<_nd;i++) {
      _l[i] /= grid->_processors[i];
      
      std::cout << GridLogMessage << "_l[i] = " << _l[i] << ", _bs[i] = " << _bs[i] << std::endl;
      assert(!(_l[i] % _bs[i])); // lattice must accommodate choice of blocksize
      
      int r = _l[i] / _l_cb[i];
      assert(!(_bs[i] % r)); // checkerboarding must accommodate choice of blocksize
      _bs_cb[i] = _bs[i] / r;
      _block_sites *= _bs_cb[i];
      _nb[i] = _l[i] / _bs[i];
      _nb_o[i] = _nb[i] / _grid->_simd_layout[i];
      if (_nb[i] % _grid->_simd_layout[i]) { // simd must accommodate choice of blocksize
	std::cout << GridLogMessage << "Problem: _nb[" << i << "] = " << _nb[i] << " _grid->_simd_layout[" << i << "] = " << _grid->_simd_layout[i] << std::endl;
	assert(0);
      }
      _blocks *= _nb[i];
      _o_blocks *= _nb_o[i];
      _cf_size *= _l[i];
    }

    _cf_size *= 12 / 2;
    _cf_block_size = _cf_size / _blocks;
    _cf_o_block_size = _cf_size / _o_blocks;
    
    std::cout << GridLogMessage << "BlockedGrid:" << std::endl;
    std::cout << GridLogMessage << " _l     = " << _l << std::endl;
    std::cout << GridLogMessage << " _l_cb     = " << _l_cb << std::endl;
    std::cout << GridLogMessage << " _l_cb_o     = " << _l_cb_o << std::endl;
    std::cout << GridLogMessage << " _bs    = " << _bs << std::endl;
    std::cout << GridLogMessage << " _bs_cb    = " << _bs_cb << std::endl;

    std::cout << GridLogMessage << " _nb    = " << _nb << std::endl;
    std::cout << GridLogMessage << " _nb_o    = " << _nb_o << std::endl;
    std::cout << GridLogMessage << " _blocks = " << _blocks << std::endl;
    std::cout << GridLogMessage << " _o_blocks = " << _o_blocks << std::endl;
    std::cout << GridLogMessage << " sizeof(vCoeff_t) = " << sizeof(vCoeff_t) << std::endl;
    std::cout << GridLogMessage << " _cf_size = " << _cf_size << std::endl;
    std::cout << GridLogMessage << " _cf_block_size = " << _cf_block_size << std::endl;
    std::cout << GridLogMessage << " _block_sites = " << _block_sites << std::endl;
    std::cout << GridLogMessage << " _grid->oSites() = " << _grid->oSites() << std::endl;
    
    //    _grid->Barrier();
    //abort();
  }
  
  void block_to_coor(int b, std::vector<int>& x0) {
    
    Coordinate bcoor;
    bcoor.resize(_nd);
    x0.resize(_nd);
    assert(b < _o_blocks);
    Lexicographic::CoorFromIndex(bcoor,b,_nb_o);
    int i;
    
    for (i=0;i<_nd;i++) {
      x0[i] = bcoor[i]*_bs_cb[i];
    }
    
    //std::cout << GridLogMessage << "Map block b -> " << x0 << std::endl;
    
  }
  
  void block_site_to_o_coor(const std::vector<int>& x0, Coordinate& coor, int i) {
    Lexicographic::CoorFromIndex(coor,i,_bs_cb);
    for (int j=0;j<_nd;j++)
      coor[j] += x0[j];
  }
  
  int block_site_to_o_site(const std::vector<int>& x0, int i) {
    Coordinate coor;  coor.resize(_nd);
    block_site_to_o_coor(x0,coor,i);
    Lexicographic::IndexFromCoor(coor,i,_l_cb_o);
    return i;
  }
  
  vCoeff_t block_sp(int b, const Field& x, const Field& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    vCoeff_t ret = 0.0;
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      ret += TensorRemove(innerProduct(x._odata[ss],y._odata[ss]));
    }
    
    return ret;
    
  }
  
  vCoeff_t block_sp(int b, const Field& x, const std::vector< ComplexD >& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    constexpr int nsimd = sizeof(vCoeff_t) / sizeof(Coeff_t);
    int lsize = _cf_o_block_size / _block_sites;
    
    std::vector< ComplexD > ret(nsimd);
    for (int i=0;i<nsimd;i++)
      ret[i] = 0.0;
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      
      int n = lsize / nsimd;
      for (int l=0;l<n;l++) {
	for (int j=0;j<nsimd;j++) {
	  int t = lsize * i + l*nsimd + j;
	  
	  ret[j] += conjugate(((Coeff_t*)&x._odata[ss]._internal)[l*nsimd + j]) * y[t];
	}
      }
    }
    
    vCoeff_t vret;
    for (int i=0;i<nsimd;i++)
      ((Coeff_t*)&vret)[i] = (Coeff_t)ret[i];
    
    return vret;
    
  }
  
  
  // vector - lane operations sometimes needed for non-performance-critical code
  template<class T>
  void vcaxpy_lane(int l, iScalar<T>& r,const vCoeff_t& a,const iScalar<T>& x,const iScalar<T>& y) {
    vcaxpy_lane(l,r._internal,a,x._internal,y._internal);
  }
  
  template<class T,int N>
  void vcaxpy_lane(int l, iVector<T,N>& r,const vCoeff_t& a,const iVector<T,N>& x,const iVector<T,N>& y) {
    for (int i=0;i<N;i++)
      vcaxpy_lane(l,r._internal[i],a,x._internal[i],y._internal[i]);
  }
  
  void vcaxpy_lane(int l, vCoeff_t& r,const vCoeff_t& a,const vCoeff_t& x,const vCoeff_t& y) {
    ((Coeff_t*)&r)[l] = ((Coeff_t*)&a)[l]*((Coeff_t*)&x)[l] + ((Coeff_t*)&y)[l];
  }
  
  void block_lane_caxpy(int b, int l, Field& ret, const vCoeff_t& a, const Field& x, const Field& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      vcaxpy_lane(l,ret._odata[ss],a,x._odata[ss],y._odata[ss]);
    }
    
  }
  
  
  template<class T>
  void vcaxpy(iScalar<T>& r,const vCoeff_t& a,const iScalar<T>& x,const iScalar<T>& y) {
    vcaxpy(r._internal,a,x._internal,y._internal);
  }
  
  template<class T,int N>
  void vcaxpy(iVector<T,N>& r,const vCoeff_t& a,const iVector<T,N>& x,const iVector<T,N>& y) {
    for (int i=0;i<N;i++)
      vcaxpy(r._internal[i],a,x._internal[i],y._internal[i]);
  }
  
  void vcaxpy(vCoeff_t& r,const vCoeff_t& a,const vCoeff_t& x,const vCoeff_t& y) {
    r = a*x + y;
  }
  
  void block_caxpy(int b, Field& ret, const vCoeff_t& a, const Field& x, const Field& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    autoView(ret_v, ret, CpuWrite);
    autoView(x_v, x, CpuRead);
    autoView(y_v, y, CpuRead);
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      vcaxpy(ret_v[ss],a,x_v[ss],y_v[ss]);
    }
    
  }
  

  
  
  template<class T>
  void vcaxpy2(iScalar<T>& r,
	       const vCoeff_t& a1,const iScalar<T>& x1,
	       const vCoeff_t& a2,const iScalar<T>& x2,
	       const iScalar<T>& y) {
    vcaxpy2(r._internal,
	    a1,x1._internal,
	    a2,x2._internal,
	    y._internal);
  }
  
  template<class T,int N>
  void vcaxpy2(iVector<T,N>& r,
	       const vCoeff_t& a1,const iVector<T,N>& x1,
	       const vCoeff_t& a2,const iVector<T,N>& x2,
	       const iVector<T,N>& y) {
    for (int i=0;i<N;i++)
      vcaxpy2(r._internal[i],
	      a1,x1._internal[i],
	      a2,x2._internal[i],
	      y._internal[i]);
  }
  
  void vcaxpy2(vCoeff_t& r,
	       const vCoeff_t& a1,const vCoeff_t& x1,
	       const vCoeff_t& a2,const vCoeff_t& x2,
	       const vCoeff_t& y) {
    r = a1*x1 + a2*x2 + y;
  }
  
  void block_caxpy2(int b, Field& ret, 
		    const vCoeff_t& a1, const Field& x1, 
		    const vCoeff_t& a2, const Field& x2, 
		    const Field& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      vcaxpy2(ret._odata[ss],a1,x1._odata[ss],a2,x2._odata[ss],y._odata[ss]);
    }
    
  }
  

  
  
  template<class T>
  void vcaxpy4(iScalar<T>& r,
	       const vCoeff_t& a1,const iScalar<T>& x1,
	       const vCoeff_t& a2,const iScalar<T>& x2,
	       const vCoeff_t& a3,const iScalar<T>& x3,
	       const vCoeff_t& a4,const iScalar<T>& x4,
	       const iScalar<T>& y) {
    vcaxpy4(r._internal,
	    a1,x1._internal,
	    a2,x2._internal,
	    a3,x3._internal,
	    a4,x4._internal,
	    y._internal);
  }
  
  template<class T,int N>
  void vcaxpy4(iVector<T,N>& r,
	       const vCoeff_t& a1,const iVector<T,N>& x1,
	       const vCoeff_t& a2,const iVector<T,N>& x2,
	       const vCoeff_t& a3,const iVector<T,N>& x3,
	       const vCoeff_t& a4,const iVector<T,N>& x4,
	       const iVector<T,N>& y) {
    for (int i=0;i<N;i++)
      vcaxpy4(r._internal[i],
	      a1,x1._internal[i],
	      a2,x2._internal[i],
	      a3,x3._internal[i],
	      a4,x4._internal[i],
	      y._internal[i]);
  }
  
  void vcaxpy4(vCoeff_t& r,
	       const vCoeff_t& a1,const vCoeff_t& x1,
	       const vCoeff_t& a2,const vCoeff_t& x2,
	       const vCoeff_t& a3,const vCoeff_t& x3,
	       const vCoeff_t& a4,const vCoeff_t& x4,
	       const vCoeff_t& y) {
    r = a1*x1 + a2*x2 + a3*x3 + a4*x4 + y;
  }
  
  void block_caxpy4(int b, Field& ret, 
		    const vCoeff_t& a1, const Field& x1, 
		    const vCoeff_t& a2, const Field& x2, 
		    const vCoeff_t& a3, const Field& x3, 
		    const vCoeff_t& a4, const Field& x4, 
		    const Field& y) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    autoView(ret_v, ret, CpuWrite);
    autoView(x1_v, x1, CpuRead);
    autoView(x2_v, x2, CpuRead);
    autoView(x3_v, x3, CpuRead);
    autoView(x4_v, x4, CpuRead);
    autoView(y_v, y, CpuRead);
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      vcaxpy4(ret_v[ss],
	      a1,x1_v[ss],
	      a2,x2_v[ss],
	      a3,x3_v[ss],
	      a4,x4_v[ss],
	      y_v[ss]);
    }
    
  }
  
  void block_caxpy(int b, std::vector< ComplexD >& ret, const vCoeff_t& a, const Field& x, const std::vector< ComplexD >& y) {
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    constexpr int nsimd = sizeof(vCoeff_t) / sizeof(Coeff_t);
    int lsize = _cf_o_block_size / _block_sites;
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      
      int n = lsize / nsimd;
      for (int l=0;l<n;l++) {
	vCoeff_t r = a* ((vCoeff_t*)&x._odata[ss]._internal)[l];
	
	for (int j=0;j<nsimd;j++) {
	  int t = lsize * i + l*nsimd + j;
	  ret[t] = y[t] + ((Coeff_t*)&r)[j];
	}
      }
    }
    
  }
  
  void block_set(int b, Field& ret, const std::vector< ComplexD >& x) {
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    int lsize = _cf_o_block_size / _block_sites;
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      
      for (int l=0;l<lsize;l++)
	((Coeff_t*)&ret._odata[ss]._internal)[l] = (Coeff_t)x[lsize * i + l]; // convert precision
    }
    
  }
  
  void block_get(int b, const Field& ret, std::vector< ComplexD >& x) {
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    int lsize = _cf_o_block_size / _block_sites;
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      
      for (int l=0;l<lsize;l++)
	x[lsize * i + l] = (ComplexD)((Coeff_t*)&ret._odata[ss]._internal)[l];
    }
    
  }
  
  template<class T>
  void vcscale(iScalar<T>& r,const vCoeff_t& a,const iScalar<T>& x) {
    vcscale(r._internal,a,x._internal);
  }
  
  template<class T,int N>
  void vcscale(iVector<T,N>& r,const vCoeff_t& a,const iVector<T,N>& x) {
    for (int i=0;i<N;i++)
      vcscale(r._internal[i],a,x._internal[i]);
  }
  
  void vcscale(vCoeff_t& r,const vCoeff_t& a,const vCoeff_t& x) {
    r = a*x;
  }
  
  void block_cscale(int b, const vCoeff_t& a, Field& ret) {
    
    std::vector<int> x0;
    block_to_coor(b,x0);
    
    for (int i=0;i<_block_sites;i++) { // only odd sites
      int ss = block_site_to_o_site(x0,i);
      vcscale(ret._odata[ss],a,ret._odata[ss]);
    }
  }
  
  void getCanonicalBlockOffset(int cb, Coordinate& x0) {
    const int ndim = 5;
    assert(_nb.size() == ndim);
    Coordinate _nbc = Coordinate({ _nb[1], _nb[2], _nb[3], _nb[4], _nb[0] });
    Coordinate _bsc = Coordinate({ _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] });
    x0.resize(ndim);
    
    assert(cb >= 0);
    assert(cb < _nbc[0]*_nbc[1]*_nbc[2]*_nbc[3]*_nbc[4]);
    
    Lexicographic::CoorFromIndex(x0,cb,_nbc);
    int i;
    
    for (i=0;i<ndim;i++) {
      x0[i] *= _bsc[i];
    }
    
    //if (cb < 2)
    //	std::cout << GridLogMessage << "Map: " << cb << " To: " << x0 << std::endl;
  }

  void pokeBlockOfVectorCanonical(int cb,Field& v,const std::vector<float>& buf) {
    Coordinate _bsc = Coordinate({ _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] });
    Coordinate ldim = v.Grid()->LocalDimensions();
    Coordinate cldim = Coordinate({ ldim[1], ldim[2], ldim[3], ldim[4], ldim[0] });
    const int _nbsc = _bs_cb[0]*_bs_cb[1]*_bs_cb[2]*_bs_cb[3]*_bs_cb[4];
    // take canonical block cb of v and put it in canonical ordering in buf
    Coordinate cx0;
    getCanonicalBlockOffset(cb,cx0);
    
  autoView(v_view, v, CpuRead);
#pragma omp parallel
    {
      Coordinate co0,cl0;
      co0=cx0; cl0=cx0;
      
#pragma omp for
      for (int i=0;i<_nbsc;i++) {
	Lexicographic::CoorFromIndex(co0,2*i,_bsc); // 2* for eo
	for (int j=0;j<(int)_bsc.size();j++)
	  cl0[j] = cx0[j] + co0[j];
	
	Coordinate l0 = Coordinate({ cl0[4], cl0[0], cl0[1], cl0[2], cl0[3] });
	int oi = v.Grid()->oIndex(l0);
	int ii = v.Grid()->iIndex(l0);
	int lti = i;
	
	//if (cb < 2 && i<2)
	//  std::cout << GridLogMessage << "Map: " << cb << ", " << i << " To: " << cl0 << ", " << cx0 << ", " << oi << ", " << ii << std::endl;
	
	for (int s=0;s<4;s++)
	  for (int c=0;c<3;c++) {
	    Coeff_t& ld = ((Coeff_t*)&v_view[oi]._internal._internal[s]._internal[c])[ii];
	    int ti = 12*lti + 3*s + c;
	    ld = Coeff_t(buf[2*ti+0], buf[2*ti+1]);
	  }
      }
    }
  }
  
//   void peekBlockOfVectorCanonical(int cb,const Field& v,std::vector<float>& buf) {
//     Coordinate _bsc = Coordinate({ _bs[1], _bs[2], _bs[3], _bs[4], _bs[0] });
//     Coordinate ldim = v.Grid()->LocalDimensions();
//     Coordinate cldim = Coordinate({ ldim[1], ldim[2], ldim[3], ldim[4], ldim[0] });
//     const int _nbsc = _bs_cb[0]*_bs_cb[1]*_bs_cb[2]*_bs_cb[3]*_bs_cb[4];
//     // take canonical block cb of v and put it in canonical ordering in buf
//     Coordinate cx0;
//     getCanonicalBlockOffset(cb,cx0);
//
//     buf.resize(_cf_block_size * 2);
//
// #pragma omp parallel
//     {
//       Coordinate co0,cl0;
//       co0=cx0; cl0=cx0;
//
// #pragma omp for
//       for (int i=0;i<_nbsc;i++) {
// 	Lexicographic::CoorFromIndex(co0,2*i,_bsc); // 2* for eo
// 	for (int j=0;j<(int)_bsc.size();j++)
// 	  cl0[j] = cx0[j] + co0[j];
//
// 	std::vector<int> l0 = { cl0[4], cl0[0], cl0[1], cl0[2], cl0[3] };
// 	int oi = v.Grid()->oIndex(l0);
// 	int ii = v.Grid()->iIndex(l0);
// 	int lti = i;
//
// 	//if (cb < 2 && i<2)
// 	//  std::cout << GridLogMessage << "Map: " << cb << ", " << i << " To: " << cl0 << ", " << cx0 << ", " << oi << ", " << ii << std::endl;
//
// 	for (int s=0;s<4;s++)
// 	  for (int c=0;c<3;c++) {
// 	    Coeff_t& ld = ((Coeff_t*)&v._odata[oi]._internal._internal[s]._internal[c])[ii];
// 	    int ti = 12*lti + 3*s + c;
// 	    buf[2*ti+0] = ld.real();
// 	    buf[2*ti+1] = ld.imag();
// 	  }
//       }
//     }
//   }
  
  int globalToLocalCanonicalBlock(int slot,const std::vector<int>& src_nodes,int nb) {
    // processor coordinate
    int _nd = (int)src_nodes.size();
    std::vector<int> _src_nodes = src_nodes;
    std::vector<int> pco(_nd);
    Lexicographic::CoorFromIndex(pco,slot,_src_nodes);
    std::vector<int> cpco = { pco[1], pco[2], pco[3], pco[4], pco[0] };
    
    // get local block
    std::vector<int> _nbc = { _nb[1], _nb[2], _nb[3], _nb[4], _nb[0] };
    assert(_nd == 5);
    std::vector<int> c_src_local_blocks(_nd);
    for (int i=0;i<_nd;i++) {
      assert(_grid->_fdimensions[i] % (src_nodes[i] * _bs[i]) == 0);
      c_src_local_blocks[(i+4) % 5] = _grid->_fdimensions[i] / src_nodes[i] / _bs[i];
    }
    std::vector<int> cbcoor(_nd); // coordinate of block in slot in canonical form
    Lexicographic::CoorFromIndex(cbcoor,nb,c_src_local_blocks);
    
    // cpco, cbcoor
    std::vector<int> clbcoor(_nd);
    for (int i=0;i<_nd;i++) {
      int cgcoor = cpco[i] * c_src_local_blocks[i] + cbcoor[i]; // global block coordinate
      int pcoor = cgcoor / _nbc[i]; // processor coordinate in my Grid
      int tpcoor = _grid->_processor_coor[(i+1)%5];
      if (pcoor != tpcoor)
	return -1;
      clbcoor[i] = cgcoor - tpcoor * _nbc[i]; // canonical local block coordinate for canonical dimension i
    }
    
    int lnb;
    Lexicographic::IndexFromCoor(clbcoor,lnb,_nbc);
    //std::cout << "Mapped slot = " << slot << " nb = " << nb << " to " << lnb << std::endl;
    return lnb;
  }
  
};



template<class Field, class Allocator = std::allocator<Field> >
class cev_BasisFieldVector {
public:
  int _Nm;
  
  typedef typename Field::scalar_type Coeff_t;
  typedef typename Field::vector_type vCoeff_t;
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  
  std::vector<Field, Allocator> _v; // _Nfull vectors
  
  void report(int n,GridBase* value) {
    
    std::cout << GridLogMessage << "BasisFieldVector allocated:\n";
    std::cout << GridLogMessage << " Delta N = " << n << "\n";
    std::cout << GridLogMessage << " Size of full vectors (size) = " << 
    ((double)n*sizeof(vobj)*value->oSites() / 1024./1024./1024.) << " GB\n"; // zyd: This is memory used by each process; I used 4 processes per node
    std::cout << GridLogMessage << " Size = " << _v.size() << " Capacity = " << _v.capacity() << std::endl;
    
    value->Barrier();
    
    // if (value->IsBoss()) {
    //   system("cat /proc/meminfo");
    // }
    
    value->Barrier();
    
  }
  
  cev_BasisFieldVector(int Nm,GridBase* value, const Allocator& a = std::allocator<Field>() ) : _Nm(Nm), _v(Nm,value,a) {
    report(Nm,value);
  }
  
  cev_BasisFieldVector(const Allocator& a = std::allocator<Field>() ) : _Nm(0), _v(a) {
  }
  
  ~cev_BasisFieldVector() {
  }
  
  Field& operator[](int i) {
    return _v[i];
  }
  
  void orthogonalize(Field& w, int k) {
    
    for(int j=0; j<k; ++j){
      Coeff_t ip = (Coeff_t)innerProduct(_v[j],w);
      w = w - ip*_v[j];
    }
    
  }

  size_t size() const {
    return _Nm;
  }
  
  void resize(int n, GridBase* value) {
    if (n > _Nm)
      _v.reserve(n);
    
    _v.resize(n,value);

    if (n < _Nm)
      _v.shrink_to_fit();

    report(n - _Nm,value);

    _Nm = n;    
  }

  void resize(int n) {
    assert(_v.size() > 0);
    resize(n,_v[0].Grid());
  }

  void deflate(const std::vector<RealD>& eval,const Field& src_orig,Field& result) {
    result = zero;
    int N = (int)_v.size();
    for (int i=0;i<N;i++) {
      Field& tmp = _v[i];
      axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
    }
  }
  
}; 

template<typename Field>
class cev_BlockProjector {
public:

  typedef typename Field::scalar_type Coeff_t;
  typedef typename Field::vector_type vCoeff_t;

  cev_BasisFieldVector<Field>& _evec;
  cev_BlockedGrid<Field>& _bgrid;

  cev_BlockProjector(cev_BasisFieldVector<Field>& evec, cev_BlockedGrid<Field>& bgrid) : _evec(evec), _bgrid(bgrid) {
  }

  void removeZeroBlocks(RealD thres = 1e-7) {
    GridStopWatch sw;
    sw.Start();

    constexpr int nsimd = sizeof(vCoeff_t) / sizeof(Coeff_t);

    int cnt = 0;

    Field tmpzero(_evec._v[0].Grid());
    tmpzero = zero;

#pragma omp parallel
    {
      int lcnt = 0;

#pragma omp for
      for (int b=0;b<_bgrid._o_blocks;b++) {

	for (int l=0;l<nsimd;l++) {

	  int j=0;
	  for (int i=0;i<_evec._Nm;i++) {
	  
	    // i == source, j == dest
	    auto nrm0 = _bgrid.block_sp(b,_evec._v[i],_evec._v[i]);
	    Coeff_t* n = (Coeff_t*)&nrm0;

	    if (n[l].real() > thres) {
	      // v[j] = v[i]
	      if (i!=j)
		_bgrid.block_lane_caxpy(b,l,_evec._v[j],1.0,_evec._v[i],tmpzero);
	      
	      j++;
	    } else {
	      lcnt++;
	    }
	  }
	  
	  // now set blocks from j to Nm to zero
	  for (int i=j;i<_evec._Nm;i++)
	    _bgrid.block_lane_caxpy(b,l,_evec._v[i],0.0,tmpzero,tmpzero);	
	}
      }

#pragma omp critical
      {
	cnt += lcnt;
      }
    }
    sw.Stop();
    std::cout << GridLogMessage << "Removed zero blocks in " << sw.Elapsed() << " (" << ((RealD)cnt / (RealD)_bgrid._blocks / (RealD)_evec._Nm) 
	      << " below threshold)" << std::endl;
  }

  template<typename CoarseField>
  void coarseToFine(const CoarseField& in, Field& out) {

    out = Zero();
    out.Checkerboard() = _evec._v[0].Checkerboard();

    autoView(in_v, in, CpuRead);
    int Nbasis = sizeof(in_v[0]._internal._internal) / sizeof(in_v[0]._internal._internal[0]);
    assert(Nbasis == _evec._Nm);

    // autoView(in_v, in, CpuRead);
#pragma omp parallel for
    for (int b=0;b<_bgrid._o_blocks;b++) {

      int nw = _evec._Nm / 4;

      for (int j=0;j<4*nw;j+=4) {
	_bgrid.block_caxpy4(b,out,
			    in_v[b]._internal._internal[j+0],_evec._v[j+0],
			    in_v[b]._internal._internal[j+1],_evec._v[j+1],
			    in_v[b]._internal._internal[j+2],_evec._v[j+2],
			    in_v[b]._internal._internal[j+3],_evec._v[j+3],
			    out);
      }

      for (int j=4*nw;j<_evec._Nm;j++) {
	_bgrid.block_caxpy(b,out,in_v[b]._internal._internal[j],_evec._v[j],out);
      }

    }

    //std::cout << GridLogMessage << "Check cTF norm2(in) = " << norm2(in) << " norm2(out) = " << norm2(out) << std::endl;
  }

  template<typename CoarseField>
  void fineToCoarse(const Field& in, CoarseField& out) {

    out = zero;

    int Nbasis = sizeof(out._odata[0]._internal._internal) / sizeof(out._odata[0]._internal._internal[0]);
    assert(Nbasis == _evec._Nm);


    Field tmp(_bgrid.Grid());
    tmp = in;
    
    GridStopWatch gsw;
    gsw.Start();
#pragma omp parallel for
    for (int b=0;b<_bgrid._o_blocks;b++) {
      for (int j=0;j<_evec._Nm;j++) {
	// |rhs> -= <j|rhs> |j>
	auto c = _bgrid.block_sp(b,_evec._v[j],tmp);
	_bgrid.block_caxpy(b,tmp,-c,_evec._v[j],tmp); // may make this more numerically stable
	out._odata[b]._internal._internal[j] = c;
      }
    }
    gsw.Stop();

    std::cout<<GridLogMessage<< "Timing fineToCoarse: " << gsw.Elapsed() << std::endl;
    //std::cout << GridLogMessage << "Check fTC norm2(in) = " << norm2(in) << " norm2(out) = " << norm2(out) << std::endl;
  }

  template<typename CoarseField>
    void deflateFine(cev_BasisFieldVector<CoarseField>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    result = zero;
    for (int i=0;i<N;i++) {
      Field tmp(result.Grid());
      coarseToFine(_coef._v[i],tmp);
      axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
    }
  }

  template<typename CoarseField,typename Allocator>
    void deflateCoarse(cev_BasisFieldVector<CoarseField,Allocator>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    CoarseField src_coarse(_coef._v[0].Grid());
    CoarseField result_coarse = src_coarse;
    result_coarse = zero;
    fineToCoarse(src_orig,src_coarse);

    //std::cout << GridLogMessage << "N: " << N << std::endl;
    //std::cout << GridLogMessage << "Eval: " << eval << std::endl;

    for (int i=0;i<N;i++) {
      //std::cout << GridLogMessage << "norm2(v[i]) = " << norm2(_coef._v[i]) << std::endl;
      axpy(result_coarse,TensorRemove(innerProduct(_coef._v[i],src_coarse)) / eval[i],_coef._v[i],result_coarse);
    }
    coarseToFine(result_coarse,result);
  }

  template<typename CoarseField,typename Allocator>
    void deflate(cev_BasisFieldVector<CoarseField,Allocator>& _coef,const std::vector<RealD>& eval,int N,const Field& src_orig,Field& result) {
    // Deflation on coarse Grid is much faster, so use it by default.  Deflation on fine Grid is kept for legacy reasons for now.
    deflateCoarse(_coef,eval,N,src_orig,result);
  }

};

#define SHRT_UMAX 65535
#define FP16_BASE 1.4142135623730950488
#define FP16_COEF_EXP_SHARE_FLOATS 10
   static float unmap_fp16_exp(unsigned short e) {
     float de = (float)((int)e - SHRT_UMAX / 2);
     return ::pow( FP16_BASE, de );
   }

static void read_floats(char* & ptr, float* out, int64_t n) {
  float* in = (float*)ptr;
  ptr += 4*n;

  for (int64_t i=0;i<n;i++)
    out[i] = in[i];
}

template<typename OPT>
static void read_floats_fp16(char* & ptr, OPT* out, int64_t n, int nsc) {

  int64_t nsites = n / nsc;
  if (n % nsc) {
    fprintf(stderr,"Invalid size in write_floats_fp16\n");
    exit(4);
  }

  unsigned short* in = (unsigned short*)ptr;
  ptr += 2*(n+nsites);

  // do for each site
  for (int64_t site = 0;site<nsites;site++) {

    OPT* ev = &out[site*nsc];

    unsigned short* bptr = &in[site*(nsc + 1)];

    unsigned short exp = *bptr++;
    OPT max = unmap_fp16_exp(exp);
    OPT min = -max;

    for (int i=0;i<nsc;i++) {
      ev[i] = fp_unmap( *bptr++, min, max, SHRT_UMAX );
    }

  }

}

   static void canonical_block_to_coarse_coordinates(GridBase* _coarsegrid,int nb,int& ii,int& oi) {
      // canonical nb needs to be mapped in a coordinate on my coarsegrid (ii,io)
      Coordinate _l = _coarsegrid->LocalDimensions();
      Coordinate _cl = Coordinate({ _l[1], _l[2], _l[3], _l[4], _l[0] });
      Coordinate _cc(_l.size());
      Lexicographic::CoorFromIndex(_cc,nb,_cl);
      Coordinate _c = Coordinate({ _cc[4], _cc[0], _cc[1], _cc[2], _cc[3] });
      ii = _coarsegrid->iIndex(_c);
      oi = _coarsegrid->oIndex(_c);
    }


template<typename vtype, int N > using cev_CoarseSiteFieldGeneral = iScalar< iVector<vtype, N> >;
template<int N> using cev_CoarseSiteFieldF = cev_CoarseSiteFieldGeneral< vComplexF, N >;
template<int N> using cev_CoarseLatticeFermionF = Lattice< cev_CoarseSiteFieldF<N> >;

static float fp_unmap(int val, float min, float max, int N) {
  return min + (float)(val + 0.5) * (max - min)  / (float)( N + 1 );
}

template<int Nbasis>
bool uncompress_evec_from_directory_basis(const char* dir, std::vector<double>& eval,
					  std::vector<LatticeFermionF>& uncompressed_evec,
					  int ngroups, GridBase* _grid, int neig,   
					  std::vector<uint32_t> & crc32, std::vector<uint32_t>& b,
					  std::vector<uint32_t> & nn, int nkeep_single, int blocks,
					  int _FP16_COEF_EXP_SHARE_FLOATS) {

  FILE* f;
  int nkeep = Nbasis;

  char hostname[1024];
  gethostname(hostname, 1024);

  if (_grid->IsBoss())
    std::cout << GridLogMessage << "Resize to " << neig << std::endl;

  eval.resize(neig);
  
  _grid -> show_decomposition();
  std::cout << "Before resizing evecs" << std::endl; print_memory();
  uncompressed_evec.resize(neig, _grid);
  std::cout << "Capacity of uncompressed_evec after resize: " << uncompressed_evec.capacity() << std::endl;
  std::cout << "After resizing evecs" << std::endl; print_memory();

  // infrastructure
  std::vector<int> block_size(5);
  assert( b.size() == block_size.size() );
  assert( b.size() == 5 );
  block_size[0] = b[0];
  block_size[1] = b[1];
  block_size[2] = b[2];
  block_size[3] = b[3];
  block_size[4] = b[4];

  std::cout << "Before resizing cev_BasisFieldVector 1" << std::endl; print_memory();
  cev_BasisFieldVector<LatticeFermionF> evec(nkeep, _grid);
  std::cout << "After resizing cev_BasisFieldVector 1" << std::endl; print_memory();

  cev_BlockedGrid<LatticeFermionF>      bgrid(_grid,block_size);
  cev_BlockProjector<LatticeFermionF>   pr(evec,bgrid);
  std::vector<int> coarseFourDimLatt;
  for (int i=0;i<4;i++) {
    coarseFourDimLatt.push_back(bgrid._nb[1+i] * bgrid._grid->_processors[1+i]);
    if (_grid->IsBoss())
      std::cout << GridLogMessage << "CoarseGrid[" << i << "] = " << coarseFourDimLatt[i] << std::endl;
  }
  assert(bgrid._grid->_processors[0] == 1);
  GridCartesian* UCoarseGridF = SpaceTimeGrid::makeFourDimGrid(coarseFourDimLatt, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridCartesian* FCoarseGridF = SpaceTimeGrid::makeFiveDimGrid(bgrid._nb[0],UCoarseGridF);
// #define Nbasis 100
  std::cout << "Before resizing cev_BasisFieldVector 2" << std::endl; print_memory();
  cev_BasisFieldVector<cev_CoarseLatticeFermionF<Nbasis> > coef(neig,FCoarseGridF);
  std::cout << "After resizing cev_BasisFieldVector 2" << std::endl; print_memory();

  // now get read geometry
  std::map<int, std::vector<int> > slots;
  std::vector<int> slot_lvol, lvol;
  int64_t slot_lsites;
  int ntotal;
  std::vector<int> _nn(nn.begin(),nn.end());
  cev_get_read_geometry(_grid,_nn,
			slots,slot_lvol,lvol,slot_lsites,
			ntotal);
  int _nd = (int)lvol.size();

  // slot layout
  int nperdir = ntotal / 32;
  if (nperdir < 1)
    nperdir=1;
  
  // add read groups
  for (int ngroup=0;ngroup<ngroups;ngroup++) { // the ngroups divides the read tasks into ngroups by mpi_rank % ngroups == groupid .
    
    bool action = _grid->ThisRank() % ngroups == ngroup;
    
    std::cout << GridLogMessage << "Reading in group " << ngroup << " / " << ngroups << std::endl;
    
    // load all necessary slots and store them appropriately
    for (auto sl=slots.begin();sl!=slots.end();sl++) {
      
      std::vector<int>& idx = sl->second;
      int slot = sl->first;
      std::vector<float> rdata;
       
      char buf[4096];
      
      if (action) {
	// load one slot vector
	sprintf(buf,"%s/%2.2d/%10.10d.compressed",dir,slot/nperdir,slot);
	f = fopen(buf,"rb");
	if (!f) {
	  fprintf(stderr,"Node %s cannot read %s\n",hostname,buf); fflush(stderr);
	  return false;
	}
      }
      
      uint32_t crc = 0x0;
      off_t size;
	   
      GridStopWatch gsw;
      _grid->Barrier();
      gsw.Start();
      
      std::vector<char> raw_in(0);
      if (action) {
	fseeko(f,0,SEEK_END);
	size = ftello(f);
	fseeko(f,0,SEEK_SET);
	
  std::cout << "Before raw_in.resize" << std::endl; print_memory();
  std::cout << "size: " << size << std::endl;  // Since 00000000.compressed is 16GB, size should about around 17179869184
  std::cout << ">>>>>>> raw_in size is " << size/1024./1024./1024. << " GB" << std::endl;
	raw_in.resize(size);
  std::cout << "After raw_in.resize" << std::endl; print_memory();

	assert(fread(&raw_in[0],size,1,f) == 1);
      }
      
      _grid->Barrier();
      gsw.Stop();
      
      RealD totalGB = (RealD)size / 1024./1024./1024 * _grid->_Nprocessors;
      RealD seconds = gsw.useconds() / 1e6;
      
      if (action) {
	std::cout << GridLogMessage << "[" << slot << "]  Read " << totalGB << " GB of compressed data at " << totalGB/seconds << " GB/s" << std::endl;
	
	uint32_t crc_comp = cev_crc32((unsigned char*)&raw_in[0],size,0);
	
	if (crc_comp != crc32[slot]) {
	  std::cout << "Node " << hostname << " found crc mismatch for file " << buf << " (" << std::hex << crc_comp << " vs " << crc32[slot] << std::dec << ")" << std::endl;
	  std::cout << "Byte size: " << size << std::endl;
	}
	
	assert(crc_comp == crc32[slot]);
      }
      
      _grid->Barrier();
      
      if (action) {
	fclose(f);
      }
      
      char* ptr = &raw_in[0];
      
      GridStopWatch gsw2;
      gsw2.Start();

      typedef typename LatticeFermionF::scalar_type Coeff_t;
      typedef typename cev_CoarseLatticeFermionF<Nbasis>::scalar_type CoeffCoarse_t;

      if (action) {
	int nsingleCap = nkeep_single;
	if (nkeep < nsingleCap)
	  nsingleCap = nkeep;
	
	int _cf_block_size = slot_lsites * 12 / 2 / blocks;
	
#define FP_16_SIZE(a,b)  (( (a) + (a/b) )*2)

	// first read single precision basis vectors
#pragma omp parallel
	{
	  std::vector<float> buf(_cf_block_size * 2);
#pragma omp for
	  for (int nb=0;nb<blocks;nb++) {
	    for (int i=0;i<nsingleCap;i++) {
	      char* lptr = ptr + buf.size()*(i + nsingleCap*nb)*4;
	      read_floats(lptr, &buf[0], buf.size() );
	      int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
	      if (mnb != -1)
		pr._bgrid.pokeBlockOfVectorCanonical(mnb,pr._evec._v[i],buf);
	    }
	  }
	  
#pragma omp barrier
#pragma omp single
	  {
	    ptr = ptr + buf.size()*nsingleCap*blocks*4;
	  }
	  
	}

	// TODO: at this point I should add a checksum test for block_sp(nb,v,v) for all blocks, then I would know that the mapping
	// to blocks is OK at this point; after that ...
	
	// then read fixed precision basis vectors
#pragma omp parallel
	{
	  std::vector<float> buf(_cf_block_size * 2);
#pragma omp for
	  for (int nb=0;nb<blocks;nb++) {
	    for (int i=nsingleCap;i<(int)pr._evec.size();i++) {
	      char* lptr = ptr + FP_16_SIZE( buf.size(), 24 )*((i-nsingleCap) + (pr._evec.size() - nsingleCap)*nb);
	      read_floats_fp16(lptr, &buf[0], buf.size(), 24);
	      int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
	      if (mnb != -1)
		pr._bgrid.pokeBlockOfVectorCanonical(mnb,pr._evec._v[i],buf);
	    }
	  }
	  
#pragma omp barrier
#pragma omp single
	  {
	    ptr = ptr + FP_16_SIZE( buf.size()*(pr._evec.size() - nsingleCap)*blocks, 24 );
	  }
	}
	
#pragma omp parallel
	{
	  std::vector<float> buf1(nkeep_single*2);
	  std::vector<float> buf2((nkeep - nkeep_single)*2);
	  
#pragma omp for
	  for (int j=0;j<(int)coef.size();j++) {
      autoView(coef_v_view, coef._v[j], CpuRead);

	    for (int nb=0;nb<blocks;nb++) {
	      // get local coordinate on coarse grid
	      int ii,oi;
	      int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
	      if (mnb != -1)
		canonical_block_to_coarse_coordinates(coef._v[0].Grid(),mnb,ii,oi);
	      
	      char* lptr = ptr + (4*buf1.size() + FP_16_SIZE(buf2.size(), _FP16_COEF_EXP_SHARE_FLOATS))*(nb + j*blocks);
	      int l;
	      read_floats(lptr, &buf1[0], buf1.size() );
	      if (mnb != -1) {
		for (l=0;l<nkeep_single;l++) {

		  ((CoeffCoarse_t*)&coef_v_view[oi]._internal._internal[l])[ii] = CoeffCoarse_t(buf1[2*l+0],buf1[2*l+1]);
		}
	      }
	      read_floats_fp16(lptr, &buf2[0], buf2.size(), _FP16_COEF_EXP_SHARE_FLOATS);
	      if (mnb != -1) {
		for (l=nkeep_single;l<nkeep;l++) {
		  ((CoeffCoarse_t*)&coef_v_view[oi]._internal._internal[l])[ii] = CoeffCoarse_t(buf2[2*(l-nkeep_single)+0],buf2[2*(l-nkeep_single)+1]);
		}
	      }
	      
	    }
    }
	}
	
	// set checkerboard
	for (int i=0;i<(int)pr._evec.size();i++)
	  pr._evec._v[i].Checkerboard() = Odd;
	
	gsw2.Stop();
	seconds=gsw2.useconds()/1e6;
	std::cout << GridLogMessage << "Processed " << totalGB << " GB of compressed data at " << totalGB/seconds << " GB/s" << std::endl;
      }
    }
  }
#undef FP_16_SIZE

  // std::cout << GridLogMessage << "coef[0] = " << coef[0] << std::endl;
  //
  // std::cout << GridLogMessage << "basis[0] = " << pr._evec._v[0] << std::endl;

  for (int i=0;i<neig;i++) {
    if(i%100 == 0) {
      std::cout << GridLogMessage << "Uncompress " << i  << "/" << neig << std::endl;
    }
    pr.coarseToFine(coef[i],uncompressed_evec[i]);
  }
  std::cout << GridLogMessage << "Finished uncompressing eigenvectors" << std::endl;
  return true;
}

  bool read_evals(GridBase* _grid, const char* fn, std::vector<RealD>& evals) {
    
    FILE* f = 0;
    uint32_t status = 0;
    if (_grid->IsBoss()) {
      f = fopen(fn,"rt");
      status = f ? 1 : 0;
    }
    _grid->GlobalSum(status);
    
    if (!status)
      return false;
    
    uint32_t N;
    if (f)
      assert(fscanf(f,"%d\n",&N)==1);
    else
      N = 0;
    _grid->GlobalSum(N);
    
    std::cout << "Reading " << N << " eigenvalues" << std::endl;
    
    evals.resize(N);
    
    for (int i=0;i<N;i++) {
      if (f)
	assert(fscanf(f,"%lf",&evals[i])==1);
      else
	evals[i] = 0;
    }
    
    _grid->GlobalSumVector(&evals[0],evals.size());
    
    if (f)
      fclose(f);
    return true;
  }


bool uncompress_evec_from_directory(const char* dir, std::vector<double>& eval,
				    std::vector<LatticeFermionF>& uncompressed_evec, int ngroups = 1) {

  assert(uncompressed_evec.size()); // need as input at least one evec from which we take the Grid

  GridBase* _grid = uncompressed_evec[0].Grid();
  //const BasisFieldVector<Field>& basis = pr._evec;

   // first read metadata
  char buf[4096];
  sprintf(buf,"%s/metadata.txt",dir);

  std::vector<uint32_t> s,b,nb,nn,crc32;
  s.resize(5); b.resize(5); nb.resize(5); nn.resize(5);
  uint32_t neig, nkeep, nkeep_single, blocks, _FP16_COEF_EXP_SHARE_FLOATS;
  uint32_t nprocessors = 1;

  FILE* f = 0;
  uint32_t status = 0;
  if (_grid->IsBoss()) {
    f = fopen(buf,"rb");
    status=f ? 1 : 0;
  }
  _grid->GlobalSum(status);
  std::cout << GridLogMessage << "Read params status " << status << std::endl;

  if (!status) {
    return false;
  }

  sprintf(buf,"%s/eigen-values.txt",dir);
  read_evals(_grid,buf,eval);

#define _IRL_READ_INT(buf,p) if (f) { assert(fscanf(f,buf,p)==1); } else { *(p) = 0; } _grid->GlobalSum(*(p));

  for (int i=0;i<5;i++) {
    sprintf(buf,"s[%d] = %%d\n",i);
    _IRL_READ_INT(buf,&s[(i+1)%5]);
  }
  for (int i=0;i<5;i++) {
    sprintf(buf,"b[%d] = %%d\n",i);
    _IRL_READ_INT(buf,&b[(i+1)%5]);
  }
  for (int i=0;i<5;i++) {
    sprintf(buf,"nb[%d] = %%d\n",i);
    _IRL_READ_INT(buf,&nb[(i+1)%5]);
  }
  _IRL_READ_INT("neig = %d\n",&neig);
  _IRL_READ_INT("nkeep = %d\n",&nkeep);
  _IRL_READ_INT("nkeep_single = %d\n",&nkeep_single);
  _IRL_READ_INT("blocks = %d\n",&blocks);
  _IRL_READ_INT("FP16_COEF_EXP_SHARE_FLOATS = %d\n",&_FP16_COEF_EXP_SHARE_FLOATS);
  
  for (int i=0;i<5;i++) {
    assert(_grid->FullDimensions()[i] % s[i] == 0);
    nn[i] = _grid->FullDimensions()[i] / s[i];
    nprocessors *= nn[i];
  }

  std::cout << GridLogMessage << "Reading data that was generated on node-layout " << nn << std::endl;
  
  crc32.resize(nprocessors);
  for (int i =0;i<nprocessors;i++) {
    sprintf(buf,"crc32[%d] = %%X\n",i);
    _IRL_READ_INT(buf,&crc32[i]);
  }
  
#undef _IRL_READ_INT

  if (f)
    fclose(f);

  if (nkeep == 1000) {  // Change this number and the template argument to whatever number of 
    return uncompress_evec_from_directory_basis<1000>(dir, eval,uncompressed_evec, ngroups, _grid, neig, crc32, b, nn, nkeep_single, blocks, _FP16_COEF_EXP_SHARE_FLOATS);
  } else {
    if (_grid->IsBoss())
      std::cout << GridLogMessage << "Basis size " << nkeep << " needs to manually be added to the source code" << std::endl;
  }

  return false;
}

void zyd_read_compressed_evecs(const std::string &path, GridRedBlackCartesian *FrbGrid_f, std::vector<LatticeFermionF> &evec, std::vector<double> &eval, int ngroups = 1) {
  // ngroups: It divides the read tasks into ngroups. At a time, only the processes in a group work; each process read the entire file (e.g. lanczos.output/02/0000000002.compressed, ~16GB) and redistribute it to the correct nodes.
  // to save memory, ngroups should be the equal to the number of ranks per node; then, only one process in a node would be reading the file (e.g. lanczos.output/02/0000000002.compressed, ~16GB) at a time.

  evec.resize(1,FrbGrid_f);  // uncompress_evec_from_directory uses the grid of the evec[0].Grid().
  eval.clear();

  uncompress_evec_from_directory(path.c_str(), eval, evec, ngroups);

  // for (int i=0;i<eval.size();i++) 
  for (int i=0; i<eval.size(); i+=200) {
    std::cout << GridLogMessage << "Test: norm2(evec[" << i << "])=" << norm2(evec[i]) << std::endl;
  }
}






}   // end of namespace Grid


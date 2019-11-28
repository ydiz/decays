#pragma once

#include <complex>
#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>

// #include <qlat/grid.h>
#include <headers/pGG.h>

namespace qlat{

// this function is deleted in new version of QLattice
inline void gt_inverse(GaugeTransform& gt, const GaugeTransform& gt0)
{
  TIMER("gt_inverse");
  gt.init(geo_resize(gt0.geo));
  const Geometry& geo = gt.geo;
  qassert(is_matching_geo_mult(gt.geo, gt0.geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    gt.get_elem(xl) = matrix_adjoint(gt0.get_elem(xl));
  }
}

}


std::vector<int> read_mpi_coor(const std::string &prefix) {

	std::ifstream f(prefix + "/geo-info.txt");
	std::string s;
	std::vector<int> mpi_coor;
	while(getline(f, s)) {
		if(s.size() >= 18 && s.substr(0, 18) == "geo.geon.size_node") {
			// cout << s<< endl;
			int i = std::stoi(s.substr(s.find("=")+2));
			mpi_coor.push_back(i);
		}
	}
	assert(mpi_coor.size()==4);
	f.close();
	return mpi_coor;
}


template <class T>
struct TypeMap{
	typedef int type;
};

template<>
struct TypeMap<Grid::ComplexF> {
	using type = Grid::vComplexF;
};

template<>
struct TypeMap<Grid::ComplexD> {
	using type = Grid::vComplexD;
};

namespace Grid {
namespace QCD {
using LatticeLoop = Lattice<iScalar<iScalar<iVector<vComplex, 4> > >>;

}}

void grid_convert(Grid::QCD::LatticeLoop& grid_loop, const qlat::FieldM<qlat::Complex, 4>& qlat_loop)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_loop.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>

    qlat::Vector<qlat::Complex> qlat_loop_site = qlat_loop.get_elems_const(xl); // qlat_loop_site is a vector of Complex; vector size is 16
		assert(qlat_loop_site.size()==4);
		
    typename LatticeLoop::vector_object::scalar_object grid_loop_site;
		// assert(sizeof(qlat_loop_site) == sizeof(grid_loop_site));

		// Complex *p_qlat = (Complex *)&qlat_loop_site; // T is either ComplexF or ComplexD
		// std::copy(p_qlat, p_qlat + 16, (Complex *)&grid_loop_site);	
		std::copy(qlat_loop_site.data(), qlat_loop_site.data() + 4, (Complex *)&grid_loop_site);	

    Coordinate c(coor);
    pokeLocalSite(grid_loop_site, grid_loop, c);
  }
}


// only for double precision PionGGElemField
void grid_convert(Grid::QCD::LatticePGG& grid_pgg, const qlat::PionGGElemField& qlat_pgg)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_pgg.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>
    auto ms = qlat_pgg.get_elems_const(xl); // ms is a vector of WilsonMatrix; vector size is 1
		// qlat::WilsonMatrix qlat_prop_site = ms[0];
		auto qlat_pgg_site = ms[0];
		assert(ms.size()==1);
		
		PGGElem grid_pgg_site;
		assert(sizeof(qlat_pgg_site) == sizeof(grid_pgg_site));

		Complex *p_qlat = (Complex *)&qlat_pgg_site; // T is either ComplexF or ComplexD
	
		std::copy(p_qlat, p_qlat + 16, (Complex *)&grid_pgg_site);	

    // pokeLocalSite(grid_pgg_site, grid_pgg, coor);
    Coordinate c(coor);
    pokeLocalSite(grid_pgg_site, grid_pgg, c);
  }
}

// For cheng's three point function 
void grid_convert(Grid::QCD::LatticePGG& grid_pgg, const qlat::FieldM<qlat::Complex, 16>& qlat_pgg)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_pgg.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>

    qlat::Vector<qlat::Complex> qlat_pgg_site = qlat_pgg.get_elems_const(xl); // qlat_pgg_site is a vector of Complex; vector size is 16
		assert(qlat_pgg_site.size()==16);
		
		PGGElem grid_pgg_site;
		// assert(sizeof(qlat_pgg_site) == sizeof(grid_pgg_site));

		// Complex *p_qlat = (Complex *)&qlat_pgg_site; // T is either ComplexF or ComplexD
		// std::copy(p_qlat, p_qlat + 16, (Complex *)&grid_pgg_site);	
		std::copy(qlat_pgg_site.data(), qlat_pgg_site.data() + 16, (Complex *)&grid_pgg_site);	

    // pokeLocalSite(grid_pgg_site, grid_pgg, coor);
    Coordinate c(coor);
    pokeLocalSite(grid_pgg_site, grid_pgg, c);
  }
}



// For luchang's three point function // file "decay" and "fission"
// Here, input is a 8x8 matrix. mu=0-7, nu=0-7; we take only the top left submatrix
// see https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/luchang/0001%20pion-gg%20status/
void grid_convert(Grid::QCD::LatticePGG& grid_pgg, const qlat::FieldM<qlat::Complex, 64>& qlat_pgg)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_pgg.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>

    qlat::Vector<qlat::Complex> qlat_pgg_site = qlat_pgg.get_elems_const(xl); // qlat_pgg_site is a vector of Complex; vector size is 16
		assert(qlat_pgg_site.size()==64);
		
		PGGElem grid_pgg_site;

    for(int mu=0; mu<4; ++mu) {
      std::copy(qlat_pgg_site.data() + 8*mu, qlat_pgg_site.data() + 8*mu + 4, ((Complex *)&grid_pgg_site)+4*mu);	
    }

    // pokeLocalSite(grid_pgg_site, grid_pgg, coor);
    Coordinate c(coor);
    pokeLocalSite(grid_pgg_site, grid_pgg, c);
  }
}





void grid_convert(Grid::LatticeColourMatrix& grid_gt, const qlat::GaugeTransform& qlat_gt)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_gt.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>
    auto ms = qlat_gt.get_elems_const(xl); // ms is a vector of WilsonMatrix; vector size is 1
		// qlat::WilsonMatrix qlat_prop_site = ms[0];
		auto qlat_gt_site = ms[0];
		assert(ms.size()==1);
		
		typename LatticeColourMatrix::vector_object::scalar_object grid_gt_site;
		assert(sizeof(qlat_gt_site) == sizeof(grid_gt_site));

		Complex *p_qlat = (Complex *)&qlat_gt_site; // T is either ComplexF or ComplexD
	
		std::copy(p_qlat, p_qlat + 9, (Complex *)&grid_gt_site);	

    // pokeLocalSite(grid_gt_site, grid_gt, coor);
    Coordinate c(coor);
    pokeLocalSite(grid_gt_site, grid_gt, c);
  }
}

template<class T> // T can be ComplexF or ComplexD
typename std::enable_if<std::is_same<T, Grid::ComplexF>::value || std::is_same<T, Grid::ComplexD>::value, void>::type 
grid_convert(Grid::Lattice<Grid::iSpinColourMatrix<typename TypeMap<T>::type >>& grid_prop, const qlat::Propagator4dT<T>& qlat_prop)
{
  using namespace Grid;
  using namespace Grid::QCD;
  const qlat::Geometry& geo = qlat_prop.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const qlat::Coordinate xl = geo.coordinate_from_index(index); // get  local coordinate
    // std::vector<int> coor = grid_convert(xl); // just copy the four components of xl to a vector<int>
    std::vector<int> coor(xl.begin(), xl.end()); // just copy the four components of xl to a vector<int>
    auto ms = qlat_prop.get_elems_const(xl); // ms is a vector of WilsonMatrix; vector size is 1
	// qlat::WilsonMatrix qlat_prop_site = ms[0];
		auto qlat_prop_site = ms[0];
		assert(ms.size()==1);
		
		Grid::iSpinColourMatrix< T > grid_prop_site;
		assert(sizeof(qlat_prop_site) == sizeof(grid_prop_site));

		T *p_qlat = (T *)&qlat_prop_site; // T is either ComplexF or ComplexD
		
		for(int row=0; row<12; ++row)
			for(int column=0; column<12; ++column) {
				int grid_spin_row = row/3;
				int grid_color_row = row%3;
				int grid_spin_column = column/3;
				int grid_color_column = column%3;
				grid_prop_site()(grid_spin_row, grid_spin_column)(grid_color_row, grid_color_column) = *(p_qlat + 12*row + column);
			}

		// pokeLocalSite(grid_prop_site, grid_prop, coor);
    Coordinate c(coor);
		pokeLocalSite(grid_prop_site, grid_prop, c);
  }
}

namespace qlat {

void print_qlat_field_site(const PionGGElemField &field, const std::vector<int> coor) {
	std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
	qlat::Coordinate site_coor(coor[0], coor[1], coor[2], coor[3]);
	auto x = qlat::field_get_elems(field, site_coor);
	for(size_t mu=0; mu<field.geo.multiplicity; ++mu)  // for gauge field, multiplicity is 4; for propagator, it is 1
	{
		std::cout << "mu = " << mu << std::endl;
		// std::cout << x[mu].em() << std::endl;
		std::cout << show_pgge(x[mu]) << std::endl;
	}
}

template<class T>
void print_qlat_field_site(const T &field, const std::vector<int> coor) {
	std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
	qlat::Coordinate site_coor(coor[0], coor[1], coor[2], coor[3]);
	auto x = qlat::field_get_elems(field, site_coor);
	for(size_t mu=0; mu<field.geo.multiplicity; ++mu)  // for gauge field, multiplicity is 4; for propagator, it is 1
	{
		std::cout << "mu = " << mu << std::endl;
		std::cout << x[mu].em() << std::endl;
	}
}

}


// only for gauge_field
template<class T>
void print_field(const T &field) {
  std::vector<int> g_size(4);
  for(size_t i=0; i<4; ++i) g_size[i] = field.geo.node_site[i] * field.geo.geon.size_node[i];
  for(size_t x0=0; x0<g_size[0]; ++x0)
	for(size_t x1=0; x1<g_size[1]; ++x1)
		for(size_t x2=0; x2<g_size[2]; ++x2)
			for(size_t x3=0; x3<g_size[3]; ++x3)
			{
				std::cout << "[ " << x3 << " " << x2 << " " << x1 << " " << x0 << " ]" << std::endl;
				qlat::Coordinate coor(x3, x2, x1, x0);
				auto x = qlat::field_get_elems(field, coor);
				for(size_t mu=0; mu<field.geo.multiplicity; ++mu)  // for gauge field, multiplicity is 4; for propagator, it is 1
				{
					std::cout << "mu = " << mu << std::endl;
					std::cout << x[mu].em() << std::endl;
				}
			}
}

namespace Grid {
namespace QCD {
// FIXME: grid_convert has not been defined
void read_loop(LatticeLoop &lat, const std::string &path) {
  qlat::FieldM<Complex, 4> qlat_loop;
  // qlat::dist_read_field_double(qlat_pgg, path);
	dist_read_field_double_from_float(qlat_loop, path);
  grid_convert(lat, qlat_loop);
}
// for both wall and point propagators
void read_qlat_propagator(LatticePropagator &lat, const std::string &path) {
	qlat::Propagator4d qlat_prop;
    // std::cout << GridLogMessage << "before reading" << std::endl;
	dist_read_field_double_from_float(qlat_prop, path);
    // std::cout << GridLogMessage << "after reading" << std::endl;
	grid_convert(lat, qlat_prop);
    // std::cout << GridLogMessage << "after grid_convert" << std::endl;
}

void read_qlat_propagator_no_dist(LatticePropagator &lat, const std::string &path) {
	qlat::Propagator4d qlat_prop;
    std::cout << GridLogMessage << "no dist before reading" << std::endl;
	read_field_double_from_float(qlat_prop, path);
    std::cout << GridLogMessage << "no dist after reading" << std::endl;
	grid_convert(lat, qlat_prop);
    // std::cout << GridLogMessage << "no dist after grid_convert" << std::endl;
}


void read_cheng_PGG(LatticePGG &lat, const std::string &path) { // read file xxxxx/decay_cheng
  qlat::FieldM<qlat::Complex, 16> qlat_pgg;
    std::cout << GridLogMessage << "before reading" << std::endl;
  qlat::read_field(qlat_pgg, path);
    std::cout << GridLogMessage << "after reading" << std::endl;
  qlat::to_from_big_endian_64(qlat::get_data(qlat_pgg));

  grid_convert(lat, qlat_pgg);
}


void read_luchang_PGG(LatticePGG &lat, const std::string &path) { // read file xxxxx/decay, and xxxx/fission
  qlat::FieldM<qlat::Complex, 64> qlat_pgg;
    std::cout << GridLogMessage << "before reading" << std::endl;
  qlat::read_field(qlat_pgg, path);
    std::cout << GridLogMessage << "after reading" << std::endl;
  qlat::to_from_big_endian_64(qlat::get_data(qlat_pgg));
  grid_convert(lat, qlat_pgg);
}


// void read_cheng_PGG(LatticePGG &lat, const std::string &path) {
//   qlat::PionGGElemField qlat_pgg;
//   dist_read_field(qlat_pgg, path);
//   grid_convert(lat, qlat_pgg);
// }
//
// void read_luchang_PGG(LatticePGG &lat, const std::string &path) {
//   // qlat::PionGGElemField qlat_pgg;
//   qlat::FieldM<Complex, 16> qlat_pgg;
//   qlat::dist_read_field_double(qlat_pgg, path);
//   grid_convert(lat, qlat_pgg);
// }

}}


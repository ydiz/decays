#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdlib>

template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
	using namespace Grid;
	std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
	typename T::vector_object::scalar_object site;
	peekSite(site, field, coor);
	std::cout << site << std::endl;
}

bool dirExists(const std::string &path){
	struct stat info;
	if( stat( path.c_str(), &info ) == 0 ) return true; // dir does exist
	else return false;
}


int get_t(const std::string &path) {
	return std::stoi(path.substr(path.find("=") + 1));
}

void get_ts(const std::string &path, std::vector<int> &ts, std::map<int, std::string> &subdirs) {

	ts.clear();
	subdirs.clear();
	DIR *dir;
	dir = opendir(path.c_str());
	assert(dir!=NULL); // make sure directory exists
	struct dirent *entry;

	std::string subdir_name;
	while ((entry = readdir (dir)) != NULL) {
		// printf ("%s\n", entry->d_name);
		subdir_name = std::string(entry->d_name);
		if(subdir_name.substr(0, 2) == "t=") {
			int t = get_t(subdir_name); 
			ts.push_back(t);
			subdirs.insert(std::pair<int, std::string>(t, path + "/" + subdir_name));
		}
	}
	closedir (dir);

	std::sort(ts.begin(), ts.end());
}


std::vector<int> get_xg(const std::string &path) {

	std::stringstream ss;
	ss.str(path.substr(path.find("(") + 1));

	std::vector<int> ret(4);
	for(int &x: ret) { 
		ss >> x; 
		ss.ignore(); // extract comma and ignore it
	}

	return ret;
}

void get_xgs(const std::string &path, std::vector<std::vector<int>> &xgs, std::map<std::vector<int>, std::string> &subdirs, char quark='l') {
  xgs.clear();
  subdirs.clear();

  std::string type;
  if(quark=='l') type = "0";
  else if(quark=='s') type = "1";
  else assert(0);

	DIR *dir;
	dir = opendir(path.c_str());
	assert(dir!=NULL); // make sure directory exists
	struct dirent *entry;

	std::string subdir_name;
	while ((entry = readdir (dir)) != NULL) {
		// printf ("%s\n", entry->d_name);
		subdir_name = std::string(entry->d_name);
		// if(subdir_name.substr(0, 3) == "xg=" && subdir_name.substr(subdir_name.find("type"), 6) == ("type="+type) && subdir_name.substr(subdir_name.find("accuracy"), 10) == "accuracy=0") {
		if(subdir_name.substr(0, 3) == "xg=" && subdir_name.substr(subdir_name.find("type"), 6) == ("type="+type) && subdir_name.substr(subdir_name.find("accuracy")) == "accuracy=0") {
			std::vector<int> xg = get_xg(subdir_name); 
			xgs.push_back(xg);
			subdirs.insert(std::pair<std::vector<int>, std::string>(xg, path + "/" + subdir_name));
		}
	}
	closedir (dir);
}




namespace Grid {
namespace QCD {

template<class T>
void writeScidac(T& field, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
  if(field.Grid()->IsBoss()) {
    std::string base_dir = filename.substr(0, filename.rfind('/'));
    std::cout << "base_dir: " << base_dir << std::endl;
    system(("mkdir -p " + base_dir).c_str());
    // mkdir(base_dir.c_str(), 0777);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  emptyUserRecord record;
  ScidacWriter WR(field.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
	WR.writeScidacFieldRecord(field, record);
  WR.close();
};

template<class T>
void readScidac(T& field, const std::string &filename){
  emptyUserRecord record;
  ScidacReader RD;
  RD.open(filename);
	RD.readScidacFieldRecord(field, record);
  RD.close();
};

// template<class T>
void writeScidac_prop_d2f(LatticePropagatorD& prop, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
  static GridCartesian * grid_f = SpaceTimeGrid::makeFourDimGrid(prop.Grid()->FullDimensions(), GridDefaultSimd(Nd,vComplexF::Nsimd()), prop.Grid()->_processors); 
  LatticePropagatorF prop_f(grid_f);
  precisionChange(prop_f, prop);

  if(prop_f.Grid()->IsBoss()) {
    std::string base_dir = filename.substr(0, filename.rfind('/'));
    std::cout << "base_dir: " << base_dir << std::endl;
    system(("mkdir -p " + base_dir).c_str());
    // mkdir(base_dir.c_str(), 0777);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  emptyUserRecord record;
  ScidacWriter WR(prop_f.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
  WR.writeScidacFieldRecord(prop_f, record);
  WR.close();
};

void readScidac_prop_f2d(LatticePropagatorD& prop, const std::string &filename){
  static GridCartesian * grid_f = SpaceTimeGrid::makeFourDimGrid(prop.Grid()->FullDimensions(), GridDefaultSimd(Nd,vComplexF::Nsimd()), prop.Grid()->_processors); 
  LatticePropagatorF prop_f(grid_f);

  assert(dirExists(filename));

  emptyUserRecord record;
  ScidacReader RD;
  RD.open(filename);
  RD.readScidacFieldRecord(prop_f, record);
  RD.close();

  precisionChange(prop, prop_f);
};






}}

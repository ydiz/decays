#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include <sys/stat.h>

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
  emptyUserRecord record;
  ScidacWriter WR(field._grid->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
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








}}

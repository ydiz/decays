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

std::vector<std::vector<int>> coor_from_file(const std::string &fname) {

  using namespace std;
  ifstream f(fname);

  vector<vector<int>> rst;

  std::string tmp;
  while(getline(f, tmp)) {
    vector<int> coor;
    stringstream ss(tmp);
    int i;
    while(ss >> i) coor.push_back(i);
    rst.push_back(coor);
  }

  f.close();
  cout << "There are " << rst.size() << " point sources" << endl;
  return rst;
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





/////////////////////////////////////
// Taken from  Hadrons/Hadrons/A2AVectors.hpp, Hadrons/Hadrons/Global.cpp 
// usage: A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
/////////////////////////////////////

#define MAX_PATH_LENGTH 512u

// recursive mkdir /////////////////////////////////////////////////////////////
int mkdir(const std::string dirName)
{
  if (!dirName.empty() and access(dirName.c_str(), R_OK|W_OK|X_OK))
  {
    mode_t mode755;
    char   tmp[MAX_PATH_LENGTH];
    char   *p = NULL;
    size_t len;

    mode755 = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;

    snprintf(tmp, sizeof(tmp), "%s", dirName.c_str());
    len = strlen(tmp);
    if(tmp[len - 1] == '/')
    {
      tmp[len - 1] = 0;
    }
    for(p = tmp + 1; *p; p++)
    {
      if(*p == '/')
      {
        *p = 0;
        ::mkdir(tmp, mode755);
        *p = '/';
      }
    }

    return ::mkdir(tmp, mode755);
  }
  else
  {
    return 0;
  }
}

std::string dirname(const std::string &s)
{
  constexpr char sep = '/';
  size_t         i   = s.rfind(sep, s.length());

  if (i != std::string::npos)
  {
    return s.substr(0, i);
  }
  else
  {
    return "";
  }
}

void makeFileDir(const std::string filename, GridBase *g)
{
  bool doIt = true;

  if (g)
  {
    doIt = g->IsBoss();
  }
  if (doIt)
  {
    std::string dir    = dirname(filename);
    int         status = mkdir(dir);

    if (status)
    {
      std::cout <<  "cannot create directory '" + dir << std::endl;
    }
  }
}






class A2AVectorsIo
{
  public:
    struct Record: Serializable
  {
    GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
        unsigned int, index);
    Record(void): index(0) {}
  };
  public:
    template <typename Field>
      static void write(const std::string fileStem, std::vector<Field> &vec, 
          const bool multiFile, const int trajectory = -1);
    template <typename Field>
      static void read(std::vector<Field> &vec, const std::string fileStem,
          const bool multiFile, const int trajectory = -1);
  private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
        const bool multiFile)
    {
      std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

      if (multiFile)
      {
        return stem + t;
      }
      else
      {
        return stem + t + ".bin";
      }
    }
};


  template <typename Field>
void A2AVectorsIo::write(const std::string fileStem, std::vector<Field> &vec, 
    const bool multiFile, const int trajectory)
{
  Record       record;
  GridBase     *grid = vec[0].Grid();
  ScidacWriter binWriter(grid->IsBoss());
  std::string  filename = vecFilename(fileStem, trajectory, multiFile);

  if (multiFile)
  {
    std::string fullFilename;

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

      std::cout << GridLogMessage << "Writing vector " << i << std::endl;
      makeFileDir(fullFilename, grid);
      binWriter.open(fullFilename);
      record.index = i;
      binWriter.writeScidacFieldRecord(vec[i], record);
      binWriter.close();
    }
  }
  else
  {
    makeFileDir(filename, grid);
    binWriter.open(filename);
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      std::cout << GridLogMessage << "Writing vector " << i << std::endl;
      record.index = i;
      binWriter.writeScidacFieldRecord(vec[i], record);
    }
    binWriter.close();
  }
}

  template <typename Field>
void A2AVectorsIo::read(std::vector<Field> &vec, const std::string fileStem, 
    const bool multiFile, const int trajectory)
{
  Record       record;
  ScidacReader binReader;
  std::string  filename = vecFilename(fileStem, trajectory, multiFile);

  if (multiFile)
  {
    std::string fullFilename;

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

      // LOG(Message) << "Reading vector " << i << std::endl;
      std::cout << "Reading vector " << i << std::endl;
      binReader.open(fullFilename);
      binReader.readScidacFieldRecord(vec[i], record);
      binReader.close();
      if (record.index != i)
      {
        std::cout << "vector index mismatch" << std::endl;
        // HADRONS_ERROR(Io, );
      }
    }
  }
  else
  {
    binReader.open(filename);
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      // LOG(Message) << "Reading vector " << i << std::endl;
      std::cout << "Reading vector " << i << std::endl;
      binReader.readScidacFieldRecord(vec[i], record);
      if (record.index != i)
      {
        std::cout << "vector index mismatch" << std::endl;
        // HADRONS_ERROR(Io, "vector index mismatch");
      }
    }
    binReader.close();
  }
}





















template<class T>
void writeScidac(T& field, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
  std::cout << "writing to" << filename << std::endl;
  makeFileDir(filename, field.Grid());
  // if(field.Grid()->IsBoss()) {
  //   std::string base_dir = filename.substr(0, filename.rfind('/'));
  //   std::cout << "base_dir: " << base_dir << std::endl;
  //   assert(dirExists(base_dir));
  //   // system(("mkdir -p " + base_dir).c_str());      //  sometimes has error
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  emptyUserRecord record;
  ScidacWriter WR(field.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
	WR.writeScidacFieldRecord(field, record);
  WR.close();
};

template<class T>
void readScidac(T& field, const std::string &filename){
  std::cout << "reading from" << filename << std::endl;
  assert(dirExists(filename));
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

  makeFileDir(filename, prop.Grid());
  // if(prop_f.Grid()->IsBoss()) {
  //   std::string base_dir = filename.substr(0, filename.rfind('/'));
  //   std::cout << "base_dir: " << base_dir << std::endl;
  //   system(("mkdir -p " + base_dir).c_str());
  //   // mkdir(base_dir.c_str(), 0777);
  // }
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
















}

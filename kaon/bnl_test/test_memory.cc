#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

void print_memory() {
  using namespace std;

  map<string, double> memInfo; // item -> memory in GB
  ifstream f("/proc/meminfo");
  while(f) {
    string item, tmp;
    long long int size;

    f >> item >> size;
    getline(f, tmp); // move to the next line

    if(!item.empty()) {
      item.pop_back(); // remove last character
      memInfo[item] = double(size) / 1024 / 1024; // kB -> GB
    }
  }

  std::cout << "MemTotal: " << memInfo["MemTotal"] << " GB " << std::endl;
  std::cout << "Non cache/buffer MemUsed: " << memInfo["MemTotal"] - memInfo["MemFree"] - memInfo["Buffers"] - memInfo["Cached"] << " GB " << std::endl;
  std::cout << "Buffers/Cached: " << memInfo["Buffers"] + memInfo["Cached"] << " GB " << std::endl;
  std::cout << "MemAvailable: " << memInfo["MemAvailable"] << " GB " << std::endl; // Memavailable is roughly memfree + buffers + cached,

}

int main() {
  MPI_Init(NULL, NULL);
  int Nnodes;
  MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0) cout << "Number of nodes: " << Nnodes << endl;
  if(rank==0) print_memory();
  if(rank==0) cout << string(40, '=') << endl;
  vector<double> vec(1024 * 1024 * 1024); // should allocate 8GB memory
  if(rank==0) print_memory();
  std::cout << "Finised" << std::endl;

  // std::cout << "=============" << std::endl;
  //
  // print_memory();
  // vector<double> vec2(1024 * 1024 * 1024);
  // print_memory();
  // std::cout << "malloc " << std::endl;
  // void *buffer = malloc(1024 * 1024 * 1024);
  // print_memory();

}


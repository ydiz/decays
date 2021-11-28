#include "kaon/kaon.h"

using namespace std;
using namespace Grid;

vector<vector<int>> select_half_point_sources(vector<vector<int>> &xgs) {

  sort(xgs.begin(), xgs.end(), [](const vector<int> &v1, const vector<int> &v2) {return v1[3]<v2[3];});

  vector<vector<int>> xgs_new;

  int i=0, j=0;
  while(i<xgs.size()) {
    while(j<xgs.size()-1 && xgs[j][3]==xgs[j+1][3]) ++j; // j points to the last element with the same time coordinate

    int n = j - i + 1;
    for(int k=i; k<=i+(n-1)/2; ++k) xgs_new.push_back(xgs[k]);

    j += 1;
    i = j;
  }
  return xgs_new;
}

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

#ifdef CUTH_FREE_FIELD
  int max_uv_sep = 3;
  int min_tsep = 2;
  int max_tsep = 4;
  Env env("FreeField_8nt8");
#else
  int max_uv_sep = 16;
  int min_tsep = 6;
  int max_tsep = 14;
  Env env("24ID");
#endif

  env.setup_traj(1234);

  // vector<vector<int>> xgs = env.xgs_s;


  vector<vector<int>> xgs_new = select_half_point_sources(env.xgs_s);

  for(auto &x: xgs_new) {
    std::cout << x << std::endl;
  }

  std::cout << GridLogMessage << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}

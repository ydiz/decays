#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int distance(int t1, int t2, int T) {
  int tmp = std::abs(t1 - t2);
  if(tmp <= T/2) return tmp;
  else return T-tmp;
}

int left_distance(int t1, int t2, int T) { // the distance from t1 to t2 if you can only move to the left 
  int tmp = t1 - t2;
  if(tmp >= 0) return tmp;
  else return T+tmp;
}

int right_distance(int t1, int t2, int T) { // the distance from t1 to t2 if you can only move to the right
  int rst = T - left_distance(t1, t2, T);
  if(rst==T) return 0;
  else return rst;
}



int main() {

	// // int xt = 6;
  // std::vector<int> x {1,1,1,6};
  // std::vector<int> xp {1,1,1,0};
  // int t_min = 10;
  //
  //
	// for(xp[3]=0; xp[3]<64; ++xp[3]){
  //   int t_wall = (rightPoint(x[3], xp[3], 64) + t_min) % 64;
  //   int t_sep = distance(x[3], t_wall, 64);
  //
	// 	cout << "xpt: " << xp[3] << " my result: " << t_wall << " " << t_sep << "\t";
  // }

  int T = 64;
  int t_wall = 10;

	for(int t=0; t<64; ++t){
		cout << "t: " << t << " left_distance: " << left_distance(t_wall, t, T)  << " right_distance: " << right_distance(t_wall, t, T) << endl;
  }

	return 0;
}

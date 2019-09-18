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

int rightPoint(int t_base, int t, int T) { // after shift t_base to 0, determine wheter t is on the right or t_base is on the right
  int tmp = t - t_base; // shift t_base to 0
  if(tmp > 0 && tmp <=T/2) return t;
  else return t_base; // tmp < 0 || tmp > T/2
}

inline int mod(const int x, const int len)
{
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline int smod(const int x, const int len)
{
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}


int main() {

	// int xt = 6;
  std::vector<int> x {1,1,1,6};
  std::vector<int> xp {1,1,1,0};
  int t_min = 10;


	for(xp[3]=0; xp[3]<64; ++xp[3]){
    int t_wall = (rightPoint(x[3], xp[3], 64) + t_min) % 64;
    int t_sep = distance(x[3], t_wall, 64);

		cout << "xpt: " << xp[3] << " my result: " << t_wall << " " << t_sep << "\t";

    // cheng's way
    // int t_wall;
    // int t_sep;
    //
    int diff = smod(x[3] - xp[3], 64);
    if (diff >= 0)
    {
      t_wall = mod(x[3] + t_min, 64);
    } else {
      t_wall = mod(xp[3] + t_min, 64);
    }
    // t_sep = mod(t_wall - x[3], 64);
    t_sep = std::abs(smod(t_wall - x[3], 64));

		cout << "cheng's "<< t_wall << " " << t_sep  << endl;
  }

	return 0;
}

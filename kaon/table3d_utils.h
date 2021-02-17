using namespace std;
using namespace Grid;

vector<vector<vector<Complex>>> initialize_table3d(int T) { // generate table3d[tK][xt][vt]; initialize it to 0
    vector<vector<vector<Complex>>> table3d(T);  
    for(int tK=0; tK<T; ++tK) {
      table3d[tK].resize(T);
      for(int xt=0; xt<T; ++xt) {
        table3d[tK][xt].resize(T, 0.);
      }
    }
    return table3d;
}

vector<vector<Complex>> table3d_to_table2d(const vector<vector<vector<Complex>>> &table3d) { // table3d[tK][xt][vt] -> table2d[tK][vt] with xt = 0
  int T = table3d.size();

  vector<vector<Complex>> table2d(T);
  for(auto &x: table2d) x.resize(T, 0.);


  for(int tK=0; tK<T; ++tK) { 
    for(int vt=0; vt<T; ++vt) {
      for(int xt=0; xt<T; ++xt) {
        table2d[tK][vt] += table3d[(tK+xt)%T][xt][(vt + xt)%T];
      }
    }
  }

  for(int tK=0; tK<T; ++tK)   
    for(int vt=0; vt<T; ++vt) 
      table2d[tK][vt] /= double(T);

  return table2d;
}





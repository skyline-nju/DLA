#include <iostream>
#include <cstdio>
#include <cmath>
#include "run.h"

using namespace std;

int main(int argc, char *argv[]) {
  int N = 1000;
  vector<Disk> cluster;
  cluster.reserve(N);
  cluster.push_back(Disk(0, 0));
  Grid grid(15000, 20);
  Ran myran(2);
  run(cluster, 100000, &myran, grid);
  cal_fractal_dim(cluster);
  output_xyz(cluster);
}
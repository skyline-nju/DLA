#include <iostream>
#include <cstdio>
#include <cmath>
#include "run.h"
#include "vec.h"

using namespace std;

int main(int argc, char *argv[]) {
  //int N = 500001;
  //vector<Disk> cluster;
  //cluster.reserve(N);
  //cluster.push_back(Disk(0, 0));
  //Grid grid(12000, 20);
  //Ran myran(14);
  //run(cluster, N, &myran, grid);
  //cal_fractal_dim(cluster);
  //Disk::output_xyz(cluster);

  int N = atoi(argv[1]);
  vector<Rect> rect;
  rect.reserve(N);
  Cell<Rect> cell(atoi(argv[2]), ceil(2 * Rect::Rab) + 1);
  int seed = atoi(argv[3]);
  run(rect, N, cell, seed, atof(argv[4]));
}                        
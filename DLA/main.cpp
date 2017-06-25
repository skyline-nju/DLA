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

  int N = 100;
  vector<Rect> rect;
  rect.reserve(N);
  Ran myran(9);
  Cell cell(500, ceil(2 * Rect::Rab) + 1);
  run(rect, N, cell, &myran);
  //run(rect, N, &myran);
  Rect::output(rect);
}                        
#include <iostream>
#include <cstdio>
#include <cmath>
#include "run.h"
#include "vec.h"

using namespace std;

int main(int argc, char *argv[]) {
  //int N = 1000;
  //vector<Disk> cluster;
  //cluster.reserve(N);
  //cluster.push_back(Disk(0, 0));
  //Grid grid(15000, 20);
  //Ran myran(2);
  //run(cluster, N, &myran, grid);
  //cal_fractal_dim(cluster);
  //output_xyz(cluster);

  clock_t beg = clock();
  int N = 60;
  vector<Rect> rect;
  rect.reserve(N);
  Ran myran(99);
  //Cell cell(500, ceil(2 * Rect::Rab) + 1);
  //run(rect, N, cell, &myran);
  run(rect, N, &myran);
  Rect::output(rect);
  //cout << "Elapsed time = " << (clock() - beg) / CLOCKS_PER_SEC << endl;
}                        
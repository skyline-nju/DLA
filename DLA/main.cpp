#include <iostream>
#include <cstdio>
#include <cmath>
#include "run.h"
#include "vec.h"
#include "rotate.h"

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

  //int N = 100;
  //vector<Rect> rect;
  //rect.reserve(N);
  //Ran myran(5);
  //run(rect, N, &myran);
  //Rect::output(rect);

  vector<Point> point_set;
  vector<Segment> line_set;

  Vec2<double> O(0, 0);
  Vec2<double> M(1, 1);
  Vec2<double> N(1, -1);
  get_increase_segment_set(O, M, N, point_set, line_set);
  for (auto i : line_set) {
    cout << i.p1.d << endl;
    cout << i.p2.d << endl;
  }
}                        
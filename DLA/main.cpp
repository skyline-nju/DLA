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

  Rect a(0, 0, 0);
  vector<Rect> b;
  b.push_back(a);
  Vec2<double> Delta(4, 5);
  a.shift(Delta);
  a.rotate(PI / 10);
  b.push_back(a);
  Rect::output(b);

  double angle;
  bool flag = false;
  b[0].collideR(b[1], true, angle, flag);
  if (flag) {
    Rect c = b[0];
    c.rotate(angle);
    //c.rotate(-acos(0.943373));
    b.push_back(c);
    Rect::output(b);
    cout << angle / PI * 180 << endl;
  } else {
    cout << "no collision" << endl;
  }


}                        
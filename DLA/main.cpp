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

  int N = 10;
  vector<Rect> rect;
  rect.reserve(N);
  Ran myran(5);
  run(rect, N, &myran);
  Rect::output(rect);

  //Rect a(0, 0, 60);
  //vector<Rect> b;
  //b.push_back(a);
  //Vec2<double> Delta(18, 0);
  //a.shift(Delta);
  //a.rotate(-PI / 4);
  //b.push_back(a);
  //Rect::output(b);

  //double angle;
  //bool CW = false;
  //RotStatus status(PI / 36);
  //vector<Vector2D> point_set;
  //vector<Segment> segment_set;
  //b[0].get_segment_set(CW, point_set, segment_set);
  //b[0].collideR(point_set, segment_set, b[1], CW, status);
  //if (status.flag) {
  //  Rect c = b[0];
  //  angle = acos(status.cos_angle);
  //  if (CW)
  //    angle = -angle;
  //    
  //  c.rotate(angle);
  //  b.push_back(c);
  //  Rect::output(b);
  //  cout << angle / PI * 180 << endl;
  //  cout << fixed << setprecision(14) << status.contact_point << endl;
  //  cout << "angle = " << angle << endl;
  //  for (int i = 0; i < 4; i++) {
  //    cout << c.vertex[i] << endl;
  //  }
  //} else {
  //  cout << status.cos_angle << endl;
  //  cout << "no collision" << endl;
  //}


  //double angle;
  //bool flag = false;
  //b[0].collideR(b[1], false, angle, flag);
  //if (flag) {
  //  Rect c = b[0];
  //  c.rotate(angle);
  //  b.push_back(c);
  //  Rect::output(b);
  //  cout << angle / PI * 180 << endl;
  //} else {
  //  cout << "no collision" << endl;
  //}

}                        
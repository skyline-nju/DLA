#include <iostream>
#include <cstdio>
#include <cmath>
#include "run.h"
#include "vec.h"
#include "rotate.h"

using namespace std;

struct Foo
{
  Foo(int x):b(x) {}
  void add(int a) const {
    cout << "a + b = " << a + b << endl;
  }
  void show(); 
  //void (*f)(int);
  int b;
  void(Foo::*f2)(int a) const;
};

void (Foo::*f)(int a) const ;
void Foo::show() {
  f2 = &Foo::add;
  (this->*f2)(2);
}
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

  int N = 50;
  vector<Rect> rect;
  rect.reserve(N);
  Ran myran(18);
  Cell cell(1500, ceil(2 * Rect::Rab) + 1);
  run(rect, N, cell, &myran);
  Rect::output(rect);
  cout << cell.get_row(0) << endl;

  //Foo A(2);
  //A.add(1);
  //f = &Foo::add;
  //(A.*f)(2);
  //A.show();
}                        
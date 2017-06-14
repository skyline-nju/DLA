#ifndef PARTICLE_H
#define PARTICLE_H

struct Disk
{
  Disk() { x = y = 0; }
  Disk(double x0, double y0) :x(x0), y(y0) {}
  double get_rr(const Disk &d2) {
    return (x - d2.x) * (x - d2.x) + (y - d2.y) * (y - d2.y);
  }
  double x;
  double y;
};

struct Rod
{
  Rod() {}

  double xc;
  double yc;
  double ux;
  double uy;
  static double a;
  static double b;
};
#endif

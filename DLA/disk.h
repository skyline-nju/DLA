#ifndef DISK_H
#define DISK_H
#include <vector>
#include <fstream>

struct Disk
{
  Disk() { x = y = 0; }
  Disk(double x0, double y0) :x(x0), y(y0) {}
  double get_rr(const Disk &d2);
  void collide(const Disk &d, double ux, double uy,
    double l, double &l_hit, bool &flag) const;
  static void output(const std::vector<Disk> &cluster);
  static void output_xyz(const std::vector<Disk> &cluster);

  double x;
  double y;

  static double sigma;
  static double rp;
};

inline double Disk::get_rr(const Disk &d2) {
  return (x - d2.x) * (x - d2.x) + (y - d2.y) * (y - d2.y);
}
#endif

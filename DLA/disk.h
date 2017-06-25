#ifndef DISK_H
#define DISK_H
#include <vector>
#include <fstream>
#include "grid.h"
#include "rand.h"
#include "comn.h"
#include "vec.h"

struct Status {
  Status(double lm) : l(lm), flag(false), count(0) {}
  double l;
  int count;
  int neighbor_tag;
  bool flag;
};

struct Disk
{
  Disk() { x = y = 0; }
  Disk(double x0, double y0) :x(x0), y(y0) {}
  double get_rr(const Disk &d2);
  void collide(const Disk &my_disk, double ux, double uy,
               int tag, Status &status) const;
  void move(const std::vector<Disk> &cluster, Ran *myran, const Grid &grid,
            Status &status);
  void slip(const std::vector<Disk> &cluster, Ran *myran, const Grid &grid,
            const Status &status);
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

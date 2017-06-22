#ifndef PARTICLE_H
#define PARTICLE_H
#include "vec.h"
#include "comn.h"
#include "rotate.h"
#include "grid.h"

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

struct Rect
{
  Rect() {}
  Rect(double xc, double yc, double theta);
  Rect(double xc, double yc, double ux, double uy);
  void cal_vertex();
  void get_mov_dir(int idx0, Vec2<double> &u) const;
  void collide_longitudinal(int idx0, const Vec2<double> &u, const Rect &rect,
                            double l, double &l_hit, bool &flag_collide) const;
  void collide_transverse(int idx0, const Vec2<double> &u, const Rect &rect,
                          double l, double &l_hit, bool &flag_collide) const;
  void translate(const std::vector<Rect> &cluster, double lm,
                 int idx0, bool &collided);
  void translate(const std::vector<Rect> &cluster, const Cell &cell,
                 double lm, int idx0, bool &collided);
  void rotate(const std::vector<Rect> &cluster, double theta_m,
              bool clockwise, bool &collided);
  void rotate(const std::vector<Rect> &cluster, const Cell &cell, 
              double theta_m, bool clockwise, bool &collided);
  void get_segment_set(bool CW, std::vector<Vector2D> &my_point_set,
                       std::vector<Segment> &my_segment_set) const;
  void shift(const Vec2<double> &delta);
  static void output(const std::vector<Rect> &rect);
  static void output(const std::vector<Rect> &rect, const char *f);

  Vec2<double> center;
  Vec2<double> orient;
  Vec2<double> vertex[4];

  static double a;
  static double b;
  static double La;
  static double Lb;
  static double Rab;
  static void (Rect::*collide_wrapper)(int, const Vec2<double> &,
                                       const Rect &, double, double &,
                                       bool &) const;
};

inline double Disk::get_rr(const Disk &d2) {
  return (x - d2.x) * (x - d2.x) + (y - d2.y) * (y - d2.y);
}

inline void Rect::cal_vertex() {
  Vec2<double> dR1(orient.x * a + orient.y * b, orient.y * a - orient.x * b);
  Vec2<double> dR2(orient.x * a - orient.y * b, orient.y * a + orient.x * b);
  vertex[0] = center + dR1;
  vertex[1] = center + dR2;
  vertex[2] = center - dR1;
  vertex[3] = center - dR2;
}

inline void Rect::shift(const Vec2<double> &Delta) {
  center += Delta;
  cal_vertex();
}

void dis_point_edge(double &d, const double *X, const double *Y,
                    int size, int im, double LY);
#endif
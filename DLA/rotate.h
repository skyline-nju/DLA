/* 参考文献：
   多边形的旋转靠接算法，胡华， 计算机研究与发展， 1998.
   An Algorithm For Rotating Polygons Contaction, Hu Hua,
   Computer Research & Development 1998.
*/
#ifndef ROTATE_H
#define ROTATE_H
#include <vector>
#include "vec.h"

struct Vector2D
{
  Vector2D() {}
  Vector2D(const Vec2<double> &AB) : vec(AB) { magnitude = vec.module(); }
  Vector2D(const Vec2<double> &AB, double l) : vec(AB), magnitude(l) {}

  Vec2<double> vec;
  double magnitude;
};

struct Segment
{
  Segment() {}
  Segment(const Vector2D &a, const Vector2D &b, bool increase):
      beg(a), end(b), increasing(increase) {}
  bool intersect(const Vector2D &OP) const;
  Vector2D beg;
  Vector2D end;
  bool increasing;
};

inline bool Segment::intersect(const Vector2D &OP) const {
  if (increasing) {
    return OP.magnitude >= beg.magnitude && OP.magnitude <= end.magnitude;
  } else {
    return OP.magnitude <= beg.magnitude && OP.magnitude >= end.magnitude;
  }
}

struct RotStatus
{
  RotStatus(double max_angle): cos_angle(cos(max_angle)), flag(false) {}
  double cos_angle;
  double sin_angle;
  Vec2<double> contact_point;
  Segment contact_line;
  bool flag;
};

void get_decrease_segment_set(const Vec2<double> &O,
                              const Vec2<double> &M,
                              const Vec2<double> &N,
                              std::vector<Vector2D> &point_set,
                              std::vector<Segment> &line_set);

void get_increase_segment_set(const Vec2<double> &O,
                              const Vec2<double> &M,
                              const Vec2<double> &N,
                              std::vector<Vector2D> &point_set,
                              std::vector<Segment> &line_set);


bool rotate_contact(const Vector2D &P, const Segment &MN,
                    bool clockwise, double &cos_angle);

void rotate_contact(const Vector2D &OP, const Segment &MN,
                    bool CW, RotStatus &status);

void check(const std::vector<Vector2D> &point_set,
           const std::vector<Segment> &segment_set, bool clockwise,
           bool &flag, double &max_cos);

void check(const std::vector<Vector2D> &point_set,
           const std::vector<Segment> &segment_set,
           bool clockwise, RotStatus & status);

void get_min_angle(const Vec2<double> &O,
                   const Vec2<double> *A, int nA,
                   const Vec2<double> *B, int nB,
                   bool clockwise, double &min_angle, bool &flag);

void get_min_angle(const Vec2<double> &O,
                   const std::vector<Vector2D> &point_set_A,
                   const std::vector<Segment> &line_set_A,
                   const Vec2<double> *B, int nB, bool CW,
                   RotStatus &status);

#endif
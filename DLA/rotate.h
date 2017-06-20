/* 参考文献：
   多边形的旋转靠接算法，胡华， 计算机研究与发展， 1998.
   An Algorithm For Rotating Polygons Contaction, Hu Hua,
   Computer Research & Development 1998.
*/
#ifndef ROTATE_H
#define ROTATE_H
#include <vector>
#include "vec.h"

struct Point
{
  Point() {}
  Point(const Vec2<double> &_OP) : vec(_OP) { l = vec.module(); }
  Point(const Vec2<double> &_OP, double lOP) : vec(_OP), l(lOP) {};
  Point(const Vec2<double> &P, const Vec2<double> &O) :vec(P - O) {
    l = vec.module(); }

  Vec2<double> vec;  // vector OP
  double l;         // length of OP
};

struct Segment
{
  Segment() {}
  Segment(const Vec2<double> &OH, double lOM, double lON) :
    H(Point(OH)), l_max(lOM), l_min(lON) {}
  Segment(const Point &pH, double lOM, double lON) :
    H(pH), l_max(lOM), l_min(lON) {}
  Vec2<double> unit_vec() const { return Vec2<double>(-H.vec.y / H.l, H.vec.x / H.l); }
  
  Point H;
  double l_max;
  double l_min;
};

void get_decrease_segment_set(const Vec2<double> &O,
                              const Vec2<double> &M,
                              const Vec2<double> &N,
                              std::vector<Point> &point_set,
                              std::vector<Segment> &line_set);

void get_increase_segment_set(const Vec2<double> &O,
                              const Vec2<double> &M,
                              const Vec2<double> &N,
                              std::vector<Point> &point_set,
                              std::vector<Segment> &line_set);


bool rotate_contact(const Point &P, const Segment &MN, bool clockwise,
                    const Vec2<double> &MN_hat, double &cos_angle);

void check(const std::vector<Point> &point_set,
           const std::vector<Segment> &segment_set, bool clockwise,
           bool reverse, bool &flag, int &idx_P, int &idx_S, double &max_cos);

void get_min_angle(const Vec2<double> &O,
                   const Vec2<double> *A, int nA,
                   const Vec2<double> *B, int nB,
                   bool clockwise, double &min_angle, bool &flag);
#endif
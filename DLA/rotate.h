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
  Point(const Vec2<double> &_pos, const Vec2<double> &O) :
    pos(_pos) { d = pos.distance(O); }

  Vec2<double> pos;
  double d;         // distance to the origin point O
};

struct Segment
{
  Segment() {}
  Segment(const Point &q1, const Point &q2) : p1(q1), p2(q2) {}
  Point p1; // starting point
  Point p2; // ending point
  double d; // distance to the origin point O
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

void get_min_angle(const Vec2<double> &O, const Vec2<double> *A, int nA,
                   const Vec2<double> *B, int nB, bool clockwise);

void get_ratate_angle(const Vec2<double> &O, const Point &P,
                      const Segment &MN, double &angle);
#endif
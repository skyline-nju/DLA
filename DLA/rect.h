#ifndef RECT_H
#define RECT_H
#include "vec.h"
#include "comn.h"
#include "rotate.h"
#include "grid.h"


struct TranStatus
{
  TranStatus(double l) : l_hit(l), flag(false) {}
  void update(double d, bool vertex_to_edge,
              int contact_vertex, int contact_edge, int tag);
  double l_hit;
  int neighbor_tag;
  int cur_neighbor;
  int idx_vertex;
  int idx_edge;
  bool vertex_edge;                 
  bool flag;
};

inline void TranStatus::update(double d, bool vertex_to_edge,
                               int contact_vertex, int contact_edge, int tag) {
  if (d <= l_hit) {
    l_hit = d;
    idx_vertex = contact_vertex;
    idx_edge = contact_edge;
    neighbor_tag = tag;
    vertex_edge = vertex_to_edge;
    flag = true;
 }
}

struct Rect
{
  Rect() {}
  Rect(int tag0): tag(tag0) {}
  Rect(double xc, double yc, double theta, int tag0);
  Rect(double xc, double yc, double ux, double uy, int tag0);
  void cal_vertex();
  void get_mov_dir(int idx0, Vec2<double> &u) const;
  void collide_longitudinal(int idx0, const Vec2<double> &u, const Rect &rect,
                            double l, double &l_hit, bool &flag_collide) const;
  void collide_longitudinal(int idx0, const Vec2<double> &u, const Rect &rect,
                            TranStatus &status) const;
  void collide_transverse(int idx0, const Vec2<double> &u, const Rect &rect,
                          double l, double &l_hit, bool &flag_collide) const;
  void collide_transverse(int idx0, const Vec2<double> &u, const Rect &rect,
                          TranStatus &status) const;
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
  void get_segment_set(bool CW, std::vector<Vector2D>& my_point_set,
                       std::vector<Segment>& my_segment_set,
                       std::vector<int>& my_point_idx,
                       std::vector<int>& my_segment_idx);
  void shift(const Vec2<double> &delta);
  static void output(const std::vector<Rect> &rect);
  static void output(const std::vector<Rect> &rect, const char *f);

  int tag;
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

void dis_point_edge(double &d, int &contact_vertex, int &contact_edge,
                    const double *X, const double *Y, int size,
                    int im, double LY);
#endif
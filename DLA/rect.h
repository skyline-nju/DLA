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
  int idx_vertex;
  int idx_edge;
  bool vertex_edge;                 
  bool flag;
};

struct Rect
{
  Rect(): next(NULL) {}
  Rect(int tag0): tag(tag0), next(NULL) {}
  Rect(double xc, double yc, double theta, int tag0);
  Rect(double xc, double yc, double ux, double uy, int tag0);
  void cal_vertex();
  void get_mov_dir(int idx0, Vec2<double> &u) const;
  void collide_longitudinal(int idx0, const Vec2<double> &u, const Rect &rect,
                            TranStatus &status) const;
  void collide_transverse(int idx0, const Vec2<double> &u, const Rect &rect,
                          TranStatus &status) const;
  void translate(const std::vector<Rect> &cluster, const Cell<Rect> &cell,
                 double lm, int idx0, bool &collided);
  void rotate(const std::vector<Rect> &cluster, const Cell<Rect> &cell, 
              double theta_m, bool clockwise, bool &collided);
  void rotate_around(const std::vector<Rect> &cluster,const Cell<Rect> &cell,
                     int idx_exc, const Vec2<double> &contact_point,
                     double angle, bool &blocked);
  template <typename T>
  void rotate_around_contact_pnt(const std::vector<Rect> &cluster,
                                 const Cell<Rect> &cell,
                                 const T& contact_status);
  template <typename T>
  void tilt(const std::vector<Rect> &cluster, const Cell<Rect> &cell,
            const T& contact_status, double tilting_angle);
  void get_segment_set(bool CW, std::vector<Vector2D> &my_point_set,
                       std::vector<Segment> &my_segment_set) const;
  void get_segment_set(bool CW, std::vector<Vector2D>& my_point_set,
                       std::vector<Segment>& my_segment_set,
                       std::vector<int>& my_point_idx,
                       std::vector<int>& my_segment_idx);
  void shift(const Vec2<double> &delta);
  void output(std::ofstream &out);
  static void output(const std::vector<Rect> &rect);
  static void output(const std::vector<Rect> &rect, const char *f);

  int tag;
  Vec2<double> center;
  Vec2<double> orient;
  Vec2<double> vertex[4];
  Rect *next;

  static double a;
  static double b;
  static double La;
  static double Lb;
  static double Rab;
  static double tilt_angle;
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

void dis_point_edge(double &d, int &contact_vertex, int &contact_edge,
                    const double *X, const double *Y, int size,
                    int im, double LY);

template <class T>
void show_contact(const T &status, std::ofstream &fout1,
                  std::ofstream &fout2, Vec2<double> *my_vertex,
                  const std::vector<Rect> &cluster) {
  if (status.flag) {
    if (status.vertex_edge) {
      fout1 << my_vertex[status.idx_vertex] << endl;
      fout2 << cluster[status.neighbor_tag].vertex[status.idx_edge] << "\t"
            << cluster[status.neighbor_tag].vertex[(status.idx_edge + 1) % 4]
            << endl;
    } else {
      fout1 << cluster[status.neighbor_tag].vertex[status.idx_vertex] << endl;
      fout2 << my_vertex[status.idx_edge] << "\t"
            << my_vertex[(status.idx_edge + 1) % 4] << endl;
    }
  }
}

template<typename T>
void Rect::rotate_around_contact_pnt(const std::vector<Rect>& cluster,
                                     const Cell<Rect> & cell,
                                     const T & contact_status) {
  if (contact_status.flag) {
    double angle_CW, angle_CCW;
    Vec2<double> u1;
    Vec2<double> normal;
    int neighbor = contact_status.neighbor_tag;
    Vec2<double> contact_point;
    int iv = contact_status.idx_vertex;
    int ie = contact_status.idx_edge;
    if (contact_status.vertex_edge) {
      get_mov_dir(iv, u1);
      cluster[neighbor].get_mov_dir(ie, normal);
      contact_point = vertex[iv];
      Vec2<double> u2(-u1.y, u1.x);
      angle_CW = acos(-u1.dot(normal));
      angle_CCW = acos(u2.dot(normal));
    } else {
      cluster[neighbor].get_mov_dir(iv, u1);
      get_mov_dir(ie, normal);
      contact_point = cluster[neighbor].vertex[iv];
      Vec2<double> u2(-u1.y, u1.x);
      angle_CCW = acos(-u1.dot(normal));
      angle_CW = acos(u2.dot(normal));
    }
    bool flag_contact;
    rotate_around(cluster, cell, neighbor, contact_point, angle_CCW, flag_contact);
  }
}

template<typename T>
void Rect::tilt(const std::vector<Rect>& cluster, const Cell<Rect> & cell,
                const T & contact_status, double tilting_angle) {
  if (contact_status.flag && (contact_status.idx_edge == 1 || contact_status.idx_edge == 3)) {
    int iv = contact_status.idx_vertex;
    int ie = contact_status.idx_edge;
    int neighbor = contact_status.neighbor_tag;
    Vec2<double> contact_point;
    Vec2<double> e1;
    Vec2<double> e2;
    if (contact_status.vertex_edge) {
      // new rod's vertex contacts with old rod's edge
      get_mov_dir(iv, e1);
      cluster[neighbor].get_mov_dir((ie + 1) % 4, e2);
      contact_point = vertex[iv];
      double cos_angle = -e1.dot(e2);
      if (cos_angle > 1)
        cos_angle = 1;
      double angle0 = acos(cos_angle);
      bool flag_blocked;
      if (iv == 0 || iv == 2) {
        if (tilting_angle < 0) {
          double d_theta = angle0 + tilting_angle;
          rotate_around(cluster, cell, neighbor, contact_point, d_theta, flag_blocked);
        } else {
          rotate_around(cluster, cell, neighbor, contact_point, angle0, flag_blocked);
          if (!flag_blocked) {
            contact_point = cluster[neighbor].vertex[(ie + 1) % 4];
            rotate_around(cluster, cell, neighbor, contact_point, tilting_angle, flag_blocked);
          }
        }
      } else {
        if (tilting_angle < 0) {
          rotate_around(cluster, cell, neighbor, contact_point, angle0 - 0.5 * PI, flag_blocked);
          if (!flag_blocked) {
            contact_point = cluster[neighbor].vertex[ie];
            rotate_around(cluster, cell, neighbor, contact_point, tilting_angle, flag_blocked);
          }
        } else {
          double d_theta = tilting_angle - (0.5 * PI - angle0);
          rotate_around(cluster, cell, neighbor, contact_point, d_theta, flag_blocked);
        }
      }
    } else {
      // new rod's edge contacts with old rod's vertex
      cluster[neighbor].get_mov_dir(iv, e1);
      get_mov_dir((ie + 1) % 4, e2);
      contact_point = cluster[neighbor].vertex[iv];
      double cos_angle = -e1.dot(e2);
      if (cos_angle > 1)
        cos_angle = 1;
      double angle0 = acos(cos_angle);
      bool flag_blocked;
      if (iv == 0 || iv == 2) {
        if (tilting_angle < 0) {
          rotate_around(cluster, cell, neighbor, contact_point, -angle0, flag_blocked);
          if (!flag_blocked) {
            contact_point = vertex[(ie+1) % 4];
            rotate_around(cluster, cell, neighbor, contact_point, tilting_angle, flag_blocked);
          }
        } else {
          double dtheta = tilting_angle - angle0;
          rotate_around(cluster, cell, neighbor, contact_point, dtheta, flag_blocked);
        }
      } else {
        if (tilting_angle < 0) {
          double d_theta = 0.5 * PI - angle0 + tilting_angle;
          rotate_around(cluster, cell, neighbor, contact_point, d_theta, flag_blocked);
        } else {
          rotate_around(cluster, cell, neighbor, contact_point, 0.5 * PI - angle0, flag_blocked);
          if (!flag_blocked) {
            contact_point = vertex[ie];
            rotate_around(cluster, cell, neighbor, contact_point, tilting_angle, flag_blocked);
          }
        }
      }
    }
  }
}


#endif
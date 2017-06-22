#include "rotate.h"

using namespace std;

void get_decrease_segment_set(const Vec2<double>& O,
                              const Vec2<double>& M,
                              const Vec2<double>& N,
                              vector<Vector2D>& point_set,
                              vector<Segment>& line_set) {
  Vec2<double> MN = N - M;
  Vec2<double> MO = O - M;
  double MN_dot_MO = MN.dot(MO);
  if (MN_dot_MO > 0) {
    Vector2D vOM(-MO);
    double MN_square = MN.square();
    point_set.push_back(vOM);
    if (MN_dot_MO >= MN_square) {
      Vector2D vON(N - O);
      line_set.push_back(Segment(vOM, vON, false));
      point_set.push_back(vON);
    } else {
      Vec2<double> OH = MN * (MN_dot_MO / MN_square) - MO;
      Vector2D vOH(OH);
      line_set.push_back(Segment(vOM, vOH, false));
    }
  }
}

void get_increase_segment_set(const Vec2<double>& O,
                              const Vec2<double>& M,
                              const Vec2<double>& N,
                              vector<Vector2D>& point_set,
                              vector<Segment>& line_set) {
  Vec2<double> NM = M - N;
  Vec2<double> NO = O - N;
  double NM_dot_NO = NM.dot(NO);
  if (NM_dot_NO > 0) {
    Vector2D vON(-NO);
    double NM_square = NM.square();
    point_set.push_back(vON);
    if (NM_dot_NO >= NM_square) {
      Vector2D vOM(M - O);
      line_set.push_back(Segment(vOM, vON, true));
      point_set.push_back(vOM);
    } else {
      Vec2<double> OH = NM * (NM_dot_NO / NM_square) - NO;
      Vector2D vOH(OH);
      line_set.push_back(Segment(vOH, vON, true));
    }
  }
}

void get_segment_set(const Vec2<double> &O, const Vec2<double> *vertex, int n,
                     vector<Vector2D> &point_set,
                     vector<Segment> &line_set, bool CW) {
  point_set.reserve(2 * n);
  line_set.reserve(n);
  if (CW) {
    for (int i = 0; i < n; i++) {
      int j = i - 1;
      if (j < 0) j += n;
      get_increase_segment_set(O, vertex[i], vertex[j], point_set, line_set);
    }
  } else {
    for (int i = 0; i < n; i++) {
      int j = i - 1;
      if (j < 0) j += n;
      get_decrease_segment_set(O, vertex[i], vertex[j], point_set, line_set);
    }
  }

}

void rotate_contact(const Vector2D &OP, const Segment &MN, bool CW,
                    RotStatus &status) {
  if (MN.intersect(OP)) {
    Vec2<double> vMN = MN.end.vec - MN.beg.vec;
    double rr = OP.magnitude * OP.magnitude;
    double a = vMN.square();
    double b = 2 * MN.beg.vec.dot(vMN);
    double c = MN.beg.magnitude * MN.beg.magnitude - rr;
    double Delta = sqrt(b * b - 4 * a * c);
    double lambda = (-b + Delta) / (2 * a);
    if (lambda > 1)
      lambda = (-b - Delta) / (2 * a);
    Vec2<double> vOQ = MN.beg.vec + lambda * vMN;
    double cross_product = OP.vec.cross(vOQ);
    if ((CW && cross_product < 0) || (!CW && cross_product > 0)) {
      double cos_angle = vOQ.dot(OP.vec) / rr;
      if (cos_angle >= status.cos_angle) {
        status.cos_angle = cos_angle;
        //status.sin_angle = cross_product / rr;
        status.contact_point = vOQ;
        status.flag = true;
      }
    }
  }
}

void for_each_pair(const vector<Vector2D>& point_set,
                   const vector<Segment>& segment_set,
                   bool clockwise, RotStatus & status) {
  for (auto line = segment_set.cbegin(); line != segment_set.cend(); ++line) {
    for (auto pnt = point_set.cbegin(); pnt != point_set.cend(); ++pnt) {
      rotate_contact(*pnt, *line, clockwise, status);
    }
  }
}

void get_min_angle(const Vec2<double> &O, const vector<Vector2D> &point_set_A,
                   const vector<Segment> &line_set_A,
                   const Vec2<double> *B, int nB, bool CW,
                   RotStatus &status) {
  vector<Vector2D> point_set_B;
  vector<Segment> line_set_B;
  get_segment_set(O, B, nB, point_set_B, line_set_B, CW);
  for_each_pair(point_set_A, line_set_B, CW, status);
  for_each_pair(point_set_B, line_set_A, !CW, status);
}


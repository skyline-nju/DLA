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

bool rotate_contact(const Vector2D & OP, const Segment & MN, bool CW,
                    double & cos_angle) {
  bool flag = false;
  if (MN.intersect(OP)) {
    Vec2<double> vMN = MN.end.vec - MN.beg.vec;
    double rr = OP.magnitude * OP.magnitude;
    double a = vMN.square();
    double b = 2 * MN.beg.vec.dot(vMN);
    double c = MN.beg.magnitude * MN.beg.magnitude - rr;
    double Delta2 = b * b - 4 * a * c;
    if (Delta2 >= 0) {
      double Delta = sqrt(Delta2);
      double t = (-b + Delta) / (2 * a);
      if (t > 1)
        t = (-b - Delta) / (2 * a);
      Vec2<double> vOQ = MN.beg.vec + t * vMN;
      double cross = OP.vec.cross(vOQ);
      if ((CW && cross < 0) || (!CW && cross > 0)) {
        flag = true;
        cos_angle = vOQ.dot(OP.vec) / rr;
      }
    }
  }
  return flag;
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
      if (cos_angle > status.cos_angle) {
        status.cos_angle = cos_angle;
        status.contact_point = vOQ;
        status.flag = true;
      }
    }
  }
}

void check(const vector<Vector2D>& point_set,
           const vector<Segment>& segment_set, bool clockwise,
           bool &flag, double &max_cos) {
  flag = false;
  for (int i = 0, nS = segment_set.size(); i < nS; i++) {
    for (int j = 0, nP = point_set.size(); j < nP; j++) {
      double cos_angle;
      if (rotate_contact(point_set[j], segment_set[i],
                         clockwise, cos_angle)) {
        if (flag && max_cos < cos_angle) {
          max_cos = cos_angle;
        } else if (!flag) {
          flag = true;
          max_cos = cos_angle;
        }
      } 
    }
  }
}

void check(const vector<Vector2D>& point_set,
           const vector<Segment>& segment_set,
           bool clockwise, RotStatus & status) {
  for (int i = 0, nS = segment_set.size(); i < nS; i++) {
    for (int j = 0, nP = point_set.size(); j < nP; j++) {
      rotate_contact(point_set[i], segment_set[j], clockwise, status);
    }
  }
}

void get_min_angle(const Vec2<double>& O,
                   const Vec2<double>* A, int nA,
                   const Vec2<double>* B, int nB,
                   bool clockwise, double &min_angle, bool &flag) {
  vector<Vector2D> point_set_A;
  vector<Vector2D> point_set_B;
  vector<Segment> line_set_A;
  vector<Segment> line_set_B;
  if (clockwise) {
    for (int i = 0; i < nA; i++) {
      int j = i - 1;
      if (j < 0) j += nA;
      get_decrease_segment_set(O, A[i], A[j], point_set_A, line_set_A);
    }
    for (int i = 0; i < nB; i++) {
      int j = i - 1;
      if (j < 0) j += nB;
      get_increase_segment_set(O, B[i], B[j], point_set_B, line_set_B);
    }
  } else {
    for (int i = 0; i < nA; i++) {
      int j = i - 1;
      if (j < 0) j += nA;
      get_increase_segment_set(O, A[i], A[j], point_set_A, line_set_A);
    }
    for (int i = 0; i < nB; i++) {
      int j = i - 1;
      if (j < 0) j += nB;
      get_decrease_segment_set(O, B[i], B[j], point_set_B, line_set_B);
    }
  }
  bool flag_A;
  double max_cos_A;
  bool flag_B;
  double max_cos_B;
  check(point_set_A, line_set_B, clockwise, flag_A, max_cos_A);
  check(point_set_B, line_set_A, !clockwise, flag_B, max_cos_B);
  if (!flag_A && !flag_B) {
    flag = false;
  } else if (flag_A && flag_B) {
    flag = true;
    if (max_cos_A > max_cos_B) {
      min_angle = acos(max_cos_A);
    } else {
      min_angle = acos(max_cos_B);
    }
  } else if (flag_A) {
    flag = true;
    min_angle = acos(max_cos_A);
  } else {
    flag = true;
    min_angle = acos(max_cos_B);
  }
  if (flag) {
    if (clockwise)
      min_angle = -min_angle;
  }
}

void get_min_angle(const Vec2<double> &O, const vector<Vector2D> &point_set_A,
                   const vector<Segment> &line_set_A,
                   const Vec2<double> *B, int nB, bool CW,
                   RotStatus &status) {
  vector<Vector2D> point_set_B;
  vector<Segment> line_set_B;
  point_set_B.reserve(2 * nB);
  line_set_B.reserve(nB);
  if (CW) {
    for (int i = 0; i < nB; i++) {
      int j = i - 1;
      if (j < 0) j += nB;
      get_increase_segment_set(O, B[i], B[j], point_set_B, line_set_B);
    }
  } else {
    for (int i = 0; i < nB; i++) {
      int j = i - 1;
      if (j < 0) j += nB;
      get_decrease_segment_set(O, B[i], B[j], point_set_B, line_set_B);
    }
  }
  check(point_set_A, line_set_B, CW, status);
  check(point_set_B, line_set_A, !CW, status);
}


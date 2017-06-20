#include "rotate.h"

using namespace std;

void get_decrease_segment_set(const Vec2<double>& O,
                              const Vec2<double>& M,
                              const Vec2<double>& N,
                              vector<Point>& point_set,
                              vector<Segment>& line_set) {
  Vec2<double> MN = N - M;
  Vec2<double> MO = O - M;
  double MN_dot_MO = MN.dot(MO);
  if (MN_dot_MO > 0) {
    Point pM(-MO);
    double MN_square = MN.square();
    Vec2<double> OH = MN * (MN_dot_MO / MN_square) - MO;
    Point pH(OH);
    point_set.push_back(pM);
    if (MN_dot_MO >= MN_square) {
      Point pN(N, O);
      line_set.push_back(Segment(pH, pM.l, pN.l));
      point_set.push_back(pN);
    } else {
      line_set.push_back(Segment(pH, pM.l, pH.l));
    }
  }
}

void get_increase_segment_set(const Vec2<double>& O,
                              const Vec2<double>& M,
                              const Vec2<double>& N,
                              vector<Point>& point_set,
                              vector<Segment>& line_set) {
  get_decrease_segment_set(O, N, M, point_set, line_set);
}

bool rotate_contact(const Point & P, const Segment & MN, bool clockwise,
                    const Vec2<double>& MN_hat, double & cos_angle) {
  bool flag = false;
  if (P.l >= MN.l_min && P.l <= MN.l_max) {
    double rr = P.l * P.l;
    double l_HQ = sqrt(rr - MN.H.l * MN.H.l);
    Vec2<double> OQ = MN.H.vec + MN_hat * l_HQ;
    cos_angle = OQ.dot(P.vec) / rr;
    double product = P.vec.cross(OQ);
    if ((clockwise && product < 0) || (!clockwise && product > 0)) {
      flag = true;
      if (clockwise) {
        cout << "CW: " << P.vec << "\t" << cos_angle << endl;
        cout << "MN: " << MN_hat << endl;
      } else {
        cout << "CCW: " << P.vec << "\t" << cos_angle << endl;
        cout << "MN: " << MN_hat << endl;
      }
      cout << endl;
    }
  } 
  return flag;
}

void check(const vector<Point>& point_set,
           const vector<Segment>& segment_set, bool clockwise, bool reverse,
           bool &flag, int &idx_P, int &idx_S, double &max_cos) {
  int nP = point_set.size();
  int nS = segment_set.size();
  flag = false;
  for (int i = 0; i < nS; i++) {
    Vec2<double> unit_vec = segment_set[i].unit_vec();
    if (reverse) unit_vec = -unit_vec;
    for (int j = 0; j < nP; j++) {
      double cos_angle;
      if (rotate_contact(point_set[j], segment_set[i],
                         clockwise, unit_vec, cos_angle)) {
        if (flag && max_cos < cos_angle) {
          max_cos = cos_angle;
          idx_P = j;
          idx_S = i;
        } else if (!flag) {
          flag = true;
          max_cos = cos_angle;
          idx_P = j;
          idx_S = i;
        }
      } 
    }
  }
  cout << "max cos = " << max_cos << endl;
}

void get_min_angle(const Vec2<double>& O,
                   const Vec2<double>* A, int nA,
                   const Vec2<double>* B, int nB,
                   bool clockwise, double &min_angle, bool &flag) {
  vector<Point> point_set_A;
  vector<Point> point_set_B;
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
  int idx_Pnt_A;
  int idx_Line_A;
  double max_cos_A;
  bool flag_B;
  int idx_Pnt_B;
  int idx_Line_B;
  double max_cos_B;
  bool reverse = clockwise == true ? false : true;
  check(point_set_A, line_set_B, clockwise, reverse,
        flag_A, idx_Pnt_A, idx_Line_B, max_cos_A);
  check(point_set_B, line_set_A, !clockwise, reverse,
        flag_B, idx_Pnt_B, idx_Line_A, max_cos_B);
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


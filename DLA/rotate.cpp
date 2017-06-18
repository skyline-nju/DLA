#include "rotate.h"

using namespace std;

void get_decrease_segment_set(const Vec2<double>& O,
                              const Vec2<double>& M,
                              const Vec2<double>& N,
                              vector<Point>& point_set,
                              vector<Segment>& line_set) {
  Vec2<double> MN = N - M;
  Vec2<double> MO = O - M;
  double res_dot = MN.dot(MO);
  if (res_dot > 0) {
    double L_MN = MN.module();
    double L_MQ = res_dot / L_MN;
    if (L_MQ > L_MN) {
      Point pM(M, O);
      Point pN(N, O);
      point_set.push_back(pM);
      point_set.push_back(pN);
      line_set.push_back(Segment(pM, pN));
    } else {
      Point pM(M, O);
      Point pQ(M + MN * L_MQ / L_MN, O);
      point_set.push_back(pM);
      line_set.push_back(Segment(pM, pQ));
      if (L_MN == L_MQ) {
        Point pN(N, O);
        point_set.push_back(pN);
      }
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

void get_min_angle(const Vec2<double>& O, const Vec2<double>* A, int nA,
                   const Vec2<double>* B, int nB, bool clockwise) {
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
      if (j < 0) j += nA;
      get_increase_segment_set(O, B[i], B[j], point_set_B, line_set_B);
    }
  } else {

  }
}

void get_ratate_angle(const Vec2<double>& O, const Point & P,
                      const Segment & MN, double & angle) {
}

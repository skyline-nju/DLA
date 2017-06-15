#include "particle.h"

using namespace std;

double Rod::a = 7;
double Rod::b = 1;

void Rod::collide(const Rod &rod, double vx, double vy, int first_vertex,
                  double l, double &lhit, bool &flag) {
  double L_parallel;  // edge parallel to the moving direction.
  double L_vertical;
  if (first_vertex == 0 || first_vertex == 2) {
    L_parallel = 2 * a;
    L_vertical = 2 * b;
  } else {
    L_parallel = 2 * b;
    L_vertical = 2 * a;
  }
  bool cross_flag = false;
  for (int i = 0; i < 4; i++) {
    double lambda = -vy * (rod.vertex_x[i] - vertex_x[first_vertex])
                    +vx * (rod.vertex_y[i] - vertex_y[first_vertex]);
    if (lambda >= 0 && lambda <= L_vertical) {
      cross_flag = true;
      break;
    }
  }
  if (!cross_flag) {
    lhit = l;
    flag = false;
  } else {
    for (int i = 0; i < 4; i++) {

    }
  }

}

void Rod::check_overlape(int first_vertex, double vx, double vy, const Rod & rod) {
  double xp[2];
  xp[0] = vx * vertex_x[first_vertex] + vy * vertex_y[first_vertex];
  int second_vertex = first_vertex == 3 ? 0 : first_vertex + 1;
  xp[1] = vx * vertex_x[second_vertex] + vy * vertex_y[second_vertex];
  double Xp[4];
  bool flag = false;
  int sign = 0;
  for (int i = 0; i < 4; i++) {
    Xp[i] = vx * rod.vertex_x[i] + vy * rod.vertex_y[i];
    if (Xp[i] >= xp[0] && Xp[i] <= xp[1]) {
      flag = true;
      break;
    } else if (sign == 0) {
      sign = Xp[i] < xp[0] ? -1 : 1;
    } else {
      if ((sign < 0 && Xp[i] >= xp[0]) || (sign > 0 && Xp[i] <= xp[1])) {
        flag = true;
        break;
      }
    }
  }
}

void Rod::cal_vertex() {
  double dx1 = ux * a + uy * b;
  double dx2 = ux * a - uy * b;
  double dy1 = uy * a - ux * b;
  double dy2 = uy * a + ux * b;
  vertex_x[0] = xc + dx1;
  vertex_x[1] = xc + dx2;
  vertex_x[2] = xc - dx1;
  vertex_x[3] = xc - dx2;
  vertex_y[0] = yc + dy1;
  vertex_y[1] = yc + dy2;
  vertex_y[2] = yc - dy1;
  vertex_y[3] = yc - dy2;
}


#include "particle.h"

using namespace std;

double Disk::sigma = 2;
double Disk::rp = 1;
double Rect::a = 7;
double Rect::b = 1;
double Rect::La = 2 * Rect::a;
double Rect::Lb = 2 * Rect::b;


void Disk::collide(const Disk & p1, double ux, double uy,
                   double l, double &l_hit, bool &flag) const {
  double dx = x - p1.x;
  double dy = y - p1.y;
  double B = 2 * (ux * dx + uy * dy);
  double C = dx * dx + dy * dy - sigma * sigma;
  double Delta = B * B - 4 * C;
  if (Delta < 0) {
    flag = false;
  } else {
    double sqrt_Delta = sqrt(Delta);
    l_hit = B + sqrt_Delta < 0 ?
      0.5 * (-B - sqrt_Delta) : 0.5 * (-B + sqrt_Delta);
    if (l_hit > l || l_hit < 0) {
      flag = false;
    } else {
      flag = true;
    }
  }
}

void Disk::output(const vector<Disk> &cluster) {
  char fname[100];
  snprintf(fname, 100, "N%d.dat", cluster.size());
  ofstream fout(fname);
  for (int i = 0; i < cluster.size(); i++)
    fout << cluster[i].x << "\t" << cluster[i].y << endl;
}

void Disk::output_xyz(const vector<Disk> &cluster) {
  char fname[100];
  snprintf(fname, 100, "N%d.xyz", cluster.size());
  ofstream fout(fname);
  fout << cluster.size() << endl;
  fout << "type\tx\ty\n";
  for (int i = 0; i < cluster.size(); i++)
    fout << "X\t" << cluster[i].x << "\t" << cluster[i].y << endl;
}

Rect::Rect(double xc, double yc, double theta): center(xc, yc) {
  double angle = theta / 180.0 * PI;
  orient.x = cos(angle);
  orient.y = sin(angle);
  cal_vertex();
}

Rect::Rect(double xc, double yc, double ux, double uy):
           center(xc, yc), orient(ux, uy) {
  cal_vertex();
}

void Rect::collide1(int idx0, const Vec2<double>& u, const Rect & rect,
                    double l, double &l_hit, bool &flag_collide) const {
  double LY = idx0 == 0 || idx0 == 2 ? Lb : La;
  double X[4];
  bool flag_all_negative = true;
  for (int i = 0; i < 4; i++) {
    X[i] = u.dot(rect.vertex[i] - vertex[idx0]);
    if (X[i] > 0)
      flag_all_negative = false;
  }
  if (flag_all_negative) {
    l_hit = l;
    flag_collide = false;
  } else {
    bool flag_all_gt_L = true;
    bool flag_all_lt_zero = true;
    Vec2<double> v(-u.y, u.x);
    double Y[4];
    for (int i = 0; i < 4; i++) {
      Y[i] = v.dot(rect.vertex[i] - vertex[idx0]);
      if (Y[i] >= 0)
        flag_all_lt_zero = false;
      if (Y[i] <= LY)
        flag_all_gt_L = false;
    }
    if (flag_all_lt_zero || flag_all_gt_L) {
      l_hit = l;
      flag_collide = false;
    } else {
      int im = get_idx_min(X, 4);
      if (Y[im] >= 0 && Y[im] <= LY) {
        l_hit = X[im];
      } else {
        dis_point_edge(l_hit, X, Y, 4, im, LY);
      }
      flag_collide = l_hit <= l ? true : false;
    }
  }
}

void Rect::collide2(int idx0, const Vec2<double>& u, const Rect & rect,
                    double l, double &l_hit, bool &flag_collide) const {
  double LY = La;
  double X[4];
  double Y[4];
  bool flag_valid[4];
  int n_valid = 0;
  Vec2<double> v(-u.y, u.x);
  for (int i = 0; i < 4; i++) {
    Vec2<double> dR(rect.vertex[i] - vertex[idx0]);
    X[i] = u.dot(dR);
    Y[i] = v.dot(dR);
    if (X[i] > 0 && Y[i] >=0 && Y[i] <= LY) {
      flag_valid[i] = true;
      n_valid++;
    } else {
      flag_valid[i] = false;
    }
  }
  if (n_valid == 0) {
    l_hit = l;
    flag_collide = false;
  } else if (n_valid >= 3) {
    int im = get_idx_min(X, 4);
    l_hit = X[im];
    flag_collide = l_hit <= l ? true : false;
  } else {
    int im = get_idx_min(X, 4);
    if (flag_valid[im]) {
      l_hit = X[im];
      flag_collide = l_hit <= l ? true : false;
    } else {
      dis_point_edge(l_hit, X, Y, 4, im, LY);
      flag_collide = l_hit <= l ? true : false;
    }
  }
}

void Rect::output(const vector<Rect> &rect, const char *filename) {
  ofstream fout(filename);
  for (int i = 0; i < rect.size(); i++) {
    double x = rect[i].vertex[3].x;
    double y = rect[i].vertex[3].y;
    double theta = atan2(rect[i].orient.y, rect[i].orient.x) / PI * 180;
    fout << x << "\t" << y << "\t" << theta << endl;
  }
  fout.close();
}

void Rect::output(const vector<Rect> &rect) {
  int size = rect.size();
  char fname[100];
  snprintf(fname, 100, "rect_%d.dat", size);
  output(rect, fname);
}

void Rect::get_mov_dir(int idx0, Vec2<double> &u) const {
  switch (idx0) {
    case 0: {
      u.x = orient.x;
      u.y = orient.y;
      break;
    } case 1: {
      u.x = -orient.y;
      u.y = orient.x;
      break;
    } case 2: {
      u.x = -orient.x;
      u.y = -orient.y;
      break;
    } case 3: {
      u.x = orient.y;
      u.y = -orient.x;
      break;
    } default: {
      break;
    }
  }
}

void dis_point_edge(double &d, const double *X, const double *Y,
                    int size, int im, double LY) {
  int ipre = im - 1;
  int inxt = im + 1;
  if (ipre == -1)
    ipre = size - 1;
  else if (inxt == size)
    inxt = 0;
  int i2;
  if (Y[im] > LY) {
    i2 = Y[ipre] <= LY ? ipre : inxt;
  } else {
    i2 = Y[ipre] >= 0 ? ipre : inxt;
  }
  double k = (X[i2] - X[im]) / (Y[i2] - Y[im]);
  if (k >= 0) {
    d = k * (-Y[im]) + X[im];
  } else {
    d = k * (LY - Y[im]) + X[im];
  }
}
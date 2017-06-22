#include "particle.h"
#include <functional>

using namespace std;
ofstream fout("vertex.txt");
double Disk::sigma = 2;
double Disk::rp = 1;
double Rect::a = 7;
double Rect::b = 1;
double Rect::La = 2 * Rect::a;
double Rect::Lb = 2 * Rect::b;
double Rect::Rab = sqrt(Rect::a * Rect::a + Rect::b * Rect::b);
void (Rect::*Rect::collide_wrapper)(int, const Vec2<double> &, const Rect &,
                                    double, double &, bool &) const;

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

void Rect::collide_longitudinal(int idx0, const Vec2<double>& u, const Rect & rect,
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

void Rect::collide_longitudinal(int idx0, const Vec2<double>& u,
                                const Rect &rect, TranStatus &status) const {
  double LY = Lb;
  double X[4];
  bool all_X_negative = true;
  for (int i = 0; i < 4; i++) {
    X[i] = u.dot(rect.vertex[i] - vertex[idx0]);
    if (X[i] > 0)
      all_X_negative = false;
  }
  if (!all_X_negative) {
    bool all_Y_gt_L = true;
    bool all_Y_lt_zero = true;
    Vec2<double> v(-u.y, u.x);
    double Y[4];
    for (int i = 0; i < 4; i++) {
      Y[i] = v.dot(rect.vertex[i] - vertex[idx0]);
      if (Y[i] >= 0)
        all_Y_lt_zero = false;
      if (Y[i] <= LY)
        all_Y_gt_L = false;
    }
    if (!all_Y_lt_zero && !all_Y_gt_L) {
      int im = get_idx_min(X, 4);
      if (Y[im] >= 0 && Y[im] <= LY ) {
        if (X[im] <= status.l_hit) {
          status.l_hit = X[im];
          status.flag = true;
          status.contact_point = rect.vertex[im];
        }
      } else {
        double d;
        int idx;
        dis_point_edge(d, idx, X, Y, 4, im, LY);
        idx += idx0;
        if (d <= status.l_hit) {
          status.l_hit = d;
          status.flag = true;
          status.contact_point = vertex[idx] + u * d;
        }
      }
    }
  }
}

void Rect::collide_transverse(int idx0, const Vec2<double>& u, const Rect & rect,
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

void Rect::collide_transverse(int idx0, const Vec2<double>& u,
                              const Rect &rect, TranStatus &status) const {
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
    if (X[i] > 0 && Y[i] >= 0 && Y[i] <= LY) {
      flag_valid[i] = true;
      n_valid++;
    } else {
      flag_valid[i] = false;
    }
  }
  if (n_valid >= 3) {
    int im = get_idx_min(X, 4);
    if (X[im] <= status.l_hit) {
      status.flag = true;
      status.contact_point = rect.vertex[im];
      status.l_hit = X[im];
      if (status.l_hit < 0)
        cout << "L = " << status.l_hit << endl;
    }
  } else if (n_valid > 0) {
    int im = get_idx_min(X, 4);
    if (flag_valid[im]) {
      if (X[im] <= status.l_hit) {
        status.flag = true;
        status.contact_point = rect.vertex[im];
        status.l_hit = X[im];
      }
    } else {
      double d;
      int delta_n;
      dis_point_edge(d, delta_n, X, Y, 4, im, LY);
      if (d <= status.l_hit) {
        status.flag = true;
        status.contact_point = vertex[(idx0 + delta_n) % 4] + u * d;
        status.l_hit = d;
      }
    }
  }
}

void Rect::translate(const std::vector<Rect>& cluster, double lm,
                     int idx0, bool & collided) {
  collided = false;
  double l = lm;
  Vec2<double> u;
  get_mov_dir(idx0, u);
  for (int i = 0, size = cluster.size(); i < size; i++) {
    bool flag = false;
    double l_hit;
    (this->*collide_wrapper)(idx0, u, cluster[i], lm, l_hit, flag);
    if (flag) {
      collided = true;
      if (l > l_hit) l = l_hit;
    }
  }
  center += u * l;
  cal_vertex();
}

void Rect::translate2(const std::vector<Rect>& cluster, double lm,
                      int idx0, bool & collided) {
  collided = false;
  TranStatus status(lm);
  Vec2<double> u;
  get_mov_dir(idx0, u);
  for (int i = 0, size = cluster.size(); i < size; i++) {
    if (idx0 == 0 || idx0 == 2) {
      collide_longitudinal(idx0, u, cluster[i], status);
    } else {
      collide_transverse(idx0, u, cluster[i], status);
    }
  }
  if (status.flag) {
    collided = true;
    fout << status.contact_point << endl;
  } else {
    collided = false;
  }
  center += u * status.l_hit;
  cal_vertex();
}

void Rect::translate(const std::vector<Rect>& cluster, const Cell & cell,
                     double lm, int idx0, bool & collided) {
  collided = false;
  double l = lm;
  Vec2<double> u;
  get_mov_dir(idx0, u);
  int my_col = cell.get_col(center.x);
  int my_row = cell.get_row(center.y);
  if (!cell.is_isolate(my_col, my_row)) {
    for (int row = my_row - 1; row <= my_row + 1; row++) {
      for (int col = my_col - 1; col <= my_col + 1; col++) {
        int idx = cell.get_idx(col, row);
        auto beg = cell.tag[idx].cbegin();
        auto end = cell.tag[idx].cend();
        for (auto iter = beg; iter != end; ++iter) {
          bool flag = false;
          double l_hit;
          (this->*collide_wrapper)(idx0, u, cluster[*iter], lm, l_hit, flag);
          if (flag) {
            collided = true;
            if (l > l_hit) l = l_hit;
          }
        }
      }
    }
  }
  center += u * l;
  cal_vertex();
}

void Rect::rotate(const std::vector<Rect>& cluster, double theta_m,
                  bool CW, bool & collided) {
  double theta;
  RotStatus status(theta_m);
  vector<Vector2D> point_set;
  vector<Segment> segment_set;
  get_segment_set(CW, point_set, segment_set);
  for (int i = 0, size = cluster.size(); i < size; i++) {
    get_min_angle(center, point_set, segment_set,
                  cluster[i].vertex, 4, CW, status);
  }
  if (status.flag) {
    collided = true;
    theta = acos(status.cos_angle);
    fout << status.contact_point << endl;
  } else {
    collided = false;
    theta = theta_m;
  }
  if (CW) theta = -theta;
  orient.rotate(theta);
  cal_vertex();
}

void Rect::rotate(const std::vector<Rect>& cluster, const Cell &cell,
                  double theta_m, bool CW, bool & collided) {
  collided = false;
  double theta = theta_m;
  int my_col = cell.get_col(center.x);
  int my_row = cell.get_row(center.y);
  if (!cell.is_isolate(my_col, my_row)) {
    RotStatus status(theta_m);
    vector<Vector2D> point_set;
    vector<Segment> segment_set;
    get_segment_set(CW, point_set, segment_set);
    auto f = bind(get_min_angle, cref(center), cref(point_set), cref(segment_set),
                  placeholders::_1, 4, CW, ref(status));
    cell.for_each_neighbor(my_col, my_row, cluster, f);
    if (status.flag) {
      collided = true;
      theta = acos(status.cos_angle);
    }
  }
  if (CW) theta = -theta;
  orient.rotate(theta);
  cal_vertex();
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
  case 0:
    u.x = orient.x;
    u.y = orient.y;
    collide_wrapper = &Rect::collide_longitudinal;
    break;
  case 1:
    u.x = -orient.y;
    u.y = orient.x;
    collide_wrapper = &Rect::collide_transverse;
    break;
  case 2:
    u.x = -orient.x;
    u.y = -orient.y;
    collide_wrapper = &Rect::collide_longitudinal;
    break;
  case 3:
    u.x = orient.y;
    u.y = -orient.x;
    collide_wrapper = &Rect::collide_transverse;
    break; 
  default:
    break;
  }
}

void Rect::get_segment_set(bool CW, vector<Vector2D> &my_point_set,
                           vector<Segment> &my_segment_set) const {
  my_point_set.reserve(4);
  my_segment_set.reserve(4);
  if (CW) {
    for (int i = 0; i < 4; i++) {
      int j = i == 0 ? 3 : i - 1;
      Vector2D OM(vertex[i] - center, Rab);
      double d = i == 0 || i == 2 ? b : a;
      Vector2D ON((vertex[i] + vertex[j]) * 0.5 - center, d);
      my_point_set.push_back(OM);
      my_segment_set.push_back(Segment(OM, ON, false));
    }
  } else {
    for (int i = 0; i < 4; i++) {
      int j = i == 0 ? 3 : i - 1;
      double d = i == 0 || i == 2 ? b : a;
      Vector2D OM((vertex[i] + vertex[j]) * 0.5 - center, d);
      Vector2D ON(vertex[j] - center, Rab);
      my_point_set.push_back(ON);
      my_segment_set.push_back(Segment(OM, ON, true));
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

void dis_point_edge(double &d, int &delta_n, const double *X, const double *Y,
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
    delta_n = 0;
  } else {
    d = k * (LY - Y[im]) + X[im];
    delta_n = 1;
  }
}

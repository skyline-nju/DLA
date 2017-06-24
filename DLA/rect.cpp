#include "rect.h"

using namespace std;
ofstream fout2("vertex.txt");
ofstream fout3("edge.txt");

double Rect::a = 7;
double Rect::b = 1;
double Rect::La = 2 * Rect::a;
double Rect::Lb = 2 * Rect::b;
double Rect::Rab = sqrt(Rect::a * Rect::a + Rect::b * Rect::b);
void (Rect::*Rect::collide_wrapper)(int, const Vec2<double> &, const Rect &,
                                    double, double &, bool &) const;

Rect::Rect(double xc, double yc, double theta, int tag0):
           center(xc, yc), tag(tag0) {
  double angle = theta / 180.0 * PI;
  orient.x = cos(angle);
  orient.y = sin(angle);
  cal_vertex();
}

Rect::Rect(double xc, double yc, double ux, double uy, int tag0):
           center(xc, yc), orient(ux, uy), tag(tag0) {
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
        status.update(X[im], false, im, idx0, rect.tag);
      } else {
        double d;
        int contact_vertex = idx0;
        int contact_edge;
        dis_point_edge(d, contact_vertex, contact_edge, X, Y, 4, im, LY);
        status.update(d, true, contact_vertex, contact_edge, rect.tag);
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
    status.update(X[im], false, im, idx0, rect.tag);
  } else if (n_valid > 0) {
    int im = get_idx_min(X, 4);
    if (flag_valid[im]) {
      status.update(X[im], false, im, idx0, rect.tag);
    } else {
      double d;
      int contact_vertex = idx0;
      int contact_edge;
      dis_point_edge(d, contact_vertex, contact_edge, X, Y, 4, im, LY);
      status.update(d, true, contact_vertex, contact_edge, rect.tag);
    }
  }
}

void Rect::translate(const std::vector<Rect>& cluster, double lm,
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
  } else {
    collided = false;
  }
  center += u * status.l_hit;
  cal_vertex();
  //if (status.flag) {
  //  cout << cluster.size() << "\t" << status.neighbor_tag << endl;
  //  if (status.vertex_edge) {
  //    fout2 << vertex[status.idx_vertex] << endl;
  //    fout3 << cluster[status.neighbor_tag].vertex[status.idx_edge] << "\t"
  //          << cluster[status.neighbor_tag].vertex[(status.idx_edge + 1) % 4] << endl;
  //  }
  //  else {
  //    fout2 << cluster[status.neighbor_tag].vertex[status.idx_vertex] << endl;
  //    fout3 << vertex[status.idx_edge] << "\t"
  //      << vertex[(status.idx_edge + 1) % 4] << endl;
  //  }
  //}
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
    //auto lambda = [&, idx0]
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
  vector<int> point_idx;
  vector<int> segment_idx;
  get_segment_set(CW, point_set, segment_set, point_idx, segment_idx);
  for (int i = 0, size = cluster.size(); i < size; i++) {
    //get_min_angle(center, point_set, segment_set,
    //  cluster[i].vertex, 4, CW, status);
    get_min_angle(point_set, segment_set, point_idx, segment_idx,
                  center, cluster[i].vertex, 4, CW, i, status);
  }
  if (status.flag) {
    collided = true;
    theta = acos(status.cos_angle);
    //fout << status.contact_point << endl;
  } else {
    collided = false;
    theta = theta_m;
  }
  if (CW) theta = -theta;
  orient.rotate(theta);
  cal_vertex();
  if (status.flag) {
    if (status.vertex_edge) {
      fout2 << vertex[status.idx_vertex] << endl;
      fout3 << cluster[status.neighbor_tag].vertex[status.idx_edge] << "\t"
        << cluster[status.neighbor_tag].vertex[(status.idx_edge + 1) % 4] << endl;
    } else {
      fout2 << cluster[status.neighbor_tag].vertex[status.idx_vertex] << endl;
      fout3 << vertex[status.idx_edge] << "\t"
        << vertex[(status.idx_edge + 1) % 4] << endl;
    }
  }
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
    auto lambda = [&, CW](const Vec2<double> *B) {
      get_min_angle(center, point_set, segment_set, B, 4, CW, status);
    };
    cell.for_each_neighbor(my_col, my_row, cluster, lambda);
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
    fout << rect[i].tag << "\t" << x << "\t" << y << "\t" << theta << endl;
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

void Rect::get_segment_set(bool CW, vector<Vector2D>& my_point_set, 
                           vector<Segment>& my_segment_set,
                           vector<int>& my_point_idx,
                           vector<int>& my_segment_idx) {
  my_point_set.reserve(4);
  my_segment_set.reserve(4);
  my_point_idx.reserve(4);
  my_segment_idx.reserve(4);
  if (CW) {
    for (int i = 0; i < 4; i++) {
      int j = i == 0 ? 3 : i - 1;
      Vector2D OM(vertex[i] - center, Rab);
      double d = i == 0 || i == 2 ? b : a;
      Vector2D ON((vertex[i] + vertex[j]) * 0.5 - center, d);
      my_point_set.push_back(OM);
      my_point_idx.push_back(i);
      my_segment_set.push_back(Segment(OM, ON, false));
      my_segment_idx.push_back(j);
    }
  } else {
    for (int i = 0; i < 4; i++) {
      int j = i == 0 ? 3 : i - 1;
      double d = i == 0 || i == 2 ? b : a;
      Vector2D OM((vertex[i] + vertex[j]) * 0.5 - center, d);
      Vector2D ON(vertex[j] - center, Rab);
      my_point_set.push_back(ON);
      my_point_idx.push_back(j);
      my_segment_set.push_back(Segment(OM, ON, true));
      my_segment_idx.push_back(j);
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

void dis_point_edge(double &d, int &contact_vertex, int &contact_edge, 
                    const double *X, const double *Y, int size,
                    int im, double LY) {
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
  contact_edge = i2 == ipre ? i2 : im;
  double k = (X[i2] - X[im]) / (Y[i2] - Y[im]);
  if (k >= 0) {
    d = k * (-Y[im]) + X[im];
  } else {
    d = k * (LY - Y[im]) + X[im];
    contact_vertex += 1;
    if (contact_vertex == 4)
      contact_vertex = 0;
  }
}
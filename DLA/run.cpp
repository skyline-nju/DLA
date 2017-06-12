#include "run.h"
#include "particle.h"

using namespace std;

double sigma = 2;
double rp = 1;
double Lmin = 1;
double Rmax = 0;
double Rmax2 = 0;
double R_release = Rmax + 10;
double Rkill = R_release * 5;

void launch(Disk &p, Ran *myran) {
  double theta = myran->doub() * 2.0 * PI;
  p.x = R_release * cos(theta);
  p.y = R_release * sin(theta);
}

void collision(Disk &p0, const Disk & p1, double ux, double uy,
               double l, double &l_hit) {
  double dx = p0.x - p1.x;
  double dy = p0.y - p1.y;
  double B = 2 * (ux * dx + uy * dy);
  double C = dx * dx + dy * dy - sigma * sigma;
  double Delta = B * B - 4 * C;
  if (Delta < 0) {
    l_hit = -1;
  } else {
    l_hit = 0.5 * (-B + sqrt(Delta));
    if (l_hit > l || l_hit < 0) {
      l_hit = -1;
    }
  }
}

bool one_step(Disk &p0, std::vector<Disk> &cluster, Ran *myran, Grid &grid) {
  bool flag;
  int col0 = int(p0.x - grid.x_left);
  int row0 = int(p0.y - grid.y_lower);
  int d_wc = grid.dis[col0 + row0 * grid.ncols];
  double theta = myran->doub() * 2.0 * PI;
  double ux = cos(theta);
  double uy = sin(theta);
  if (d_wc <= 4) {
    double l_hit_min = Lmin;
    for (int row = row0 - 4; row <= row0 + 4; row++) {
      int tmp = row * grid.ncols;
      for (int col = col0 - 4; col <= col0 + 4; col++) {
        int idx = grid.cluster[col + tmp];
        if (idx > 0) {
          double l_hit;
          collision(p0, cluster[idx], ux, uy, Lmin, l_hit); 
          if (l_hit >= 0 && l_hit < l_hit_min) {
            l_hit_min = l_hit;
          }
        }
      }
    }
    p0.x += l_hit_min * ux;
    p0.y += l_hit_min * uy;
    if (l_hit_min < Lmin) {
      flag = true;
      cluster.push_back(p0);
      grid.update(col0, row0);
    } else {
      flag = false;
    }
  } else {
    double l = d_wc - 4;
    p0.x += l * ux;
    p0.y += l * uy;
    flag = false;
  }
  return flag;
}

void run(vector<Disk> &cluster, int nPar, Ran *myran, Grid &grid) {
  while (cluster.size() < nPar) {
    Disk p0;
    launch(p0, myran);
    while (true) {
      double r2 = p0.x * p0.x + p0.y * p0.y;
      if (one_step(p0, cluster, myran, grid)) {
        if (r2 > Rmax2) {
          Rmax2 = r2;
          Rmax = sqrt(Rmax2);
          R_release = Rmax + 5;
          Rkill = R_release * 5;
        }
        break;
      } else if (r2 > Rkill * Rkill){
        break;
      }
    }
  }
}

void cal_fractal_dim(const std::vector<Disk>& cluster) {
  double r = 1;
  int size = cluster.size();
  while (r < Rmax) {
    int count = 0;
    double r2 = r * r;
    for (int i = 0; i < size; i++) {
      if (cluster[i].x * cluster[i].x + cluster[i].y * cluster[i].y <= r2) {
        count++;
       }
    }
    cout << r << "\t" << count << endl;
    r *= 1.2;
  }
}

void output(const vector<Disk> &cluster) {
  char fname[100];
  snprintf(fname, 100, "N%d.dat", cluster.size());
  ofstream fout(fname);
  for (int i = 0; i < cluster.size(); i++)
    fout << cluster[i].x << "\t" << cluster[i].y << endl;
}

void output_xyz(const vector<Disk> &cluster) {
  char fname[100];
  snprintf(fname, 100, "N%d.xyz", cluster.size());
  ofstream fout(fname);
  fout << cluster.size() << endl;
  fout << "type\tx\ty\n";
  for (int i = 0; i < cluster.size(); i++)
    fout << "X\t" << cluster[i].x << "\t" << cluster[i].y << endl;
}

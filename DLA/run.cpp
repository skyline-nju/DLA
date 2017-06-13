#include "run.h"
#include "particle.h"

using namespace std;

double sigma = 2;
double rp = 1;
double Lmin = 1;
double Rmax = 0;
double Rmax2 = 0;
double R_release = Rmax + 5;
double Rkill = R_release * 5;

void launch(Disk &p, Ran *myran) {
  double theta = myran->doub() * 2.0 * PI;
  p.x = R_release * cos(theta);
  p.y = R_release * sin(theta);
  //cout << p.x << "\t" << p.y << "\t" << R_release << endl;
}

void collision(Disk &p0, const Disk & p1, double ux, double uy,
               double l, double &l_hit, bool &flag) {
  double dx = p0.x - p1.x;
  double dy = p0.y - p1.y;
  double B = 2 * (ux * dx + uy * dy);
  double C = dx * dx + dy * dy - sigma * sigma;
  double Delta = B * B - 4 * C;
  if (Delta < 0) {
    flag = false;
    //cout << "Delta = " << Delta << endl;
  } else {
    double sqrt_Delta = sqrt(Delta);
    l_hit = B + sqrt_Delta < 0 ? 
            0.5 * (-B - sqrt_Delta) : 0.5 * (-B + sqrt_Delta);
    if (l_hit > l || l_hit < 0) {
      flag = false;
    } else {
      flag = true;
    }
    //cout << "l_hit = " << l_hit << "\t" << p0.x << "\t" << p0.y << endl;
  }
}

bool one_step(Disk &p0, vector<Disk> &cluster, Ran *myran, Grid &grid) {
  bool flag_stick;
  int col0 = int(p0.x - grid.x_left);
  int row0 = int(p0.y - grid.y_lower);
  int d_wc = grid.get_dis(col0, row0);
  double theta = myran->doub() * 2.0 * PI;
  double ux = cos(theta);
  double uy = sin(theta);
  double step_size;
  if (d_wc <= 4) {
    int count_collision = 0;
    double l_hit_min;
    for (int row = row0 - 4; row <= row0 + 4; row++) {
      for (int col = col0 - 4; col <= col0 + 4; col++) {
        int tag = grid.get_tag(col, row);
        if (tag > 0) {
          bool flag;
          double l_hit;
          collision(p0, cluster[tag - 1], ux, uy, Lmin, l_hit, flag);
          if (flag) {
            if (count_collision && l_hit_min > l_hit) {
              l_hit_min = l_hit;
            } else {
              l_hit_min = l_hit;
            }
            count_collision++;
          }
        }
      }
    }
    if (count_collision > 0) {
      //cout << "count = " << count_collision << endl;
      step_size = l_hit_min;
      flag_stick = true;
    } else {
      step_size = Lmin;
      flag_stick = false;
    }
  } else {
    step_size = d_wc - 4;
    flag_stick = false;
  }
  //cout << d_wc << "\tstep size = " << step_size << "\t" << p0.get_rr(cluster[0]);
  p0.x += step_size * ux;
  p0.y += step_size * uy;
  //cout << "\t" << p0.get_rr(cluster[0]) << endl;
  if (flag_stick) {
    cluster.push_back(p0);
    grid.update(p0.x, p0.y);
  }
  return flag_stick;
}

void run(vector<Disk> &cluster, int nPar, Ran *myran, Grid &grid) {
  while (cluster.size() < nPar) {
    Disk p0;
    launch(p0, myran);
    while (true) {
      bool flag = one_step(p0, cluster, myran, grid);
      double r2 = p0.x * p0.x + p0.y * p0.y;
      if (flag) {
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

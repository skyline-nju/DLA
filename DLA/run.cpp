#include "run.h"
#include "particle.h"

using namespace std;

double Lmin = 1;
double Rmax = 0;
double Rmax2 = 0;
double R_release = Rmax + 18;
double Rkill = R_release * 5;

void launch(Disk &p, Ran *myran) {
  double theta = myran->doub() * 2.0 * PI;
  p.x = R_release * cos(theta);
  p.y = R_release * sin(theta);
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
          p0.collide(cluster[tag - 1], ux, uy, Lmin, l_hit, flag);
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
  p0.x += step_size * ux;
  p0.y += step_size * uy;
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

void launch(Rect &p, Ran *myran) {
  double theta = myran->doub() * 2.0 * PI;
  double angle = myran->doub() * 360;
  p = Rect(R_release * cos(theta), R_release * sin(theta), angle);
}

bool one_step(Rect & p, std::vector<Rect>& cluster, Ran * myran) {
  bool flag_stick;
  int idx0 = int(myran->doub() * 4);
  Vec2<double> u;
  p.get_mov_dir(idx0, u);
  double l = Lmin;
  int count = 0;
  for (int i = 0; i < cluster.size(); i++) {
    bool flag_collide = false;
    double l_hit;
    //p.collide(idx0, u, cluster[i], Lmin, l_hit, flag_collide);
    if (idx0 == 0 || idx0 == 2)
      p.collide1(idx0, u, cluster[i], Lmin, l_hit, flag_collide);
    else
      p.collide2(idx0, u, cluster[i], Lmin, l_hit, flag_collide);
    if (flag_collide) {
      if (l > l_hit) {
        l = l_hit;
      }
      count++;
    }
  }
  p.center += l * u;
  p.cal_vertex();
  if (count > 0) {
    flag_stick = true;
    cluster.push_back(p);
  } else {
    flag_stick = false;
  }
  return flag_stick;
}

void run(std::vector<Rect>& cluster, int nPar, Ran * myran) {
  cluster.push_back(Rect(0, 0, 0));
  vector<Rect> traj;
  while (cluster.size() < nPar) {
    cout << cluster.size() << endl;
    Rect p0;
    launch(p0, myran);
    while (true) {
      bool flag = one_step(p0, cluster, myran);
      double r2 = p0.center.x * p0.center.x + p0.center.y * p0.center.y;
      if (flag) {
        if (r2 > Rmax2) {
          Rmax2 = r2;
          Rmax = sqrt(Rmax2);
          R_release = Rmax + 18;
          Rkill = R_release * 5;
        }
        break;
      } else if (r2 > Rkill * Rkill) {
        break;
      }
    }
  }
  Rect::output(traj, "traj.dat");
}

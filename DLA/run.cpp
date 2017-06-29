#include "run.h"
#include "disk.h"
#include "rect.h"
#include "rotate.h"

using namespace std;

double Lmin = 1;
double Rmax = 0;
double Rmax2 = 0;
double R_release = Rmax + 18;
double Rkill = R_release * 4;
double alpha = 7;
double theta_m = 6 / 49.0 * log(0.5 * alpha) / log(alpha);
// Da / Dr = 6 / La * La * log(alpha / 2) / log(alpha)

void launch(Disk &p, Ran *myran) {
  double theta = myran->doub() * 2.0 * PI;
  p.x = R_release * cos(theta);
  p.y = R_release * sin(theta);
}

bool one_step(Disk &p0, vector<Disk> &cluster, Ran *myran, Grid &grid) {
  Status status(Lmin);
  p0.move(cluster, myran, grid, status);
  if (status.flag)
    p0.slip(cluster, myran, grid, status);
  return status.flag;
}

void run(vector<Disk> &cluster, int nPar, Ran *myran, Grid &grid) {
  while (cluster.size() < nPar) {
    Disk p0;
    launch(p0, myran);
    while (true) {
      if (one_step(p0, cluster, myran, grid)) {
        cluster.push_back(p0);
        grid.update(p0.x, p0.y);
        double r2 = p0.x * p0.x + p0.y * p0.y;
        if (cluster.size() % 50000 == 0)
          cout << "N = " << cluster.size() << endl;
        if (r2 > Rmax2) {
          Rmax2 = r2;
          Rmax = sqrt(Rmax2);
          R_release = Rmax + 4;
          Rkill = R_release * 4;
        }
        break;
      } else if (p0.x * p0.x + p0.y * p0.y > Rkill * Rkill){
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
  double angle = myran->doub() * 2.0 * PI;
  p.center.x = R_release * cos(theta);
  p.center.y = R_release * sin(theta);
  p.orient.x = cos(angle);
  p.orient.y = sin(angle);
  p.cal_vertex();
}


bool one_step(Rect &p, vector<Rect> &cluster, Cell<Rect> &cell, Ran *myran) {
  int idx0 = int(myran->doub() * 6);
  bool collided = false;
  if (idx0 == 0 || idx0 == 2) {
    p.translate(cluster, cell, Lmin, idx0, collided);
  } else if (idx0 == 1 || idx0 == 3){
    p.translate(cluster, cell, Lmin * 0.5, idx0, collided);
  } else if (idx0 == 4) {
    p.rotate(cluster, cell, theta_m, true, collided);
  } else {
    p.rotate(cluster, cell, theta_m, false, collided);
  }
  return collided;
}

void run(vector<Rect>& cluster, int nPar, Cell<Rect> & cell, Ran * myran, 
         ofstream &fout) {
  cluster.push_back(Rect(0, 0, 0, 0));
  cell.update(&cluster[0]);
  cout << "initialized" << endl;
  while (cluster.size() < nPar) {
    Rect p0(cluster.size());
    launch(p0, myran);
    while (true) {
      if (one_step(p0, cluster, cell, myran)) {
        cluster.push_back(p0);
        cell.update(&cluster[cluster.size() - 1]);
        p0.output(fout);
        cout << "n = " << cluster.size() << endl;
        double rr = p0.center.square();
        if (rr > Rmax2) {
          Rmax2 = rr;
          Rmax = sqrt(Rmax2);
          R_release = Rmax + 18;
          Rkill = R_release * 4;
        }
        break;
      } else if (p0.center.square() > Rkill * Rkill) {
        break;
      }
    }
  }
}

void run(vector<Rect>& cluster, int nPar, Cell<Rect> & cell, int seed,
         double tilt_theta) {
  char filename[100];
  snprintf(filename, 100, "%d_%g_%d.dat", nPar, tilt_theta, seed);
  Rect::output(cluster, filename);
  ofstream fout(filename);
  Rect::tilt_angle = tilt_theta / 180 * PI;
  Ran myran(seed);
  run(cluster, nPar, cell, &myran, fout);
}


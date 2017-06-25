#include "disk.h"
#include <cmath>


using namespace std;

double Disk::sigma = 2;
double Disk::rp = 1;

void Disk::collide(const Disk & disk, double ux, double uy,
                   int tag, Status & status) const {
  double dx = x - disk.x;
  double dy = y - disk.y;
  double B = 2 * (ux * dx + uy * dy);
  double C = dx * dx + dy * dy - sigma * sigma;
  double Delta = B * B - 4 * C;
  if (Delta >= 0) {
    double sqrt_Delta = sqrt(Delta);
    double l_hit = B + sqrt_Delta < 0 ? 0.5 * (-B - sqrt_Delta) :
                                        0.5 * (-B + sqrt_Delta);
    if (l_hit <= status.l && l_hit > 0 ) {
      status.flag = true;
      status.l = l_hit;
      status.neighbor_tag = tag;
      status.count++;
    }
  }
}

void Disk::move(const vector<Disk> &cluster, Ran * myran, const Grid & grid,
                Status &status) {
  int col0 = int(x - grid.x_left);
  int row0 = int(y - grid.y_lower);
  double theta = myran->doub() * 2.0 * PI;
  double ux = cos(theta);
  double uy = sin(theta);
  if (col0 < 0 || col0 >= grid.ncols || row0 < 0 || row0 >= grid.nrows) {
    status.l = 20;
  } else {
    int d_wc = grid.get_dis(col0, row0);
    if (d_wc <= 4) {
      for (int row = row0 - 4; row <= row0 + 4; row++) {
        for (int col = col0 - 4; col <= col0 + 4; col++) {
          int tag = grid.get_tag(col, row);
          if (tag >= 0) {
            collide(cluster[tag], ux, uy, tag, status);
          }
        }
      }
    } else {
      status.l = d_wc - 4;
    }
  }
  x += status.l * ux;
  y += status.l * uy;
}

void Disk::slip(const vector<Disk>& cluster, Ran * myran, const Grid & grid,
                const Status & status) {
  int i0 = status.neighbor_tag;
  double x0 = cluster[i0].x;
  double y0 = cluster[i0].y;
  double theta0 = atan2(y - y0, x - x0);
  int col0 = grid.get_col(cluster[i0].x);
  int row0 = grid.get_row(cluster[i0].y);
  double flag_CW = false;
  double flag_CCW = false;
  double min_angle_CW;
  double min_angle_CCW;
  for (int row = row0 - 4; row <= row0 + 4; row++) {
    for (int col = col0 - 4; col <= col0 + 4; col++) {
      int i = grid.get_tag(col, row);
      if (i >= 0 && i != i0) {
        double dx = cluster[i].x - x0;
        double dy = cluster[i].y - y0;
        double d_square = dx * dx + dy * dy;
        if (d_square <= 4 * sigma * sigma) {
          double phi = atan2(dy, dx);
          double dphi = acos(sqrt(d_square) / (2 * sigma));

          double angle_CW = theta0 - (phi + dphi);
          if (angle_CW < 0)
            angle_CW += 2 * PI;
          if (flag_CW && min_angle_CW > angle_CW) {
            min_angle_CW = angle_CW;
          } else if (!flag_CW) {
            flag_CW = true;
            min_angle_CW = angle_CW;
          }

          double angle_CCW = (phi - dphi) - theta0;
          if (angle_CCW < 0)
            angle_CCW += 2 * PI;
          if (flag_CCW && min_angle_CCW > angle_CCW) {
            min_angle_CCW = angle_CCW;
          } else if (!flag_CCW) {
            flag_CCW = true;
            min_angle_CCW = angle_CCW;
          }
        }
      }
    }
  }
  //if (flag_CW) {
  //  Vec2<double> OP(x - x0, y - y0);
  //  OP.rotate(-min_angle_CW);
  //  x = x0 + OP.x;
  //  y = y0 + OP.y;
  //}
  if (flag_CCW) {
    Vec2<double> OP(x - x0, y - y0);
    OP.rotate(min_angle_CCW);
    x = x0 + OP.x;
    y = y0 + OP.y;
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
#include "disk.h"
#include <cmath>

using namespace std;

double Disk::sigma = 2;
double Disk::rp = 1;

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
    l_hit = B + sqrt_Delta < 0 ? 0.5 * (-B - sqrt_Delta) : 
                                 0.5 * (-B + sqrt_Delta);
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
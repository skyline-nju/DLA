#ifndef PARTICLE_H
#define PARTICLE_H

struct Disk
{
  Disk() { x = y = 0; }
  Disk(double x0, double y0) :x(x0), y(y0) {}
  double get_rr(const Disk &d2) {
    return (x - d2.x) * (x - d2.x) + (y - d2.y) * (y - d2.y);
  }
  double x;
  double y;
};

struct Rod
{
  Rod() {}
  void collide(const Rod &rod, double vx, double vy, int first_vertex,
               double l, double &lhit, bool &flag);

  void cal_vertex();
  void check_overlape(int first_vertex, double vx, double vy, const Rod &rod);
  double xc;
  double yc;
  double ux;
  double uy;
  double vertex_x[4];
  double vertex_y[4];
  static double a;
  static double b;
};

inline double projection(double x0, double y0, double u1, double u2,
                         double x, double y) {
  return u1 * (x - x0) + u2 * (y - y0);
}


#endif

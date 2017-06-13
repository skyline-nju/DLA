#ifndef GRID_H
#define GRID_H

class Grid
{
public:
  Grid(int n, int dmax);
  ~Grid();
  void update(int coli, int rowi);
  void update(double x, double y) {
    update(int(x - x_left), int(y - y_lower));
  }
  int get_dis(int col, int row) { return dis[col + row * ncols]; }
  int get_tag(int col, int row) { return cluster[col + row * ncols]; }
  void show();
  void show(int m);
  int *cluster;
  int *dis;
  int *vicinity;
  int ncols;
  int nrows;
  int Dmax;
  int cluster_size;
  double x_left;
  double y_lower;
};

void show_grid(int *grid, int n);

void show_grid(int *grid, int n, int m);
#endif

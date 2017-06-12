#ifndef GRID_H
#define GRID_H

class Grid
{
public:
  Grid(int n, int dmax);
  ~Grid();
  void update(int coli, int rowi);
  void show();
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
#endif

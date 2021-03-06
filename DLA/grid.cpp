#include "grid.h"
using namespace std;

void set_vicinity(int *grid, int m) {
  int ncols = 2 * m + 1;
  int nrows = ncols;
  for (int row = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      int idx = col + row * ncols;
      int dx = fabs(col - m);
      int dy = fabs(row - m);
      if (dx == 0 && dy == 0) {
        grid[idx] = 0;
      } else if (dx == 0) {
        grid[idx] = dy;
      } else if (dy == 0) {
        grid[idx] = dx;
      } else {
        grid[idx] = floor(sqrt(dx * dx + dy * dy));
      }
      if (grid[idx] > m)
        grid[idx] = m;
    }
  }
}

Grid::Grid(int n, int dmax): ncols(n), nrows(n), Dmax(dmax) {
  cluster = new int[n * n];
  cluster_size = 0;
  dis = new unsigned char[n * n];
  vicinity = new int[(2 * Dmax + 1) * (2 * Dmax + 1)];
  cout << sizeof(cluster[0]) * n * n  / 1024 / 1024<< endl;
  cout << sizeof(dis[0]) * n * n / 1024 / 1024 << endl;
  set_vicinity(vicinity, Dmax);
  for (int i = 0; i < ncols * nrows; i++) {
    cluster[i] = -1;
    dis[i] = Dmax;
  }
  int col0, row0;
  if (n % 2 == 0) {
    col0 = row0 = n / 2;
    x_left = y_lower = -n / 2;
  } else {
    col0 = row0 =  (n - 1) / 2;
    x_left = y_lower = -(n - 1) / 2;
  }
  dis[col0 + ncols * row0] = 0;
  update(col0, row0);
}


Grid::~Grid() {
  delete[] cluster;
  delete[] dis;
  delete[] vicinity;
}

void Grid::update(int coli, int rowi) {
  int idxi = coli + rowi * ncols;
  cluster_size++;
  //cout << "[" << coli << ", " << rowi << "] = " << cluster_size << endl;
  cluster[idxi] = cluster_size - 1;

  int k = 0;
  for (int row = rowi - Dmax; row <= rowi + Dmax; row++) {
    int tmp = row * ncols;
    for (int col = coli - Dmax; col <= coli + Dmax; col++) {
      int idx = col + tmp;
      if (vicinity[k] < dis[idx])
        dis[idx] = vicinity[k];
      k++;
    }
  }
}

void Grid::show() {
  cout << "vicinity: " << endl;
  show_grid(vicinity, 2 * Dmax + 1);
  cout << "distance: " << endl;
  show_grid(dis, ncols);
  cout << "cluster: " << endl;
  show_grid(cluster, ncols);
}

void Grid::show(int m) {
  //cout << "vicinity: " << endl;
  //show_grid(vicinity, 2 * Dmax + 1);
  cout << "distance: " << endl;
  show_grid(dis, ncols, m);
  cout << "cluster: " << endl;
  show_grid(cluster, ncols, m);
}




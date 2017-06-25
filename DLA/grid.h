#ifndef GRID_H
#define GRID_H
#include <list>
#include <vector>

class Grid
{
public:
  Grid(int n, int dmax);
  ~Grid();
  void update(int coli, int rowi);
  void update(double x, double y);
  int get_dis(int col, int row) const { return dis[col + row * ncols]; }
  int get_tag(int col, int row) const { return cluster[col + row * ncols]; }
  int get_col(double x) const { return int(x - x_left); }
  int get_row(double y) const { return int(y - y_lower); }
  void show();
  void show(int m);

  int *cluster;
  unsigned char *dis;
  int *vicinity;
  int ncols;
  int nrows;
  int Dmax;
  int cluster_size;
  double x_left;
  double y_lower;
};

inline void Grid::update(double x, double y) {
  update(int(x - x_left), int(y - y_lower));
}

class Cell
{
public:
  Cell(int n, double l0);
  int get_col(double x) const { return int((x - x_left) * one_over_l); }
  int get_row(double y) const { return int((y - y_lower) * one_over_l); }
  int get_idx(int col, int row) const { return col + ncols * row; }
  int get_idx(double x, double y) const;
  bool is_isolate(int col, int row) const;
  bool is_isolate(double x, double y) const;
  void update(int coli, int rowi);
  void update(double x, double y);
  template <class T, class UnaryFunc>
  void for_each_neighbor(int colc, int rowc, const std::vector<T> &vec,
                         UnaryFunc f) const;


  std::vector<std::list<int>> tag;
  std::vector<bool> isolate;
  int ncols;
  int nrows;
  int ncell;
  int particle_num;
  double x_left;
  double y_lower;
  double l;
  double one_over_l;
};

inline int Cell::get_idx(double x, double y) const {
  return get_col(x) + ncols * get_row(y);
}

inline bool Cell::is_isolate(int col, int row) const {
  return isolate[col + ncols * row];
}

inline bool Cell::is_isolate(double x, double y) const {
  return isolate[get_idx(x, y)];
}

inline void Cell::update(double x, double y) {
  update(get_col(x), get_row(y));
}

template <class T, class UnaryFunc>
void Cell::for_each_neighbor(int colc, int rowc, const std::vector<T> &vec,
                             UnaryFunc f) const {
  for (int row = rowc - 1; row <= rowc + 1; row++) {
    int tmp = row * ncols;
    for (int col = colc - 1; col <= colc + 1; col++) {
      int idx = col + tmp;
      for (auto iter = tag[idx].cbegin(); iter != tag[idx].cend(); ++iter) {
        f(vec[*iter]);
      }
    }
  }
}

template <typename T>
void show_grid(T *grid, int n) {
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      int idx = col + row * n;
      cout << grid[idx] << " ";
    }
    cout << endl;
  }
}

template <typename T>
void show_grid(T *grid, int n, int m) {
  int col_c = n % 2 == 0 ? n / 2 : (n - 1) / 2;
  int row_c = col_c;
  for (int row = row_c - m; row <= row_c + m; row++) {
    for (int col = col_c - m; col <= col_c + m; col++) {
      int idx = col + row * n;
      cout << grid[idx] << " ";
    }
    cout << endl;
  }
}

#endif

#ifndef GRID_H
#define GRID_H
#include <list>
#include <vector>

void show_grid(int *grid, int n);

void show_grid(int *grid, int n, int m);

class Grid
{
public:
  Grid(int n, int dmax);
  ~Grid();
  void update(int coli, int rowi);
  void update(double x, double y);
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
  void for_each_neighbor(int colc, int rowc, const T &vec, UnaryFunc f) const;


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
void Cell::for_each_neighbor(int colc, int rowc, const T &vec, UnaryFunc f) const {
  for (int row = rowc - 1; row <= rowc + 1; row++) {
    int tmp = row * ncols;
    for (int col = colc - 1; col <= colc + 1; col++) {
      int idx = col + tmp;
      for (auto iter = tag[idx].cbegin(); iter != tag[idx].cend(); ++iter) {
        f(vec[*iter].vertex);
      }
    }
  }
}

#endif

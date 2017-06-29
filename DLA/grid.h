#ifndef GRID_H
#define GRID_H
#include <vector>
#include <iostream>
#include <cmath>

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

template <typename T>
class Cell
{
public:
  Cell(int n, double l0);
  int get_col(double x) const { return int((x - x_left) * one_over_l); }
  int get_row(double y) const { return int((y - y_lower) * one_over_l); }
  int get_idx(int col, int row) const { return col + ncols * row; }
  int get_idx(double x, double y) const;
  bool out_of_range(int col, int row) const;
  bool is_isolate(int col, int row) const { return isolate[get_idx(col, row)]; }
  bool is_isolate(double x, double y) const;
  void update(T *p);
  template <class UnaryFunc>
  void for_each_neighbor(int colc, int rowc, UnaryFunc f) const;
  template <class UnaryFunc>
  void for_each_neighbor(int colc, int rowc, int idx_exclude, UnaryFunc f) const;

  std::vector<T *> cell_list;
  bool *isolate;
  int ncols;
  int nrows;
  int ncell;
  int particle_num;
  double x_left;
  double y_lower;
  double l;
  double one_over_l;
};

template <typename T>
Cell<T>::Cell(int n, double l0) : ncols(n), nrows(n), l(l0), ncell(n*n) {
  cell_list.reserve(ncell);
  isolate = new bool[ncell];
  for (int i = 0; i < ncell; i++) {
    T* p = NULL;
    cell_list.push_back(p);
    isolate[i] = true;
  }
  int col0, row0;
  if (n % 2 == 0) {
    col0 = row0 = n / 2;
    x_left = y_lower = (-n / 2) * l;
  } else {
    col0 = row0 = (n - 1) / 2;
    x_left = y_lower = -(n - 1) / 2 * l;
  }
  particle_num = 0;
  one_over_l = 1 / l;
  std::cout << "xl = " << x_left << std::endl;
  std::cout << "yl = " << y_lower << std::endl;
  std::cout << "l = " << l << "\t1/l = " << one_over_l << std::endl;
  std::cout << "ncols = " << ncols << std::endl;
}

template <class T>
void Cell<T>::update(T *p) {
  int coli = get_col(p->center.x);
  int rowi = get_row(p->center.y);
  int idx = get_idx(coli, rowi);
  for (int row = rowi - 1; row <= rowi + 1; row++) {
    int tmp = row * ncols;
    for (int col = coli - 1; col <= coli + 1; col++) {
      isolate[col + tmp] = false;
    }
  }
  p->next = cell_list[idx];
  cell_list[idx] = p;
  particle_num++;
}

template <typename T>
inline int Cell<T>::get_idx(double x, double y) const {
  return get_col(x) + ncols * get_row(y);
}

template <typename T>
inline bool Cell<T>::out_of_range(int col, int row) const {
  return col < 0 || col >= ncols || row < 0 || row >= nrows;
}

template <typename T>
inline bool Cell<T>::is_isolate(double x, double y) const {
  int col = get_col(x);
  int row = get_row(y);
  return is_isolate(col, row);
}

template <class T>
template <class UnaryFunc>
void Cell<T>::for_each_neighbor(int colc, int rowc, UnaryFunc f) const {
  for (int row = rowc - 1; row <= rowc + 1; row++) {
    int tmp = row * ncols;
    for (int col = colc - 1; col <= colc + 1; col++) {
      int idx = col + tmp;
      if (cell_list[idx]) {
        T* curNode = cell_list[idx];
        do {
          f(*curNode);
          curNode = curNode->next;
        } while (curNode);
      }
    }
  }
}

template <class T>
template<class UnaryFunc>
void Cell<T>::for_each_neighbor(int colc, int rowc, int idx_exclude,
                                UnaryFunc f) const {
  for (int row = rowc - 1; row <= rowc + 1; row++) {
    int tmp = row * ncols;
    for (int col = colc - 1; col <= colc + 1; col++) {
      int idx = col + tmp;
      if (cell_list[idx]) {
        T* curNode = cell_list[idx];
        do {
          if (curNode->tag != idx_exclude)
            f(*curNode);
          curNode = curNode->next;
        } while (curNode);
      }
    }
  }
}

template <typename T>
void show_grid(T *grid, int n) {
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      int idx = col + row * n;
      std::cout << grid[idx] << " ";
    }
    std::cout << "\n";
  }
}

template <typename T>
void show_grid(T *grid, int n, int m) {
  int col_c = n % 2 == 0 ? n / 2 : (n - 1) / 2;
  int row_c = col_c;
  for (int row = row_c - m; row <= row_c + m; row++) {
    for (int col = col_c - m; col <= col_c + m; col++) {
      int idx = col + row * n;
      std::cout << grid[idx] << " ";
    }
    std::cout << "\n";
  }
}

#endif

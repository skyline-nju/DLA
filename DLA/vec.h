#ifndef VEC_H
#define VEC_H
#include <cmath>
#include <iostream>

template <typename T>
struct Vec2
{
  T x;
  T y;

  Vec2() { x = y = 0; }
  Vec2(T a) { x = y = a; }
  Vec2(T a, T b) : x(a), y(b) {}
  Vec2(const Vec2<T> &a) :x(a.x), y(a.y) {}
  
  bool operator == (const Vec2<T> &a) const { return x == a.x && y == a.y; }
  bool operator != (const Vec2<T> &a) const { return x != a.x || y != a.y; }
  Vec2<T> operator -() const { return Vec2<T>(-x, -y); }

  template <typename T2>
  Vec2<T> operator +(const Vec2<T2> &a) const { return Vec2<T>(x + a.x, y + a.y); }
  Vec2<T> operator +(double a) const { return Vec2<T>(x + a, y + a); }
  Vec2<T> operator +(int a) const { return Vec2<T>(x + a, y + a); }
  friend Vec2<T> operator +(double a, const Vec2<T> &b) { return b + a; }
  friend Vec2<T> operator +(int a, const Vec2<T> &b) { return b + a; }

  template <typename T2>
  Vec2<T> operator -(const Vec2<T2> &a) const { return Vec2<T>(x - a.x, y - a.y); }
  Vec2<T> operator -(double a) const { return Vec2<T>(x - a, y - a); }
  Vec2<T> operator -(int a) const { return Vec2<T>(x - a, y - a); }
  friend Vec2<T> operator -(double a, const Vec2<T> &b) { return b - a; }
  friend Vec2<T> operator -(int a, const Vec2<T> &b) { return b - a; }


  template <typename T2>
  Vec2<T> operator *(const Vec2<T2> &a) const { return Vec2<T>(x * a.x, y * a.y); }
  Vec2<T> operator *(double a) const { return Vec2<T>(x * a, y * a); }
  Vec2<T> operator *(int a) const { return Vec2<T>(x * a, y * a); }
  friend Vec2<T> operator *(double a, const Vec2<T> &b) { return b * a; }
  friend Vec2<T> operator *(int a, const Vec2<T> &b) { return b * a; }

  Vec2<T> operator /(double a) const { return Vec2<T>(x / a, y / a); }
  Vec2<T> operator /(int a) const { return Vec2<T>(x / a, y / a); }


  template <typename T2>
  void operator +=(const Vec2<T2> &a) { x += a.x; y += a.y; }
  template <typename T2>
  void operator +=(T2 a) { x += a; y += a; }
  template <typename T2>
  void operator -=(const Vec2<T2> &a) { x -= a.x; y -= a.y; }
  template <typename T2>
  void operator -=(T2 a) { x -= a; y -= a; }
  template <typename T2>
  void operator *=(T2 a) { x *= a; y *= a; }
  template <typename T2>
  void operator /=(T2 a) { x /= a; y /= a; }

  template <typename T2>
  double dot(const Vec2<T2> &a) const { return x * a.x + y * a.y; }
  template <typename T2>
  double cross(const Vec2<T2> &a) const { return x * a.y - y * a.x; }
  template <typename T2>
  double distance(const Vec2<T2> &a) const { return sqrt(square_distance(a)); }
  template <typename T2> 
  double square_distance(const Vec2<T2> &a) const;
  double square() const { return x * x + y * y; }
  double module() const { return sqrt(x * x + y * y); }
  void rotate(double theta);
  Vec2<T> get_rotated_vec(double theta) const;

  friend std::ostream& operator <<(std::ostream &output, const Vec2<T> &a) {
    output << a.x << "\t" << a.y;
    return output;
  }
};

template <class T>
template <class T2>
inline double Vec2<T>::square_distance(const Vec2<T2> &a) const {
  double dx = x - a.x;
  double dy = y - a.y;
  return dx * dx + dy * dy;
}

template <class T>
inline void Vec2<T>::rotate(double theta) {
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double u = x * cos_theta - y * sin_theta;
  double v = x * sin_theta + y * cos_theta;
  x = u;
  y = v;
}

template <class T>
inline  Vec2<T> Vec2<T>::get_rotated_vec(double theta) const {
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double u = x * cos_theta - y * sin_theta;
  double v = x * sin_theta + y * cos_theta;
  return Vec2<T>(u, v);
}
#endif
